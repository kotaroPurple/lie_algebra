
import numpy as np
import matplotlib.pyplot as plt

SHOW_TRAJECTORY = True


def make_data(number: int, r: float = 1.0):
    """
    極座標で3次元データを作る
    x = r.sin_theta.cos_phi
    y = r.sin_theta.sin_phi
    z = r.cos_theta
    w = 1 (同次座標)
    """
    # prepare
    number_sq = int(np.sqrt(number))
    theta_min, theta_max = np.deg2rad(-20), np.deg2rad(20)
    phi_min, phi_max = np.deg2rad(-180), np.deg2rad(180)
    thetas = np.linspace(theta_min, theta_max, number_sq)
    phis = np.linspace(phi_min, phi_max, number_sq)
    thetas, phis = np.meshgrid(thetas, phis)
    thetas = thetas.flatten()
    phis = phis.flatten()
    # make x,y,z
    xs = r * np.sin(thetas) * np.cos(phis)
    ys = r * np.sin(thetas) * np.sin(phis)
    zs = r * np.cos(thetas)
    ws = np.ones_like(xs)
    xyz = np.c_[xs, ys, zs, ws]
    return xyz


def move_data(data: np.ndarray, transformation_matrix: np.ndarray):
    """
    data: (N,4) # 同次座標
    transformation_matrix: (4,4)
    """
    result = np.dot(data, transformation_matrix.T)
    return result


def transformation_matrix_to_vector(trans_matrix: np.ndarray):
    """
    変換行列をベクトルに変換する
    """
    # 回転ベクトル
    R = trans_matrix[:3, :3]
    rot_vector, theta = _rotation_matrix_to_vector(R)
    # 並進ベクトル
    trans = trans_matrix[:3, 3]
    rho = translation_vec_to_lie_translation(trans, rot_vector, theta)
    # connect
    result = np.r_[theta * rot_vector, rho]
    return result


def translation_vec_to_lie_translation(in_trans: np.ndarray, rot_vector: np.ndarray, theta: float):
    """
    いわゆる並進ベクトルを lie algebra での並進ベクトルに変換する
    """
    inv_J = calculate_inverse_J_matrix(rot_vector, theta)
    rho = np.dot(inv_J, in_trans)
    return rho


def calculate_J_matrix(vec: np.ndarray, theta: float):
    """
    J = (sin_theta / theta).I + (1-sin_theta / theta).l.l.T + (1-cos_theta) / theta l^
    """
    J = np.sin(theta) / theta * np.eye(3)
    J += (1. - np.sin(theta)) / theta * np.dot(vec[:, np.newaxis], vec[np.newaxis, :])
    J += (1. - np.cos(theta)) / theta * rot_vector_to_skew_symmetric_matrix(vec)
    return J


def calculate_inverse_J_matrix(vec: np.ndarray, theta: float):
    """
    J_inv = (t/2 / tan(t/2)).I + (1- (t/2 / tan(t/2))).l.l.T - t/2.l^
    """
    J_inv = (theta / 2) / np.tan(theta/2) * np.eye(3)
    J_inv += (1-theta/2/np.tan(theta/2)) * np.dot(vec[:, np.newaxis], vec[np.newaxis, :])
    J_inv -= theta/2 * rot_vector_to_skew_symmetric_matrix(vec)
    return J_inv


def rotation_matrix_to_vector(rot_matrix: np.ndarray):
    """
    回転行列を回転ベクトルに変換する
    """
    rot_vector, theta = _rotation_matrix_to_vector(rot_matrix)
    return theta * rot_vector


def _rotation_matrix_to_vector(rot_matrix: np.ndarray):
    rot_vector = -np.array([rot_matrix[1, 2] - rot_matrix[2, 1],
                            rot_matrix[2, 0] - rot_matrix[0, 2],
                            rot_matrix[0, 1] - rot_matrix[1, 0]])
    rot_vector = normalize_vector(rot_vector)
    trace_r = np.trace(rot_matrix)
    theta = np.arccos((trace_r - 1) / 2.)
    return rot_vector, theta


def vector_to_transformation_matrix(vec: np.ndarray):
    """
    回転ベクトル, 並進ベクトル (6次元) から変換行列を求める
    """
    # prepare
    rot_vector, trans_vector = vec[:3], vec[3:]
    normal_vector = normalize_vector(rot_vector)
    theta = np.sqrt(np.sum(rot_vector**2))
    # rotation matrix
    rot_matrix = np.cos(theta) * np.eye(3)
    rot_matrix += (1 - np.cos(theta)) \
        * (np.dot(normal_vector[:, np.newaxis], normal_vector[np.newaxis, :]))
    rot_matrix += np.sin(theta) * rot_vector_to_skew_symmetric_matrix(normal_vector)
    # translation
    J = calculate_J_matrix(normal_vector, theta)
    trans_for_T = np.dot(J, trans_vector)
    # transfromation matrix
    T = np.eye(4)
    T[:3, :3] = rot_matrix
    T[:3, 3] = trans_for_T
    return T


def rot_vector_to_skew_symmetric_matrix(vec: np.ndarray):
    """
    回転ベクトルを反対称行列にする (-A.T = A)
    """
    result = np.zeros((3, 3))
    result[0, 1] = -vec[2]
    result[0, 2] = vec[1]
    result[1, 0] = vec[2]
    result[1, 2] = -vec[0]
    result[2, 0] = -vec[1]
    result[2, 1] = vec[0]
    return result


def normalize_vector(vec: np.ndarray):
    """
    ベクトルを正規化する
    """
    normal_vector = vec / np.linalg.norm(vec)
    return normal_vector


def calculate_cost(dataA: np.ndarray, dataB: np.ndarray, t_matrix: np.ndarray):
    """
    二乗誤差を求める
    cost = 1/N * Sum_i |T.Ai - Bi|^2
    """
    predicted_B = move_data(dataA, t_matrix)
    cost = np.mean((predicted_B - dataB)**2)
    return cost


def calculate_cost_gradient(dataA: np.ndarray, dataB: np.ndarray, transform_a2b: np.ndarray):
    """
    二乗誤差の微分を求める
    cost = 1/N * Sum_i |T.Ai - Bi|^2
    grad_theta = -Sum_i (R.Ai + t) x (Bi)
    grad_trans =  Sum_i (R.Ai + t - Bi)
    """
    predicted_B = move_data(dataA, transform_a2b)  # (N,4)
    cross_result = np.cross(predicted_B[:, :3], dataB[:, :3])  # (N,3)
    grad_theta = -np.sum(cross_result, axis=0)  # (3,)
    grad_trans = np.sum(predicted_B[:, :3] - dataB[:, :3], axis=0)  # (3,)
    grad = np.r_[grad_theta, grad_trans]
    return grad


def calculate_center(data: np.ndarray):
    center = np.mean(data[:, :3], axis=0)
    return center


def main():
    # make data
    radius = 1.0
    number = 100
    xyz_original = make_data(number, r=radius)
    # correct transformation matrix
    rot_vector = normalize_vector(np.array([1, 1, 0]))
    rot_angle = np.deg2rad(90)
    trans_vector = np.array([0.2, 0.2, 0.2])
    norm_trans = np.linalg.norm(trans_vector)
    trans_vector_lie = translation_vec_to_lie_translation(trans_vector, rot_vector, rot_angle)
    correct_t_vector = np.r_[rot_angle * rot_vector, trans_vector_lie]
    correct_T = vector_to_transformation_matrix(correct_t_vector)
    xyz_moved = move_data(xyz_original, correct_T)

    # predict transformation matrix
    init_rot_vector = normalize_vector(np.array([0, 0, 1]))
    init_angle = np.deg2rad(5)
    init_trans_vector = np.array([0.3, 0.3, 0.3])
    init_trans_lie = translation_vec_to_lie_translation(
        init_trans_vector, init_rot_vector, init_angle)
    init_vector = np.r_[init_angle * init_rot_vector, init_trans_lie]
    init_matrix = vector_to_transformation_matrix(init_vector)
    current_T = init_matrix.copy()
    updated_t_vec = init_vector.copy()
    lr = 0.001
    scores = list()
    store_t_vector = list()
    store_center = list()

    # main
    iterations = 2000
    for _ in range(iterations):
        # calculate score, gradient
        _score = calculate_cost(xyz_original, xyz_moved, current_T)
        grad = calculate_cost_gradient(xyz_original, xyz_moved, current_T)
        # update transformation matrix
        updated_t_vec -= lr * grad
        current_T = vector_to_transformation_matrix(updated_t_vec)
        # store data
        scores.append(_score)
        store_t_vector.append(updated_t_vec.copy())
        # calculate center
        tmp_moved = move_data(xyz_original, current_T)
        center = calculate_center(tmp_moved)
        store_center.append(center)

    # predict data
    predicted_xyz = move_data(xyz_original, current_T)

    # to np.ndarray
    store_t_vector = np.array(store_t_vector)
    store_center = np.array(store_center)

    #
    # show result
    #

    # show score vs iteration
    fig = plt.figure(figsize=(12, 8))
    ax1 = fig.add_subplot(221)
    ax1.plot(scores)
    ax1.set_title('score')

    # show 3d data (original, correct data, predicted data)
    ax2 = fig.add_subplot(222, projection='3d')
    ax2.scatter(
        xyz_original[:, 0], xyz_original[:, 1], xyz_original[:, 2],
        color='purple', label='original', s=5)
    ax2.scatter(
        xyz_moved[:, 0], xyz_moved[:, 1], xyz_moved[:, 2],
        color='blue', label='moved', s=5)
    ax2.scatter(
        predicted_xyz[:, 0], predicted_xyz[:, 1], predicted_xyz[:, 2],
        color='red', label='predicted', s=5)
    ax2.plot(store_center[:, 0], store_center[:, 1], store_center[:, 2], color='orange')
    ax2.set_xlim3d(-radius-norm_trans, radius+norm_trans)
    ax2.set_ylim3d(-radius-norm_trans, radius+norm_trans)
    ax2.set_zlim3d(-radius-norm_trans, radius+norm_trans)
    ax2.set_title('xyz data')
    plt.legend()

    # rotation vectors, translation vectors (3d)
    ax3 = fig.add_subplot(223, projection='3d')
    if SHOW_TRAJECTORY:
        _p = ax3.scatter(
            store_t_vector[:, 0], store_t_vector[:, 1], store_t_vector[:, 2],
            c=np.arange(iterations) / iterations, cmap='rainbow', s=5)
    else:
        _p = ax3.scatter(
            store_t_vector[:, 0], store_t_vector[:, 1], store_t_vector[:, 2],
            c=scores, cmap='gist_rainbow', s=5)
    ax3.scatter(correct_t_vector[0], correct_t_vector[1], correct_t_vector[2], color='red')
    ax3.set_title('rotation vector')
    fig.colorbar(_p)

    ax4 = fig.add_subplot(224, projection='3d')
    if SHOW_TRAJECTORY:
        _p = ax4.scatter(
            store_t_vector[:, 3], store_t_vector[:, 4], store_t_vector[:, 5],
            c=np.arange(iterations) / iterations, cmap='rainbow', s=5)
    else:
        _p = ax4.scatter(
            store_t_vector[:, 3], store_t_vector[:, 4], store_t_vector[:, 5],
            c=scores, cmap='gist_rainbow', s=5)
    ax4.scatter(correct_t_vector[3], correct_t_vector[4], correct_t_vector[5], color='red')
    ax4.set_title('translation vector')
    fig.colorbar(_p)

    print('---- correct transformation vector ----')
    print(correct_t_vector)
    print('---- predicted transformation vector ----')
    print(updated_t_vec)

    print()
    print('---- correct transformation matrix ----')
    print(correct_T)
    print('---- predicted transformation matrix ----')
    print(vector_to_transformation_matrix(updated_t_vec))

    plt.show()


if __name__ == '__main__':
    main()
