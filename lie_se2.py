
import numpy as np
import matplotlib.pyplot as plt

SHOW_TRAJECTORY = True


def make_data(number: int, r: float = 1.0):
    """
    極座標で2次元データを作る
    x = r.cos_phi
    y = r.sin_phi
    w = 1 (同次座標)
    """
    # prepare
    phi_min, phi_max = np.deg2rad(-30), np.deg2rad(30)
    phis = np.linspace(phi_min, phi_max, number)
    # make x,y,w
    xs = r * np.cos(phis)
    ys = r * np.sin(phis)
    ws = np.ones_like(xs)
    xyw = np.c_[xs, ys, ws]
    return xyw


def move_data(data: np.ndarray, transformation_matrix: np.ndarray):
    """
    data: (N,3) # 同次座標
    transformation_matrix: (3,3)
    """
    result = np.dot(data, transformation_matrix.T)
    return result


def transformation_matrix_to_vector(trans_matrix: np.ndarray):
    """
    変換行列をベクトルに変換する
    """
    # 回転ベクトル
    R = trans_matrix[:2, :2]
    theta = _rotation_matrix_to_vector(R)
    # 並進ベクトル
    trans = trans_matrix[:2, 2]
    rho = translation_vec_to_lie_translation(trans, theta)
    # connect
    result = np.r_[theta, rho]
    return result


def translation_vec_to_lie_translation(in_trans: np.ndarray, theta: float):
    """
    いわゆる並進ベクトルを lie algebra での並進ベクトルに変換する
    """
    inv_J = calculate_inverse_J_matrix(theta)
    rho = np.dot(inv_J, in_trans)
    return rho


def calculate_J_matrix(theta: float):
    """
    J = (sin_theta / theta).I + (1-cos_theta) / theta . A
    A = ((0, -1), (1, 0))
    """
    J = np.sin(theta) / theta * np.eye(2)
    J += (1. - np.cos(theta)) / theta * get_skew_symmetric_matrix()
    return J


def calculate_inverse_J_matrix(theta: float):
    """
    J_inv = (t/2 / tan(t/2)).I - (t/2).A
    """
    theta_half = theta / 2
    J_inv = theta_half / np.tan(theta_half) * np.eye(2)
    J_inv -= theta_half * get_skew_symmetric_matrix()
    return J_inv


def rotation_matrix_to_vector(rot_matrix: np.ndarray):
    """
    回転行列を回転ベクトルに変換する
    """
    theta = _rotation_matrix_to_vector(rot_matrix)
    return theta


def _rotation_matrix_to_vector(rot_matrix: np.ndarray):
    trace_r = np.trace(rot_matrix)
    theta = np.arccos(trace_r / 2.)
    return theta


def vector_to_transformation_matrix(vec: np.ndarray):
    """
    回転ベクトル, 並進ベクトル (6次元) から変換行列を求める
    """
    # prepare
    rot_vector, trans_vector = vec[:1], vec[1:]
    theta = np.abs(rot_vector)
    # rotation matrix
    rot_matrix = np.cos(theta) * np.eye(2)
    rot_matrix += get_skew_symmetric_matrix() * np.sin(theta)
    # translation
    J = calculate_J_matrix(theta)
    trans_for_T = np.dot(J, trans_vector)
    # transfromation matrix
    T = np.eye(3)
    T[:2, :2] = rot_matrix
    T[:2, 2] = trans_for_T
    return T


def get_skew_symmetric_matrix():
    result = np.zeros((2, 2))
    result[0, 1] = -1
    result[1, 0] = 1
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
    grad_theta = Sum_i (R.Ai + t)^T S.Bi
    grad_trans = Sum_i (R.Ai + t - Bi)
    S = [[0, -1], [1, 0]]
    """
    predicted_B = move_data(dataA, transform_a2b)  # (N,3)
    skew_matrix = get_skew_symmetric_matrix()
    s_dot_b = np.dot(dataB[:, :2], skew_matrix.T)  # (N,2)
    grad_theta = np.sum(predicted_B[:, :2] * s_dot_b)
    grad_trans = np.sum(predicted_B[:, :2] - dataB[:, :2], axis=0)  # (2,)
    grad = np.r_[grad_theta, grad_trans]
    return grad


def calculate_center(data: np.ndarray):
    center = np.mean(data[:, :2], axis=0)
    return center


def main():
    # make data
    radius = 1.0
    number = 50
    xy_original = make_data(number, r=radius)
    # correct transformation matrix
    rot_angle = np.deg2rad(90)
    trans_vector = np.array([0.1, 0.1])
    norm_trans = np.linalg.norm(trans_vector)
    trans_vector_lie = translation_vec_to_lie_translation(trans_vector, rot_angle)
    correct_t_vector = np.r_[rot_angle, trans_vector_lie]
    correct_T = vector_to_transformation_matrix(correct_t_vector)
    xy_moved = move_data(xy_original, correct_T)

    # predict transformation matrix
    init_angle = np.deg2rad(45)
    init_trans_vector = np.array([0.2, -0.2])
    init_trans_lie = translation_vec_to_lie_translation(init_trans_vector, init_angle)
    init_vector = np.r_[init_angle, init_trans_lie]
    init_matrix = vector_to_transformation_matrix(init_vector)
    current_T = init_matrix.copy()
    updated_t_vec = init_vector.copy()
    lr = 0.002
    scores = list()
    store_t_vector = list()
    store_center = list()

    # main
    iterations = 1000
    for _ in range(iterations):
        # calculate score, gradient
        _score = calculate_cost(xy_original, xy_moved, current_T)
        grad = calculate_cost_gradient(xy_original, xy_moved, current_T)
        # update transformation matrix
        updated_t_vec -= lr * grad
        current_T = vector_to_transformation_matrix(updated_t_vec)
        # store data
        scores.append(_score)
        store_t_vector.append(updated_t_vec.copy())
        # calculate center
        tmp_moved = move_data(xy_original, current_T)
        center = calculate_center(tmp_moved)
        store_center.append(center)

    # predict data
    predicted_xy = move_data(xy_original, current_T)

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
    ax2 = fig.add_subplot(222)
    ax2.scatter(xy_original[:, 0], xy_original[:, 1], color='purple', label='original', s=5)
    ax2.scatter(xy_moved[:, 0], xy_moved[:, 1], color='blue', label='moved', s=5)
    ax2.scatter(predicted_xy[:, 0], predicted_xy[:, 1], color='red', label='predicted', s=5)
    ax2.plot(store_center[:, 0], store_center[:, 1], color='orange')
    ax2.set_xlim(-radius-norm_trans, radius+norm_trans)
    ax2.set_ylim(-radius-norm_trans, radius+norm_trans)
    ax2.set_title('xy data')
    plt.legend()

    # rotation vectors, translation vectors (2d)
    ax3 = fig.add_subplot(223)
    if SHOW_TRAJECTORY:
        _p = ax3.scatter(
            np.arange(iterations), store_t_vector[:, 0],
            c=np.arange(iterations) / iterations, cmap='rainbow', s=5)
    else:
        _p = ax3.scatter(
            np.arange(iterations), store_t_vector[:, 0],
            c=scores, cmap='gist_rainbow', s=5)
    ax3.scatter(iterations, correct_t_vector[0], color='red')
    ax3.set_title('rotation vector')
    fig.colorbar(_p)

    ax4 = fig.add_subplot(224)
    if SHOW_TRAJECTORY:
        _p = ax4.scatter(
            store_t_vector[:, 1], store_t_vector[:, 2],
            c=np.arange(iterations) / iterations, cmap='rainbow', s=5)
    else:
        _p = ax4.scatter(
            store_t_vector[:, 1], store_t_vector[:, 2],
            c=scores, cmap='gist_rainbow', s=5)
    ax4.scatter(correct_t_vector[1], correct_t_vector[2], color='red')
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
