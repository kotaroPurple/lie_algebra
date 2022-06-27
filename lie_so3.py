
import numpy as np
import matplotlib.pyplot as plt


def make_data(number: int, r: float = 1.0):
    """
    極座標で3次元データを作る
    x = r.sin_theta.cos_phi
    y = r.sin_theta.sin_phi
    z = r.cos_theta
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
    xyz = np.c_[xs, ys, zs]
    return xyz


def rotate_data(data: np.ndarray, rot_matrix: np.ndarray):
    """
    data: (N,3)
    rot_matrix: (3,3)
    """
    result = np.dot(data, rot_matrix.T)
    return result


def rotation_matrix_to_vector(rot_matrix: np.ndarray):
    """
    回転行列を回転ベクトルに変換する
    """
    rot_vector = -np.array([rot_matrix[1, 2] - rot_matrix[2, 1],
                            rot_matrix[2, 0] - rot_matrix[0, 2],
                            rot_matrix[0, 1] - rot_matrix[1, 0]])
    rot_vector = normalize_vector(rot_vector)
    trace_r = np.trace(rot_matrix)
    theta = np.arccos((trace_r - 1) / 2.)
    result = theta * rot_vector
    return result


def vector_to_rotation_matrix(rot_vector: np.ndarray):
    """
    回転ベクトルから回転行列を求める
    """
    normal_vector = normalize_vector(rot_vector)
    theta = np.sqrt(np.sum(rot_vector**2))
    rot_matrix = np.cos(theta) * np.eye(3)
    rot_matrix += (1 - np.cos(theta)) \
        * (np.dot(normal_vector[:, np.newaxis], normal_vector[np.newaxis, :]))
    rot_matrix += np.sin(theta) * rot_vector_to_skew_symmetric_matrix(normal_vector)
    return rot_matrix


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


def calculate_cost(dataA: np.ndarray, dataB: np.ndarray, rot_a2b: np.ndarray):
    """
    二乗誤差を求める
    cost = 1/N * Sum_i |R.Ai - Bi|^2
    """
    predicted_B = rotate_data(dataA, rot_a2b)
    cost = np.mean((predicted_B - dataB)**2)
    return cost


def calculate_cost_gradient(dataA: np.ndarray, dataB: np.ndarray, rot_a2b: np.ndarray):
    """
    二乗誤差の微分を求める
    cost = 1/N * Sum_i |R.Ai - Bi|^2
    grad = -Sum_i (R.Ai) x (Bi - R.Ai)
         = -Sum_i (R.Ai) * (Bi)
    """
    predicted_B = rotate_data(dataA, rot_a2b)
    cross_result = np.cross(predicted_B, dataB)
    grad = -np.sum(cross_result, axis=0)  # (3,)
    return grad


def calculate_center(data: np.ndarray):
    center = np.mean(data[:, :3], axis=0)
    return center


def main():
    # correct rotation matrix
    rot_vector = normalize_vector(np.array([1, 0, 0]))
    rot_angle = np.deg2rad(150)
    correct_R = vector_to_rotation_matrix(rot_angle * rot_vector)
    # make data
    radius = 1.0
    xys_original = make_data(100, r=radius)
    # rotate data
    xys_moved = rotate_data(xys_original, correct_R)  # original (dataA) -> moved (dataB)

    # predict rotation matrix
    init_vector = normalize_vector(np.array([1, 1, 0]))
    init_angle = np.deg2rad(30)
    init_matrix = vector_to_rotation_matrix(init_angle * init_vector)
    current_rot = init_matrix.copy()
    updated_rot_vec = init_angle * init_vector
    lr = 0.001
    scores = list()
    store_rot_vector = list()
    store_center = list()

    for _ in range(1000):
        # calculate score, gradient
        _score = calculate_cost(xys_original, xys_moved, current_rot)
        grad = calculate_cost_gradient(xys_original, xys_moved, current_rot)
        # update rotation matrix
        # _rot_vec = rotation_matrix_to_vector(current_rot)
        updated_rot_vec -= lr * grad
        current_rot = vector_to_rotation_matrix(updated_rot_vec)
        # store data
        scores.append(_score)
        store_rot_vector.append(updated_rot_vec.copy())
        # calculate center
        tmp_moved = rotate_data(xys_original, current_rot)
        center = calculate_center(tmp_moved)
        store_center.append(center)

    # predict data
    predicted_xys = rotate_data(xys_original, current_rot)

    # to np.ndarray
    store_rot_vector = np.array(store_rot_vector)
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
        xys_original[:, 0], xys_original[:, 1], xys_original[:, 2],
        color='purple', label='original')
    ax2.scatter(
        xys_moved[:, 0], xys_moved[:, 1], xys_moved[:, 2],
        color='blue', label='moved')
    ax2.scatter(
        predicted_xys[:, 0], predicted_xys[:, 1], predicted_xys[:, 2],
        color='red', label='predicted')
    ax2.plot(store_center[:, 0], store_center[:, 1], store_center[:, 2], color='orange')
    ax2.set_xlim3d(-radius, radius)
    ax2.set_ylim3d(-radius, radius)
    ax2.set_zlim3d(-radius, radius)
    ax2.set_title('xyz data')
    plt.legend()

    # rotation vectors (3d)
    correct_rot_vector = rotation_matrix_to_vector(correct_R)
    ax3 = fig.add_subplot(223, projection='3d')
    _p = ax3.scatter(
        store_rot_vector[:, 0], store_rot_vector[:, 1], store_rot_vector[:, 2], c=scores)
    ax3.scatter(correct_rot_vector[0], correct_rot_vector[1], correct_rot_vector[2], color='red')
    ax3.set_title('rotation vector with colorized score')
    fig.colorbar(_p)
    plt.show()


if __name__ == '__main__':
    main()
