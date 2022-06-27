
import numpy as np
import matplotlib.pyplot as plt
from plyfile import PlyData  # https://github.com/dranjan/python-plyfile
from lie_se3 import (
    move_data,
    # transformation_matrix_to_vector,
    vector_to_transformation_matrix,
    calculate_cost,
    calculate_cost_gradient,
    normalize_vector,
    translation_vec_to_lie_translation
)

# 同階層の data にある ply ファイルを読み込む
import os
_dirpath = os.path.dirname(__file__)
PLY_FILEPATH = os.path.join(_dirpath, 'data/bun000.ply')

SHOW_TRAJECTORY = True


def read_ply_file(filepath: str):
    """
    ply ファイルを読み込み np.ndarray (N,4) にする
    """
    pcd = PlyData.read(filepath)
    vertex = pcd['vertex']
    (x, y, z) = (vertex[t] for t in ('x', 'y', 'z'))
    n = x.size
    xyz = np.c_[x, y, z, np.ones(n)]  # (N,4)
    return xyz


def choice_data(data: np.ndarray, ratio: float = 0.5):
    """
    データを間引く
    """
    number = data.shape[0]
    use_number = int(number * ratio)
    perm = np.random.permutation(number)
    data = data[perm[:use_number], :]
    print(f'choice data: {number} -> {use_number} (ratio= {ratio})')
    return data


def main():
    # read data (ply)
    print('=> loading PLY file..')
    xyz_original = read_ply_file(PLY_FILEPATH)
    xyz_original = choice_data(xyz_original, ratio=0.1)
    print(f'=> loaded {PLY_FILEPATH}')

    # enlarge data
    xyz_original[:, :3] *= 2

    # move to center
    center = np.mean(xyz_original[:, :3], axis=0)
    xyz_original[:, :3] -= center

    # data range
    xmin, xmax = xyz_original[:, 0].min(), xyz_original[:, 0].max()
    ymin, ymax = xyz_original[:, 1].min(), xyz_original[:, 1].max()
    zmin, zmax = xyz_original[:, 2].min(), xyz_original[:, 2].max()
    print('=> data range')
    print(f'xmin, xmax: {xmin}, {xmax}')
    print(f'ymin, ymax: {ymin}, {ymax}')
    print(f'zmin, zmax: {zmin}, {zmax}')

    # correct transformation matrix
    rot_vector = normalize_vector(np.array([1, 1, 0]))
    rot_angle = np.deg2rad(50)
    trans_vector = np.array([0.1, 0.1, 0.1])
    # trans_vector = np.array([0.0, 0.0, 0.0])
    norm_trans = np.linalg.norm(trans_vector)
    trans_vector_lie = translation_vec_to_lie_translation(trans_vector, rot_vector, rot_angle)
    correct_t_vector = np.r_[rot_angle * rot_vector, trans_vector_lie]
    correct_T = vector_to_transformation_matrix(correct_t_vector)
    xyz_moved = move_data(xyz_original, correct_T)

    # xyz_moved info
    # moved_center = np.mean(xyz_moved[:, :3], axis=0)  # xyz

    # predict transformation matrix
    init_rot_vector = normalize_vector(np.array([1, 0, 0]))
    init_angle = np.deg2rad(5)
    init_trans_vector = np.array([0.00, 0.00, 0.00])
    # init_trans_vector = moved_center
    init_trans_lie = translation_vec_to_lie_translation(
        init_trans_vector, init_rot_vector, init_angle)
    init_vector = np.r_[init_angle * init_rot_vector, init_trans_lie]
    init_matrix = vector_to_transformation_matrix(init_vector)
    current_T = init_matrix.copy()
    updated_t_vec = init_vector.copy()
    lr = 0.0001
    scores = list()
    store_t_vector = list()

    iterations = 1000

    print(f'=> predict : iteration= {iterations}')
    for i in range(iterations):
        # calculate score, gradient
        _score = calculate_cost(xyz_original, xyz_moved, current_T)
        grad = calculate_cost_gradient(xyz_original, xyz_moved, current_T)
        # update transformation matrix
        updated_t_vec -= lr * grad
        current_T = vector_to_transformation_matrix(updated_t_vec)
        # store data
        scores.append(_score)
        store_t_vector.append(updated_t_vec.copy())
        # update learning rate
        if i == iterations // 2:
            lr /= 2.

    # predict data
    predicted_xyz = move_data(xyz_original, current_T)

    #
    # show result
    #

    # show score vs iteration
    fig = plt.figure(figsize=(12, 8))
    ax1 = fig.add_subplot(221)
    ax1.plot(scores)
    ax1.set_title('score')

    # rotation vectors, translation vectors (3d)
    store_t_vector = np.array(store_t_vector)
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

    # show 3d data (original, correct data, predicted data)
    fig2 = plt.figure(figsize=(8, 8))
    ax2 = fig2.add_subplot(111, projection='3d')
    ax2.scatter(
        xyz_original[:, 0], xyz_original[:, 1], xyz_original[:, 2],
        color='purple', label='original', s=5)
    ax2.scatter(
        xyz_moved[:, 0], xyz_moved[:, 1], xyz_moved[:, 2],
        color='blue', label='moved', s=5)
    ax2.scatter(
        predicted_xyz[:, 0], predicted_xyz[:, 1], predicted_xyz[:, 2],
        color='red', label='predicted', s=5)
    ax2.set_xlim3d(xmin-norm_trans, xmax+norm_trans)
    ax2.set_ylim3d(ymin-norm_trans, ymax+norm_trans)
    ax2.set_zlim3d(zmin-norm_trans, zmax+norm_trans)
    ax2.set_title('xyz data')
    plt.legend()

    print()
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
