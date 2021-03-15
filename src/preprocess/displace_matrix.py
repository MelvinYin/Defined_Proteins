"""
Displace a matrix file to the left or right, filling remaining spots with the
composition frequencies.
"""
import matplotlib.pyplot as plt
import os

from config import paths
from utils import seq_logo
from preprocess import crop_matrix


def displace(displacement, matrix_file, output_file,
             composition_file=paths.COMPOSITION):
    matrix = []
    with open(matrix_file, 'r') as file:
        for line in file:
            matrix.append((line.strip().split(" ")))
    composition = []
    with open(composition_file, 'r') as file:
        for line in file:
            split_line = line.strip().split(" ")
            if len(split_line) != 2:
                continue
            composition.append(split_line[1])

    output_matrix = []
    assert len(composition) == len(matrix[0])
    if displacement > 0:
        # shift right
        for i in range(displacement):
            output_matrix.append(composition)
        for i in range(len(matrix) - displacement):
            output_matrix.append(matrix[i])
    else:
        # shift left
        displacement *= -1
        for i in range(displacement, len(matrix)):
            output_matrix.append(matrix[i])
        for i in range(displacement):
            output_matrix.append(composition)
    assert len(output_matrix) == len(matrix)
    with open(output_file, 'w') as file:
        for line in output_matrix:
            file.write(" ".join(line) + "\n")


def test_displace():
    direct_from_nbdb = os.path.join(paths.USER_INPUT, "GxxGxG_pssm.txt")
    cropped = os.path.join(paths.USER_INPUT, "GxxGxG_pssm_cropped.txt")
    test_output = os.path.join(paths.TEST, 'test.txt')
    crop_matrix.crop(direct_from_nbdb, cropped)

    displace(-2, cropped, test_output, paths.COMPOSITION)
    seq_logo.build_logo_nbdb(test_output)
    os.remove(cropped)
    os.remove(test_output)
    plt.show()

