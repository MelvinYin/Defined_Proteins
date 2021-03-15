import os

from config import paths
from preprocess.converge import wrapper, conv_to_meme


def encode_proteome(proteome_fname, output):
    conv_folder = os.path.join(paths.CONV_FOLDER, "binaries")
    bash_exec = paths.BASH_EXEC
    return wrapper.encode_proteome(proteome_fname, output, conv_folder,
                                   bash_exec)


def encode_blosum(output):
    conv_folder = os.path.join(paths.CONV_FOLDER, "binaries")
    bash_exec = paths.BASH_EXEC
    return wrapper.encode_blosum(output, conv_folder,
                                   bash_exec)


def encode_matrix(input_matrix, output):
    conv_folder = os.path.join(paths.CONV_FOLDER, "binaries")
    bash_exec = paths.BASH_EXEC
    return wrapper.encode_matrix(input_matrix, output, conv_folder,
                                 bash_exec)

from preprocess.converge import meme_to_conv

def convert_meme_to_conv(meme, composition, matrix, minimal=False):
    if minimal:
        ret_code = meme_to_conv.convert_minimal(meme, composition, matrix)
    else:
        ret_code = meme_to_conv.convert_full(meme, composition, matrix)
    return ret_code

def encode_input_matrix(matrix, filename):
    with open(filename, 'w') as file:
        for line in matrix:
            for term in line[:-1]:
                file.write(str(term)+",")
            file.write(str(line[-1]))
        file.write("\n")

import shutil
def run(input_matrix, profile_length, kmatches, output_matrix,
        debug_meme=None, iteration=5):
    shutil.copyfile(input_matrix, paths.CONV_INPUT_MATRIX)

    wrapper.converge_calculate(profile_length, kmatches, output_matrix,
                               iteration=iteration)
    os.remove(paths.CONV_INPUT_MATRIX)
    if debug_meme is not None:
        conv_to_meme.convert(output_matrix, paths.COMPOSITION, debug_meme)