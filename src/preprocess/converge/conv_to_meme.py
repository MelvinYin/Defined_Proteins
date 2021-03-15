import re

import numpy as np

from utils import generic


def convert(input_conv_path, composition_path, output_path):
    # Converge output to minimal
    # Converts converge motif format to minimal meme format
    # see http://meme-suite.org/doc/examples/sample-protein-motif.meme

    # input_conv=''output.4.matrix.0''
    # composition='composition.txt'
    # output="meme_format.txt"
    num_matches, pssm = _parse_converge_output(input_conv_path)
    composition_map = _parse_converge_composition(composition_path)
    _format_minimal_from_conv(num_matches, pssm, composition_map, output_path)
    return


def _parse_converge_output(filename):
    pssm = []
    num_matches = None
    with open(filename, "r") as file:
        for line in file:
            if num_matches is None:
                num_matches = int(line)
                continue
            if not line.strip():
                break
            if re.match("[A-Za-z]", line.strip()):
                continue
            values = re.findall("[0-9]+", line)
            matrix_line = np.array(list(map(int, values)), dtype=int)
            pssm.append(matrix_line)
    pssm = np.array(pssm, dtype=int)
    return num_matches, pssm


def _parse_converge_composition(filename):
    composition_map = dict()
    with open(filename, "r") as file:
        for line in file:
            if re.match("[A-Z]", line):
                alphabet = line[0]
                composition = line[2:]
                composition_map[alphabet] = float(composition)
                continue
    summed_composition = sum(composition_map.values())
    for key, value in composition_map.items():
        composition_map[key] = value / summed_composition
    return composition_map


def _format_minimal_from_conv(num_matches, pssm, composition_map, output):
    alphabet_str = "".join(generic.AA_values)
    with open(output, 'w') as file:
        file.write("MEME version 4\n")
        file.write("\n")
        file.write("ALPHABET= " + alphabet_str + "\n")
        file.write("\n")
        file.write("Background letter frequencies\n")
        for i, alphabet in enumerate(generic.AA_values):
            composition = round(composition_map[alphabet], 4)
            file.write(f"{alphabet} {composition} ")
            if (i != 0) and (i % 9 == 0):
                file.write("\n")
        file.write("\n")
        file.write("\n")
        m_count = 0

        m_count += 1
        file.write(f"MOTIF MEME-1\n")
        file.write(f"letter-probability matrix: alength= 20 w= {len(pssm)} "
                   f"nsites= {num_matches} E= 0.000\n")
        # E is just some random number, filled in by subsequent eval calc.
        for line in pssm:
            total_counts = sum(line)
            for count in line:
                file.write("{:.6f} ".format(count/total_counts))
            file.write("\n")
        file.write("\n")
    return


# if __name__ == "__main__":
#     import os
#     from config import paths
#
#     input_conv_path = os.path.join(paths.ROOT, "output_mine")
#     composition_path = os.path.join(paths.ROOT, "composition.txt")
#     output_path = os.path.join(paths.ROOT, "output_meme.txt")
#     convert(input_conv_path, composition_path, output_path)
