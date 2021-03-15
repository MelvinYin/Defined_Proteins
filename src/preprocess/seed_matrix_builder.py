import numpy as np
import os

from config import paths
from preprocess.search import search_converters
from utils import generic



def from_aligned(selected_seqs, aligned_seqs, output):
    # Length of actual motif
    FINAL_MOTIF_LEN = 30

    search_input_fname = os.path.join(paths.DEBUG, "search_input.txt")

    matrix_len = search_converters.seq_to_search_input(aligned_seqs,
                                                       search_input_fname)
    aligned_sequences = search_converters.search_get_seq_direct(
        search_input_fname, selected_seqs, matrix_len)
    output_matrix = np.zeros((FINAL_MOTIF_LEN, 20), dtype=int)
    for seq in aligned_sequences:
        for i, char in enumerate(seq):
            if char in generic.AA1_to_index:
                output_matrix[i][generic.AA1_to_index[char]] += 1

    with open(output, 'w') as file:
        for line in output_matrix:
            for i, value in enumerate(line):
                to_write = str(value)
                if i != len(line) - 1:
                    to_write += ","
                file.write(to_write)
            file.write("\n")
    return len(aligned_sequences)

def from_nbdb(pssm_file, num_seqs, output):
    """
    For nbdb pssm files
    GxGGxG_pssm: 4000
    """
    output_matrix = convert_nbdb_matrix_to_conv_encoder(pssm_file, num_seqs)
    with open(output, 'w') as file:
        for line in output_matrix:
            for i, value in enumerate(line):
                to_write = str(value)
                if i != len(line) - 1:
                    to_write += ","
                file.write(to_write)
            file.write("\n")


def convert_nbdb_matrix_to_conv_encoder(nbdb_file, num_seqs):
    matrix = [[] for __ in range(30)]
    with open(nbdb_file, 'r') as file:
        for i, line in enumerate(file):
            probs = line.split(" ")
            assert len(probs) == 20
            for prob in probs:
                matrix[i].append(round(float(prob) * num_seqs))
    for i in matrix:
        if sum(i) > num_seqs:
            if i[0] > 0:
                i[0] -= sum(i) - num_seqs
            else:
                i[1] -= sum(i) - num_seqs
        elif sum(i) < num_seqs:
            i[0] += num_seqs - sum(i)
    for i in matrix:
        assert sum(i) == num_seqs
    return matrix
