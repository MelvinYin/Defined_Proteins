from collections import defaultdict, OrderedDict
import logging
import numpy as np
import os
import re
import subprocess

from config import paths
from utils import generic


def search_run(input_file, search_seqs):
    search_input = os.path.join(paths.SEARCH_DIR, "extracted.matrix")
    search_output = os.path.join(paths.SEARCH_DIR, "output.txt")
    if os.path.isfile(search_input):
        os.remove(search_input)
    if os.path.isfile(search_output):
        os.remove(search_output)
    # output_meme_to_search_input(output_meme, search_input, length)
    length = conv_output_to_search_input(input_file, search_input)
    assert os.path.isfile(search_input)
    assert os.path.isfile(search_seqs)
    command = f"{paths.SEARCH_EXEC} {search_input} {search_seqs} " \
              f"{search_output} {length}"
    subprocess.run(command, shell=True)
    assert os.path.isfile(search_output)
    motif_pos = search_output_to_motif_pos_scan(search_output)
    return motif_pos

def conv_output_to_search_input(input_file, search_input):
    assert os.path.isfile(input_file)
    nsites = None
    matrix = []
    with open(input_file, 'r') as file:
        for line in file:
            if nsites is None:
                nsites = int(line.strip())
                continue
            if line.strip().startswith("A,C,D"):
                continue
            if nsites is not None and line.strip():
                line_arr = []
                for term in line.strip().split(","):
                    # line_arr.append(term)
                    line_arr.append(str(round(int(term)/nsites, 6)))
                matrix.append(line_arr)
                assert len(matrix[-1]) == 20, matrix

    with open(search_input, 'w') as file:
        file.write("BEGIN\n")
        file.write(f"MATRIX K={nsites} L={len(matrix)}\n")
        for line in matrix:
            for i, term in enumerate(line):
                if i != len(line) - 1:
                    file.write(term + " ")
                else:
                    file.write(term + "\n")
        file.write("END")
    return len(matrix)

def search_get_seq_direct(search_input, seq_file, length=30):
    search_output = os.path.join(paths.SEARCH_DIR, "output.txt")
    if os.path.isfile(search_output):
        os.remove(search_output)
    assert os.path.isfile(paths.SEARCH_EXEC), paths.SEARCH_EXEC
    assert os.path.isfile(search_input)
    assert os.path.isfile(seq_file)
    command = f"{paths.SEARCH_EXEC} {search_input} {seq_file} " \
              f"{search_output} {length}"
    subprocess.run(command, shell=True)
    assert os.path.isfile(search_output)
    seqs = search_output_to_motif_seq(search_output)
    os.remove(search_output)
    return seqs

def get_map():
    pdb_cid_map = dict()
    cluster_map = defaultdict(list)
    cluster_i = 0
    with open(paths.PDB_CLUSTER, 'r') as file:
        for line in file:
            for term in line.strip().split(" "):
                pdb_id, cid = term.strip().split("_")
                pdb_id = pdb_id.lower()
                pdb_cid_map[(pdb_id, cid)] = cluster_i
                cluster_map[cluster_i].append((pdb_id, cid))
            cluster_i += 1
    return pdb_cid_map, cluster_map

def search_output_to_motif_pos(search_output):
    motif_positions = dict()
    cid = None
    pdb_id = None

    with open(search_output, 'r') as file:
        for line in file:
            if not line.strip():
                continue
            if line.startswith(">"):
                pdb_id, cid = line[1:].strip().split("_")[:2]
                cid = cid[0]
                continue

            if not re.match("[0-9]+", line.strip()):
                continue
            assert pdb_id is not None
            assert cid is not None
            motif_pos = int(line.strip())
            if pdb_id not in motif_positions:
                motif_positions[pdb_id] = dict()
                motif_positions[pdb_id]['sno_markers'] = [motif_pos]
                motif_positions[pdb_id]['cid'] = [cid]
            else:
                motif_positions[pdb_id]['sno_markers'].append(motif_pos)
                motif_positions[pdb_id]['cid'].append(cid)
            pdb_id = None
            cid = None
    motif_positions = OrderedDict(sorted(motif_positions.items()))
    return motif_positions


def search_output_to_motif_pos_scan(search_output):
    motif_positions = dict()
    cid = None
    pdb_id = None
    seen = set()
    pdb_cid_map, cluster_map = get_map()
    with open(search_output, 'r') as file:
        for line in file:
            if not line.strip():
                continue
            if line.startswith(">"):
                pdb_id, cid = line[1:].strip().split("_")[:2]
                cid = cid[0]
                continue

            if not re.match("[0-9]+", line.strip()):
                continue
            assert pdb_id is not None
            assert cid is not None
            motif_pos = int(line.strip())
            if (pdb_id, cid) in seen:
                pdb_id = None
                cid = None
                continue
            if (pdb_id, cid) in pdb_cid_map:
                cluster_i = pdb_cid_map[(pdb_id, cid)]
                for term in cluster_map[cluster_i]:
                    seen.add(term)
            else:
                print((pdb_id, cid))
            if pdb_id not in motif_positions:
                motif_positions[pdb_id] = dict()
                motif_positions[pdb_id]['sno_markers'] = [motif_pos]
                motif_positions[pdb_id]['cid'] = [cid]
            else:
                motif_positions[pdb_id]['sno_markers'].append(motif_pos)
                motif_positions[pdb_id]['cid'].append(cid)
            pdb_id = None
            cid = None
    motif_positions = OrderedDict(sorted(motif_positions.items()))
    return motif_positions

def search_output_to_motif_seq(search_output, length=30):
    output = []
    cid = None
    pdb_id = None
    with open(search_output, 'r') as file:
        for line in file:
            if line.startswith(">"):
                pdb_id, cid = line[1:].strip().split("_")[:2]
                cid = cid[0]
                continue
            assert pdb_id is not None
            assert cid is not None
            if re.match("[0-9]+", line.strip()):
                continue
            seq = line.strip()
            if len(seq) != length:
                continue
            output.append(seq)
    return output



def _matrix_builder(aligned_path):
    alphabets = set(generic.AA3_to_AA1.values())
    AA_to_index = {AA: i for i, AA in enumerate(sorted(alphabets))}
    matrix_counter = None
    print(aligned_path)
    with open(aligned_path, 'r') as file:
        for line in file:
            if line.startswith(">"):
                continue
            line = line.strip().upper()
            if matrix_counter is None:
                matrix_counter = np.zeros((len(line), len(alphabets)),
                                          dtype=int)
            for i, char in enumerate(line):
                try:
                    AA_index = AA_to_index[char]
                except KeyError:
                    # Key not found for some reason
                    continue
                matrix_counter[i, AA_index] += 1
    return matrix_counter

def _matrix_to_search_input(matrix_counter, filename):
    length = len(matrix_counter)
    nsites = sum(matrix_counter[0])
    with open(filename, 'w') as file:
        file.write("BEGIN\n")
        file.write(f"MATRIX K={nsites} L={length}\n")
        for row in matrix_counter:
            output = []
            for i, term in enumerate(row):
                output.append(str(round(term / nsites, 7)))
                if i != len(row) - 1:
                    output.append(" ")
            file.write("".join(output) + "\n")
        file.write("END")


def seq_to_search_input(sequences, output):
    matrix_counter = _matrix_builder(sequences)
    _matrix_to_search_input(matrix_counter, output)
    return len(matrix_counter)