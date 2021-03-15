"""
Converts a meme output with a single motif to a converge matrix and
composition file.

Single, because that's what we need for now. We need two separate parsers,
one for the full meme, and one for minimal.
Biopython's meme parser (as of v1.73) is rather buggy, so skipping it for now.
"""

from collections import OrderedDict
import os
import re

from config import paths

ALPHABETS = "ACDEFGHIKLMNPQRSTVWY"

def convert_full(input_meme_path, output_composition, output_matrix):
    # Converge output to minimal
    # Converts converge motif format to minimal meme format
    # see http://meme-suite.org/doc/examples/sample-protein-motif.meme

    # input_conv=''output.4.matrix.0''
    # composition='composition.txt'
    # output="meme_format.txt"

    composition_map, matrices, n_sites = _parse_full_meme_input(input_meme_path)
    _write_composition(composition_map, output_composition)
    _write_matrix(matrices, n_sites, output_matrix)
    return


def convert_minimal(input_meme_path, output_composition, output_matrix):
    composition_map, n_sites, matrices = _parse_minimal_meme_input(
        input_meme_path)
    _write_composition(composition_map, output_composition)
    _write_matrix(matrices, n_sites, output_matrix)
    return

def _parse_minimal_meme_input(input_path):
    with open(input_path, 'r') as file:
        lines = file.readlines()

    # Extract matrix
    all_counts = []
    start_matrix = False
    for line in lines:
        if not start_matrix:
            if not line.startswith("letter-probability matrix: "):
                continue
            else:
                start_matrix = True
                continue
        if not line.strip():
            break
        counts_in_str = line.strip().split(" ")
        assert len(counts_in_str) == 20, counts_in_str
        all_counts.append(list(float(i) for i in counts_in_str))
    assert len(all_counts) == 30

    # Extract composition
    alphabet_probs = dict()
    start_composition = False
    for line in lines:
        if not start_composition:
            if not line.startswith("Background letter frequencies"):
                continue
            else:
                start_composition = True
                continue
        if not line.strip():
            break
        split_line = line.strip().split(" ")
        assert len(split_line) % 2 == 0
        for i in range(len(split_line) // 2):
            alphabet, prob = split_line[i * 2], split_line[i * 2 + 1]
            prob = float(prob)
            alphabet_probs[alphabet] = prob
    print(alphabet_probs)
    assert len(alphabet_probs) == len(ALPHABETS)

    # Extract n_sites
    n_sites = None
    for line in lines:
        if not line.startswith("letter-probability matrix"):
            continue
        nsite_match = re.search(r"nsites= ([0-9]+)", line)
        assert nsite_match is not None
        n_sites = int(nsite_match[1])
        break
    assert n_sites is not None

    return alphabet_probs, n_sites, all_counts

def _write_matrix(matrices, n_sites, output_path):
    formatted_lines = []
    formatted_lines.append("PROFILE 4\n")
    formatted_lines.append("BEGIN\n")
    formatted_lines.append("ORIGIN\n")
    formatted_lines.append(f"MATRIX ID=0 K={n_sites} L=50\n")
    formatted_lines.append("50        " + "        ".join(list(ALPHABETS)) +
                           "\n")
    for i, probs in enumerate(matrices):
        line = ""
        if i <= 9:
            line += f" {i}"
        else:
            line += str(i)
        for prob in probs:
            prob_in_str = "{:.6f}".format(prob)
            line += f" {prob_in_str}"
        line += "\n"
        formatted_lines.append(line)
    formatted_lines.append("END")
    with open(output_path, 'w') as file:
        file.writelines(formatted_lines)

def _parse_full_meme_input(path):
    with open(path, 'r') as file:
        lines = file.readlines()
    # Extract matrix
    all_counts = []
    start_matrix = False
    for line in lines:
        if not start_matrix:
            if not line.startswith("letter-probability matrix: "):
                continue
            else:
                start_matrix = True
                continue
        if line.startswith("------------------------"):
            break
        counts_in_str = line.split("  ")
        assert len(counts_in_str) == 20
        all_counts.append(list(float(i) for i in counts_in_str))
    assert len(all_counts) == 50
    # Extract composition
    alphabet_probs = dict()
    start_composition = False
    for line in lines:
        if not start_composition:
            if not line.startswith("Letter frequencies in dataset"):
                continue
            else:
                start_composition = True
                continue
        if line.startswith("Background letter frequencies"):
            break
        split_line = line.strip().split(" ")
        assert len(split_line) % 2 == 0
        for i in range(len(split_line) // 2):
            alphabet, prob = split_line[i*2], split_line[i*2+1]
            prob = float(prob)
            alphabet_probs[alphabet] = prob
    assert len(alphabet_probs) == len(ALPHABETS)
    # Extract n_sites
    n_sites = None
    for line in lines:
        if not re.search("sites =", line):
            continue
        n_sites_match = re.search("sites = ([0-9]+)", line)
        assert n_sites_match is not None
        n_sites = int(n_sites_match[1])
        break
    assert n_sites is not None
    return alphabet_probs, all_counts, n_sites

def _write_composition(composition_counts, output_path):
    composition_counts = OrderedDict(sorted(composition_counts.items(),
                                            reverse=True))
    formatted_lines = []
    for alphabet, prob in composition_counts.items():
        prob_in_str = "{:.6f}".format(prob*100)
        line = f"{alphabet},{prob_in_str}\n"
        formatted_lines.append(line)
    with open(output_path, 'w') as file:
        file.writelines(formatted_lines)

if __name__ == "__main__":
    input_meme_path = os.path.join(paths.ROOT, 'meme.txt')
    output_composition = os.path.join(paths.ROOT, 'composition.txt')
    output_matrix = os.path.join(paths.ROOT, 'output_mine.txt')
    convert_full(input_meme_path, output_composition, output_matrix)

# if __name__ == "__main__":
#     input_meme_path = os.path.join(paths.ROOT, 'conv_meme.txt')
#     output_composition = os.path.join(paths.ROOT, 'comp.txt')
#     output_matrix = os.path.join(paths.ROOT, 'out.txt')
#     # convert_full(input_meme_path, output_composition, output_matrix)
#     convert_minimal(input_meme_path, output_composition, output_matrix)




# path = os.path.join(paths.ROOT, 'meme.txt')
# # a = meme.read(path)
# with open(path, 'r') as file:
#     # print(a)
#     a = meme.read(file)
#     # print(a)
# for motif in a:
#     # print(motif.motif_name, motif.sequence_name, motif.strand,
#     #       motif.pvalue)
#     print(motif.counts)


#
# def convert(input_meme_path, output_composition, output_matrix):
#     # Converge output to minimal
#     # Converts converge motif format to minimal meme format
#     # see http://meme-suite.org/doc/examples/sample-protein-motif.meme
#
#     # input_conv=''output.4.matrix.0''
#     # composition='composition.txt'
#     # output="meme_format.txt"
#
#     composition_map, matrices = _parse_meme_input(input_meme_path)
#     _write_composition(composition_map, output_composition)
#     _write_matrix(matrices, output_matrix)
#     return
#
#
#
#
#
# def _parse_meme_input(filename):
#     alphabets = ""
#     length = 30
#     matrices = OrderedDict()
#     matrix = []
#     nsite = 0
#     matrix_count = 0
#     with open(filename, "r") as file:
#         for line in file:
#             if line.startswith("BEGIN") and matrix_count != 0:
#                 assert len(matrix) == length, len(matrix)
#                 motif_name = f"MEME-{matrix_count}"
#                 matrices[motif_name] = (nsite, matrix)
#                 assert nsite != 0
#                 matrix = []
#                 nsite = 0
#                 continue
#             if line.startswith("MATRIX"):
#                 matrix_count += 1
#                 match = re.search(r"K=([0-9]+)", line)
#                 if match is None:
#                     raise AssertionError
#                 nsite = int(match[1])
#                 continue
#             if (line.startswith("50") or line.startswith("30")):
#                 if not alphabets:
#                     matched_alphabets = re.findall("[A-Z]", line)
#                     alphabets = "".join(matched_alphabets)
#                 continue
#             if re.match(" [0-9]", line) or re.match("[0-9]+", line):
#                 probs = re.findall(r"[0-1]\.[0-9]+", line)
#                 assert len(probs) == len(alphabets)
#                 matrix.append(probs)
#                 continue
#         if matrix:
#             motif_name = f"MEME-{matrix_count}"
#             matrices[motif_name] = (nsite, matrix)
#     return alphabets, matrices
#
#
# def _parse_converge_composition(filename):
#     composition_map = dict()
#     with open(filename, "r") as file:
#         for line in file:
#             if re.match("[A-Z]", line):
#                 alphabet = line[0]
#                 composition = line[2:]
#                 composition_map[alphabet] = float(composition)
#                 continue
#     summed_composition = sum(composition_map.values())
#     for key, value in composition_map.items():
#         composition_map[key] = value / summed_composition
#     return composition_map
#
#
# def _format_minimal_from_conv(alphabets, composition_map, matrices, output):
#     m_to_write = list(range(len(matrices)))
#     with open(output, 'w') as file:
#         file.write("MEME version 4\n")
#         file.write("\n")
#         file.write("ALPHABET= " + alphabets + "\n")
#         file.write("\n")
#         file.write("Background letter frequencies\n")
#         for i, alphabet in enumerate(alphabets):
#             composition = round(composition_map[alphabet], 4)
#             file.write(f"{alphabet} {composition} ")
#             if (i != 0) and (i % 9 == 0):
#                 file.write("\n")
#         file.write("\n")
#         file.write("\n")
#         m_count = 0
#         while matrices:
#             motif_name, (nsite, matrix) = matrices.popitem(last=False)
#             if m_count not in m_to_write:
#                 m_count += 1
#                 continue
#             m_count += 1
#             file.write(f"MOTIF {motif_name}")
#             file.write("\n")
#             file.write(f"letter-probability matrix: alength= 20 w= 30 "
#                        f"nsites= {nsite} E= 0.000")  # alength = len(alphabets)
#             # E is just some random number, filled in by subsequent eval calc.
#             # w = width of motif
#             file.write("\n")
#             for line in matrix:
#                 to_write = ""
#                 for prob in line:
#                     to_write += prob + " "
#                 file.write(to_write)
#                 file.write("\n")
#             file.write("\n")
#     return
#
#
# if __name__ == "__main__":
#     import os
#     from config import paths
#
#     input_conv_path = os.path.join(paths.ROOT, "output.4.matrix.0")
#     composition_path = os.path.join(paths.ROOT, "composition.txt")
#     output_path = os.path.join(paths.ROOT, "output_meme.txt")
#     convert(input_conv_path, composition_path, output_path)