import os

from config import paths
from utils import generic

def screen(fasta_path=None, output_path=None):
    if not fasta_path:
        fasta_fname = "mg_50.fasta"
        fasta_path = os.path.join(paths.ROOT, 'data', 'user', 'input', fasta_fname)
    output = []
    alphabets = set(generic.AA3_to_AA1.values())
    to_keep = True
    with open(fasta_path, 'r') as file:
        title = next(file)
        current_seq = ""
        for line in file:
            if line.startswith(">"):
                if to_keep:
                    output.append(title)
                    output.append(current_seq)
                    title = line
                    current_seq = ""
                    continue
                else:
                    title = line
                    current_seq = ""
                    to_keep = True
                    continue
            if not to_keep:
                continue
            for char in line.strip():
                if char not in alphabets:
                    to_keep = False
                    break
            current_seq += line
        if to_keep:
            output.append(title)
            output.append(current_seq)
    if output_path is None:
        output_path = fasta_path
    with open(output_path, 'w') as file:
        file.writelines(output)

if __name__ == "__main__":

    screen(os.path.join(paths.ROOT, "mg_100.fasta"), os.path.join(paths.ROOT,
                                                                  "mg_100_screened.fasta"))

