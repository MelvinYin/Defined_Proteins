import os

from config import paths

def remove(fasta_path=None, output_path=None):
    if fasta_path is None:
        fasta_fname = "mg_50.fasta"
        fasta_path = os.path.join(paths.ROOT, 'data', 'user', 'input', fasta_fname)
    output = []
    with open(fasta_path, 'r') as file:
        title = next(file)
        current_seq = ""
        for line in file:
            if line.startswith(">"):
                if len(current_seq) > 50:
                    output.append(title)
                    output.append(current_seq)
                    title = line
                    current_seq = ""
                    continue
                else:
                    title = line
                    current_seq = ""
                    continue
            current_seq += line
        if len(current_seq) > 50:
            output.append(title)
            output.append(current_seq)
    if output_path is None:
        output_path = fasta_path
    with open(output_path, 'w') as file:
        file.writelines(output)

if __name__ == "__main__":
    remove(os.path.join(paths.ROOT, "mg_100_screened.fasta"),
           os.path.join(paths.ROOT, "mg_100_filtered.fasta"))

