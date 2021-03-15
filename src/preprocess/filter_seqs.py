def delete_short_seqs(fasta_fname, threshold):
    """
    Delete sequences with len < threshold.
    """
    output_str = []
    current_seq = []
    header = None
    with open(fasta_fname, 'r') as file:
        for line in file:
            if line.startswith(">"):
                if current_seq:
                    if header is None:
                        raise Exception
                    current_seq = "".join(current_seq)
                    if len(current_seq) >= threshold:
                        output_str.append(header)
                        output_str.append(current_seq)
                    current_seq = []
                header = line
                continue
            current_seq.append(line)
    output_str = ''.join(output_str)
    with open(fasta_fname, 'w') as file:
        file.write(output_str)