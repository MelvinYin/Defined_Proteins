import contextlib
import datetime
import logging
import os
from urllib import request


def warn_if_exist(path, filetype='file'):
    assert filetype in ('file', 'folder')
    if path is not None:
        if filetype == 'file':
            if os.path.isfile(path):
                print(f"File in <{path}> exists. Replacing.")
        else:
            if os.path.isdir(path):
                print(f"Folder in <{path}> exists. Replacing.")


def quit_if_missing(path, filetype='file'):
    assert filetype in ('file', 'folder')
    if filetype == 'file':
        if not os.path.isfile(path):
            logging.error(f"File in <{path}> missing, exiting.")
            raise Exception
    else:
        if not os.path.isdir(path):
            logging.error(f"Folder in <{path}> missing, exiting.")
            raise Exception


def setup_debug_folder(global_debug_folder):
    quit_if_missing(global_debug_folder, filetype="folder")
    timestamp = datetime.datetime.now().isoformat()
    debug_folder = os.path.join(global_debug_folder, timestamp)
    warn_if_exist(debug_folder, filetype="folder")
    os.mkdir(debug_folder)
    return debug_folder


def download_pdb_files(seq_cid_map,
                       output_folder,
                       file_suffix='.pdb',
                       url_template='https://files.rcsb.org/view/{}.pdb'):
    # This downloads the .pdb files listed in pdb_list, from rcsb server.
    if not os.path.isdir(output_folder):
        os.mkdir(output_folder)
    stored_pdb_files = set(os.listdir(output_folder))
    for pname in seq_cid_map.keys():
        pname = pname.lower()
        url = url_template.format(pname.strip())
        output_path = os.path.join(output_folder, pname+file_suffix)
        if pname+file_suffix not in stored_pdb_files:
            with contextlib.closing(request.urlopen(url)) as contents:
                with open(output_path, 'w') as output_file:
                    output_file.write(contents.read().decode("utf-8"))


def read_fasta(seq_file):
    with open(seq_file, 'r') as file:
        seq_lines = file.readlines()
    header_seq_map = dict()
    header = None
    current_seq = []
    for line in seq_lines:
        if line.startswith(">"):
            if current_seq:
                seq = "".join(current_seq)
                header_seq_map[header] = seq
                header = line[1:].strip()
                current_seq = []
            else:
                header = line[1:].strip()
        else:
            current_seq.append(line.strip())
    if header:
        seq = "".join(current_seq)
        header_seq_map[header] = seq
    return header_seq_map


def write_fasta(header_seq_map, output, line_len=60):
    with open(output, 'w') as file:
        for header, seq in header_seq_map.items():
            file.write(">" + header+"\n")
            while len(seq) > line_len:
                file.write(seq[:60]+"\n")
                seq = seq[60:]
            else:
                file.write(seq+"\n")


_aa_index = [('ALA', 'A'),
             ('CYS', 'C'),
             ('ASP', 'D'),
             ('GLU', 'E'),
             ('PHE', 'F'),
             ('GLY', 'G'),
             ('HIS', 'H'),
             ('HSE', 'H'),
             ('HSD', 'H'),
             ('ILE', 'I'),
             ('LYS', 'K'),
             ('LEU', 'L'),
             ('MET', 'M'),
             ('MSE', 'M'),
             ('ASN', 'N'),
             ('PRO', 'P'),
             ('GLN', 'Q'),
             ('ARG', 'R'),
             ('SER', 'S'),
             ('THR', 'T'),
             ('VAL', 'V'),
             ('TRP', 'W'),
             ('TYR', 'Y')]

AA3_to_AA1 = dict(_aa_index)

AA1_to_index = dict()
AA_values = sorted(set(AA3_to_AA1.values()))
assert len(AA_values) == 20
for i, AA in enumerate(AA_values):
    AA1_to_index[AA] = i
