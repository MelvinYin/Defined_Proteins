from collections import OrderedDict
import logging
import matplotlib.pyplot as plt
import os
import pickle
import traceback

from config import paths
from preprocess.converge import conv_interface
from pdb_component import pdb_interface
from preprocess import filter_seqs, seed_matrix_builder, clean_fasta_alphabet, \
    crop_matrix, displace_matrix
from preprocess.search import search_converters
from utils import plots


def from_aligned(input_seqs, aligned_seqs, seq_file, output, key=""):
    seed_matrix = os.path.join(paths.DEBUG, f'{key}_seed_matrix.txt')
    numseq = seed_matrix_builder.from_aligned(input_seqs, aligned_seqs,
                                              seed_matrix)
    debug_logo = os.path.join(paths.DEBUG, f"{key}_aligned_logo.png")
    plots.plot_logo_from_matrix(seed_matrix)
    plt.savefig(debug_logo)
    find(seed_matrix, numseq, seq_file, output)


def from_nbdb(input_nbdb, output, seq_file, displacement, key=""):
    """
    num_seqs don't really matter here because of normalisation
    we cannot use numseq from NBDB file, cos it's too big. So just use a
    random number, it doesn't really matter.
    """

    pssm_shortened = os.path.join(paths.DEBUG, f'{key}_pssm_shortened.txt')
    pssm_aligned = os.path.join(paths.DEBUG, f'{key}_pssm_aligned.txt')
    seed_matrix = os.path.join(paths.DEBUG, f'{key}_seed_matrix.txt')

    GxGGxG_numseq = 2000

    crop_matrix.crop(input_nbdb, pssm_shortened)
    if displacement != 0:
        displace_matrix.displace(displacement, pssm_shortened, pssm_aligned)
    else:
        pssm_aligned = pssm_shortened

    seed_matrix_builder.from_nbdb(pssm_aligned, GxGGxG_numseq, seed_matrix)
    debug_logo = os.path.join(paths.DEBUG, f"{key}_nbdb_logo.png")
    plots.plot_logo_from_matrix(seed_matrix)
    plt.savefig(debug_logo)
    find(seed_matrix, GxGGxG_numseq, seq_file, output)

def find(matrix_file, num_seqs, pdb_seq_file, output, motif_len=30):
    conv_output = paths.CONV_OUTPUT
    conv_interface.run(matrix_file, motif_len, num_seqs, conv_output)

    pdb_cid_motif_raw = search_converters.search_run(conv_output,
                                                     pdb_seq_file)
    pdb_cids = []

    for pdb_id, values in pdb_cid_motif_raw.items():
        for cid in values['cid']:
            pdb_cids.append((pdb_id, cid))

    pdb_cid_seq = dict()
    print(len(pdb_cids))
    if os.path.isfile(paths.RCSB_SEQS):
        with open(paths.RCSB_SEQS, 'rb') as file:
            rcsb_seqs = pickle.load(file)
    else:
        rcsb_seqs = dict()

    for i, (pdb_id, cid) in enumerate(pdb_cids):
        if not i % 10:
            print(i)
        if (pdb_id.upper(), cid.upper()) in rcsb_seqs:
            pdb_cid_seq[(pdb_id, cid)] = rcsb_seqs[(pdb_id.upper(), cid.upper())]
        else:
            try:
                ATOM = pdb_interface.get_info_for(pdb_id)[0]
                ATOM_cid = ATOM[ATOM.cid == cid]
                if ATOM_cid is None:
                    continue
                seq = pdb_interface._extract_seq_from_df(ATOM_cid)
                if seq is None:
                    continue
            except Exception as e:
                print(f"get_seq_for() fails for pdb_id/cid {pdb_id}/{cid}. "
                              f"Skipping.")
                print(f"Traceback: <{traceback.format_exc()}>")
                print(f"Error_msg: <{e}>")
                continue
            pdb_cid_seq[(pdb_id, cid)] = seq
            rcsb_seqs[(pdb_id.upper(), cid.upper())] = seq
    with open(paths.RCSB_SEQS, 'wb') as file:
        pickle.dump(rcsb_seqs, file, -1)
    pdb_cid_seq = OrderedDict(sorted(pdb_cid_seq.items()))
    pdb_structure_seqs = os.path.join(paths.DEBUG, "pdb_structure_seqs.txt")
    with open(pdb_structure_seqs, 'w') as file:
        for (pdb_id, cid), seq in pdb_cid_seq.items():
            file.write(f">{pdb_id}_{cid}\n")
            file.write(seq + "\n")
    clean_fasta_alphabet.screen(pdb_structure_seqs, pdb_structure_seqs)
    filter_seqs.delete_short_seqs(pdb_structure_seqs, motif_len)

    motif_positions = search_converters.search_run(conv_output,
                                                   pdb_structure_seqs)
    with open(output, 'wb') as file:
        pickle.dump(motif_positions, file, -1)