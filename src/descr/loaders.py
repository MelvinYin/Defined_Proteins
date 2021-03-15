import os
import logging
import traceback

# import numpy as np

from config import paths
from parsers import loader
# from descr.parsers import pdb_list_parser

# pylint: disable=invalid-name
def load_pdb_info(motif_map, pdb_dir=paths.PDB_FILES):
    """
    Original form have these lines: "ATOM", "ANISOU", "HETATM", "TER", but
    ANISOU and TER are not used, so only keeping ATOM and HETATM.
    """
    pdb_info_data_map = dict()
    for i, (p_name, properties) in enumerate(motif_map.items()):
        logging.info(f"{i}: {p_name}")
        filepath = os.path.join(pdb_dir, p_name + ".pdb")
        sno_markers, cid = properties['sno_markers'], properties['cid']
        try:
            file_data = _load_data(filepath)
        except Exception as e:
            print(f"_load_data() fails for file {filepath}. Skipping.")
            print(f"Traceback: <{traceback.format_exc()}>")
            print(f"Error_msg: <{e}>")
            continue
        ATOM, HETATM, hb = file_data
        for marker in sno_markers:
            pdb_info_data_map[(p_name, marker, cid)] = (ATOM, HETATM, hb)
    return pdb_info_data_map


def _inplace_AA3_substitution(res, sno, AA3_TO_AA1):
    """
    Not used, not deprecated.
    """
    assert len(res) == len(sno)

    for i in range(len(res)):
        try:
            res.iloc[i] = AA3_TO_AA1[res.iloc[i]]
        except KeyError:
            res.iloc[i] = "X"
            sno.iloc[i] = 0
    return res, sno

def _MODRES_sub(main_res, MODRES_res, MODRES_std_res_name):
    for ATOM_i, ATOM_res in enumerate(main_res):
        for MODRES_i, res in enumerate(MODRES_res):
            if ATOM_res == res:
                main_res[ATOM_i] = MODRES_std_res_name[MODRES_i]
    return main_res

def _load_data(file_path_no_suffix):
    """
    For AA3_to_AA1, copy makes it much faster, because otherwise df will
    update itself with every change.
    """
    pdf_files = loader.Loader(file_path_no_suffix)
    # hb_files = Loader(file_path_no_suffix + ".hb2")

    ATOM = pdf_files.parse_with('ATOMParser')
    MODRES = pdf_files.parse_with('MODRESParser')
    HETATM = pdf_files.parse_with('HETATMParser')

    if not MODRES.empty:
        ATOM_res = ATOM.res.copy()
        ATOM.res = _MODRES_sub(ATOM_res, MODRES.res, MODRES.std_res_name)
    if HETATM.empty:
        HETATM = None
    try:
        hb = pdf_files.parse_with('HbondParser')
    except:
        hb = None
    return ATOM, HETATM, hb

# pylint: enable=invalid-name
