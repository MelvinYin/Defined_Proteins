import os
import logging
import traceback

# import numpy as np

from pdb_component.parsers import loader
import pickle
from config import paths
# from descr.parsers import pdb_list_parser

# pylint: disable=invalid-name
# def load_pdb_info(motif_map, pdb_dir=paths.PDB_FOLDER):
#     """
#     Original form have these lines: "ATOM", "ANISOU", "HETATM", "TER", but
#     ANISOU and TER are not used, so only keeping ATOM and HETATM.
#     """
#     pdb_info_data_map = dict()
#     preloaded_pdb_files = set(os.listdir(paths.PRELOADED_PDB_FOLDER))
#     for i, (p_name, properties) in enumerate(motif_map.items()):
#         logging.info(f"{i}: {p_name}")
#         pdb_id = p_name.lower()
#         sno_markers, cid = properties['sno_markers'], properties['cid']
#         if pdb_id+".pkl" in preloaded_pdb_files:
#             preloaded_filepath = os.path.join(paths.PRELOADED_PDB_FOLDER,
#                                               pdb_id + ".pkl")
#             with open(preloaded_filepath, 'rb') as file:
#                 ATOM, HETATM, hb = pickle.load(file)
#         else:
#             filepath = os.path.join(pdb_dir, p_name + ".pdb")
#             try:
#                 file_data = _load_data(filepath)
#             except Exception as e:
#                 logging.error(f"_load_data() fails for file {filepath}. Skipping.")
#                 logging.error(f"Traceback: <{traceback.format_exc()}>")
#                 logging.error(f"Error_msg: <{e}>\n\n")
#                 continue
#             ATOM, HETATM, hb = file_data
#         for marker in sno_markers:
#             pdb_info_data_map[(p_name, marker, cid)] = (ATOM, HETATM, hb)
#     return pdb_info_data_map


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


def load_pdb_info(pdb_code):
    pdb_code = pdb_code.lower().strip()
    filepath = os.path.join(paths.PDB_FILES, pdb_code + ".pdb")
    try:
        file_data = _load_data(filepath)
    except Exception as e:
        print(f"_load_pdb_info() fails for file {pdb_code}. Skipping.")
        print(f"Traceback: <{traceback.format_exc()}>")
        print(f"Error_msg: <{e}>\n\n")
        return False
    output_path = os.path.join(paths.PDB_PARSED, pdb_code + ".pkl")
    with open(output_path, 'wb') as file:
        pickle.dump(file_data, file, -1)
    return True

def _load_data(file_path):
    """
    For AA3_to_AA1, copy makes it much faster, because otherwise df will
    update itself with every change.
    """
    pdf_files = loader.Loader(file_path)

    ATOM = pdf_files.parse_with('ATOMParser')
    MODRES = pdf_files.parse_with('MODRESParser')
    HETATM = pdf_files.parse_with('HETATMParser')
    try:
        hb = pdf_files.parse_with('HbondParser')
    except:
        hb = None
    if not MODRES.empty:
        ATOM_res = ATOM.res.copy()
        ATOM.res = _MODRES_sub(ATOM_res, MODRES.res, MODRES.std_res_name)
    if HETATM.empty:
        HETATM = None
    return ATOM, HETATM, hb

# pylint: enable=invalid-name
