import contextlib
import logging
import os
import pickle
import traceback
from urllib import request
import urllib.error
import boto3
import botocore

from config import paths
from pdb_component import pdb_utils, loaders


def get_seq_for(pdb_code, cid=None):
    # logging.debug(f"{pdb_code}: {cid}")
    filedata = get_info_for(pdb_code)
    if filedata is None:
        return None
    ATOM, HETATM, hb = filedata
    del HETATM
    del hb
    if cid:
        ATOM_cid = ATOM[ATOM.cid == cid]
        seq = _extract_seq_from_df(ATOM_cid)
        return seq
    cid_seq_map = dict()
    for current_cid, ATOM_cid in ATOM.groupby("cid"):
        seq = _extract_seq_from_df(ATOM_cid)
        cid_seq_map[current_cid] = seq
    return cid_seq_map


def get_info_for(pdb_code, reset=False):
    pdb_suffix = pdb_code.lower().strip() + ".pkl"
    if reset:
        get_success = download(pdb_code, silent=False)
        if not get_success:
            print(f"download(pdb_code) failed for {pdb_code}")
            return None
        get_success = loaders.load_pdb_info(pdb_code)
        if not get_success:
            print(f"load_pdb_info(pdb_code) failed for {pdb_code}")
            return None
    else:
        if pdb_suffix not in paths.PDB_PARSED_SET:
            s3 = boto3.client('s3', aws_access_key_id="AKIAY6UR252SQUQ3OSWZ",
                              aws_secret_access_key="08LQj"
                                                    "+ryk9SMojG18vERXKKzhNSYk5pLhAjrIAVX")
            output_path = os.path.join(paths.PDB_PARSED, pdb_suffix)
            print(f"S3: {pdb_suffix}")
            with open(output_path, 'wb') as f:
                try:
                    s3.download_fileobj('definedproteins', pdb_suffix, f)
                    paths.PDB_PARSED_SET = set(os.listdir(paths.PDB_PARSED))
                    get_success = True
                    print(f"S3 Success")
                except:
                    print(f"S3 Fail")
                    if pdb_code.lower().strip() + ".pdb" in paths.PDB_FILES_SET:
                        get_success = loaders.load_pdb_info(pdb_code)
                    else:
                        get_success = download(pdb_code, silent=False)
                        if get_success:
                            paths.PDB_FILES_SET = set(os.listdir(paths.PDB_FILES))
                            get_success = loaders.load_pdb_info(pdb_code)
                        else:
                            print(f"download(pdb_code) failed for {pdb_code}")
            if not get_success:
                print(f"get_info_for(pdb_code) failed for {pdb_code}")
                return None
            paths.PDB_PARSED_SET = set(os.listdir(paths.PDB_PARSED))
    filepath = os.path.join(paths.PDB_PARSED, pdb_suffix)
    with open(filepath, 'rb') as file:
        output = pickle.load(file)
    return output

def preload_all():
    for filename in paths.PDB_FILES_SET:
        pdb_code = filename.split(".")[0]
        loaders.load_pdb_info(pdb_code)


def download_new(pdb_code, silent=False):
    pdb_code = pdb_code.lower().strip()

    url = f"https://files.rcsb.org/download/{pdb_code}.pdb"

    output_path = os.path.join(paths.PDB_FILES, pdb_code + ".pdb")
    try:
        with contextlib.closing(request.urlopen(url)) as contents:
            with open(output_path, 'w') as output_file:
                output_file.write(contents.read().decode("utf-8"))
    except urllib.error.HTTPError as e:
        if not silent:
            logging.info(f"download() fails for file {output_path}. Probably "
                         f"invalid pdb_code.")
            logging.info(f"Traceback: <{traceback.format_exc()}>")
            logging.info(f"Error_msg: <{e}>\n")
        return False
    assert os.path.isfile(output_path)
    return True

def download(pdb_code, silent=False):
    pdb_code = pdb_code.lower().strip()

    url = pdb_utils.PDB_URL_TEMPLATE.format(pdb_code)
    output_path = os.path.join(paths.PDB_FILES, pdb_code+".pdb")
    try:
        with contextlib.closing(request.urlopen(url)) as contents:
            with open(output_path, 'w') as output_file:
                output_file.write(contents.read().decode("utf-8"))
    except urllib.error.HTTPError as e:
        print(f"download() fails for file {output_path}. Probably "
             f"invalid pdb_code.")
        print(f"Traceback: <{traceback.format_exc()}>")
        print(f"Error_msg: <{e}>\n")
        if not silent:
            logging.info(f"download() fails for file {output_path}. Probably "
                         f"invalid pdb_code.")
            logging.info(f"Traceback: <{traceback.format_exc()}>")
            logging.info(f"Error_msg: <{e}>\n")
        return False
    assert os.path.isfile(output_path)
    return True


def _extract_seq_from_df(df):
    # Assumption that df is screened for cid already, so res is unique
    # Assumption that df is sorted by sno
    seq = []
    for index, row in df.iterrows():
        curr_sno = len(seq) + 1
        row_sno = row.sno
        if row_sno < curr_sno:
            continue
        elif row_sno == curr_sno:
            AA3 = row.res
            assert isinstance(AA3, str)
            try:
                AA1 = pdb_utils.AA3_to_AA1[AA3]
            except IndexError:
                AA1 = "X"
        else:
            diff = row_sno - curr_sno
            seq.extend(["X"] * diff)
            AA3 = row.res
            assert isinstance(AA3, str)
            try:
                AA1 = pdb_utils.AA3_to_AA1[AA3]
            except IndexError:
                AA1 = "X"
        seq.append(AA1)
    seq = "".join(seq)
    return seq


