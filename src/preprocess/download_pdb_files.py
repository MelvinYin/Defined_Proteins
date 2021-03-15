import contextlib
import os
from urllib import request

def download(pname_cid_map,
             pdb_folder,
             file_suffix='.pdb',
             url_template='https://files.rcsb.org/view/{}.pdb'):
    stored_pdb_files = set(os.listdir(pdb_folder))
    for pname in pname_cid_map.keys():
        pname = pname.lower()
        url = url_template.format(pname.strip())
        output_path = os.path.join(pdb_folder, pname + file_suffix)
        if pname + file_suffix not in stored_pdb_files:
            with contextlib.closing(request.urlopen(url)) as contents:
                with open(output_path, 'w') as output_file:
                    output_file.write(contents.read().decode("utf-8"))

def trim_pname_cid(pname_cid_map, pdb_folder):
    fnames = []
    for fname in os.listdir(pdb_folder):
        split_fname = fname.split(".")
        assert len(split_fname) == 2
        pname = split_fname[0]
        fnames.append(pname)
    fnames = set(fnames)
    pname_list = list(pname for pname in pname_cid_map.keys())
    for pname in pname_list:
        if pname not in fnames:
            del pname_cid_map[pname]
