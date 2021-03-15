_aa_index = [('ALA', 'A'), ('CYS', 'C'), ('ASP', 'D'), ('GLU', 'E'),
             ('PHE', 'F'), ('GLY', 'G'), ('HIS', 'H'), ('HSE', 'H'),
             ('HSD', 'H'), ('ILE', 'I'), ('LYS', 'K'), ('LEU', 'L'),
             ('MET', 'M'), ('MSE', 'M'), ('ASN', 'N'), ('PRO', 'P'),
             ('GLN', 'Q'), ('ARG', 'R'), ('SER', 'S'), ('THR', 'T'),
             ('VAL', 'V'), ('TRP', 'W'), ('TYR', 'Y')]

AA3_to_AA1 = dict(_aa_index)

PDB_URL_TEMPLATE = 'https://files.rcsb.org/view/{}.pdb'