import os

ROOT = "/".join(os.path.dirname(__file__).split("/")[:-2])
RAW_ROOT = ROOT
BASH_EXEC = "/bin/bash"

DATA = os.path.join(ROOT, 'data')
SRC = os.path.join(ROOT, 'src')

DEBUG = os.path.join(DATA, 'debug')
INTERNAL = os.path.join(DATA, 'internal')
USER_INPUT = os.path.join(DATA, 'input')
USER_OUTPUT = os.path.join(DATA, 'output')

# for folder in (DEBUG, USER_OUTPUT):
#     if not os.path.isdir(folder):
#         os.mkdir(folder)

CEQLOGO_INPUT = os.path.join(INTERNAL, "ceqlogo")
DESCRS = os.path.join(INTERNAL, "descrs")
MATCHERS = os.path.join(INTERNAL, 'matchers')
MOTIF_POSITIONS = os.path.join(INTERNAL, 'motif_positions')
RCSB_SEQS = os.path.join(INTERNAL, 'rcsb_seqs.pkl')
PDB_CLUSTER = os.path.join(INTERNAL, 'pdb_cluster.txt')
PDB_FILES = os.path.join(INTERNAL, "pdb_files")
PDB_PARSED = os.path.join(INTERNAL, "pdb_files_parsed")

# for folder in (PDB_FILES, PDB_PARSED):
#     if not os.path.isdir(folder):
#         os.mkdir(folder)
#
# PDB_FILES_SET = set(os.listdir(PDB_FILES))
# PDB_PARSED_SET = set(os.listdir(PDB_PARSED))

COMPOSITION = os.path.join(INTERNAL, "composition.txt")
RCSB_SEQS_FASTA = os.path.join(USER_INPUT, "rcsb_seqs_full.txt")

PREPROCESS = os.path.join(SRC, "preprocess")
SEARCH_DIR = os.path.join(PREPROCESS, 'search')
# SEARCH_EXEC = os.path.join(SEARCH_DIR, 'search')
SEARCH_EXEC = os.path.join(SRC, "preprocess", "search", "source", "search")

CONV_FOLDER = os.path.join(PREPROCESS, 'converge')
CONV_EXEC = os.path.join(CONV_FOLDER, 'calculator')
CONV_INPUT_MATRIX = os.path.join(CONV_FOLDER, 'input_matrix.txt')
CONV_OUTPUT = os.path.join(CONV_FOLDER, 'converged_matrix.txt')
PROTEOME_BINARY = os.path.join(INTERNAL, 'proteome_binary')

HB_EXEC = os.path.join(SRC, "pdb_component", "parsers", "hb", "hb_calculator")
def initialise(base_dir):
    ROOT = "/".join(base_dir.split("/")[:-2])
    assert os.path.isdir(ROOT)
    BASH_EXEC = "/bin/bash"

    DATA = os.path.join(ROOT, 'data')
    SRC = os.path.join(ROOT, 'src')
    assert os.path.isdir(DATA)
    assert os.path.isdir(SRC)

    DEBUG = os.path.join(DATA, 'debug')
    INTERNAL = os.path.join(DATA, 'internal')
    USER_INPUT = os.path.join(DATA, 'input')
    USER_OUTPUT = os.path.join(DATA, 'output')

    CEQLOGO_INPUT = os.path.join(INTERNAL, "ceqlogo")
    DESCRS = os.path.join(INTERNAL, "descrs")
    MATCHERS = os.path.join(INTERNAL, 'matchers')
    MOTIF_POSITIONS = os.path.join(INTERNAL, 'motif_positions')
    PDB_CLUSTER = os.path.join(INTERNAL, 'pdb_cluster.txt')
    PDB_FILES = os.path.join(INTERNAL, "pdb_files")
    PDB_PARSED = os.path.join(INTERNAL, "pdb_files_parsed")

    PDB_FILES_SET = set(os.listdir(PDB_FILES))
    PDB_PARSED_SET = set(os.listdir(PDB_PARSED))

    COMPOSITION = os.path.join(INTERNAL, "composition.txt")
    RCSB_SEQS_FASTA = os.path.join(USER_INPUT, "rcsb_seqs_full.txt")

    PREPROCESS = os.path.join(SRC, "preprocess")
    SEARCH_DIR = os.path.join(PREPROCESS, 'search')
    # SEARCH_EXEC = os.path.join(SEARCH_DIR, 'search')
    SEARCH_EXEC = os.path.join(SRC, "preprocess", "search", "search")

    CONV_FOLDER = os.path.join(PREPROCESS, 'converge')
    CONV_EXEC = os.path.join(CONV_FOLDER, 'calculator')
    CONV_INPUT_MATRIX = os.path.join(CONV_FOLDER, 'input_matrix.txt')
    CONV_OUTPUT = os.path.join(CONV_FOLDER, 'converged_matrix.txt')
    PROTEOME_BINARY = os.path.join(INTERNAL, 'proteome_binary')

    HB_EXEC = os.path.join(SRC, "pdb_component", "parsers", "hb",
                           "hb_calculator")
    globals().update(locals())
    assert os.path.isfile(CONV_EXEC)
    assert os.path.isfile(SEARCH_EXEC)
