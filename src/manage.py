#!/usr/bin/env python
"""Django's command-line utility for administrative tasks."""
import os
import sys
from config import paths

def main(*args, **kwargs):
    paths.initialise(os.path.join(os.path.dirname(os.path.realpath(__file__)), "__main__"))
    assert os.path.isfile(paths.SEARCH_EXEC), paths.SEARCH_EXEC
    # import subprocess
    # from config import paths
    # search_input = os.path.join(paths.SEARCH_DIR, "extracted.matrix")
    # search_output = os.path.join(paths.SEARCH_DIR, "output.txt")
    # pdb_seq_file = paths.RCSB_SEQS_FASTA
    # # assert os.path.isfile(desired)
    # assert os.path.isfile(search_input)
    # assert os.path.isfile(pdb_seq_file)
    # command = f"{desired} {search_input} {pdb_seq_file} " \
    #           f"{search_output} 30"
    # ret_code = subprocess.run(command, shell=True)
    # print(ret_code)
    # import sys
    # sys.exit()

    os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'site_main.settings')
    try:
        from django.core.management import execute_from_command_line
    except ImportError as exc:
        raise ImportError(
            "Couldn't import Django. Are you sure it's installed and "
            "available on your PYTHONPATH environment variable? Did you "
            "forget to activate a virtual environment?"
        ) from exc
    execute_from_command_line(sys.argv)


if __name__ == '__main__':
    main()
