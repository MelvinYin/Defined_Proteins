import logging
import os

from preprocess.parsers import pdb_parser


class Loader:
    def __init__(self, filepath):
        assert isinstance(filepath, str) and os.path.isfile(filepath)
        self.filepath = filepath

    def _load_file(self, filepath):
        with open(filepath) as file:
            lines = [line for line in file]
        return lines

    def _load_parser(self, parser_name):
        assert isinstance(parser_name, str)
        if hasattr(pdb_parser, parser_name):
            parser = getattr(pdb_parser, parser_name)
            assert issubclass(parser, pdb_parser.BaseParser), \
                f"{parser_name} is not a parser in module pdb_parser"
        else:
            logging.error(f"{parser_name} is not in module pdb_parser")
            raise Exception
        return parser

    def parse_with(self, parsername):
        parser = self._load_parser(parsername)
        if parsername == 'HbondParser':
            loaded_parser = parser(self.filepath)
        else:
            filelines = self._load_file(self.filepath)
            loaded_parser = parser(filelines)
        parsed = loaded_parser.parsed
        return parsed
