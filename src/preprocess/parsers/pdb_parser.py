import abc
import os
import re
import subprocess

import numpy as np
import pandas as pd

from config import paths

# pylint: disable=anomalous-backslash-in-string, arguments-differ
class BaseParser(abc.ABC):
    @abc.abstractmethod
    def split_line(self, arg):
        """
        :return: Dict[str, str]
        """
        return

    @abc.abstractmethod
    def convert_type(self, arg):
        """
        :return: Dict[str, Union[int, float, str]]
        """
        return

    @abc.abstractmethod
    def check_validity(self, arg):
        """
        :return: Bool
        """
        return

    @abc.abstractmethod
    def remove_tmp_keys(self, arg):
        return

    def to_df(self, data):
        """
        :return: pd.DataFrame
        """
        if data:
            df_index = data[0].keys()
            df = pd.DataFrame(data, columns=df_index)
        else:
            # empty df allowing merging to work properly later
            df = pd.DataFrame()
        return df

    def parse_filedata(self, file_data, record_name):
        """
        :return: pd.DataFrame
        """
        data_list = []
        for line in file_data:
            if line.startswith(record_name):
                splitted = self.split_line(line)
                # if not self.check_validity(splitted):
                #     continue
                screened = self.remove_tmp_keys(splitted)
                converted = self.convert_type(screened)
                data_list.append(converted)
        parsed = self.to_df(data_list)
        return parsed

###############################################################################

class ATOMParser(BaseParser):
    def __init__(self, file_data):
        record_name = "ATOM  "
        self._check_for_duplicate = []
        self.parsed = self.parse_filedata(file_data, record_name)

    def split_line(self, line):
        """
        Splitting via other ways (e.g. re, split) that depend on there being a
        whitespace in between entries does not work for some files.
        :param line: str
        :return: Dict[str, str]
        """
        splitted = dict(_record_name=line[:6],
                        ano=line[6:11],
                        aname=line[12:16],
                        _altLoc=line[16],
                        res=line[17:20],
                        cid=line[21],
                        sno=line[22:26],
                        _AChar=line[27],
                        coord=tuple([line[31:38], line[39:46], line[47:54]]),
                        occupancy=line[55:60],
                        tempfactor=line[61:66],
                        elementsymbol=line[77:78],
                        _charge=line[79:80])

        for key, value in splitted.items():
            if key != "coord":
                splitted[key] = value.strip()
            else:
                splitted[key] = tuple(x.strip() for x in splitted["coord"])
        return splitted

    def convert_type(self, line):
        """
        :param line: Dict[str, str]
        :return: Dict[str, Union[int, float, str]]
        """
        line["ano"] = int(line["ano"])
        line["sno"] = int(line["sno"])
        line["coord"] = np.array([float(i) for i in line["coord"][:3]])
        line["occupancy"] = float(line["occupancy"])
        line["tempfactor"] = float(line["tempfactor"])
        return line

    def check_validity(self, line):
        """
        :param line: Dict[str, str]
        :return: Bool
        """
        amino_acids = frozenset(["ALA", "CYS", "ASP", "GLU", "PHE", "GLY",
                                 "HIS", "ILE", "LYS", "LEU", "MET", "ASN",
                                 "PRO", "GLN", "ARG", "SER", "THR", "VAL",
                                 "TRP", "TYR"])

        # strip() changes " " to ""
        assert line["_AChar"] == ""

        # terminating atom, but may have other cid below.
        if line["aname"] == "OXT":
            return False

        if line["_altLoc"] != " ":
            values = (line['aname'], line['res'], line['cid'], line['sno'])
            if values not in self._check_for_duplicate:
                self._check_for_duplicate.append(values)
            else:
                return False

        assert re.fullmatch("[0-9]+", line["ano"]) \
               and line["ano"] != "0"
        assert re.fullmatch("[NCOSHP][ABGDEZH]?[XT]?[T]?([0-3][0-3]?)?",
                            line["aname"]), line
        assert line["res"] in amino_acids
        assert re.fullmatch("[A-Z0-9]", line["cid"]), line['cid']
        assert re.fullmatch("-?[0-9]+", line["sno"]), line["sno"]
        assert re.fullmatch("-?[0-9]{1,3}\.[0-9]{3}", line["coord"][0])
        assert re.fullmatch("-?[0-9]{1,3}\.[0-9]{3}", line["coord"][1])
        assert re.fullmatch("-?[0-9]{1,3}\.[0-9]{3}", line["coord"][2])
        assert re.fullmatch("[01]\.[0-9]{2}", line["occupancy"])
        assert re.fullmatch("[0-9]{1,2}\.[0-9]{2}", line["tempfactor"])
        assert re.fullmatch("[NCOSH]", line["elementsymbol"])
        assert re.fullmatch("[1-9+\-]?", line["_charge"])
        return True

    def remove_tmp_keys(self, line):
        del line["_record_name"]
        del line["_altLoc"]
        del line["_AChar"]
        del line["_charge"]
        assert len(line) == 9
        return line

###############################################################################

class HETATMParser(BaseParser):
    def __init__(self, file_data):
        record_name = "HETATM"
        self._check_for_duplicate = []
        self.parsed = self.parse_filedata(file_data, record_name)

    # same as ATOM
    def split_line(self, line):
        """
        :param line: str
        :return: Dict[str, str]
        """
        splitted = dict(_record_name=line[:6],
                        ano=line[6:11],
                        aname=line[12:16],
                        _altLoc=line[16],
                        res=line[17:20],
                        cid=line[21],
                        sno=line[22:26],
                        _AChar=line[27],
                        coord=tuple([line[31:38], line[39:46], line[47:54]]),
                        occupancy=line[55:60],
                        tempfactor=line[61:66],
                        elementsymbol=line[77:78],
                        _charge=line[79:80])

        for key, value in splitted.items():
            if key != "coord":
                splitted[key] = value.strip()
            else:
                splitted[key] = tuple(x.strip() for x in splitted["coord"])
        return splitted

    # same as ATOM
    def convert_type(self, line):
        """
        :param line: Dict[str, str]
        :return: Dict[str, Union[int, float, str]]
        """
        line["ano"] = int(line["ano"])
        line["sno"] = int(line["sno"])
        line["coord"] = np.array([float(i) for i in line["coord"][:3]])
        line["occupancy"] = float(line["occupancy"])
        line["tempfactor"] = float(line["tempfactor"])
        return line

    def check_validity(self, line):
        """
        :param line: Dict[str, str]
        :return: Bool
        """
        assert line["_AChar"] == ""

        if line["_altLoc"] != " ":
            values = (line['aname'], line['res'],
                      line['cid'], line['sno'])
            if values not in self._check_for_duplicate:
                self._check_for_duplicate.append(values)
            else:
                return False

        assert re.fullmatch("[0-9]+", line["ano"])
        assert line["ano"] != "0"
        # No aname check.
        # No residue check.
        assert re.fullmatch("[A-Z0-9]", line["cid"])
        assert re.fullmatch("[0-9]+", line["sno"]), line["sno"]
        assert re.fullmatch("-?[0-9]{1,3}\.[0-9]{3}", line["coord"][0])
        assert re.fullmatch("-?[0-9]{1,3}\.[0-9]{3}", line["coord"][1])
        assert re.fullmatch("-?[0-9]{1,3}\.[0-9]{3}", line["coord"][2])
        assert re.fullmatch("[01]\.[0-9]{2}", line["occupancy"])
        assert re.fullmatch("[0-9]{1,2}\.[0-9]{2}", line["tempfactor"])
        assert re.fullmatch("[1-9+\-]?", line["_charge"])
        # No elementsymbol check.
        return True

    def remove_tmp_keys(self, line):
        del line["_record_name"]
        del line["_altLoc"]
        del line["_AChar"]
        del line["_charge"]
        assert len(line) == 9
        return line

###############################################################################

class HbParser(BaseParser):
    def __init__(self, file_data):
        """
        Check self.record_name with file data manually to make sure it tally,
        otherwise lines will get skipped.
        """
        record_name = ("A", "B", "P", "I")
        self.parsed = self.parse_filedata(file_data, record_name)

    def split_line(self, line):
        """
        d = Donor, a = Acceptor
        :param line: str
        :return: Dict[str, str]
        """
        splitted = dict(d_cid=line[0],
                        d_sno=line[1:5],
                        _d_insertion_code=line[5],
                        d_res=line[6:9],
                        d_aname=line[10:14],

                        a_cid=line[14],
                        a_sno=line[15:19],
                        _a_insertion_code=line[19],
                        a_res=line[20:23],
                        a_aname=line[24:28],

                        d_a_dist=line[28:33],
                        atom_category=line[33:36],
                        d_a_gap=line[36:39],
                        d_a_CA_dist=line[40:45],
                        d_a_angle_at_H=line[46:51],
                        H_a_dist=line[52:57],
                        H_a_antecedent_angle_at_a=line[58:63],
                        d_a_antecedent_angle_at_a=line[64:69],
                        hb_count=line[70:75])

        for key, value in splitted.items():
            splitted[key] = value.strip()
        return splitted

    def convert_type(self, line):
        """
        :param line: Dict[str, str]
        :return: Dict[str, Union[int, float, str]]
        """
        line["d_sno"] = int(line["d_sno"])
        line["a_sno"] = int(line["a_sno"])
        line["d_a_dist"] = float(line["d_a_dist"])
        line["d_a_gap"] = int(line["d_a_gap"])
        line["d_a_CA_dist"] = float(line["d_a_CA_dist"])
        line["d_a_angle_at_H"] = float(line["d_a_angle_at_H"])
        line["H_a_dist"] = float(line["H_a_dist"])
        line["H_a_antecedent_angle_at_a"] \
            = float(line["H_a_antecedent_angle_at_a"])
        line["d_a_antecedent_angle_at_a"] \
            = float(line["d_a_antecedent_angle_at_a"])
        line["hb_count"] = int(line["hb_count"])
        return line

    def check_validity(self, line):
        """
        d_aname = {OXT|O2A|O2B|N4B|O3A|O3C|O5C|O6C|O3C|O1S|O1P}
        I in cid from 1kap
        Amino acid guesses:

        NAG = N-acetylglucosamine
        NDG = 2-(ACETYLAMINO)-2-DEOXY-A-D-GLUCOPYRANOSE
        https://www.rcsb.org/ligand/NDG
        What's ACG? http://www.rcsb.org/structure/2YMV

        Pretty sure O2' in aname is an error

        :param line: Dict[str, str]
        :return: Bool
         """
        # ZN from 1h71, which has P as starting
        # pylint: disable=invalid-name
        amino_acids_and_H2O = frozenset(["ALA", "CYS", "ASP", "GLU", "PHE",
                                         "GLY", "HIS", "ILE", "LYS", "LEU",
                                         "MET", "ASN", "PRO", "GLN", "ARG",
                                         "SER", "THR", "VAL", "TRP", "TYR",
                                         "HOH", "NAG", "NDG", "ACG", "GDP",
                                         "SO4", "MES", 'MIS', 'MPD', 'BMA',
                                         'NO3', 'GOL', 'ZN'])
        # pylint: enable=invalid-name

        assert line["_d_insertion_code"] == "-"
        assert line["_a_insertion_code"] == "-"

        assert re.fullmatch("[ABPI]", line["d_cid"])
        assert re.fullmatch("[ABPI]", line["a_cid"])
        assert re.fullmatch("[0-9]{4}", line["d_sno"])
        assert re.fullmatch("[0-9]{4}", line["a_sno"])
        assert line["d_res"] in amino_acids_and_H2O
        assert line["a_res"] in amino_acids_and_H2O
        assert re.fullmatch(
            "OXT|O2'|[NCOSH][0-7]?[NCOSHP]?[ABCDE]?[ABGDEZH]?([0-7][0-7]?)?"
            , line["d_aname"])
        assert re.fullmatch(
            "OXT|O2'|[NCOSH][0-7]?[NCOSHP]?[ABCDE]?[ABGDEZH]?([0-7][0-7]?)?"
            , line["a_aname"])

        assert re.fullmatch("[0-9]{1}\.[0-9]{2}", line["d_a_dist"])
        assert re.fullmatch("[MSH][MSH]", line["atom_category"])
        assert re.fullmatch("-1.00|[0-9]{1,2}\.[0-9]{2}"
                            , line["d_a_CA_dist"])
        assert re.fullmatch("-1.0|[0-9]{2,3}\.[0-9]{1}"
                            , line["d_a_angle_at_H"])
        assert re.fullmatch("-1.00|[0-9]{1}\.[0-9]{2}", line["H_a_dist"])
        assert re.fullmatch("-1.0|[0-9]{2,3}\.[0-9]{1}"
                            , line["H_a_antecedent_angle_at_a"])
        assert re.fullmatch("-1.0|[0-9]{2,3}\.[0-9]{1}"
                            , line["d_a_antecedent_angle_at_a"])
        assert re.fullmatch("[0-9]{1,4}", line["hb_count"])
        return True

    def remove_tmp_keys(self, line):
        del line["_d_insertion_code"]
        del line["_a_insertion_code"]
        assert len(line) == 17
        return line


###############################################################################

class HbondParser(BaseParser):
    def __init__(self, filepath):
        assert isinstance(filepath, str)
        record_name = ("query")
        file_data = self.get_hbond_data(filepath)
        self.parsed = self.parse_filedata(file_data, record_name)

    def get_hbond_data(self, filepath):
        filepath = os.path.abspath(filepath)
        hbonds_fname = "query.hbonds"
        pdb_fname = os.path.abspath(filepath)
        hb_exec = os.path.join(paths.ROOT, 'src', 'parsers', 'hb',
                               'hb_calculator')
        process = subprocess.run([hb_exec, pdb_fname, hbonds_fname])
        if process.returncode != 0:
            msg = "Could not prepare hydrogen bonds results for %s" % pdb_fname
            raise RuntimeError(msg)
        with open(hbonds_fname, 'r') as f:
            file_data = [line for line in f]
        os.remove(hbonds_fname)
        return file_data

    def split_line(self, line):
        splitted = line.split("\t")
        columns = ['pdb_id', 'd_cid', 'd_res', 'd_sno', '_icode1',
                   'd_aname', 'a_cid', 'a_res', 'a_sno', '_icode2',
                   'a_aname', 'd_a_dist', 'a_d_dd', 'd_a_aa', 'planar1',
                   'planar2', 'atom_category']

        splitted = dict((a, b) for a, b in zip(columns, splitted))
        for key, value in splitted.items():
            splitted[key] = value.strip()
        return splitted

    def convert_type(self, line):
        line["d_sno"] = int(line["d_sno"])
        line["a_sno"] = int(line["a_sno"])
        line["d_a_dist"] = float(line["d_a_dist"])
        line["d_a_aa"] = float(line["d_a_aa"])

        line["a_d_dd"] = float(line["a_d_dd"])
        line["planar1"] = float(line["planar1"])
        line["planar2"] = float(line["planar2"])
        return line

    def check_validity(self, line):
        del line
        return True

    def remove_tmp_keys(self, line):
        for key in list(line.keys()):
            if key.startswith("_"):
                del line[key]
        assert len(line) == 15
        return line


###############################################################################

class MODRESParser(BaseParser):
    def __init__(self, file_data):
        record_name = "MODRESParser"
        self.parsed = self.parse_filedata(file_data, record_name)

    def split_line(self, line):
        splitted = dict(_record_name=line[:6],
                        id_code=line[7:11],
                        res=line[12:15],
                        cid=line[16],
                        sno=line[18:22],
                        _AChar=line[22],
                        std_res_name=line[24:27],
                        description=line[29:70])
        for key, value in splitted.items():
            splitted[key] = value.strip()
        return splitted

    def convert_type(self, line):
        line["sno"] = int(line["sno"])
        return line

    def check_validity(self, line):
        amino_acids = frozenset(["ALA", "CYS", "ASP", "GLU", "PHE", "GLY",
                                 "HIS", "ILE", "LYS", "LEU", "MET", "ASN",
                                 "PRO", "GLN", "ARG", "SER", "THR", "VAL",
                                 "TRP", "TYR", "MIS"])
        assert line["_AChar"] == " "
        assert re.fullmatch("[0-9A-Z]+", line["id_code"])
        assert line["res"] in amino_acids
        assert re.fullmatch("[AB]", line["cid"])
        assert re.fullmatch("[0-9]+", line["sno"]), line["sno"]
        assert line["std_res_name"] in amino_acids
        assert re.fullmatch("[A-Z ]+", line["description"])
        return True

    def remove_tmp_keys(self, line):
        del line["_record_name"]
        del line["_AChar"]
        assert len(line) == 6
        return line

# pylint: enable=anomalous-backslash-in-string, arguments-differ
