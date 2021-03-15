import copy
import numpy as np
import pandas as pd

from matchers.structure.s_ind_matchers.s_hbond import HbMatcher
from matchers.structure.s_ind_matchers.s_covalent import CovalentMatcher
from matchers.structure.s_ind_matchers.s_contact import ContactMatcher
from matchers.structure.s_ind_matchers.s_phipsi import PhipsiMatcher
from matchers.structure.s_ind_matchers.s_signature import SignatureMatcher
from matchers.structure.s_ind_matchers.s_vdw_new import VdwOuterMatcher
# from matchers.structure.s_ind_matchers.s_hb_new import HbondNewMatcher
from matchers.structure.s_ind_matchers.s_hbond_new2 import HbMatcher
from matchers.structure.s_ind_matchers.s_vdw_new_multiple import \
    VdwCircleFull, VdwCircleShort
from matchers.matcher_exception import Skip

class Matcher:
    def  __init__(self, cropped=True):
        self.matchers = dict()
        self.cropped = cropped

    def load(self, df):
        df_grouped = df.groupby('relative_sno')
        for relative_sno, df_per_sno in df_grouped:
            matcher = MatcherPerSno()
            try:
                matcher.load(df_per_sno)
            except Skip:
                continue
            self.matchers[relative_sno] = matcher
        return

    def query(self, queried):
        df_grouped = queried.groupby('relative_sno')
        p_all_sno = dict()
        for relative_sno, df_per_sno in df_grouped:
            try:
                matcher = self.matchers[relative_sno]
            except KeyError:
                continue
            p_per_sno = matcher.query(df_per_sno)
            if self.cropped:
                p_per_sno = p_per_sno[:3]
            p_all_sno[relative_sno] = p_per_sno
        return p_all_sno

class MatcherPerSno:
    def __init__(self, matchers=(PhipsiMatcher, SignatureMatcher, HbMatcher,
                                 VdwCircleFull, VdwCircleShort,
                                 ContactMatcher)):
        self.matchers = [x() for x in matchers]
        self.matcher_weights = dict()

    def load(self, df):
        self.identities = list(df[['filename', 'seq_marker', 'cid']].values)
        for i, identity in enumerate(self.identities):
            # cid converted to str as convention, easier to match keys later
            self.identities[i] = tuple(map(str, identity))
        self.identities = np.array(self.identities)
        for matcher in self.matchers:
            # print(f"Loading matcher {matcher.__class__.__name__}")
            matcher.load(df)
        return

    def query(self, datapoint):
        # scorers cannot change order of identity
        weight_p_raw = []
        for matcher in self.matchers:
            weight, p_raw = matcher.query(datapoint)
            assert isinstance(p_raw, np.ndarray)
            assert isinstance(weight, float)
            assert 0 <= weight <= 1
            weight_p_raw.append([weight, p_raw])
        p_all_seq = self._normalise_prob(weight_p_raw)
        order = np.argsort(-p_all_seq)
        p_selected = p_all_seq[order]
        identities = self.identities[order]
        p_matched_to_seq = list([(_id, p) for _id, p in
                                  zip(identities, p_selected)])
        return p_matched_to_seq

    def _normalise_prob(self, weight_p_raw):
        sum_weights = sum([weight_p[0] for weight_p in weight_p_raw])
        p_all_seq = None
        ref_length = len(weight_p_raw[0][1])
        for weight, p in weight_p_raw:
            assert len(p) == ref_length
            if sum_weights < 0.01:
                weight_norm = 0.99
            else:
                weight_norm = weight / sum_weights
            p_norm = p * weight_norm
            if p_all_seq is None:
                p_all_seq = p_norm
            else:
                p_all_seq = p_all_seq + p_norm
        assert max(p_all_seq) <= 1.000001   # float comparison fail otherwise
        assert min(p_all_seq) >= -0.000001  # for 1.
        return p_all_seq



