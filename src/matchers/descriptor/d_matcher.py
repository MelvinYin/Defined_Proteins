import numpy as np

from matchers.descriptor.d_ind_matchers.d_phipsi import PhipsiMatcher
from matchers.descriptor.d_ind_matchers.d_signature import SignatureMatcher
from matchers.matcher_exception import Skip


class Matcher:
    def __init__(self):
        self.matchers = dict()

        # cheap hack to tune signature_matcher to work for different descr
        self.matcher_weights = dict()

    def load(self, df):
        df_grouped = df.groupby('relative_sno')
        for relative_sno, df_per_sno in df_grouped:
            matcher = MatcherPerSno()
            try:
                matcher.load(df_per_sno)
            except Skip:
                continue
            self.matchers[relative_sno] = matcher

        # cheap hack
        self.matcher_weights['phipsi'] = 1
        self.matcher_weights['signature'] = 1
        for relative_sno, df_per_sno in df_grouped:
            del df_per_sno
            sig_weight = self.matchers[relative_sno].get_signature_weight()
            self.matcher_weights['signature'] += sig_weight
        self.matcher_weights['signature'] /= 5
        # self.matcher_weights['signature'] = min(
        #     self.matcher_weights['signature'], 1.)
        for relative_sno, df_per_sno in df_grouped:
            del df_per_sno
            self.matchers[relative_sno].set_matcher_weight(self.matcher_weights)
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
            p_all_sno[relative_sno] = p_per_sno
        return p_all_sno


class MatcherPerSno:
    # def __init__(self, matchers=(PhipsiMatcher, SignatureMatcher, HbMatcher,gxggxg
    #                              CovalentMatcher, ContactMatcher)):
    def __init__(self, matchers=(PhipsiMatcher, SignatureMatcher)):
        self.matchers = dict()
        self.matchers['phipsi'] = PhipsiMatcher()
        self.matchers['signature'] = SignatureMatcher()
        self.matcher_weights = None

    def set_matcher_weight(self, matcher_weights):
        self.matcher_weights = matcher_weights

    def get_signature_weight(self):
        return self.matchers['signature'].get_weight()

    def load(self, df):
        self.identities = list(df[['filename', 'seq_marker', 'cid']].values)
        for i, identity in enumerate(self.identities):
            # cid converted to str as convention, easier to match keys later
            self.identities[i] = tuple(map(str, identity))
        self.identities = np.array(self.identities)
        for matcher in self.matchers.values():
            matcher.load(df)
        return

    def query(self, datapoint):
        # scorers cannot change order of identity
        weight_p_raw = []
        for matcher_code, matcher in self.matchers.items():
            weight, p_raw = matcher.query(datapoint)
            assert isinstance(weight, float)
            assert isinstance(p_raw, float)
            assert 0 <= weight <= 1
            assert self.matcher_weights is not None, "matcher_weights not set"
            weight /= self.matcher_weights[matcher_code]
            weight_p_raw.append([weight, p_raw])
        p_all_seq = self._sum_prob(weight_p_raw)
        return p_all_seq

    def _sum_prob(self, weight_p_raw):
        p_full = 0
        for weight, p in weight_p_raw:
            assert weight <= 1.001
            assert weight >= -0.0001
            assert p <= 1.001
            assert p >= -0.0001
            p_full += p * weight
        p_full /= len(weight_p_raw)
        return p_full