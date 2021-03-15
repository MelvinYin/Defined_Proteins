import numpy as np

from matchers.blosum import blosum
from utils.generic import AA3_to_AA1


class SignatureMatcher:
    def __init__(self):
        self.matcher = _SignatureMatcher()

    def load(self, df):
        res = list(df.res.values)
        for i in range(len(res)):
            res[i] = AA3_to_AA1[res[i]]
        self.matcher.load(res)
        return

    def query(self, df):
        res = df.res.values
        assert len(res) == 1
        res = AA3_to_AA1[res[0]]
        weight, p = self.matcher.query(res)
        return weight, p

class _SignatureMatcher:
    def __init__(self):
        self.res_distribution = None

    def load(self, res):
        self.res = res
        __, counts = np.unique(res, return_counts=True)
        counts = np.sort(counts)[::-1]
        self.weight = sum(counts[:2]) / sum(counts)
        return

    def query(self, res):
        # Consider substitutability of res
        p_all_seqs = []
        for loaded_res in self.res:
            if res == loaded_res:
                p = 1
            else:
                key = frozenset([res, loaded_res])
                p = blosum[key]
                assert p < 1
            p_all_seqs.append(p)
        p_all_seqs = np.array(p_all_seqs)
        return self.weight, p_all_seqs

