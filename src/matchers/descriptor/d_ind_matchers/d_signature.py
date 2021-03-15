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

    def get_weight(self):
        return self.matcher.weight

class _SignatureMatcher:
    def __init__(self):
        self.res_distribution = None

    def load(self, res):
        self.res = res
        __, counts = np.unique(res, return_counts=True)
        counts = np.sort(counts)[::-1]
        self.weight = (counts[0] / sum(counts)) ** 2
        if len(counts) > 1:
            self.weight += (counts[1] / sum(counts)) ** 2 * 0.3
        return

    def query(self, res):
        # Consider substitutability of res
        terms, counts = np.unique(self.res, return_counts=True)
        score = 0
        for term, count in zip(terms, counts):
            p = 1 if term == res else blosum[frozenset([res, term])]
            score += p * count / np.sum(counts)
        assert score <= 1
        return self.weight, score

