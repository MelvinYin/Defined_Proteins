import numpy as np


class CovalentMatcher:
    def __init__(self):
        self.matcher = _CovalentMatcher()

    def load(self, df):
        covalent = df.covalent.values
        self.matcher.load(covalent)
        return

    def query(self, df):
        covalent = df.covalent.values
        assert len(covalent) == 1
        covalent = covalent[0]
        return self.matcher.query(covalent)

class _CovalentMatcher:
    def __init__(self):
        pass

    def load(self, covalent):
        self.length = len(covalent)
        self.covalent = covalent.astype(bool)
        self.weight = 0.01
        return

    def query(self, covalent):
        assert isinstance(covalent, np.integer)
        if covalent:
            p = self.covalent.astype(int)
        else:
            # parentheses needed, otherwise it'll be ~ on the converted array,
            # which gives nonsensical results.
            p = (~self.covalent).astype(int)
        return self.weight, p
