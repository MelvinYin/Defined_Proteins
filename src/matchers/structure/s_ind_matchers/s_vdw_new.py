import numpy as np
import statistics

class VdwOuterMatcher:
    def __init__(self):
        self.counts = np.array([])
        self.weight = 0

    def load(self, df):
        vdw_outer = df['vdw_outer'].values
        self.counts = np.array(vdw_outer, dtype=float)
        mean = statistics.mean(self.counts)
        stdev = statistics.stdev(self.counts)
        self.stdev = stdev
        self.values = vdw_outer
        self.weight = mean / stdev
        self.weight = min(self.weight, 0.99)
        self.weight *= 0.1

    def query(self, df):
        value = df['vdw_outer'].values
        assert len(value) == 1
        value = value[0]
        return self.weight, np.reciprocal(1 + np.abs(self.counts - value))