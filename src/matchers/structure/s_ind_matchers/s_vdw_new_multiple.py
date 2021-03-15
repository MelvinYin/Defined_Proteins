import numpy as np
import statistics


class VdwCircleShort:
    def __init__(self):
        keys = ['vdw_short_c1_count', 'vdw_short_c2_count',
                'vdw_short_c3_count']
        self.matcher_map = dict()
        for key in keys:
            self.matcher_map[key] = VdwPrototype(key)

    def load(self, df):
        for key, matcher in self.matcher_map.items():
            matcher.load(df)

    def query(self, df):
        weights, scores = 0, None
        for key, matcher in self.matcher_map.items():
            weight, score = matcher.query(df)
            weights += weight
            if scores is None:
                scores = score
            else:
                scores += score
        weights /= len(self.matcher_map)
        scores /= len(self.matcher_map)
        return weights, scores


class VdwCircleFull:
    def __init__(self):
        keys = ['vdw_full_c1_count', 'vdw_full_c2_count', 'vdw_full_c3_count']
        self.matcher_map = dict()
        for key in keys:
            self.matcher_map[key] = VdwPrototype(key)

    def load(self, df):
        for key, matcher in self.matcher_map.items():
            matcher.load(df)

    def query(self, df):
        weights, scores = 0, None
        for key, matcher in self.matcher_map.items():
            weight, score = matcher.query(df)
            weights += weight
            if scores is None:
                scores = score
            else:
                scores += score

        weights /= len(self.matcher_map)
        scores /= len(self.matcher_map)
        return weights, scores

class VdwPrototype:
    def __init__(self, key):
        self.key = key
        self.counts = np.array([])
        self.weight = 0

    def load(self, df):
        vdw_outer = df[self.key].values
        self.counts = np.array(vdw_outer, dtype=float)
        mean = statistics.mean(self.counts)
        stdev = statistics.stdev(self.counts)
        self.stdev = stdev
        self.values = vdw_outer
        if stdev < 0.01:
            self.weight = 0.1
        else:
            self.weight = mean / stdev
            self.weight = min(self.weight, 0.99)
            self.weight *= 0.1

    def query(self, df):
        value = df[self.key].values
        assert len(value) == 1
        value = value[0]
        return self.weight, np.reciprocal(1 + np.abs(self.counts - value))