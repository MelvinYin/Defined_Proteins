import numpy as np

class MockMatcher:
    # Purely as a template. Can change to abc in the future.
    def __init__(self):
        self.matcher = _MockMatcher()

    def load(self, df):
        filename = df.filename.values
        return self.matcher.load(filename)

    def query(self, df):
        filename = df.filename.values
        return self.matcher.query(filename)

class _MockMatcher:
    def __init__(self):
        pass

    def load(self, filename):
        self.length = len(filename)
        return 0.5

    def query(self, filename):
        p = np.zeros(self.length)
        return p

