import numpy as np

# Identical to CovalentMatcher

class ContactMatcher:
    def __init__(self):
        self.matcher = _ContactMatcher()

    def load(self, df):
        contact = df.contact.values
        self.matcher.load(contact)
        return

    def query(self, df):
        contact = df.contact.values
        assert len(contact) == 1
        contact = contact[0]
        return self.matcher.query(contact)

class _ContactMatcher:
    def __init__(self):
        pass

    def load(self, contact):
        self.length = len(contact)
        self.contact = contact.astype(bool)
        self.weight = 0.5
        return

    def query(self, contact):
        contact = int(contact)
        # assert isinstance(contact, np.integer)
        if contact:
            p = self.contact.astype(int)
        else:
            # parentheses needed, otherwise it'll be ~ on the converted array,
            # which gives nonsensical results.
            p = (~self.contact).astype(int)
        return self.weight, p
