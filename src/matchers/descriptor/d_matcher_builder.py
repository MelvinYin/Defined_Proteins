import os
import pickle

from config import paths
from matchers.descriptor.d_matcher import Matcher


def load_matchers():
    input_df = os.path.join(paths.DESCRS, "mg_descr.pkl")
    output_matcher = os.path.join(paths.MATCHERS,
                                  "mg_matcher_descr.pkl")
    with open(input_df, 'rb', -1) as file:
        df = pickle.load(file)
    matcher = Matcher()
    matcher.load(df)
    with open(output_matcher, 'wb') as file:
        pickle.dump(matcher, file, -1)

    input_df = os.path.join(paths.DESCRS, "efhand_descr.pkl")
    output_matcher = os.path.join(paths.MATCHERS, "efhand_matcher_descr.pkl")
    with open(input_df, 'rb', -1) as file:
        df = pickle.load(file)
    matcher = Matcher()
    matcher.load(df)
    with open(output_matcher, 'wb') as file:
        pickle.dump(matcher, file, -1)

    input_df = os.path.join(paths.DESCRS, "GxGGxG_descr.pkl")
    output_matcher = os.path.join(paths.MATCHERS, "GxGGxG_matcher_descr.pkl")
    with open(input_df, 'rb', -1) as file:
        df = pickle.load(file)
    matcher = Matcher()
    matcher.load(df)
    with open(output_matcher, 'wb') as file:
        pickle.dump(matcher, file, -1)

    input_df = os.path.join(paths.DESCRS, "GxxGxG_descr.pkl")
    output_matcher = os.path.join(paths.MATCHERS, "GxxGxG_matcher_descr.pkl")
    with open(input_df, 'rb', -1) as file:
        df = pickle.load(file)
    matcher = Matcher()
    matcher.load(df)
    with open(output_matcher, 'wb') as file:
        pickle.dump(matcher, file, -1)

    input_df = os.path.join(paths.DESCRS, "GxGxxG_descr.pkl")
    output_matcher = os.path.join(paths.MATCHERS, "GxGxxG_matcher_descr.pkl")
    with open(input_df, 'rb', -1) as file:
        df = pickle.load(file)
    matcher = Matcher()
    matcher.load(df)
    with open(output_matcher, 'wb') as file:
        pickle.dump(matcher, file, -1)
    pass


if __name__ == "__main__":
    load_matchers()

