"""
Matchers compare between a provided sequence data and those in dataset,
and return potential candidates with an associated match score.

class Matcher provide the main entry point. Call Matcher().load(df) where df
comes from descr, to load dataset. To query a particular sequence,
call query(df_of_seq). This will return, for each node/relative_sno position,
information about the sequences in the dataset that best matches the
properties of the provided sequence.

Implementation:
Matcher has 2 public methods, load and query. Load takes df of entire
dataset, query takes df of a particular sequence (so there should only be one
node per relative_sno). Requires MatcherPerSno.

MatcherPerSno.load takes in df of each node, for multiple sequences, so it is
for a particular relative_sno position. query takes in df of a node, of a seq.
Requires individual matchers.

Individual matchers should have load and query methods, which takes in single
node, multiple seqs for load, and single node, single seq for query. For load,
they should return a confidence score, between 0 and 1, that indicate how much
confidence should be placed on their prediction, based on the dataset loaded.
For query, they should return a np.array of match scores, each between 0 and 1,
for every sequence in the loaded dataset. The order of the scores should match
the order of the sequences in the loaded dataset. ###They must not change the
order of the sequences in the loaded dataset.###

PhipsiMatcher and SignatureMatcher have a public class whose load and query
extracts the desired values (phipsi and res respectively), without changing the
order in the case of load, and send it to a private class.

Todo:
Hbond as a matcher is quite tricky, because there are quite a few variables
to consider (category, length, angle, etc). Current idea is to consider only
hbonds with the same category, and disregard the rest (not yet implemented).

Will probably need to get more variables to build matchers from.
"""