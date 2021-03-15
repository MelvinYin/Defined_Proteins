import numpy as np
import math
from sklearn.mixture import BayesianGaussianMixture as GM
import sys

np.random.seed(1)
import random
random.seed(1)
import statistics

class _PhipsiMatcher:
    def __init__(self):
        # only accept from a single relative_sno, only values.
        self.to_skip = False
        self.weight_scaling_factor = 10 # so self.weight is not too low.
        self.q_scaling_factor = 1
        self.weight_accom_factor = 0.2

    def load(self, phipsis):
        self.length = len(phipsis)

        num_component = min(10, self.length)
        gm_ = GM(n_components=num_component)
        gm_.fit(X=phipsis)
        weights = gm_.weights_
        to_keep = weights > 0.05
        num_component = sum(to_keep)

        gm = GM(n_components=num_component)
        gm.fit(X=phipsis)
        precisions = gm.precisions_cholesky_

        # self.means = gm.means_
        self.phipsis = phipsis
        weight = np.mean(precisions[:, 0, 0]) \
                 + np.mean(precisions[:, 1, 1])
        weight = weight * self.weight_scaling_factor  # for matcher weight
        self.weight = min(weight, 1)
        self.weight *= self.weight_accom_factor
        covs = gm.covariances_
        cov_invs = np.array([np.linalg.inv(cov) for cov in covs])
        cluster_dist = gm.predict_proba(phipsis)
        self.cov_dist = np.einsum("ijk, li->ljk", cov_invs, cluster_dist)
        self.gm = gm   # for matcher weight
        # matcher_weight should be a product of the precision/clustering
        # behaviour of the distribution, and the posterior probability of the
        #  queried point. So, higher clustering but point does not belong in
        # distribution => other pressures acting on queried point => should
        # assign lower weight. Lower clustering and point belong => low
        # clustering means low pressure on point, so it shouldn't matter that
        #  much.
        return

    def query(self, q_phipsi):
        # get matcher weight
        q_fit_in_dist_log = self.gm.score_samples(np.array([q_phipsi]))

        assert len(q_fit_in_dist_log) == 1
        q_fit_in_dist_log = q_fit_in_dist_log[0]
        # 9.69 benchmark, for which max will be 1
        q_fit_in_dist = np.exp(q_fit_in_dist_log + 9.69)
        q_fit_in_dist = q_fit_in_dist * self.q_scaling_factor
        q_fit_in_dist = min(q_fit_in_dist, 1.)
        # scaling is a bit harsh
        matcher_weight = q_fit_in_dist * self.weight

        to_each_phipsis = self.phipsis - q_phipsi
        erf_arg = np.sqrt(np.einsum('ki, kij, kj -> k', to_each_phipsis,
                              self.cov_dist, to_each_phipsis))
        ps = np.array([min(1., 1 - math.erf(i)) for i in erf_arg])
        assert len(ps) == self.length
        return matcher_weight, ps




# class HbondNewMatcher:
#     def __init__(self):
#         self.index = []
#         self.num_structures = 0
#         self.matcher = _PhipsiMatcher()
#
#     def load(self, df):
#         points = []
#         dists = []
#         self.num_structures = len(df)
#         for i, (index, row) in enumerate(df.iterrows()):
#             # Each row is one structure. Not all rows have hbonds
#             for p1, p2, dist \
#                     in zip(row['planar1'], row['planar2'],
#                            row['d_a_dist']):
#                 self.index.append(i)
#                 points.append([p1, p2])
#                 dists.append(dist)
#         assert max(self.index) < self.num_structures
#         points = np.array(points)
#         self.matcher.load(points)
#         self.dists = np.array(dists, dtype=float)
#         mean = statistics.mean(self.dists)
#         stdev = statistics.stdev(self.dists)
#         self.stdev = stdev
#         self.weight = mean / stdev
#         self.weight = min(self.weight, 0.99)
#         self.weight *= 0.1
#
#     def query(self, df):
#         summed_scores = []
#         num_bonds = len(df['planar1'].values)
#         for p1, p2, dist in zip(df['planar1'].values, df['planar2'].values,
#                                 df['d_a_dist'].values):
#             p1 = p1[0]
#             p2 = p2[0]
#             dist = dist[0]
#             m_weight, matcher_score = self.matcher.query([p1, p2])
#             linear_score = np.reciprocal(1 + np.abs(self.dists - dist))
#             summed_score = (matcher_score + linear_score) / 2
#             weight = (self.weight + m_weight) / 2
#             summed_scores.append([weight, summed_score])
#         # this should be done for every bond in df['planar1']
#         structure_scores = [[] for __ in range(self.num_structures)]
#
#         for current_scores in summed_scores:
#             current_scores = current_scores[1]
#             assert len(self.index) == len(current_scores)
#             prev = None
#             current = []
#             for i in self.index:
#                 if prev is None:
#                     prev = i
#                 if prev == i:
#                     current.append(current_scores[i])
#                 else:
#                     structure_scores[i].append(current)
#                     current = [current_scores[i]]
#                     prev = i
#         # Selection of max bond score
#         structure_summed = np.zeros(self.num_structures, dtype=float)
#         for i, scores in enumerate(structure_scores):
#             if not scores:
#                 structure_summed[i] = 1 / ((1 + math.log(1 + num_bonds)) ** 2)
#             else:
#                 discarded = set()
#                 while len(discarded) != num_bonds:
#                     discar
#         print(num_bonds)
#         print(structure_scores)
#
#         import sys
#         sys.exit()
#         return self.matcher.query(phipsi)