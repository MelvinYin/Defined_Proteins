import numpy as np
import math
from sklearn.mixture import BayesianGaussianMixture as GM
import sys

np.random.seed(1)
import random
random.seed(1)

class _PhipsiMatcher:
    def __init__(self):
        # only accept from a single relative_sno, only values.
        self.to_skip = False
        self.weight_scaling_factor = 0.2 # so self.weight is not too low.
        self.q_scaling_factor = 1

    def load(self, phipsis):
        self.length = len(phipsis)
        if np.allclose(phipsis, np.full(phipsis.shape, 360)):
            self.to_skip = True
            return
        i_to_ignore = np.array(phipsis == np.array([360., 360.]))[:, 0]
        self.ignored_i = i_to_ignore
        phipsis = phipsis[~i_to_ignore]

        phipsi_median = np.median(phipsis, axis=0)
        phipsis = phipsis - phipsi_median
        phipsis[phipsis > 180] -= 360.
        phipsis[phipsis < -180] += 360.
        gm_ = GM(n_components=min(30, len(phipsis)))
        gm_.fit(X=phipsis)
        weights = gm_.weights_
        to_keep = weights > 0.05
        num_component = sum(to_keep)

        gm = GM(n_components=num_component)
        gm.fit(X=phipsis)
        precisions = gm.precisions_cholesky_

        # self.means = gm.means_
        self.phipsis = phipsis
        self.medians = phipsi_median
        weight = np.mean(precisions[:, 0, 0]) \
                 + np.mean(precisions[:, 1, 1])
        weight = weight * self.weight_scaling_factor  # for matcher weight
        self.weight = float(min(weight, 1.))


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
        if self.to_skip:
            return 0., 0.
        if np.allclose(q_phipsi, np.full(q_phipsi.shape, 360)):
            return 0., 0.

        q_phipsi = q_phipsi - self.medians
        q_phipsi[q_phipsi > 180] -= 360.
        q_phipsi[q_phipsi < -180] += 360.
        # get matcher weight
        q_fit_in_dist_log = self.gm.score_samples(np.array([q_phipsi]))
        assert len(q_fit_in_dist_log) == 1
        q_fit_in_dist_log = q_fit_in_dist_log[0]
        # 9.69 benchmark, for which max will be 1
        q_fit_in_dist = np.exp(q_fit_in_dist_log + 9.69)
        q_fit_in_dist = q_fit_in_dist * self.q_scaling_factor
        q_fit_in_dist = min(q_fit_in_dist, 1.)
        # matcher_weight = q_fit_in_dist * self.weight
        # return matcher_weight, ps
        return self.weight, q_fit_in_dist


class PhipsiMatcher:
    def __init__(self):
        self.matcher = _PhipsiMatcher()

    def load(self, df):
        phipsi = df[['phi', 'psi']].values
        self.matcher.load(phipsi)
        return

    def query(self, df):
        phipsi = df[['phi', 'psi']].values
        assert len(phipsi) == 1, len(phipsi)
        phipsi = phipsi[0]
        return self.matcher.query(phipsi)