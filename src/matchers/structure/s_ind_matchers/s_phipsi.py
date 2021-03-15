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
        self.weight_scaling_factor = 10 # so self.weight is not too low.
        self.q_scaling_factor = 1
        self.weight_accom_factor = 0.2

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
        num_component = min(10, self.length)
        gm_ = GM(n_components=num_component, max_iter=10000)
        gm_.fit(X=phipsis)
        weights = gm_.weights_
        to_keep = weights > 0.05
        num_component = sum(to_keep)

        gm = GM(n_components=num_component, max_iter=10000)
        gm.fit(X=phipsis)
        precisions = gm.precisions_cholesky_

        # self.means = gm.means_
        self.phipsis = phipsis
        self.medians = phipsi_median
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
        if self.to_skip:
            return 0., np.zeros(self.length)
        if np.allclose(q_phipsi, np.full(q_phipsi.shape, 360)):
            return 0., np.zeros(self.length)

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
        # scaling is a bit harsh
        matcher_weight = q_fit_in_dist * self.weight

        to_each_phipsis = self.phipsis - q_phipsi
        erf_arg = np.sqrt(np.einsum('ki, kij, kj -> k', to_each_phipsis,
                              self.cov_dist, to_each_phipsis))
        ps = np.array([1 - math.erf(i) for i in erf_arg])
        if len(ps) != self.length:
            assert any(self.ignored_i)
            i_s = np.argwhere(self.ignored_i)
            for i in i_s[::-1]:
                i = i[0]
                ps = np.concatenate([ps[:i], [0.], ps[i:]])
        assert len(ps) == self.length
        return matcher_weight, ps


class PhipsiMatcher:
    def __init__(self):
        self.matcher = _PhipsiMatcher()

    def load(self, df):
        phipsi = df[['phi', 'psi']].values
        self.matcher.load(phipsi)
        return

    def query(self, df):
        phipsi = df[['phi', 'psi']].values
        assert len(phipsi) == 1
        phipsi = phipsi[0]
        return self.matcher.query(phipsi)