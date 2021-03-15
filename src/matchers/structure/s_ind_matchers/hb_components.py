import copy
import numpy as np
from abc import ABC, abstractmethod
from collections import defaultdict, OrderedDict

np.seterr(invalid='raise')


class HbComponentMock:
    def __init__(self):
        pass

    def load(self, df):
        self.scorer_weight = 0.8    # for now
        self.num_points = len(df)
        return self.scorer_weight

    def query(self, point_df):
        probs = np.full(self.num_points, 0.5)
        return probs


class BaseComponent(ABC):
    def __init__(self):
        self.var_name = None

    def load(self, hb_df_all_seq):
        self.var_all_seq = self._get_var(hb_df_all_seq)
        return 0.8

    def _get_var(self, df):
        var_all_seq = []
        for ind_seq_data in df:
            grouped = ind_seq_data.groupby("role")
            per_seq = dict()
            for role, df_per_role in grouped:
                per_seq[role] = df_per_role[self.var_name].values
            var_all_seq.append(per_seq)
        return var_all_seq

    def query(self, point_df):
        p_all_seqs = []
        grouped = point_df.groupby("role")
        var_selected_seq = dict()
        for role, df_per_role in grouped:
            var_selected_seq[role] = df_per_role[self.var_name].values
        for var_current_seq in self.var_all_seq:
            p = self._get_score(var_selected_seq, var_current_seq)
            p_all_seqs.append(p)
        p_all_seqs = np.array(p_all_seqs)
        return p_all_seqs

    @abstractmethod
    def _get_score(self, selected, current):
        # Both are dict
        return 0.5


class CategoryComponent(BaseComponent):
    def __init__(self):
        self.var_name = 'category'

    def _get_score(self, cat_selected, cat_current):
        # cats should be dicts, with key as role and values as category
        p = 0
        roles = ('A', 'D')
        for role in roles:
            if role not in cat_selected and role not in cat_current:
                p += 1 / len(roles)
            elif role in cat_selected and role in cat_current:
                cat1 = cat_selected[role]
                cat2 = cat_current[role]
                if len(cat1) >= len(cat2):
                    _p = self._get_itemcount_between(cat1, cat2)
                    p += (1 / len(roles)) * _p
                else:
                    _p = self._get_itemcount_between(cat2, cat1)
                    p += (1 / len(roles)) * _p
            else:
                continue
        return p

    def _get_itemcount_between(self, cat1, cat2):
        """
        Number of items that exist in both lists, divided by total number of
        items in longer (cat1) list. cat2 will be altered inplace.
        """
        _cat1 = list(cat1)
        full_len = len(_cat1)
        for item in cat1:
            if item in cat2:
                del _cat1[_cat1.index(item)]
        p = (full_len - len(_cat1)) / full_len
        return p


class TempComponent:
    def __init__(self):
        self.var_name = 'd_a_dist'
        self.bins = OrderedDict()
        self.bins[2] = []
        self.bins[3] = []

    def load(self, hb_df_all_seq):
        self.var_all_seq = self._get_var(hb_df_all_seq)
        return 0.1

    # def _get_var(self, df):
    #     var_all_seq = []
    #     for ind_seq_data in df:
    #         grouped = ind_seq_data.groupby("role")
    #         per_seq = dict()
    #         for role, df_per_role in grouped:
    #             per_seq[role] = df_per_role[self.var_name].values
    #         var_all_seq.append(per_seq)
    #     return var_all_seq

    def _get_var(self, df):
        # Here, we make a hack. We are supposed to track 1 hbond from within
        # and 1 from without loop, but for convenience for now we only track
        # hbonds in general, that's all.
        # Also, the scoring part is completely messed up.
        var_all_seq = []
        for ind_seq_data in df:
            grouped = ind_seq_data.groupby("role")
            per_seq = dict()
            for role, df_per_role in grouped:
                dists = df_per_role[self.var_name].values
                bins = defaultdict(list)
                for dist in dists:
                    if dist < 3:
                        bins[2].append(dist)
                    else:
                        bins[3].append(dist)
                per_seq[role] = bins
            var_all_seq.append(per_seq)
        # print(len(var_all_seq))
        # print(var_all_seq)
        return var_all_seq

    def _get_score(self, cat_selected, cat_current):
        # cats should be dicts, with key as role and values as category
        # selected is query, with {'A': array([2.94])}
        # current is trained, with {'A': defaultdict(<class 'list'>, {2: [
        # 2.94]})}
        p = 0
        roles = ('A', 'D')
        for role in roles:
            if role not in cat_selected and role not in cat_current:
                p += 1 / len(roles)
            elif role in cat_selected and role in cat_current:
                cat1 = cat_selected[role]
                bins = defaultdict(list)
                for term in cat1:
                    if term < 3:
                        bins[2].append(term)
                    else:
                        bins[3].append(term)
                cat2 = cat_current[role]
                for key in (2, 3):
                    try:
                        _p = min(len(bins[key]), len(cat2[key])) / max(
                            len(bins[key]), len(cat2[key]))
                    except ZeroDivisionError:
                        _p = 0.99
                    p += (1 / len(roles)) * _p / 2
            else:
                continue
        return p

    def _get_itemcount_between(self, cat1, cat2):
        """
        Number of items that exist in both lists, divided by total number of
        items in longer (cat1) list. cat2 will be altered inplace.
        """
        _cat1 = list(cat1)
        full_len = len(_cat1)
        for item in cat1:
            if item in cat2:
                del _cat1[_cat1.index(item)]
        p = (full_len - len(_cat1)) / full_len
        return p

    def query(self, point_df):
        p_all_seqs = []
        grouped = point_df.groupby("role")
        var_selected_seq = dict()
        for role, df_per_role in grouped:
            var_selected_seq[role] = df_per_role[self.var_name].values
        for var_current_seq in self.var_all_seq:
            p = self._get_score(var_selected_seq, var_current_seq)
            p_all_seqs.append(p)
        p_all_seqs = np.array(p_all_seqs)
        return p_all_seqs


class BondDistComponent(BaseComponent):
    def __init__(self):
        self.var_name = 'd_a_dist'

    def _get_score(self, dist_selected, dist_current):
        # logging.debug(dist_selected)
        p = 0
        # roles = ('A', 'D')
        # for role in roles:
        #     if role not in dist_selected and role not in dist_current:
        #         p += 1 / len(roles)
        #     elif role in dist_selected and role in dist_current:
        #         cat1 = cat_selected[role]
        #         cat2 = cat_current[role]
        #         if len(cat1) >= len(cat2):
        #             _p = self._get_itemcount_between(cat1, cat2)
        #             p += (1 / len(roles)) * _p
        #         else:
        #             _p = self._get_itemcount_between(cat2, cat1)
        #             p += (1 / len(roles)) * _p
        #     else:
        #         continue
        return p
