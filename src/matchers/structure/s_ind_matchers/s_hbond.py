import copy
import pandas as pd

from .hb_components import CategoryComponent, BondDistComponent

# BondDistComponent
class _HbMatcher:
    def __init__(self, components=(CategoryComponent,)):
        self.components = [x() for x in components]

    def load(self, hb_df):
        self.component_weights = dict()
        for component in self.components:
            name = component.__class__.__name__
            weight = component.load(hb_df)
            self.component_weights[name] = weight
        self.matcher_weight = 0.5
        return

    def query(self, hb_datapoint):
        p_all_components = None
        for component in self.components:
            raw_p = component.query(hb_datapoint)
            name = component.__class__.__name__
            normalised_p = self._normalise_p(name, raw_p)
            if p_all_components is None:
                p_all_components = normalised_p
            else:
                p_all_components += normalised_p
        return self.matcher_weight, p_all_components

    def _normalise_p(self, name, p):
        current_weight = self.component_weights[name]
        total_weight = sum(self.component_weights.values())
        weight = current_weight / total_weight
        p = p * weight
        return p


class HbMatcher:
    def __init__(self):
        self.matcher = _HbMatcher()

    def load(self, df):
        hb_df = self._convert_to_hb_df(df)
        self.matcher.load(hb_df)
        return

    def query(self, df):
        hb_datapoint = self._convert_to_hb_df(df)
        assert len(hb_datapoint) == 1
        hb_datapoint = hb_datapoint[0]
        return self.matcher.query(hb_datapoint)

    def _convert_to_hb_df(self, df):
        individual_filedata = []
        for _index, row in df.iterrows():
            hb_dict = dict()
            hb_dict['role'] = row.role
            hb_dict['d_cid'] = row.d_cid
            hb_dict['d_res'] = row.d_res
            hb_dict['d_sno'] = row.d_sno
            hb_dict['d_aname'] = row.d_aname
            hb_dict['a_cid'] = row.a_cid
            hb_dict['a_res'] = row.a_res
            hb_dict['a_sno'] = row.a_sno
            hb_dict['a_aname'] = row.a_aname
            hb_dict['d_a_dist'] = row.d_a_dist
            hb_dict['a_d_dd'] = row.a_d_dd
            hb_dict['d_a_aa'] = row.d_a_aa
            hb_dict['planar1'] = row.planar1
            hb_dict['planar2'] = row.planar2
            hb_dict['category'] = row.category
            hb_df = pd.DataFrame.from_dict(hb_dict)
            individual_filedata.append(hb_df)
        return individual_filedata
