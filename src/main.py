from config import paths
import os
import pickle
import matplotlib.pyplot as plt
from collections import defaultdict
import numpy as np
from utils import plots
# with open(os.path.join(paths.INTERNAL, 'descrs', "efhand_descr.pkl"),
#           'rb') as file:
#     data = pickle.load(file)
# unique_sno = sorted(list(set(data['relative_sno'].values)))
# print(data.columns)
#
# for sno in unique_sno:
#     print(sno)
#     # print(data[data['relative_sno'] == sno])
#     for i, row in data[data['relative_sno'] == sno].iterrows():
#         print(row['planar1'])
#         print(row['planar2'])
#         print(row['d_a_dist'])
#         print("")
    # print(data[data['relative_sno'] == sno]['planar1'])
    # print(data[data['relative_sno'] == sno]['planar2'])
    # print(data[data['relative_sno'] == sno]['d_a_dist'])

# scoring function: (1/log(diff) + sum(all_scores)) / log(max_bond_count)
from preprocess import find_motif_positions
from descr import descr_main
paths.initialise(os.path.join(os.path.dirname(__file__), "random"))
pdb_seq_file = paths.RCSB_SEQS_FASTA
# matrix_file = os.path.join(paths.USER_INPUT, "GxGxxG_pssm.txt")
matrix_file = os.path.join(paths.USER_INPUT, "GxxGxG_pssm.txt")
aligned_file = os.path.join(paths.USER_INPUT, "efhand_aligned.txt")
motif_file = os.path.join(paths.INTERNAL, "motif_efhand.txt")
assert os.path.isfile(matrix_file)
assert os.path.isfile(pdb_seq_file)
# find_motif_positions.from_nbdb(matrix_file, motif_file, pdb_seq_file, 0,
#                                key="custom")
# assert os.path.isfile(pdb_seq_file)
# assert os.path.isfile(aligned_file)
# find_motif_positions.from_aligned(pdb_seq_file, aligned_file, pdb_seq_file,
#                                   motif_file, key='custom')


# import matplotlib.pyplot as plt
# assert os.path.isfile(motif_file)
# with open(motif_file, 'rb') as file:
#     motif_pos_map = pickle.load(file)
# descrs = descr_main.calculate(motif_pos_map)
#
# with open("GxxGxG_pssm.pkl", 'wb') as file:
#     pickle.dump(descrs, file, -1)
with open("GxxGxG_pssm.pkl", 'rb') as file:
    descrs = pickle.load(file)

plots.plot_hbonds_vdw_large_font(descrs)
plt.show()



# plots.plot_hbond_full_sub(descrs, "Efhand", None)
#
# with open("GxGxxG_pssm.pkl", 'rb') as file:
#     descrs = pickle.load(file)
# plots.plot_hbond_full_sub(descrs, "GxGxxG", None)
#
# with open("GxxGxG_pssm.pkl", 'rb') as file:
#     descrs = pickle.load(file)
# plots.plot_hbond_full_sub(descrs, "GxxGxG", None)
#
#
# import sys
# sys.exit()
# bins = defaultdict(int)
# for term in set(descrs['relative_sno'].values):
#     descr = descrs[descrs['relative_sno'] == term]
#     points = []
#     for term in descr['d_a_dist']:
#         for value in term:
#             points.append(value)
#     # bins = defaultdict(int)
#     for value in points:
#         bins[round(value*20)] += 1
#     pairs = defaultdict(list)
#     for p1, p2, dist in zip(descr['planar1'], descr['planar2'],
#                             descr['d_a_dist']):
#         for i, j, k in zip(p1, p2, dist):
#             if k > 3:
#                 pairs[3].append([i, j])
#             else:
#                 pairs[2].append([i, j])
#
# plt.figure()
# # dist_sorted = sorted(pairs.keys())
# # for dist in dist_sorted:
# #     values = np.array(pairs[dist]).T
# plt.scatter(bins.keys(), bins.values(), marker='x')
# plt.show()

# print(descrs.columns)
# import sys
# sys.exit()

# # plots.plot_contacts_heatmap_large_font_key(descrs, 'vdw_full_c1')
# # plots.plot_contacts_heatmap_large_font_key(descrs, 'vdw_full_c2')
# # plots.plot_contacts_heatmap_large_font_key(descrs, 'vdw_full_c3')
# title = "Van der Waals Interactions at \n5 Angstrom Threshold"
# plots.plot_contacts_heatmap_large_font_key(descrs, 'vdw_short_c1', title)
# plt.savefig(os.path.join(paths.ROOT, "output", title))
# title = "Van der Waals Interactions at \n6 Angstrom Threshold"
# plots.plot_contacts_heatmap_large_font_key(descrs, 'vdw_short_c2', title)
# plt.savefig(os.path.join(paths.ROOT, "output", title))
# title = "Van der Waals Interactions at \n8 Angstrom Threshold"
# plots.plot_contacts_heatmap_large_font_key(descrs, 'vdw_short_c3', title)
# plt.savefig(os.path.join(paths.ROOT, "output", title))
# plots.plot_signature_logo(descrs)
# plt.savefig(os.path.join(paths.ROOT, "output", "plot_signature_logo"))
#
# plots.plot_contacts_key(descrs, 'vdw_full_c1_count')
# title = "Full Van der Waals Interactions at 5 Angstrom Threshold"
# plt.title(title)
# plt.savefig(os.path.join(paths.ROOT, "output", "vdw_full_c1_count"))
# plots.plot_contacts_key(descrs, 'vdw_full_c2_count')
# title = "Full Van der Waals Interactions at 6 Angstrom Threshold"
# plt.title(title)
# plt.savefig(os.path.join(paths.ROOT, "output", "vdw_full_c2_count"))
# plots.plot_contacts_key(descrs, 'vdw_full_c3_count')
# title = "Full Van der Waals Interactions at 8 Angstrom Threshold"
# plt.title(title)
# plt.savefig(os.path.join(paths.ROOT, "output", "vdw_full_c3_count"))
#
# plots.plot_contacts_key(descrs, 'vdw_short_c1_count')
# title = "Internal Van der Waals Interactions at 5 Angstrom Threshold"
# plt.title(title)
# plt.savefig(os.path.join(paths.ROOT, "output", "vdw_short_c1_count"))
# plots.plot_contacts_key(descrs, 'vdw_short_c2_count')
# title = "Internal Van der Waals Interactions at 6 Angstrom Threshold"
# plt.title(title)
# plt.savefig(os.path.join(paths.ROOT, "output", "vdw_short_c2_count"))
# plots.plot_contacts_key(descrs, 'vdw_short_c3_count')
# title = "Internal Van der Waals Interactions at 8 Angstrom Threshold"
# plt.title(title)
# plt.savefig(os.path.join(paths.ROOT, "output", "vdw_short_c3_count"))
#
#
# plots.plot_contacts_key_sub(descrs, 'vdw_full_c1_count', 'vdw_short_c1_count')
# title = "External Van der Waals Interactions at 5 Angstrom Threshold"
# plt.title(title)
# plt.savefig(os.path.join(paths.ROOT, "output", "vdw_external_c1_count"))
# plots.plot_contacts_key_sub(descrs, 'vdw_full_c2_count', 'vdw_short_c2_count')
# title = "External Van der Waals Interactions at 6 Angstrom Threshold"
# plt.title(title)
# plt.savefig(os.path.join(paths.ROOT, "output", "vdw_external_c2_count"))
# plots.plot_contacts_key_sub(descrs, 'vdw_full_c3_count', 'vdw_short_c3_count')
# title = "External Van der Waals Interactions at 8 Angstrom Threshold"
# plt.title(title)
# plt.savefig(os.path.join(paths.ROOT, "output", "vdw_external_c3_count"))
# plt.show()
# from matchers.structure.s_matcher import Matcher
# s_matcher = Matcher(cropped=False)
# s_matcher.load(descrs)
# query = descrs[descrs['filename'] == '121p']
# # print(descrs['filename'])
# print(s_matcher.query(query))
# print(len(s_matcher.query(query)[0]))
"""
Things to talk about:
1. if you look at plot, highest close contacts occur at alpha helix, 
corresponding to lowest far contacts
2. Also, the spikes in far contacts correspond to troughs in near contacts, 
probably because the residue is exposed, and vice versa
3. Otherwise all the matchers work as planned. 

# todo:
1. generate plots for all remaining descr
2. write about how secondary structure, etc, show with figures
3. add as feature?

"""