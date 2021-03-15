from collections import defaultdict, Counter
import math
import matplotlib.pyplot as plt
import numpy as np
from operator import itemgetter
import os

from utils import seq_logo, generic


def plot_dihedral_for_diff_res(descr_full, sno_position):
    descr_full = descr_full.groupby('relative_sno')
    descr_for_sno = descr_full.get_group(sno_position)
    descr_for_sno_grouped = descr_for_sno.groupby('res')
    num_plots = len(descr_for_sno_grouped)
    num_col = math.ceil(math.sqrt(num_plots))
    num_row = math.ceil(num_plots / num_col)
    for ax_count, (res, descr_for_res) in enumerate(descr_for_sno_grouped,
                                                    start=1):
        phi = descr_for_res.phi.values
        psi = descr_for_res.psi.values
        ax = plt.subplot(num_row, num_col, ax_count)
        ax.set_title(str(res))
        ax.scatter(phi, psi, marker='.', s=1)
        ax.set_xlabel("phi")
        ax.set_ylabel('psi')
        ax.set_xlim([-180, 180])
        ax.set_ylim([-180, 180])

#pylint: disable=invalid-name
def plot_CA(descr_full):
    descr_full_grouped = descr_full.groupby('filename')
    for filename, df_per_file in descr_full_grouped:
        plt.figure()
        plt.suptitle(filename)
        df_per_file = df_per_file.groupby(['cid', 'seq_marker'])
        num_plots = len(df_per_file)
        num_col = math.ceil(math.sqrt(num_plots))
        num_row = math.ceil(num_plots / num_col)
        for ax_count, ((cid, seq_marker), df_per_set) in \
                enumerate(df_per_file, start=1):
            # subplot indices start from 1
            ax = plt.subplot(num_row, num_col, ax_count, projection='3d')
            ax.set_title(str(seq_marker) + " : " + str(cid))

            snos, CA_coords = df_per_set.sno, df_per_set.CA
            for i in CA_coords:
                assert len(i) == 3
            xs, ys, zs = list(zip(*CA_coords))

            # c for color distribution
            ax.scatter3D(xs, ys, zs, c=range(len(xs)), marker='o',
                         cmap='autumn')
            ax.plot3D(xs, ys, zs, 'k', linewidth=0.5)

            # Mark out 0th position on plot
            for i, sno in enumerate(snos):
                if sno == seq_marker:
                    ax.text(xs[i], ys[i], zs[i], sno)
# pylint: enable=invalid-name

def _for_sort(x):
    return -x[1], x[0]


def plot_logo_from_nbdb(matrix_file, num_res_considered=4, title=None):
    matrix = []
    with open(matrix_file, 'r') as file:
        for line in file:
            if not line.strip():
                continue
            matrix.append([float(i) for i in line.strip().split(" ")])
    to_logo = []
    alphabets = list(generic.AA3_to_AA1.keys())
    for matrix_row in matrix:
        output_row = []
        sorted_i = np.argsort(matrix_row)[::-1]
        for i in range(num_res_considered):
            selected_i = sorted_i[i]
            percent = matrix_row[selected_i]
            res = alphabets[selected_i]
            output_row.append((res, percent))
        to_logo.append(output_row[::-1])
    seq_logo.Logo(to_logo, 0, title=title)

def plot_logo_from_matrix(matrix_file, num_res_considered=4, title=None):
    matrix = []
    with open(matrix_file, 'r') as file:
        for line in file:
            if not line.strip():
                continue
            matrix.append([int(i) for i in line.strip().split(",")])
    to_logo = []
    alphabets = list(generic.AA3_to_AA1.keys())
    for matrix_row in matrix:
        output_row = []
        total = sum(matrix_row)
        sorted_i = np.argsort(matrix_row)[::-1]
        for i in range(num_res_considered):
            selected_i = sorted_i[i]
            percent = matrix_row[selected_i] / total
            res = alphabets[selected_i]
            output_row.append((res, percent))
        to_logo.append(output_row[::-1])
    seq_logo.Logo(to_logo, 0, title=title)

def plot_signature_logo(descr_full, num_res_considered=4, title=None):
    descr_filtered = descr_full.filter(['relative_sno', 'res'])
    min_sno = min(descr_filtered.relative_sno)
    descr_sno_grouped = descr_filtered.groupby('relative_sno')
    to_logo = [[] for _ in range(len(descr_sno_grouped))]
    for i, (sno, descr_per_sno) in enumerate(descr_sno_grouped):
        del sno
        res_unsorted = descr_per_sno.res.values
        res_unique = list(np.unique(res_unsorted, return_counts=True))
        num_seqs = sum(res_unique[1])
        res_sorted = list(zip(*sorted(zip(*res_unique), key=_for_sort)))
        res_names, res_counts = res_sorted[0][:num_res_considered], \
                                res_sorted[1][:num_res_considered]
        res_names = res_names[::-1]
        res_counts = res_counts[::-1]
        res_percents = list([i/num_seqs for i in res_counts])
        for name, percent in zip(res_names, res_percents):
            to_logo[i].append((name, percent))
    seq_logo.Logo(to_logo, min_sno, title=title)

# pylint: disable=invalid-name, dangerous-default-value
def plot_signature_bar(descr, num_res_considered=2,
                       AA3_to_AA1=generic.AA3_to_AA1):
    """
    :param num_res_considered: Number of top-count res to show per relative
    sno position.
    """
    descr = descr.filter(['relative_sno', 'res'])
    plt.figure()
    plt.suptitle("Signatures")
    ax = plt.subplot(111)

    # Generate cmap for plots
    unique_res = np.unique(descr.res)
    num_plots = len(unique_res)
    colormap = dict()
    cmap = plt.get_cmap('gist_ncar')
    for i, res in enumerate(unique_res):
        colormap[res] = cmap(i/num_plots)

    descr = descr.groupby('relative_sno')
    bottom_val = dict()
    for sno, descr_per_sno in descr:
        res_unsorted = descr_per_sno.res.values
        res_unique = list(np.unique(res_unsorted, return_counts=True))
        res_sorted = list(zip(*sorted(zip(*res_unique), key=_for_sort)))
        res_names, res_counts = res_sorted[0][:num_res_considered], \
                                res_sorted[1][:num_res_considered]
        for name, count in zip(res_names, res_counts):
            if sno in bottom_val:
                ax.bar(sno, count, bottom=bottom_val[sno], label=name,
                       color=colormap[name])
                bottom_val[sno] += count
            else:
                ax.bar(sno, count, label=name, color=colormap[name])
                bottom_val[sno] = count

    handles, labels = ax.get_legend_handles_labels()

    # Convert to single-letter code (AA3=>AA1)
    single_letter_labels = []
    for label in labels:
        single_letter_labels.append(AA3_to_AA1[label] + " / " + label)
    labels = single_letter_labels

    # Sort labels
    labels, handles = zip(*sorted(zip(labels, handles), key=itemgetter(0)))
    labels = list(labels)
    handles = list(handles)

    # Remove duplicates
    history = []
    assert len(labels) == len(handles)
    for i in range(len(labels))[::-1]:
        if labels[i] not in history:
            history.append(labels[i])
        else:
            del labels[i]
            del handles[i]
    ax.legend(handles, labels)
# pylint: enable=invalid-name, dangerous-default-value

def plot_dihedral(descr_full, remove_labels=True, add_filename=False):
    descr_grouped = descr_full.groupby('relative_sno')

    num_plots = len(descr_grouped) - 2
    num_col = math.ceil(math.sqrt(num_plots))
    num_row = math.ceil(num_plots / num_col)

    plt.figure()
    plt.suptitle("Dihedrals across different relative sno position")
    for ax_count, (relative_sno, df_per_sno) in enumerate(descr_grouped):
        if ax_count in (0, len(descr_grouped) - 1):
            # Screen off first, last position, invalid dihedral values
            continue
        phis = df_per_sno.phi.values
        psis = df_per_sno.psi.values
        filenames = df_per_sno.filename.values
        seq_markers = df_per_sno.seq_marker.values

        ax = plt.subplot(num_row, num_col, ax_count)
        ax.set_title(str(relative_sno))
        ax.scatter(phis, psis, marker='x', s=1, color='k')
        ax.set_xlabel("phi")
        ax.set_ylabel('psi')
        ax.set_xlim([-180, 180])
        ax.set_ylim([-180, 180])

        if remove_labels:
            ax.set_xlabel("")
            ax.set_ylabel('')
            ax.set_xticklabels([])
            ax.set_yticklabels([])
            ax.set_xticks([])
            ax.set_yticks([])

        if add_filename:
            for phi, psi, filename, seq_marker \
                in zip(phis, psis, filenames, seq_markers):
                filename_str = str(filename) + " " + str(seq_marker)
                ax.text(phi, psi, filename_str)


def plot_dihedral_allpos(descr_full, remove_labels=True, add_filename=False):
    descr_grouped = descr_full.groupby('relative_sno')

    num_plots = len(descr_grouped) - 2
    num_col = math.ceil(math.sqrt(num_plots))
    num_row = math.ceil(num_plots / num_col)

    plt.figure()
    plt.suptitle("Dihedrals across different relative sno position")
    for ax_count, (relative_sno, df_per_sno) in enumerate(descr_grouped):
        if ax_count in (0, len(descr_grouped) - 1):
            # Screen off first, last position, invalid dihedral values
            continue
        phis = df_per_sno.phi.values
        psis = df_per_sno.psi.values
        filenames = df_per_sno.filename.values
        seq_markers = df_per_sno.seq_marker.values
        selected_1yru = None
        selected_1a29 = None
        for i, filename in enumerate(filenames):
            if filename == "1yru":
                selected_1yru = i
            if filename == '1a29':
                selected_1a29 = i
        assert selected_1yru is not None
        assert selected_1a29 is not None


        ax = plt.subplot(num_row, num_col, ax_count)
        ax.set_title(str(relative_sno))
        ax.scatter(phis, psis, marker='x', s=1)
        ax.set_xlabel("phi")
        ax.set_ylabel('psi')
        ax.set_xlim([-180, 180])
        ax.set_ylim([-180, 180])

        ax.scatter([phis[selected_1yru]], [psis[[selected_1yru]]], marker='x',
                   s=20, color='k')
        ax.scatter([phis[selected_1a29]], [psis[[selected_1a29]]], marker='x',
                   s=20, color='r')

        if remove_labels:
            ax.set_xlabel("")
            ax.set_ylabel('')
            ax.set_xticklabels([])
            ax.set_yticklabels([])
            ax.set_xticks([])
            ax.set_yticks([])

        if add_filename:
            for phi, psi, filename, seq_marker in zip(phis, psis, filenames,
                                                      seq_markers):
                filename_str = str(filename) + " " + str(seq_marker)
                ax.text(phi, psi, filename_str)


def plot_dihedral_4pos(descr_full, remove_labels=True, add_filename=False):
    descr_grouped = descr_full.groupby('relative_sno')

    num_plots = len(descr_grouped) - 2
    num_col = math.ceil(math.sqrt(num_plots))
    num_row = math.ceil(num_plots / num_col)

    plt.figure()
    plt.suptitle("Dihedrals across different relative sno position")
    count = 0
    for ax_count, (relative_sno, df_per_sno) in enumerate(descr_grouped):
        if ax_count in (0, len(descr_grouped) - 1):
            # Screen off first, last position, invalid dihedral values
            continue
        if ax_count not in (7, 8, 9, 10):
            continue

        phis = df_per_sno.phi.values
        psis = df_per_sno.psi.values
        filenames = df_per_sno.filename.values
        seq_markers = df_per_sno.seq_marker.values
        selected_1yru = None
        selected_1a29 = None
        for i, filename in enumerate(filenames):
            if filename == "1yru":
                selected_1yru = i
            if filename == '1a29':
                selected_1a29 = i
        assert selected_1yru is not None
        assert selected_1a29 is not None
        count += 1
        ax = plt.subplot(2, 2, count)
        ax.set_title(str(relative_sno))
        ax.scatter(phis, psis, marker='x', s=1)
        ax.set_xlabel("phi")
        ax.set_ylabel('psi')
        ax.set_xlim([-180, 180])
        ax.set_ylim([-180, 180])

        ax.scatter([phis[selected_1yru]], [psis[[selected_1yru]]], marker='x',
                   s=20, color='k')
        ax.scatter([phis[selected_1a29]], [psis[[selected_1a29]]], marker='x',
                   s=20, color='r')

        if remove_labels:
            ax.set_xlabel("")
            ax.set_ylabel('')
            ax.set_xticklabels([])
            ax.set_yticklabels([])
            ax.set_xticks([])
            ax.set_yticks([])

        if add_filename:
            for phi, psi, filename, seq_marker in zip(phis, psis, filenames,
                                                      seq_markers):
                filename_str = str(filename) + " " + str(seq_marker)
                ax.text(phi, psi, filename_str)


def plot_hbonds_vdw_merged(descr):
    descr2 = descr.groupby('relative_sno')
    to_plots = dict()
    for relative_sno, df_per_sno in descr2:
        donor = df_per_sno.donor.values
        counts = [len(i) for i in donor]
        to_plots[relative_sno] = (np.mean(counts), np.std(counts))
    relative_snos = list(to_plots.keys())
    means = list(zip(*to_plots.values()))[0]
    stderrs = list(zip(*to_plots.values()))[1]
    fig, ax1 = plt.subplots()
    ax1.scatter(relative_snos, means, color='r', marker='x')
    ax1.set_ylabel("Van der Waals interactions", fontsize=40)
    ax1.spines['left'].set_color('red')
    # plt.suptitle("Num Hbonds")
    # plt.scatter(relative_snos, means)
    ax1.set_xlabel("Relative Residue Position", fontsize=40)
    ax1.set_title("VdW interactions and Hydrogen Bonds", fontsize = 40)

    ax2 = ax1.twinx()
    descr = descr.groupby(['filename', 'cid', 'seq_marker'])
    contacts_by_pos = Counter()
    for _, descr_per_set in descr:
        contacts = descr_per_set.contact.values
        relative_sno = descr_per_set.relative_sno.values
        assert len(contacts) == len(relative_sno)
        if not np.isnan(contacts[0]):
            for contact, sno in zip(contacts, relative_sno):
                contacts_by_pos[sno] += contact
    for sno in contacts_by_pos.keys():
        contacts_by_pos[sno] /= len(descr)
    # plt.suptitle("Contacts summed counts")
    ax2.scatter(contacts_by_pos.keys(), contacts_by_pos.values(), color='k',
                marker='x')
    ax2.set_ylabel("Hydrogen bonds", fontsize=40)
    ax2.spines['right'].set_color('red')
    ax2.tick_params(axis='y', colors='red', labelsize=40)
    ax2.yaxis.label.set_color('red')
    ax2.set_xlabel("Relative Residue Position", fontsize=40)


def plot_hbonds_vdw_large_font(descr):
    descr2 = descr.groupby('relative_sno')
    to_plots = dict()
    for relative_sno, df_per_sno in descr2:
        donor = df_per_sno.donor.values
        counts = [len(i) for i in donor]
        to_plots[relative_sno] = (np.mean(counts), np.std(counts))
    relative_snos = list(to_plots.keys())
    means = list(zip(*to_plots.values()))[0]
    fig, ax1 = plt.subplots(figsize=(10, 8))
    ax1.scatter(relative_snos, means, color='r', marker='x', s=70)
    ax1.set_ylabel("Van der Waals interactions", fontsize=33)
    ax1.spines['left'].set_color('red')
    # plt.suptitle("Num Hbonds")
    # plt.scatter(relative_snos, means)
    ax1.set_xlabel("Residue Position", fontsize=35)
    # ax1.set_title("VdW interactions and Hydrogen Bonds", fontsize=35)
    ax1.tick_params(axis='x', labelsize=35)
    ax1.tick_params(axis='y', labelsize=35)

    ax2 = ax1.twinx()
    descr = descr.groupby(['filename', 'cid', 'seq_marker'])
    contacts_by_pos = Counter()
    for _, descr_per_set in descr:
        contacts = descr_per_set.contact.values
        relative_sno = descr_per_set.relative_sno.values
        assert len(contacts) == len(relative_sno)
        if not np.isnan(contacts[0]):
            for contact, sno in zip(contacts, relative_sno):
                contacts_by_pos[sno] += contact
    for sno in contacts_by_pos.keys():
        contacts_by_pos[sno] /= len(descr)
    # plt.suptitle("Contacts summed counts")
    ax2.scatter(contacts_by_pos.keys(), contacts_by_pos.values(), color='k',
                marker='x', s=70)
    ax2.set_ylabel("Hydrogen bonds", fontsize=35)
    ax2.spines['right'].set_color('red')
    ax2.tick_params(axis='y', colors='red', labelsize=35)
    ax2.yaxis.label.set_color('red')
    ax2.tick_params(axis='x', labelsize=35)
    ax2.set_xlabel("Residue Position", fontsize=35)

def plot_hbonds_raw(descr):
    """
    Convert to heatmap eventually, maybe.
    """
    descr = descr.groupby('relative_sno')
    to_plots = dict()
    print(len(descr))
    for relative_sno, df_per_sno in descr:
        donor = df_per_sno.donor.values
        # print(len(donor))
        # print(len(donor[0]))
        counts = [len(i) for i in donor]
        to_plots[relative_sno] = (np.mean(counts), np.std(counts))
    relative_snos = list(to_plots.keys())
    means = list(zip(*to_plots.values()))[0]
    stderrs = list(zip(*to_plots.values()))[1]
    plt.figure()
    plt.suptitle("Num Hbonds")
    plt.errorbar(relative_snos, means, yerr=stderrs)


def plot_hbonds_heatmap_large_font(descr):
    descr = descr.groupby(['filename', 'cid', 'seq_marker'])
    heatmap_donor = np.zeros((30, 30))
    heatmap_acceptor = np.zeros((30, 30))
    # h_contacts_by_pos = Counter()
    for _, descr_per_set in descr:
        min_sno = min(descr_per_set.sno.values)
        for i, row in enumerate(descr_per_set.d_sno.values):
            for term in row:
                if term >= min_sno + 30 or term < min_sno:
                    continue
                heatmap_donor[i][term - min_sno] += 1
        for i, row in enumerate(descr_per_set.a_sno.values):
            for term in row:
                if term >= min_sno + 30 or term < min_sno:
                    continue
                heatmap_donor[i][term - min_sno] += 1
    for i in range(30):
        for j in range(30):
            heatmap_donor[i][j] /= len(descr)
            heatmap_acceptor[i][j] /= len(descr)
    fig, ax = plt.subplots(figsize=(10, 8))
    # im = ax.imshow(heatmap_donor, cmap='YlOrRd')
    im = ax.imshow(heatmap_donor, cmap='tab20')
    ax.figure.colorbar(im, format="%.1f")
    plt.title("Hydrogen Bonds")
    plt.xlabel("Donor Residue Position")
    plt.ylabel("Acceptor Residue Position")
    for i in range(30):
        plt.axhline(i + 0.5, linestyle='solid', color='lightgrey')
        plt.axvline(i + 0.5, linestyle='solid', color='lightgrey')
    plt.plot([0, 29], [0, 29], 'k--')
    plt.gca().invert_yaxis()

def plot_hbonds_heatmap(descr):
    descr = descr.groupby(['filename', 'cid', 'seq_marker'])
    heatmap_donor = np.zeros((30, 30))
    heatmap_acceptor = np.zeros((30, 30))
    # h_contacts_by_pos = Counter()
    for _, descr_per_set in descr:
        min_sno = min(descr_per_set.sno.values)
        for i, row in enumerate(descr_per_set.d_sno.values):
            for term in row:
                if term >= min_sno + 30 or term < min_sno:
                    continue
                heatmap_donor[i][term - min_sno] += 1
        for i, row in enumerate(descr_per_set.a_sno.values):
            for term in row:
                if term >= min_sno + 30 or term < min_sno:
                    continue
                heatmap_donor[i][term - min_sno] += 1
    for i in range(30):
        for j in range(30):
            heatmap_donor[i][j] /= len(descr)
            heatmap_acceptor[i][j] /= len(descr)
    fig, ax = plt.subplots(figsize=(10, 8))
    # im = ax.imshow(heatmap_donor, cmap='YlOrRd')
    im = ax.imshow(heatmap_donor, cmap='tab20')
    a = ax.figure.colorbar(im, format="%.1f")
    a.ax.tick_params(labelsize=40)
    plt.tick_params(axis='x', labelsize=40)
    plt.tick_params(axis='y', labelsize=40)
    plt.title("Hydrogen Bonds", fontsize=40)
    plt.xlabel("Residue Position", fontsize=40)
    plt.ylabel("Residue Position", fontsize=40)
    for i in range(30):
        plt.axhline(i + 0.5, linestyle='solid', color='lightgrey')
        plt.axvline(i +0.5, linestyle='solid', color='lightgrey')
    plt.plot([0, 29], [0, 29], 'k--')
    plt.gca().invert_yaxis()
    # fig, ax = plt.subplots()
    # im = ax.imshow(heatmap_acceptor, cmap='YlOrRd')
    # ax.figure.colorbar(im, format="%.1f")
    # plt.title("Hbond Acceptors")
    # plt.xlabel("Source Residue Position")
    # plt.ylabel("Destination Residue Position")
    # for i in range(30):
    #     plt.axhline(i + 0.5, linestyle='solid', color='lightgrey')
    #     plt.axvline(i + 0.5, linestyle='solid', color='lightgrey')
    # plt.plot([0, 29], [0, 29], 'k--')
    # plt.gca().invert_yaxis()

    # fig, ax = plt.subplots()
    # ax.imshow(heatmap_acceptor)
    # plt.title("acceptor")
    # plt.gca().invert_yaxis()
    # plt.show()

    # plt.figure()
    # plt.suptitle("H-Contact summed counts")
    # plt.scatter(h_contacts_by_pos.keys(), h_contacts_by_pos.values())

def plot_hbonds_percentage(descr):
    descr = descr.groupby(['filename', 'cid', 'seq_marker'])
    to_plots = defaultdict(list)
    for _, descr_per_set in descr:
        relative_sno = descr_per_set.relative_sno.values
        counts = [len(i) for i in descr_per_set.donor.values]
        total = sum(counts)
        percent = [count/total*100 for count in counts]
        for i, sno in enumerate(relative_sno):
            to_plots[sno].append(percent[i])
    relative_snos = list(to_plots.keys())
    percents = list(to_plots.values())
    means = list([np.mean(i) for i in percents])

    plt.figure()
    plt.suptitle("Num Hbonds as Percent across Fragment")
    plt.scatter(relative_snos, means)


def plot_contacts_heatmap_large_font_key(descr, key, title=None):
    descr = descr.groupby(['filename', 'cid', 'seq_marker'])
    h_contacts_by_pos = np.zeros((30, 30), dtype=float)
    for _, descr_per_set in descr:
        relative_sno = descr_per_set.relative_sno.values
        h_contacts = descr_per_set[key].values
        assert len(h_contacts) == len(relative_sno)
        # if not np.isnan(h_contacts[0]):
        for h_contact, sno in zip(h_contacts, relative_sno):
            if len(h_contact) != 30:
                continue
            h_contacts_by_pos[sno] += h_contact

    fig, ax = plt.subplots(figsize=(10, 8))
    # im = ax.imshow(h_contacts_by_pos, cmap='YlOrRd')
    im = ax.imshow(h_contacts_by_pos, cmap='tab20')

    a = ax.figure.colorbar(im, format="%.1f")
    a.ax.tick_params(labelsize=29)
    ax.tick_params(labelsize=29)
    if title is None:
        plt.title("Van der Waals Bonds", fontsize=29)
    else:
        plt.title(title, fontsize=29)
    plt.xlabel("Residue Position", fontsize=29)
    plt.ylabel("Residue Position", fontsize=29)
    for i in range(30):
        plt.axhline(i + 0.5, linestyle='solid', color='lightgrey')
        plt.axvline(i + 0.5, linestyle='solid', color='lightgrey')
    plt.plot([0, 29], [0, 29], 'k--')
    plt.gca().invert_yaxis()

def plot_contacts_heatmap_large_font_new(descr):
    descr = descr.groupby(['filename', 'cid', 'seq_marker'])
    h_contacts_by_pos = np.zeros((30, 30), dtype=float)
    for _, descr_per_set in descr:
        relative_sno = descr_per_set.relative_sno.values
        h_contacts = descr_per_set.vdw_inter_mask.values
        assert len(h_contacts) == len(relative_sno)
        # if not np.isnan(h_contacts[0]):
        for h_contact, sno in zip(h_contacts, relative_sno):
            if len(h_contact) != 30:
                continue
            h_contacts_by_pos[sno] += h_contact

    fig, ax = plt.subplots(figsize=(10, 8))
    # im = ax.imshow(h_contacts_by_pos, cmap='YlOrRd')
    im = ax.imshow(h_contacts_by_pos, cmap='tab20')

    a = ax.figure.colorbar(im, format="%.1f")
    a.ax.tick_params(labelsize=29)
    ax.tick_params(labelsize=29)
    plt.title("Van der Waals Bonds", fontsize=29)
    plt.xlabel("Residue Position", fontsize=29)
    plt.ylabel("Residue Position", fontsize=29)
    for i in range(30):
        plt.axhline(i + 0.5, linestyle='solid', color='lightgrey')
        plt.axvline(i + 0.5, linestyle='solid', color='lightgrey')
    plt.plot([0, 29], [0, 29], 'k--')
    plt.gca().invert_yaxis()

def plot_contacts_heatmap_large_font(descr):
    descr = descr.groupby(['filename', 'cid', 'seq_marker'])
    h_contacts_by_pos = np.zeros((30, 30), dtype=float)
    for _, descr_per_set in descr:
        relative_sno = descr_per_set.relative_sno.values
        h_contacts = descr_per_set.h_contacts.values
        assert len(h_contacts) == len(relative_sno)
        # if not np.isnan(h_contacts[0]):
        for h_contact, sno in zip(h_contacts, relative_sno):
            if len(h_contact) != 30:
                continue
            h_contacts_by_pos[sno] += h_contact

    fig, ax = plt.subplots(figsize=(12, 14))
    # im = ax.imshow(h_contacts_by_pos, cmap='YlOrRd')
    im = ax.imshow(h_contacts_by_pos, cmap='tab20')

    a = ax.figure.colorbar(im, format="%.1f")
    a.ax.tick_params(labelsize=40)
    ax.tick_params(labelsize=40)
    plt.title("Van der Waals Interactions", fontsize=40)
    plt.xlabel("Residue Position", fontsize=40)
    plt.ylabel("Residue Position", fontsize=40)
    for i in range(30):
        plt.axhline(i + 0.5, linestyle='solid', color='lightgrey')
        plt.axvline(i + 0.5, linestyle='solid', color='lightgrey')
    plt.plot([0, 29], [0, 29], 'k--')
    plt.gca().invert_yaxis()

def plot_contacts_heatmap(descr):
    descr = descr.groupby(['filename', 'cid', 'seq_marker'])
    h_contacts_by_pos = np.zeros((30, 30), dtype=float)
    for _, descr_per_set in descr:
        relative_sno = descr_per_set.relative_sno.values
        h_contacts = descr_per_set.h_contacts.values
        assert len(h_contacts) == len(relative_sno)
        # if not np.isnan(h_contacts[0]):
        for h_contact, sno in zip(h_contacts, relative_sno):
            if len(h_contact) != 30:
                continue
            h_contacts_by_pos[sno] += h_contact

    fig, ax = plt.subplots()
    # im = ax.imshow(h_contacts_by_pos, cmap='YlOrRd')
    im = ax.imshow(h_contacts_by_pos, cmap='tab20')

    ax.figure.colorbar(im, format="%.1f")
    plt.title("Van der Waals Bonds")
    plt.xlabel("Residue Position")
    plt.ylabel("Residue Position")
    for i in range(30):
        plt.axhline(i + 0.5, linestyle='solid', color='lightgrey')
        plt.axvline(i + 0.5, linestyle='solid', color='lightgrey')
    plt.plot([0, 29], [0, 29], 'k--')
    plt.gca().invert_yaxis()


    # fig, ax = plt.subplots()
    # ax.imshow(h_contacts_by_pos)
    # plt.gca().invert_yaxis()
    # plt.title("Van der Waals Bonds")
    # plt.xlabel("Source Residue Position")
    # plt.ylabel("Destination Residue Position")

    # for key, value in h_contacts_by_pos.items():  #     print(key)  #
    # print(value)  # # print(h_contacts_by_pos)  # plt.figure()  #
    # plt.suptitle("H-Contact summed counts")  # plt.scatter(
    # h_contacts_by_pos.keys(), h_contacts_by_pos.values())


def plot_contacts(descr):
    descr = descr.groupby(['filename', 'cid', 'seq_marker'])
    contacts_by_pos = Counter()
    for _, descr_per_set in descr:
        contacts = descr_per_set.contact.values
        relative_sno = descr_per_set.relative_sno.values
        assert len(contacts) == len(relative_sno)
        if not np.isnan(contacts[0]):
            for contact, sno in zip(contacts, relative_sno):
                contacts_by_pos[sno] += contact
    plt.figure()
    plt.suptitle("Contacts summed counts")
    plt.scatter(contacts_by_pos.keys(), contacts_by_pos.values())

def plot_hbond_full_sub(descr, key_term, sub_key):
    descr = descr.groupby(['filename', 'cid', 'seq_marker'])
    hbond_full = Counter()
    hbond_inner = Counter()
    hbond_outer = Counter()
    for _, descr_per_set in descr:
        min_sno = min(descr_per_set.sno.values)
        for i, row in enumerate(descr_per_set.d_sno.values):
            inner_count = 0
            outer_count = 0
            total_count = 0
            for term in row:
                if term >= min_sno + 30 or term < min_sno:
                    outer_count += 1
                    hbond_outer[i] += 1
                else:
                    inner_count += 1
                    hbond_inner[i] += 1
                total_count += 1
                hbond_full[i] += 1
            hbond_outer[i] += outer_count // 5
            hbond_inner[i] += outer_count // 5
            hbond_full[i] += outer_count // 5

    for term in (hbond_full, hbond_outer, hbond_inner):
        for key in term.keys():
            term[key] = term[key] // 5


    f = plt.figure()
    plt.title(f"{key_term} External Donor Hydrogen Bonds")
    plt.scatter(hbond_outer.keys(), hbond_outer.values())
    plt.ylim(bottom=0)
    plt.savefig(f"{key_term}4.png")
    f = plt.figure()
    plt.title(f"{key_term} Internal Donor Hydrogen Bonds")
    plt.scatter(hbond_inner.keys(), hbond_inner.values())
    plt.ylim(bottom=0)
    plt.savefig(f"{key_term}5.png")
    f = plt.figure()
    plt.title(f"{key_term} Total Donor Hydrogen Bonds")
    plt.scatter(hbond_full.keys(), hbond_full.values())
    plt.ylim(bottom=0)
    plt.savefig(f"{key_term}6.png")

    hbond_full = Counter()
    hbond_inner = Counter()
    hbond_outer = Counter()
    for _, descr_per_set in descr:
        min_sno = min(descr_per_set.sno.values)
        for i, row in enumerate(descr_per_set.a_sno.values):
            for term in row:
                if term >= min_sno + 30 or term < min_sno:
                    hbond_outer[i] += 1
                else:
                    hbond_inner[i] += 1
                hbond_full[i] += 1
    for term in (hbond_full, hbond_outer, hbond_inner):
        for key in term.keys():
            term[key] = term[key] // 10

    # f = plt.figure()
    # plt.title("External Acceptor Hydrogen Bonds at 3.5 Angstrom Threshold")
    # plt.scatter(hbond_outer.keys(), hbond_outer.values())
    # f = plt.figure()
    # plt.title("Internal Acceptor Hydrogen Bonds at 3.5 Angstrom Threshold")
    # plt.scatter(hbond_inner.keys(), hbond_inner.values())
    # f = plt.figure()
    # plt.title("Total Acceptor Hydrogen Bonds at 3.5 Angstrom Threshold")
    # plt.scatter(hbond_full.keys(), hbond_full.values())
    f = plt.figure()
    plt.title(f"{key_term} External Acceptor Hydrogen Bonds")
    plt.scatter(hbond_outer.keys(), hbond_outer.values())
    plt.ylim(bottom=0)
    plt.savefig(f"{key_term}1.png")
    f = plt.figure()
    plt.title(f"{key_term} Internal Acceptor Hydrogen Bonds")
    plt.scatter(hbond_inner.keys(), hbond_inner.values())
    plt.ylim(bottom=0)
    plt.savefig(f"{key_term}2.png")
    f = plt.figure()
    plt.title(f"{key_term} Total Acceptor Hydrogen Bonds")
    plt.scatter(hbond_full.keys(), hbond_full.values())
    plt.ylim(bottom=0)
    plt.savefig(f"{key_term}3.png")


def plot_contacts_key_sub(descr, key, sub_key):
    descr = descr.groupby(['filename', 'cid', 'seq_marker'])
    contacts_by_pos = Counter()
    for _, descr_per_set in descr:
        contacts = descr_per_set[key].values - descr_per_set[sub_key].values
        relative_sno = descr_per_set.relative_sno.values
        assert len(contacts) == len(relative_sno)
        if not np.isnan(contacts[0]):
            for contact, sno in zip(contacts, relative_sno):
                contacts_by_pos[sno] += contact
    plt.figure()
    # plt.suptitle("Contacts summed counts")
    plt.scatter(contacts_by_pos.keys(), contacts_by_pos.values())

def plot_contacts_key(descr, key):
    descr = descr.groupby(['filename', 'cid', 'seq_marker'])
    contacts_by_pos = Counter()
    for _, descr_per_set in descr:
        contacts = descr_per_set[key].values
        relative_sno = descr_per_set.relative_sno.values
        assert len(contacts) == len(relative_sno)
        if not np.isnan(contacts[0]):
            for contact, sno in zip(contacts, relative_sno):
                contacts_by_pos[sno] += contact
    plt.figure()
    # plt.suptitle("Contacts summed counts")
    plt.scatter(contacts_by_pos.keys(), contacts_by_pos.values())

def plot_covalents(descr):
    descr = descr.groupby(['filename', 'cid', 'seq_marker'])
    covalents_by_pos = Counter()
    for _, descr_per_set in descr:
        covalents = descr_per_set.covalent.values
        relative_sno = descr_per_set.relative_sno.values
        assert len(covalents) == len(relative_sno)
        if not np.isnan(covalents[0]):
            for covalent, sno in zip(covalents, relative_sno):
                covalents_by_pos[sno] += covalent
    plt.figure()
    plt.suptitle("Covalent summed counts")
    plt.scatter(covalents_by_pos.keys(), covalents_by_pos.values())


def plot_h_contacts(descr):
    descr = descr.groupby(['filename', 'cid', 'seq_marker'])
    h_contacts_by_pos = np.zeros((30, 30), dtype=float)
    count = 0
    for _, descr_per_set in descr:
        relative_sno = descr_per_set.relative_sno.values
        h_contacts = descr_per_set.h_contacts.values
        assert len(h_contacts) == len(relative_sno)
        # if not np.isnan(h_contacts[0]):
        for h_contact, sno in zip(h_contacts, relative_sno):
            if len(h_contact) != 30:
                continue
            h_contacts_by_pos[sno] += h_contact
        # if not np.isnan(h_contacts[0]):
        #     for h_contact, sno in zip(h_contacts, relative_sno):
        #         h_contacts_by_pos[sno] += h_contact

    # print()
    # merged = np.array(list(h_contacts_by_pos.values()))
    # merged = np.array(list(h_contacts_by_pos.values()), dtype=float)
    h_contacts_by_pos /= np.max(h_contacts_by_pos)
    fig, ax = plt.subplots()
    ax.imshow(h_contacts_by_pos)
    plt.gca().invert_yaxis()
    plt.title("Van der Waals Bonds")
    plt.xlabel("Source Residue Position")
    plt.ylabel("Destination Residue Position")

    # for key, value in h_contacts_by_pos.items():
    #     print(key)
    #     print(value)
    # # print(h_contacts_by_pos)
    # plt.figure()
    # plt.suptitle("H-Contact summed counts")
    # plt.scatter(h_contacts_by_pos.keys(), h_contacts_by_pos.values())


def plot_all(df, directory, key, save=True):
    """
    Input df should have these keys:
    ['sno', 'contact', 'covalent', 'phi', 'psi', 'region', 'ss', 'ext', 'role',
     'category', 'donor', 'acc', 'res', 'CA', 'filename', 'seq_marker', 'cid']
    """
    seqlogo_path = os.path.join(directory, f"{key}_seqlogo_prob.png")
    plot_signature_logo(df)
    if save:
        plt.savefig(seqlogo_path)

    dihedral_path = os.path.join(directory, f"{key}_dihedral.png")
    plot_dihedral(df)
    if save:
        plt.savefig(dihedral_path)

    contacts_path = os.path.join(directory, f"{key}_contacts.png")
    plot_contacts(df)
    if save:
        plt.savefig(contacts_path)

    covalents_path = os.path.join(directory, f"{key}_covalents.png")
    plot_covalents(df)
    if save:
        plt.savefig(covalents_path)
    plt.show()

    # hbond_path = os.path.join(directory, f"{key}_hbonds.png")
    # plot_hbonds_raw(df)
    # plt.savefig(hbond_path)
