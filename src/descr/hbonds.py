from collections import defaultdict
import decimal
import numpy as np
import pandas as pd
import warnings

from descr import geometry

warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)

def get_descr_hb(df_hb, df_ATOM, df_HETATM, dsr_snos):
    if df_hb is None:
        hbond_descr = _build_empty_hb_descr(dsr_snos)
        return hbond_descr

    assert not df_hb.empty
    df_ATOM = df_ATOM.filter(items=["cid", "sno", 'aname', 'coord'])
    if df_HETATM is None:
        ATOM_HETATM = df_ATOM
    else:
        df_HETATM = df_HETATM.filter(items=["cid", "sno", 'aname', 'coord'])
        ATOM_HETATM = df_ATOM.append(df_HETATM, ignore_index=True)
    ATOM_HETATM = ATOM_HETATM.set_index(['cid', 'sno', 'aname'])
    df_ATOM = df_ATOM.set_index(['cid', 'sno'])
    df_ATOM = df_ATOM[df_ATOM.aname.isin(("N", "C", "CA"))]

    hb_descr_raw = defaultdict(list)

    for _, hb_line in df_hb.iterrows():
        if hb_line.d_sno not in dsr_snos and hb_line.a_sno not in dsr_snos:
            continue
        d_info = (hb_line.d_cid, hb_line.d_sno, hb_line.d_aname)
        a_info = (hb_line.a_cid, hb_line.a_sno, hb_line.a_aname)

        other_hb_info = (
            hb_line.pdb_id, hb_line.d_cid, hb_line.d_res, hb_line.d_sno,
            hb_line.d_aname, hb_line.a_cid, hb_line.a_res, hb_line.a_sno,
            hb_line.a_aname, hb_line.d_a_dist, hb_line.a_d_dd, hb_line.d_a_aa,
            hb_line.planar1, hb_line.planar2, hb_line.atom_category)

        # Determine vector will point from doner => acceptor
        hbond_vector = _get_hbond_vector(d_info, a_info, ATOM_HETATM)

        if d_info[1] in dsr_snos:
            role = "D"
            sno, addition = _set_hb_descr(df_ATOM, hbond_vector, dsr_snos,
                                          d_info, a_info, role, other_hb_info)
            hb_descr_raw[sno].append(addition)

        if a_info[1] in dsr_snos:
            role = "A"
            sno, addition = _set_hb_descr(df_ATOM, hbond_vector, dsr_snos,
                                          a_info, d_info, role, other_hb_info)
            hb_descr_raw[sno].append(addition)

    hb_descr = _screen_duplicate(hb_descr_raw, dsr_snos)
    return hb_descr

def _build_empty_hb_descr(snos):
    hb_descr = defaultdict(list)
    for sno in snos:
        hb_descr = _fill_with_empty(hb_descr, sno)
    hb_descr = dict(hb_descr)
    return hb_descr

def _fill_with_empty(hb_descr, sno):
    hb_attributes = ['role', 'donor', 'acc', 'ext', 'pdb_id', 'd_cid', 'd_res',
                     'd_sno', 'd_aname', 'a_cid', 'a_res', 'a_sno', 'a_aname',
                     'd_a_dist', 'a_d_dd', 'd_a_aa', 'planar1', 'planar2',
                     'category']
    hb_descr['sno'].append(sno)
    for attr in hb_attributes:
        hb_descr[attr].append([])
    return hb_descr

def _screen_duplicate(hb_descr_raw, dsr_snos):
    """
    np.allclose or any float array-related comparison is really slow, hence
    the string representation. Hackish, but >100x faster.
    """
    role_swap = dict(A="D", D="A")
    hb_descr = defaultdict(list)

    for sno in dsr_snos:
        if sno not in hb_descr_raw:
            hb_descr = _fill_with_empty(hb_descr, sno)
            continue

        values = hb_descr_raw[sno]
        history = []
        roles, donors, accs, others = list(map(list, zip(*values)))
        assert len(roles) == len(donors) == len(accs)
        for i in range(len(roles))[::-1]:
            role, donor, acc = roles[i], donors[i], accs[i]
            role_reversed = role_swap[role]
            string = str(np.array((donor, acc)))
            if any([(string == item[1]) or
                    (string == item[2] and role_reversed == item[0])
                    for item in history]):
                del others[i]
                del roles[i]
                del donors[i]
                del accs[i]
            else:
                history.append([role, str(np.array((donor, acc))),
                                str(np.array((acc, donor)))])

        ext, pdb_id, d_cid, d_res, d_sno, d_aname, a_cid, a_res, a_sno, \
        a_aname, d_a_dist, a_d_dd, d_a_aa, planar1, planar2, atom_category \
            = list(zip(*others))

        hb_descr['sno'].append(sno)
        hb_descr['role'].append(roles)
        hb_descr['donor'].append(donors)
        hb_descr['acc'].append(accs)
        hb_descr['ext'].append(ext)
        hb_descr['pdb_id'].append(pdb_id)
        hb_descr['d_cid'].append(d_cid)
        hb_descr['d_res'].append(d_res)
        hb_descr['d_sno'].append(d_sno)
        hb_descr['d_aname'].append(d_aname)
        hb_descr['a_cid'].append(a_cid)
        hb_descr['a_res'].append(a_res)
        hb_descr['a_sno'].append(a_sno)
        hb_descr['a_aname'].append(a_aname)
        hb_descr['d_a_dist'].append(d_a_dist)
        hb_descr['a_d_dd'].append(a_d_dd)
        hb_descr['d_a_aa'].append(d_a_aa)
        hb_descr['planar1'].append(planar1)
        hb_descr['planar2'].append(planar2)
        hb_descr['category'].append(atom_category)

    if __debug__:
        for category in hb_descr.values():
            assert len(category) == len(hb_descr['sno'])
        categories = list(hb_descr.values())[1:]
        by_sno = list(zip(*categories))
        for sno, per_sno in enumerate(by_sno):
            for per_hbond in per_sno:
                assert len(per_hbond) == len(per_sno[0])
    hb_descr = dict(hb_descr)
    return hb_descr

def _set_hb_descr(df_ATOM, hbond_vector, dsr_snos, d_info, a_info, role,
                  other_hb_info):
    d_cid, d_sno = d_info[:2]
    a_cid, a_sno = a_info[:2]
    residue = df_ATOM.loc[d_cid, d_sno]
    doner_trans, acc_trans = _transform_hbond_vector(hbond_vector, residue)

    if d_cid == a_cid and a_sno in dsr_snos:
        ext = "INT"
    else:
        ext = "EXT"
    return d_sno, [role, doner_trans, acc_trans, [ext, *other_hb_info]]

def _transform_hbond_vector(hbond_vector, residue):
    """
    Transform hbond vector to a standardised coordinate frame, described by
    N, C, CA in ATOM, for res specified in hb. First translate, then rotate.
    """
    N_coord, C_coord, CA_coord = _get_coordinates(residue)

    rotation_matrix = _get_rotation_matrix(N_coord, C_coord, CA_coord)

    doner_translated = np.array(hbond_vector[0]) - CA_coord
    acc_translated = np.array(hbond_vector[1]) - CA_coord

    doner_transformed = np.dot(doner_translated, rotation_matrix)
    acc_transformed = np.dot(acc_translated, rotation_matrix)

    doner_transformed = np.array([round(decimal.Decimal(val), 2)
                                  for val in doner_transformed])
    acc_transformed = np.array([round(decimal.Decimal(val), 2)
                                for val in acc_transformed])

    return doner_transformed, acc_transformed

def _get_rotation_matrix(N_coord, C_coord, CA_coord):
    CA_to_N = N_coord - CA_coord
    CA_to_C = C_coord - CA_coord
    normal = geometry.crossProduct(CA_to_N, CA_to_C)
    rotation_axis = geometry.crossProduct(normal, [0., 0., 1.])
    theta = geometry.findAngle(normal, [0., 0., 1.])
    rotation_matrix = geometry.genRotMatrix(rotation_axis, theta)
    return rotation_matrix

def _get_coordinates(residue):
    """
    Dangerous: Checking aname required, as the line
    N_coord, C_coord, CA_coord = residue.coord.values
    depends on the atom order being N, CA, C,
    which depends on a sort() earlier.

    Original being:
    # N_coord = residue[residue.aname == "N"].coord.values[0]
    # C_coord = residue[residue.aname == "C"].coord.values[0]
    # CA_coord = residue[residue.aname == "CA"].coord.values[0]
    Change made for performance.
    """
    N_coord = residue[residue['aname'] == "N"].coord.values[0]
    CA_coord = residue[residue['aname'] == "CA"].coord.values[0]
    C_coord = residue[residue['aname'] == "C"].coord.values[0]
    assert len(N_coord) == 3
    assert len(C_coord) == 3
    assert len(CA_coord) == 3
    return N_coord, C_coord, CA_coord

def _get_hbond_vector(d_info, a_info, ATOM_HETATM):
    """
    Note that reversing this (d_info => a_info) causes the hbond calculated to
    point from a => d. Can be tested, shown in output file.

    Setting a nested index (d_info) for performance reasons.
    """
    assert len(d_info) == 3
    d_coord = ATOM_HETATM.loc[d_info].coord
    a_coord = ATOM_HETATM.loc[a_info].coord
    if len(d_coord) != 3:
        d_coord = d_coord.values[0]
        a_coord = a_coord.values[0]
    assert len(d_coord) == 3
    assert len(a_coord) == 3
    vector = (d_coord, a_coord)
    return vector

##############################################################################

#
# def get_descr_hb_old(df_hb, df_ATOM, df_HETATM, dsr_snos):
#     """
#     Deprecated.
#
#     For old .hb files.
#     """
#     df_ATOM = df_ATOM.filter(items=["cid", "sno", 'aname', 'coord'])
#     df_HETATM = df_HETATM.filter(items=["cid", "sno", 'aname', 'coord'])
#     df_hb = df_hb.filter(items=["d_cid", "d_sno", 'd_aname',
#                                   "a_cid", "a_sno", 'a_aname',
#                                   'atom_category', ''])
#
#     ATOM_HETATM = df_ATOM.append(df_HETATM, ignore_index=True)
#     ATOM_HETATM = ATOM_HETATM.set_index(['cid', 'sno', 'aname'])
#     df_ATOM = df_ATOM.set_index(['cid', 'sno'])
#     df_ATOM = df_ATOM[df_ATOM.aname.isin(("N", "C", "CA"))]
#
#     hb_descr_raw = defaultdict(list)
#
#     for __, hb_line in df_hb.iterrows():
#         if hb_line.d_sno not in dsr_snos and hb_line.a_sno not in dsr_snos:
#             continue
#         d_info = (hb_line.d_cid, hb_line.d_sno, hb_line.d_aname)
#         a_info = (hb_line.a_cid, hb_line.a_sno, hb_line.a_aname)
#
#         atom_category = hb_line.atom_category
#
#         # Determine vector will point from doner => acceptor
#         hbond_vector = _get_hbond_vector(d_info, a_info, ATOM_HETATM)
#
#         if d_info[1] in dsr_snos:
#             role = "D"
#             sno, addition = _set_hb_descr_old(df_ATOM, hbond_vector, dsr_snos,
#                                      d_info, a_info, role, atom_category)
#             hb_descr_raw[sno].append(addition)
#
#         if a_info[1] in dsr_snos:
#             role = "A"
#             sno, addition = _set_hb_descr_old(df_ATOM, hbond_vector, dsr_snos,
#                                      a_info, d_info, role, atom_category[::-1])
#             hb_descr_raw[sno].append(addition)
#
#     hb_descr = _screen_duplicate(hb_descr_raw, dsr_snos)
#
#     return hb_descr
#
#
# def _set_hb_descr_old(df_ATOM, hbond_vector, dsr_snos, d_info, a_info, role,
#                  atom_category):
#     """
#     Deprecated.
#
#     For old .hb files.
#     """
#     d_cid, d_sno = d_info[:2]
#     a_cid, a_sno = a_info[:2]
#     residue = df_ATOM.loc[d_cid, d_sno]
#     doner_trans, acc_trans = _transform_hbond_vector(hbond_vector, residue)
#
#     if d_cid == a_cid and a_sno in dsr_snos:
#         ext = "INT"
#     else:
#         ext = "EXT"
#     return d_sno, [ext, role, atom_category, doner_trans, acc_trans]
