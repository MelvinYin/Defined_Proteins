from __future__ import with_statement
import math, sys, os, os.path
import Bio.PDB
#import Numeric
from utils import *
from pdb_data import prepareHbonds
from collections import defaultdict
# Output DEBUG msgs
DEBUG = 0

import csv
@handlify(mode='w')
def writeHbonds(fh, hbs):
    writer = csv.writer(fh, dialect='excel-tab')
    for pdb_id in sorted(hbs):
        hbs = hbs[pdb_id]
        for chain1 in sorted(hbs):
            hbs_ch1 = hbs[chain1]
            for resno1, icode1 in hbs_ch1:
                hbs_ch1_r1 = hbs_ch1[(resno1, icode1)]
                for donor_atom in sorted(hbs_ch1_r1):
                    hbs_ch1_r1_a1 = hbs_ch1_r1[donor_atom]
                    for chain2 in sorted(hbs_ch1_r1_a1):
                        hbs_ch1_r1_a1_ch2 = hbs_ch1_r1_a1[chain2]
                        for resno2, icode2 in hbs_ch1_r1_a1_ch2:
                            hbs_ch1_r1_a1_ch2_r2 = hbs_ch1_r1_a1_ch2[(resno2, icode2)]
                            for acceptor_atom in sorted(hbs_ch1_r1_a1_ch2_r2):
                                hb = hbs_ch1_r1_a1_ch2_r2[ acceptor_atom ]
                                writer.writerow([pdb_id, chain1, hb.resname1,
                                                 resno1, icode1, donor_atom,
                                                 chain2, hb.resname2, resno2,
                                                 icode2, acceptor_atom,
                                                 "%.2f" % hb.distance,
                                                 "%.2f" % hb.a_d_dd,
                                                 "%.2f" % hb.d_a_aa,
                                                 "%.2f" % hb.planar1,
                                                 "%.2f" % hb.planar2,
                                                 hb.type])

def tofloat(x):
    try:
        return float(x)
    except:
        return None
    
@handlify(mode='r')
def readHbonds(fh):
    r = defaultdict( innerDefaultdict( innerDefaultdict( innerDefaultdict( innerDefaultdict( innerDefaultdict( dict ) ) ) ) ) )
    reader = csv.reader(fh, dialect='excel-tab')
    for row in reader:
        pdb_id, chain1, resname1, resno1, icode1, donor_atom, chain2,  resname2, resno2, icode2, acceptor_atom, distance, a_d_dd, d_a_aa, planar1, planar2, htype = row
        x = StructObject(chain1 = chain1, resname1 = resname1, resno1 = int(resno1), icode1 = icode1, donor_atom = donor_atom,
                         chain2 = chain2, resno2 = int(resno2), resname2 = resname2, icode2 = icode2, acceptor_atom = acceptor_atom,
                         distance = tofloat(distance), a_d_dd = tofloat(a_d_dd), d_a_aa = tofloat(d_a_aa), planar1 = tofloat(planar1), planar2 = tofloat(planar2), type=htype )
        r[pdb_id][ chain1 ][ (x.resno1, x.icode1) ][ x.donor_atom ][ chain2 ][ (x.resno2, x.icode2) ][ x.acceptor_atom ] = x 
    return r

def hbonds(pdb_fname):
    pdb_fname = os.path.abspath(pdb_fname)
    with tempDir() as tmp_dir:
        hb_fname = prepareHbonds(pdb_fname)
        return readHbonds(hb_fname)

################################################################
# ROSE FORMAT PARSING

@handlify()
def parseRoseData(fh):
    lc = 0
    r = defaultdict( innerDefaultdict( innerDefaultdict( dict ) ) )
    for line in fh:
        lc+=1
        if lc==1:
            continue
        resname1, resno1, atom1, resname2, resno2, atom2, distance, d_a_aa, a_d_dd, planar1, planar2 = ( line[1:4], int(line[4:8]), line[13:16],
                                                                                                         line[14:22], int(line[22:26]), line[31:34],
                                                                                                         float(line[34:41]),
                                                                                                         float(line[43:49]), float(line[49:55]),
                                                                                                         float(line[55:61]), float(line[61:67]) )
        atom1 = atom1.strip()
        atom2 = atom2.strip()
        x = StructObject(resname1 = resname1, resno1 = resno1, atom1 = atom1, resname2 = resname2, resno2 = resno2, atom2 = atom2, distance = distance, d_a_aa = d_a_aa, a_d_dd = a_d_dd,
                         planar1 = planar1, planar2 = planar2)
        r[resno1][resno2][atom1][atom2] = x
    return r

@handlify()
def parseIgorData(fh):
    r = defaultdict( innerDefaultdict( innerDefaultdict( dict ) ) )        
    for line in fh:
        atom1, resno1, atom2, resno2, d_a_aa, a_d_dd, planar1, planar2, distance = line.strip().split()
        resno1, resno2 = int(resno1), int(resno2)
        d_a_aa, a_d_dd, planar1, planar2, distance = float(d_a_aa), float(a_d_dd), float(planar1), float(planar2), float(distance)
        x = StructObject(resno1 = resno1, atom1=atom1, resno2=resno2, atom2=atom2, distance=distance, d_a_aa = d_a_aa, a_d_dd = a_d_dd, planar1=planar1, planar2=planar2)
        r[resno1][resno2][atom1][atom2] = x
    return r

if __name__=='__main__':
    pdb_dir = 'test/'
    for file in os.listdir(pdb_dir):
        if file.endswith(".pdb"):
            pdb_file = 'test/' +  file
            my_bonds = hbonds(pdb_file)
            if not my_bonds:
                print(file)
            hb_file = 'test/results_' +  file[:-4] + ".txt"
            writeHbonds(hb_file, my_bonds)

    # #x = readHbonds(hb_file)
    # #igor_fn = sys.argv[2]
    # #igor_bonds = parseIgorData(igor_fn)
    # for chain1 in sorted(hbs):
    #     if chain1!=the_chain:
    #         continue
    #     hbs_ch1 = hbs[chain1]
    #     for value in hbs_ch1:
    #         hbs_ch1_r1 = hbs_ch1[value]
    #         for donor_atom in sorted(hbs_ch1_r1):
    #             hbs_ch1_r1_a1 = hbs_ch1_r1[donor_atom]
    #             for chain2 in sorted(hbs_ch1_r1_a1):
    #                 if chain2!=chain1:
    #                     continue
    #                 hbs_ch1_r1_a1_ch2 = hbs_ch1_r1_a1[chain2]
    #                 for val2 in hbs_ch1_r1_a1_ch2:
    #                     hbs_ch1_r1_a1_ch2_r2 = hbs_ch1_r1_a1_ch2[val2]
    #                     for acceptor_atom in sorted(hbs_ch1_r1_a1_ch2_r2):
    #                         hb = hbs_ch1_r1_a1_ch2_r2[ acceptor_atom ]