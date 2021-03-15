#####################################################################################################################################################################
#
# Grzegorz M Koczyk (2007-)
# 
# Preparation of hydrogen atoms on an existing data structure.
#
#####################################################################################################################################################################
import sys, os, os.path, math
import Bio.PDB
from Bio.PDB.Atom import Atom
from Bio.PDB.Polypeptide import three_to_one
#####################################################################################################################################################################
# Constants&procedures for cleaning up residues and atoms....
# Hydrogen atoms to prune from structures
HYDRO_ATOMS = frozenset( [ "H", "D", "1H", "2H", "3H", "1D", "2D", 
"3D", "HA", "1HA", "2HA", "3HA", "1HB", "2HB", "3HB", "HB", 
"1HG", "2HG", "HG", "1HD", "2HD", "1HE", "2HE", "3HE", "HZ", 
"HH2", "1HG1", "2HG1", "3HG1", "1HG2", "2HG2", "3HG2", "HD1", "1HD1", 
"2HD1", "3HD1", "1HD2", "2HD2", "3HD2", "HD2", "HE1", "HE2", "HE3", 
"HZ3", "HH2", "HZ2", "HD1", "1HD2", "2HD2", "1DD2", "2DD2", "HE", 
"DE", "HE1", "HE2", "1HE2", "2HE2", "1DE2", "2DE2", "1HH1", "2HH1", 
"1DH1", "2DH1", "1HH2", "2HH2", "1DH2", "2DH2", "1HZ", "2HZ", "3HZ", 
"1DZ", "2DZ", "3DZ", "HD2", "HE2", "HG", "HG1", "DG1", "HH", 
"DH", "HG" ])

def isAtomToCleanup(atom):
    # Strips hydrogen and deuterium
    name = atom.id
    if name in HYDRO_ATOMS:
        return True
    else:
        return False

#from Bio.PDB import to_one_letter_code
def isResidueToCleanup(residue):
    # Strips water and hetero residues, as well as inserted residues
    if (residue.resname in to_one_letter_code):
        return False
    if residue.id[0]!=' ' or residue.id[0]=='W' or (residue.resname not in to_one_letter_code):
        return True
    #elif residue.id[-1].strip():
    #    return True
    else:
        return False        

#####################################################################################################################################################################
# CALCULATIONS ON ATOMS
def subatom(a1, a2):
    ''' Difference in atomic coordinates '''
    return a1.coord - a2.coord

def length(v):
    ''' Length of a vector '''
    return math.sqrt( sum(v_*v_ for v_ in v) )

def vsum(*args):
    ''' Sum all corresponding positions in same length vectors. '''
    return [ sum(ar) for ar in zip(*args) ]

def mult(v1, v2):
    ''' Cross-product of two 3D vectors '''
    return [ v1[1]*v2[2]-v1[2]*v2[1],             
             v1[2]*v2[0]-v1[0]*v2[2],
             v1[0]*v2[1]-v1[1]*v2[0] ]

def fix0(atom0, atom1, atom2, *args):
    ''' NH group on the peptide backbone. '''
    v1 = subatom(atom0, atom1)
    v2 = subatom(atom0, atom2)
    f1=length(v1)
    f2=length(v2)
    # Scale contents of distance vectors
    v1 = [ v_/f1 for v_ in v1 ]
    v2 = [ v_/f2 for v_ in v2 ]
    r = vsum(v1,v2)
    f=length(r)
    return [ [ r[i]/f*1.08+atom0.coord[i] for i in range(0,3) ] ]

def fix1(atom0, atom1, atom2, *args):
    ''' Sulfhydryl group on cysteine residue, hydroxyl groups... (sulphur, oxygen basically). Strangely enough - 1 of the ARG nitrogens and carbonyl group hydrogens also... '''
    v1 = subatom(atom0, atom1)
    v2 = subatom(atom0, atom2)
    f1=length(v1);
    v1 = [ v_/f1*0.5*1.08 for v_ in v1 ]
    v3 = mult(v2,v1)
    v2 = mult(v3,v1)
    f2=length(v2);
    v2 = [ v_/f2*0.866*1.08 for v_ in v2 ]
    return [ [ v2[i] + v1[i] + atom0.coord[i] for i in range(0,3) ] ]

fix2 = fix0
def fix3(atom0, atom1, atom2, *args):
    ''' NZ of lysine, ND2, NE2 of Asn and Gln '''
    v1 = subatom(atom0, atom1)
    v2 = subatom(atom0, atom2)
    f1=length(v1);
    v1 = [ v_/f1*0.5*1.08 for v_ in v1 ]
    v3 = mult(v2,v1)
    v2 = mult(v3,v1)
    f2=length(v2)
    v2 = [ v_/f2*0.866*1.08 for v_ in v2 ]
    return [ [ v1[i] + v2[i] + atom0.coord[i]
               for i in range(0,3) ],
             [ v1[i] - v2[i] + atom0.coord[i]
               for i in range(0,3) ],
             ]

def fix4(atom0, atom1, atom2, atom3, *args):
    ''' Hydrogen atom in the cases where all other positions are already taken (excluding aromatic etc.) '''
    v1 = subatom(atom0, atom1)
    v2 = subatom(atom0, atom2)
    v3 = subatom(atom0, atom3)
    f1=length(v1)
    f2=length(v2)
    f3=length(v3)
    v1 = [ v_/f1 for v_ in v1 ]
    v2 = [ v_/f2 for v_ in v2 ]
    v3 = [ v_/f3 for v_ in v3 ]
    r = vsum(v1, v2, v3)
    f=length(r)
    return [ [ r[i]/f*1.08+atom0.coord[i] for i in range(0,3) ] ]

def fix5(atom0, atom1, atom2, *args):
    ''' 1HA, 2HA on GLY and in the CH2 (methylene) groups '''
    v1 = subatom(atom0, atom1)
    v2 = subatom(atom0, atom2)
    v4 = mult(v1, v2)
    v3 = vsum(v1, v2)
    f1=length(v3)
    f2=length(v4)
    v3 = [ v_/f1*0.5*1.08 for v_ in v3 ]
    v4 = [ v_/f2*0.866*1.08 for v_ in v4 ]
    return [ [ v3[i] + v4[i] + atom0.coord[i]
               for i in range(0,3) ],
             [ v3[i] - v4[i] + atom0.coord[i]
               for i in range(0,3) ],
             ]

def _createRotationMatrix(pivot, angle0):
    angle = angle0*3.14159/180.*0.5
    lambda_ = pivot[0]*math.sin(angle)
    mu = pivot[1]*math.sin(angle)
    nu = pivot[2]*math.sin(angle)
    ro = math.cos(angle)
    rm = [ range(0,3) for i in range(0,3) ]
    rm[0][0] = lambda_*lambda_-mu*mu-nu*nu+ro*ro
    rm[1][0] = 2.*(lambda_*mu+nu*ro)
    rm[2][0] = 2.*(lambda_*nu-mu*ro)
    rm[0][1] = 2.*(lambda_*mu-nu*ro)
    rm[1][1] = mu*mu-nu*nu-lambda_*lambda_+ro*ro
    rm[2][1] = 2.*(mu*nu+lambda_*ro)
    rm[0][2] = 2.*(lambda_*nu+mu*ro)
    rm[1][2] = 2.*(mu*nu-lambda_*ro)
    rm[2][2] = nu*nu-mu*mu-lambda_*lambda_+ro*ro
    return rm
def _rotateVector(v, rm):
    return [ sum(v[i]*rm[j][i] for i in range(0,3)) for j in range(0,3) ]

def fix6(atom0, atom1, atom2, *args):
    ''' Adding three hydrogen atoms in those -CH3 groups (based on the two preceding atoms?)'''
    v1 = subatom(atom0, atom1)
    v2 = subatom(atom0, atom2)
    f1=length(v1)
    v1 = [ v_/f1*0.5*1.08 for v_ in v1 ]
    v3 = mult(v2,v1)
    v2 = mult(v1,v3)
    f2 = length(v2)
    f1_= length(v1)
    r0 = [ v2_/f2*0.866*1.08 + v1_ for v1_,v2_ in zip(v1, v2) ]
    v1 = [ v_/f1_ for v_ in v1 ]
    rm = _createRotationMatrix(v1,120.)
    r1 = _rotateVector(r0,rm)
    r2 = _rotateVector(r1,rm)
    return [ [r_[i]+atom0.coord[i] for i in range(0,3)]  for r_ in (r0, r1, r2)]

class StructObject(object):
    def __init__(self, **dictargs):    
        for k, v in dictargs.items():
            object.__setattr__(self,k, v)
    def __setattr__(self, k, v):
        object.__setattr__(self, k, v)
    def __str__(self):
        return "( " + ", ".join(["%s: %s" % (k, self.__dict__[k]) for k in sorted(self.__dict__)]) + " )"
    def __repr__(self):
        return self.__str__()

# Description of how fixes of hydrogen atoms should be performed (on which atom, using which of the geometric procedures, and in relation to what atom)
fix_dispatch = { 'GLY': { 'CA': StructObject(fix=fix5, delta1='N', delta2='C',),                          
                          },
                 'ALA': { 'CA': StructObject(fix=fix4, delta1='N', delta2='C', delta3='CB'),
                          'CB': StructObject(fix=fix6, delta1='CA', delta2='C',),
                          },
                 'VAL': { 'CA': StructObject(fix=fix4, delta1='N', delta2='C', delta3='CB'),
                          'CB': StructObject(fix=fix4, delta1='CA', delta2='CG1', delta3='CG2'),
                          'CG1': StructObject(fix=fix6, delta1='CB', delta2='CG2',),
                          'CG2': StructObject(fix=fix6, delta1='CB', delta2='CG1',),
                          },
                 'LEU': { 'CA': StructObject(fix=fix4, delta1='N', delta2='C', delta3='CB',),
                          'CB': StructObject(fix=fix5, delta1='CA', delta2='CG' ),
                          'CG': StructObject(fix=fix4, delta1='CB', delta2='CD1', delta3='CD2' ),
                          'CD1': StructObject(fix=fix6, delta1='CG', delta2='CD2' ),
                          'CD2': StructObject(fix=fix6, delta1='CG', delta2='CD1' ),
                          },
                 'ILE': { 'CA': StructObject(fix=fix4, delta1='N', delta2='C', delta3='CB' ),
                          'CB': StructObject(fix=fix4, delta1='CA', delta2='CG1', delta3='CG2' ),
                          'CG1': StructObject(fix=fix5, delta1='CB', delta2='CD1' ),
                          'CG2': StructObject(fix=fix6, delta1='CB', delta2='CG1',),
                          'CD1': StructObject(fix=fix6, delta1='CG1', delta2='CB',),
                          },
                 'MET':  { 'CA': StructObject(fix=fix4, delta1='N', delta2='C', delta3='CB' ),
                           'CB': StructObject(fix=fix5, delta1='CA', delta2='CG' ),
                           'CG': StructObject(fix=fix5, delta1='CB', delta2='SD' ),
                           'CE': StructObject(fix=fix6, delta1='SD', delta2='CG' )
                           },
                 'CYS':  { 'CA': StructObject(fix=fix4, delta1='N', delta2='C', delta3='CB' ),
                           'CB': StructObject(fix=fix5, delta1='CA', delta2='SG' ),
                           'SG': StructObject(fix=fix1, delta1='CB', delta2='CA' ),
                           },
                 'THR':  { 'CA': StructObject(fix=fix4, delta1='N', delta2='C', delta3='CB' ),
                           'CB': StructObject(fix=fix4, delta1='CA', delta2='OG1', delta3='CG2' ),
                           'OG1': StructObject(fix=fix1, delta1='CB', delta2='CA' ),
                           'CG2': StructObject(fix=fix6, delta1='CB', delta2='OG1' ),
                           },
                 'SER':  { 'CA': StructObject(fix=fix4, delta1='N', delta2='C', delta3='CB' ),
                           'CB': StructObject(fix=fix5, delta1='CA', delta2='OG' ),
                           'OG': StructObject(fix=fix1, delta1='CB', delta2='CA' ),
                           },
                 'LYS': { 'CA': StructObject(fix=fix4, delta1='N', delta2='C', delta3='CB' ),
                          'CB': StructObject(fix=fix5, delta1='CA', delta2='CG' ),
                          'CG': StructObject(fix=fix5, delta1='CB', delta2='CD' ),
                          'CD': StructObject(fix=fix5, delta1='CG', delta2='CE' ),
                          'CE': StructObject(fix=fix5, delta1='CD', delta2='NZ' ),
                          'NZ': StructObject(fix=fix3, delta1='CE', delta2='CD' ),
                          },
                 'ARG': { 'CA': StructObject(fix=fix4, delta1='N', delta2='C', delta3='CB' ),
                          'CB': StructObject(fix=fix5, delta1='CA', delta2='CG' ),
                          'CG': StructObject(fix=fix5, delta1='CB', delta2='CD' ),
                          'CD': StructObject(fix=fix5, delta1='CG', delta2='NE' ),
                          'NE': StructObject(fix=fix0, delta1='CD', delta2='CZ' ),
                          'NH1': StructObject(fix=fix3, delta1='CZ', delta2='NH2' ),
                          'NH2': StructObject(fix=fix1, delta1='CZ', delta2='NH1' ),
                          },
                 'ASP': { 'CA': StructObject(fix=fix4, delta1='N', delta2='C', delta3='CB' ),
                          'CB': StructObject(fix=fix5, delta1='CA', delta2='CG'),
                          'OD2': StructObject(fix=fix1, delta1='CG', delta2='OD1' ),
                          },
                 'ASN': { 'CA': StructObject(fix=fix4, delta1='N', delta2='C', delta3='CB' ),
                          'CB': StructObject(fix=fix5, delta1='CA', delta2='CG'),
                          'ND2': StructObject(fix=fix3, delta1='CG', delta2='OD1' ),
                          },

                 'GLU': { 'CA': StructObject(fix=fix4, delta1='N', delta2='C', delta3='CB' ),
                          'CB': StructObject(fix=fix5, delta1='CA', delta2='CG'),
                          'CG': StructObject(fix=fix5, delta1='CB', delta2='CD'),
                          'OE2': StructObject(fix=fix1, delta1='CD', delta2='OE1' ),
                          },
                 'GLN': { 'CA': StructObject(fix=fix4, delta1='N', delta2='C', delta3='CB' ),
                          'CB': StructObject(fix=fix5, delta1='CA', delta2='CG'),
                          'CG': StructObject(fix=fix5, delta1='CB', delta2='CD'),
                          'NE2': StructObject(fix=fix3, delta1='CD', delta2='OE1' ),
                          },
                 'PHE': { 'CA': StructObject(fix=fix4, delta1='N', delta2='C', delta3='CB' ),
                          'CB': StructObject(fix=fix5, delta1='CA', delta2='CG'),
                          'CD1': StructObject(fix=fix2, delta1='CG', delta2='CE1'),
                          'CD2': StructObject(fix=fix2, delta1='CD1', delta2='CE2'),
                          'CE1': StructObject(fix=fix2, delta1='CD1', delta2='CZ'),
                          'CE2': StructObject(fix=fix2, delta1='CD2', delta2='CZ'),
                          'CZ': StructObject(fix=fix2, delta1='CE1', delta2='CE2'),
                          },
                 'TYR': { 'CA': StructObject(fix=fix4, delta1='N', delta2='C', delta3='CB' ),
                          'CB': StructObject(fix=fix5, delta1='CA', delta2='CG'),
                          'CD1': StructObject(fix=fix2, delta1='CG', delta2='CE1'),
                          'CD2': StructObject(fix=fix2, delta1='CD1', delta2='CE2'),
                          'CE1': StructObject(fix=fix2, delta1='CD1', delta2='CZ'),
                          'CE2': StructObject(fix=fix2, delta1='CD2', delta2='CZ'),
                          'OH': StructObject(fix=fix1, delta1='CZ', delta2='CE2'),
                          },
                 'TRP': { 'CA': StructObject(fix=fix4, delta1='N', delta2='C', delta3='CB' ),
                          'CB': StructObject(fix=fix5, delta1='CA', delta2='CG'),
                          'CD1': StructObject(fix=fix2, delta1='CG', delta2='NE1'),
                          'NE1': StructObject(fix=fix2, delta1='CD1', delta2='CE2'),
                          'CE3': StructObject(fix=fix2, delta1='CD2', delta2='CZ3'),
                          'CZ2': StructObject(fix=fix2, delta1='CE2', delta2='CH2'),
                          'CZ3': StructObject(fix=fix2, delta1='CE3', delta2='CH2'),
                          'CH2': StructObject(fix=fix2, delta1='CZ3', delta2='CZ2'),
                          },
                 'HIS': { 'CA': StructObject(fix=fix4, delta1='N', delta2='C', delta3='CB' ),
                          'CB': StructObject(fix=fix5, delta1='CA', delta2='CG'),
                          'ND1': StructObject(fix=fix2, delta1='CG', delta2='CE1'),
                          'CD2': StructObject(fix=fix2, delta1='CG', delta2='NE2'),
                          'CE1': StructObject(fix=fix2, delta1='ND1', delta2='NE2'),
                          },
                 'PRO': { 'CA': StructObject(fix=fix4, delta1='N', delta2='C', delta3='CB' ),
                          'CB': StructObject(fix=fix5, delta1='CA', delta2='CG'),
                          'CG': StructObject(fix=fix5, delta1='CB', delta2='CD'),
                          'CD': StructObject(fix=fix5, delta1='CG', delta2='N'),
                          },
                 }

fix_dispatch['PTR']=fix_dispatch['TYR']

def addHydrogens(residue, last_c=None, last_no=None):
    ''' Add hydrogens to a given residue, renumbering atoms appropriately. Will also convert disordered atoms into their respective first representatives.'''
    rname = residue.resname
    new_atoms = []
    if last_no is None:
        last_no = 0
    if rname in fix_dispatch:
        fr = fix_dispatch[rname]
        for atom in list(residue):
            name = atom.id
            atom.serial_number = last_no            
            new_atoms.append( (name, atom) )
            last_no+=1
            # Special case is the hydrogen on peptide backbone - dispatch manually depending on last_c
            # Proline skips this step
            if name=='N' and rname!='PRO' and last_c is not None and 'CA' in residue.child_dict and 'N' in residue.child_dict:
                nname = ' H  '
                b_factor = 0.0
                occupancy = 1.0
                altloc = ' '
                fullname = name
                coords = fix0(residue['N'], last_c, residue['CA'])[0]
                new_atom = Atom(nname.strip(), coords, b_factor, occupancy, altloc, nname, last_no) 
                last_no+=1
                new_atoms.append( (new_atom.id, new_atom) )
                residue.add(new_atom)
            else:
                if name in fr:
                    f = fr[name]
                    # Sometimes a third atom is needed
                    f_delta3 = getattr(f, 'delta3', None)
                    # Check all atoms needed for adding are in place
                    if f.delta1 in residue.child_dict and f.delta2 in residue.child_dict and (f_delta3 is None or f_delta3 in residue.child_dict):
                        fixed = f.fix(atom, residue[f.delta1], residue[f.delta2], (residue[f_delta3] if f_delta3 else None))
                        # Create names according to standard - prefix is only added if more than one hydrogen is introduced
                        if len(fixed)>1:
                            names = ['%dH%s' % (i+1, atom.fullname[2:]) for i in range(0,len(fixed)) ]
                        else:
                            names = [' H%s' % atom.fullname[2:] ]
                        for nname, coords in zip(names, fixed):
                            b_factor = 0.0
                            occupancy = 1.0
                            altloc = ' '
                            fullname = name
                            new_atom = Atom(nname.strip(), coords, b_factor, occupancy, altloc, nname, last_no) 
                            last_no+=1
                            new_atoms.append( (new_atom.id, new_atom) )
                            residue.add(new_atom)
        # Restructure child lists to preserve new atom order
        residue.child_list = [v[1] for v in new_atoms]
        residue.child_dict = dict(new_atoms)
    if new_atoms:
        return new_atoms, last_no
    else:
        return None, last_no

def fixDisorderedResidue(residue):
    ''' Fixes disordered residues by returning their first respective representatives (same for individual disordered atoms). '''
    # Replacing disordered residue by first representative
    if isinstance(residue, Bio.PDB.Residue.DisorderedResidue):
        residue = residue.child_dict[ sorted(residue.child_dict.keys())[0] ]
    # Replacing disordered atoms by first representative
    atoms = []
    for atom in residue:
        if isinstance(atom, Bio.PDB.Atom.DisorderedAtom):
            atom = atom.child_dict[ sorted(atom.child_dict.keys() )[0] ]
            atom.disordered_flag=0            
        atoms.append(atom)
    # Replacing the residue children lists
    residue.child_list = atoms
    residue.child_dict = dict( (atom.id, atom) for atom in atoms )
    residue.disordered = False
    return residue
    
def remakeHydrogens(struct, cleanup_res=True):
  ''' Cleans up hydrogens, heteroatoms/non-aminoacids, fixes disordered atoms/residues and reintroduces hydrogens to a structure (will modify the structure inplace)'''
  for model in struct:
      for chain in model:
          last_c = None
          if cleanup_res:
              chain.child_list = [ fixDisorderedResidue(residue) for residue in chain if not isResidueToCleanup(residue)]
              chain.child_dict = dict( (residue.id, residue) for residue in chain )
          last_no = None
          for residue in chain:
              # In the case of discontinuities, in the chain - do not attempt to add hydrogen to peptide backbone
              if last_c is not None and last_c.get_parent().id[1]!=residue.id[1]-1:
                  last_c=None
              for atom in list(residue):
                  if isAtomToCleanup(atom):
                      try:
                          residue.detach_child(atom.id)
                      except:
                          print >> sys.stderr, atom
                          print >> sys.stderr, residue, type(residue)
                          print >> sys.stderr, residue.child_list
                          raise
              new_atoms, last_no = addHydrogens(residue, last_c, last_no)
              # TODO: Add rebuilding of residue
              try:
                  last_c = residue['C']
              except:
                  last_c = None
                    
if __name__=='__main__':
    pdb_fn = sys.argv[1]
    parser = Bio.PDB.PDBParser()
    struct = parser.get_structure('query', pdb_fn)
    remakeHydrogens(struct)
    io=Bio.PDB.PDBIO()
    io.set_structure(struct)
    io.save(sys.stdout)
       
