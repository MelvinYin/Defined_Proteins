#####################################################################################################################################################################
#
# Grzegorz M Koczyk (2007-)
# 
# Preparation of PDB files (cleaning up, vdW computation)
#
#####################################################################################################################################################################
from __future__ import with_statement
import sys, os, os.path
import logging
from contextlib import closing
import gzip, io, shutil, subprocess
from utils import *
from hprep import remakeHydrogens
import Bio.PDB

######################################################################################################################################################################
# Constants
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
"DH", "HG" ]);

# Path to executables
if 'DHCL_EXEC_ROOT' in os.environ:
    EXEC_ROOT = os.environ['DHCL_EXEC_ROOT']
else:
    EXEC_ROOT = os.path.dirname( os.path.abspath(sys.argv[0]))
# Executable for preparing hydrogen atoms
PREP23_EXECUTABLE = os.path.join(EXEC_ROOT, ("prep23.exe" if 'win' in sys.platform else 'prep23') )
# Data for above exec
PREP23_DATA = os.path.join(EXEC_ROOT, "Hi.dat")
# VdW calculation executable
VDW_EXEC = os.path.join(EXEC_ROOT, "vdw")
# Contacts calculation executable
CONTACTS_EXEC = os.path.join(EXEC_ROOT, "contacts")
# Contacts calculation executable
HBONDS_EXEC = os.path.join(EXEC_ROOT, "bin", "hbonds")

######################################################################################################################################################################
# Classes
class NoHydroSelect(Bio.PDB.Select):
    ''' Selector class for pruning hydrogens, inserted residues and heteroatoms '''
    def accept_residue(self, residue):
        # Strips water and hetero residues, as well as inserted residues
        if residue.id[0]!=' ' or residue.id[0]=='W':
            return False
        elif residue.id[-1].strip():
            return False
        else:
            return True
    def accept_atom(self, atom):
        #if atom.get_name().strip() in hydro_atoms:
        name = atom.get_name()
        # Strips hydrogen and deuterium
        if name in HYDRO_ATOMS:
            return False
        else:
            return True

class BasicSelect(Bio.PDB.Select):
    ''' Selector class for pruning inserted residues and heteroatoms '''
    def accept_residue(self, residue):
        # Strips water and hetero residues, as well as inserted residues
        if residue.id[0]!=' ' or residue.id[0]=='W':
            return False
        elif residue.id[-1].strip():
            return False
        else:
            return True
######################################################################################################################################################################
# Functions - preparation, computation
def _preparePrepExec():
    ''' Prepare prep23 executable and data in the current directory '''
    # Copy to current directory prep23 executable
    shutil.copy(PREP23_EXECUTABLE, '.')
    shutil.copy(PREP23_DATA, '.')
    if 'win32' in sys.platform:
        return 'prep23'
    else:
        return './prep23'

def preparePdb(pdb_fname, out_pdb_fname):
    ''' Prepare the PDB file with only first model and redundancies cut out '''
    # 'Absolutize' the path names - rest is done in the temporary dir
    pdb_fname = os.path.abspath(pdb_fname)
    if not os.path.exists(pdb_fname):
        raise IOError('%s does not exist' % pdb_fname)
    out_pdb_fname = os.path.abspath(out_pdb_fname)
    # Inside the temporary dir
    with tempDir() as tmp_dir:
        # Temporary names for curated input and output files
        new_pdb_fname = 'query.pdb'
        out_tmp_fname = 'out.pdb'
        # If the original PDB is packed with gzip - unpack it into a new file
        if pdb_fname.endswith('.gz'):
            rfh = gzip.open(pdb_fname,'r')
        else:
            rfh = open(pdb_fname, 'r')
        try:
            with open(new_pdb_fname, 'w') as wfh:
                wfh.write(rfh.read())
        finally:
            rfh.close()
        # Parse structure
        # Redirect standard output/error to a cStringIO,
        #so that PDBParser stops messing the output
        parser = Bio.PDB.PDBParser()
        err_fh = io.StringIO()
        sys.stdout = err_fh
        sys.stderr = err_fh
        struct = parser.get_structure('query', new_pdb_fname)
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        # Output formatted info about PDBParser's work to a log
        s=err_fh.getvalue()
        if s.strip():
            logging.info("Structure parsing generated following error message(s): \n%s\n%s\n%s" % ('-'*120, s, '-' * 120) )
        # By default use only first model
        model = struct[0]
        del struct.child_list[1:]
        # Check for discontinuities greater than 5 residues - warn about this _specifically_
        for chain in model:
            chid = chain.id
            last_rid = None
            for residue in chain:
                if last_rid is not None and rid>last_rid+5:
                    rid = residue.id[1]
                    logging.warn("Residues %s:%s-%s:%s. Results might be inaccurate, as the break in a protein chain numbers more than 5 residues." % (last_id, chain, rid, chain) )
                    last_rid = rid
        # Save structure without hydrogens
        io=Bio.PDB.PDBIO()
        io.set_structure(struct)
        io.save(new_pdb_fname)
        shutil.move(new_pdb_fname, out_pdb_fname)
        return out_pdb_fname

# TODO: add cleaning up disordered atoms/residues (prep leaves them)
def prepareWithHydrogensPrep23(pdb_fname, out_pdb_fname="wth_hydro.pdb"):
    ''' Prepare the PDB file with hydrogen data (clean up and create a new one). '''
    # 'Absolutize' the path names - rest is done in the temporary dir
    pdb_fname = os.path.abspath(pdb_fname)
    if not os.path.exists(pdb_fname) or not os.path.isfile(pdb_fname):
        raise IOError('%s does not exist or is not a file.' % pdb_fname)
    out_pdb_fname = os.path.abspath(out_pdb_fname)
    # Inside the temporary dir
    with tempDir() as tmp_dir:
        # Temporary names for curated input and output files
        new_pdb_fname = 'query.pdb'
        out_tmp_fname = 'out.pdb'
        # Prepare the sources
        prep_exec = _preparePrepExec()
        # Copy the original file into our temporary directory
        # If the original PDB is packed with gzip - unpack it into a new file
        if pdb_fname.endswith('.gz'):
            rfh = gzip.open(pdb_fname,'r')
        else:
            rfh = open(pdb_fname, 'r')
        try:
            with open(new_pdb_fname, 'w') as wfh:
                wfh.write(rfh.read())
        finally:
            rfh.close()
        # Parse structure
        # Redirect standard output/error to a cStringIO,
        #so that PDBParser stops messing the output
        parser = Bio.PDB.PDBParser()
        err_fh = io.StringIO()
        sys.stdout = err_fh
        sys.stderr = err_fh
        with open(new_pdb_fname, 'r') as rfh:
            struct = parser.get_structure('query', rfh)
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        # Output formatted info about PDBParser's work to a log
        s=err_fh.getvalue()
        if s.strip():
            logging.info("Structure parsing generated following error message(s): \n%s\n%s\n%s" % ('-'*120, s, '-' * 120) )
        # By default use only first model
        model = struct[0]
        del struct.child_list[1:]
        # Check for discontinuities greater than 5 residues - warn about this _specifically_
        for chain in model:
            chid = chain.id
            last_rid = None
            # Curate disordered residues keeping only the last
            chain.child_list = [residue for residue in chain]
            chain.child_dict = dict( (residue.id, residue) for residue in chain )
            for residue in chain:
                # Curate disordered atoms keeeping only the last
                residue.child_list = [a for a in residue]
                residue.child_dict = dict( (a.id, a) for a in residue )
                if last_rid is not None and rid>last_rid+5:
                    rid = residue.id[1]
                    logging.warn("Residues %s:%s-%s:%s. Results might be inaccurate, as the break in a protein chain numbers more than 5 residues." % (last_id, chain, rid, chain) )
                    last_rid = rid
        # Save structure without hydrogens
        io=Bio.PDB.PDBIO()
        io.set_structure(struct)
        io.save(new_pdb_fname,NoHydroSelect())
        # Run the preparation executable on the newly created PDB file
        if ( subprocess.call( "%s %s %s 1>tmp.out 2>tmp.err" % (prep_exec, new_pdb_fname, out_tmp_fname),  shell=True) != 0):
            raise RuntimeError('Could not prepare corrected structure file for %s' % pdb_fname)
        # Fix the occupancies (creating the last and final temporary PDB file)
        final_fn = "final.pdb"
        #raw_input('WAITING...')
        with open(out_tmp_fname, 'r') as rfh:
            with open(final_fn, 'w') as wfh:
                for line in rfh:
                    if line.startswith('ATOM'):
                        print >> wfh, line[:-1] + "  0.00  0.00           C"
                    else:
                        print >> wfh, line,
        # Move the output file to the desired location
        shutil.move(final_fn, out_pdb_fname)
        return out_pdb_fname

def prepareWithHydrogens(pdb_fname, out_pdb_fname="wth_hydro.pdb"):
    ''' Prepare the PDB file with hydrogen data (clean up and create a new one). '''
    # 'Absolutize' the path names - rest is done in the temporary dir
    pdb_fname = os.path.abspath(pdb_fname)
    if not os.path.exists(pdb_fname) or not os.path.isfile(pdb_fname):
        raise IOError('%s does not exist or is not a file.' % pdb_fname)
    out_pdb_fname = os.path.abspath(out_pdb_fname)
    if pdb_fname.endswith('.gz'):
        rfh = gzip.open(pdb_fname,'r')
        #print pdb_fname
    else:
        rfh = open(pdb_fname, 'r')
    try:
        # Parse structure
        parser = Bio.PDB.PDBParser()
        # Redirect standard output/error to a cStringIO,
        # so that PDBParser stops messing the output
        err_fh = io.StringIO()
        sys.stdout = err_fh
        sys.stderr = err_fh
        struct = parser.get_structure('query', rfh)
    finally:
        # Restore streams
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        # ... and close up
        rfh.close()
    # Output formatted info about PDBParser's work to a logger
    s=err_fh.getvalue()
    if s.strip():
        logging.info("Structure parsing generated following error message(s): \n%s\n%s\n%s" % ('-'*120, s, '-' * 120) )
    # By default use only first model
    # ... delete the rest
    model = struct[0]
    del struct.child_list[1:]
    # Check for discontinuities greater than 5 residues - warn about this _specifically_ (into the logger, again)
    for chain in model:
        chid = chain.id
        last_rid = None
        for residue in chain:
            if last_rid is not None and rid>last_rid+5:
                rid = residue.id[1]
                logging.warn("Residues %s:%s-%s:%s. Results might be inaccurate, as the break in a protein chain numbers more than 5 residues." % (last_id, chain, rid, chain) )
                last_rid = rid
    # Prepare the remade hydrogens
    remakeHydrogens(struct)
    # Save structure
    if out_pdb_fname.endswith('.gz'):
        with closing( gzip.open(out_pdb_fname,'w')) as wfh:
            io=Bio.PDB.PDBIO()
            io.set_structure(struct)
            io.save(wfh)
    else:
        io=Bio.PDB.PDBIO()
        io.set_structure(struct)
        io.save(out_pdb_fname)
    return out_pdb_fname

def prepareVdws(pdb_fname, atom_cs_fname="query.atom.cns", residue_cs_fname="query.residue.cns", bounds_fname="query.bds"):
    ''' Calculate the van der Waals data. '''
    # 'Absolutize' the path names - rest is done in the temporary dir
    pdb_fname = os.path.abspath(pdb_fname)
    curdir = os.getcwd()
    # Inside the temporary dir
    with tempDir() as tmp_dir:
        # Run the vdw executable on the newly created PDB file
        if ( subprocess.call( "%s %s %s %s %s" % (VDW_EXEC, pdb_fname, atom_cs_fname, residue_cs_fname, bounds_fname),  shell=True) != 0):
            raise RuntimeError("Could not prepare Lennard-Jones' results for %s" % pdb_fname)
        shutil.move(atom_cs_fname, os.path.join(curdir, atom_cs_fname))
        shutil.move(residue_cs_fname, os.path.join(curdir, residue_cs_fname))
        shutil.move(bounds_fname, os.path.join(curdir, bounds_fname))
        return atom_cs_fname, residue_cs_fname, bounds_fname

def prepareContacts(pdb_fname, atom_cs_fname="query.atom_atom.cns", bounds_fname="query.bds"):
    ''' Calculate the van der Waals data. '''
    # 'Absolutize' the path names - rest is done in the temporary dir
    pdb_fname = os.path.abspath(pdb_fname)
    atom_cs_fname = os.path.abspath(atom_cs_fname)
    bounds_fname = os.path.abspath(bounds_fname)
    curdir = os.getcwd()
    # Inside the temporary dir
    with tempDir() as tmp_dir:
        # Run the vdw executable on the newly created PDB file
        if ( subprocess.call( "%s %s %s %s" % (CONTACTS_EXEC, pdb_fname, atom_cs_fname, bounds_fname),  shell=True) != 0):
            raise RuntimeError("Could not prepare contacts' results for %s" % pdb_fname)
        #raw_input('WAIT...')
        #shutil.move(atom_cs_fname, curdir)
        #shutil.move(bounds_fname, curdir)
        return atom_cs_fname, bounds_fname

def prepareHbonds(pdb_fname, hbonds_fname="query.hbonds"):
    ''' Calculate the hydrogen bonds data. '''
    # 'Absolutize' the path names - rest is done in the temporary dir
    pdb_fname = os.path.abspath(pdb_fname)
    curdir = os.getcwd()
    # Inside the temporary dir
    with tempDir() as tmp_dir:
        # Run the vdw executable on the newly created PDB file
        if ( subprocess.call( "%s %s %s" % (HBONDS_EXEC, pdb_fname, hbonds_fname),  shell=True) != 0):
            raise RuntimeError("Could not prepare hydrogen bonds results for %s" % pdb_fname)
        shutil.move(hbonds_fname, os.path.join(curdir, hbonds_fname))
        return hbonds_fname

@handlify(is_method=False, mode='r')
def readChainBounds(fh):
    ''' Utility function to read in residues for each chain and renumber them sequentially
    Returns the following:
    - dictionary of lists of sequential numbering for each chain
    - dictionary mapping (chain, resno, icode)->(chain, sequential number)
    - dictionary mapping (chain, sequential number)->(chain, resno, icode)'''
    tmp_d = defaultdict( innerDefaultdict(list) )
    seqno_d = defaultdict( list )
    seqno2triple = defaultdict(dict)
    triple2seqno = defaultdict(dict)
    # Create temporary dictionary holding all residues by chain, resno and icode
    for line in fh:
        if line:
            chain, resno, icode = line.strip('\n').split('\t')
            resno = int(resno)
            #print chain, resno, icode
            tmp_d[chain][resno].append(icode)
    # Create the sequential dictionaries
    for chain in sorted(tmp_d):
        seq_no = last_resno = None
        tmpc_d = tmp_d[chain]
        for resno in sorted(tmpc_d):
            # If we are just starting
            if seq_no is None:
                seq_no = resno
            # If there is a break in the protein chain 
            elif resno is not None and last_resno is not None and resno>last_resno+1:
                seq_no += (resno-last_resno) - 1
            last_resno = resno
            icodes = sorted(tmpc_d[resno])
            for icode in icodes:
                #print seq_no, chain, resno, icode
                seqno2triple[chain][ seq_no ] = (resno, icode)
                triple2seqno[ chain][ (resno, icode) ] = seq_no
                seq_no += 1
    return seqno2triple, triple2seqno

def renumberBack(chain, seqno, smap):
    ''' Maps a sequential number used by our modules back to the original icode-aware numbering '''
    return smap[chain][seqno]
