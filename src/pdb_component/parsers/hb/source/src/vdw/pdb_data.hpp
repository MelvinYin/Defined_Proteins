/*
  Module: pdb_data
  Authors: Grzegorz M. Koczyk (2007), Igor Berezovsky (1990-2007)
  ***************************************************************
  Data structures and subroutines for handling PDB data.
 */

#ifndef PDB_DATA_HPP__
#define PDB_DATA_HPP__
#define IMAX 30000


#include <vector>
#include <string>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <utility>
#include <map>
#include <set>

#define DEBUG(x) 
using namespace std;

class Atom
{
public:
  // Atom name
  unsigned long m_number;
  string m_name; 
  char m_alt_loc;
  string m_res_name;
  char m_chain;
  long m_res_no;
  char m_icode;
  // Space coordinates
  double m_x;
  double m_y;
  double m_z;
public:
  Atom(unsigned long number=0, const string& name="", char alt_loc=' ', const string& res_name = "", char chain= ' ', long res_no=0, char icode=' ', double x=0., double y=0., double z=0.): 
    m_number(number),
    m_name(name),
    m_alt_loc(alt_loc),
    m_res_name(res_name),
    m_chain(chain),
    m_res_no(res_no),
    m_icode(icode),
    m_x(x), m_y(y), m_z(z) {};
};

typedef vector<Atom> AtomVector;


// Compute distance between two atoms
inline double atomDistance(Atom& atom1, Atom& atom2)
{
  DEBUG( printf("Coordinates of the first atom: %f %f %f\n", atom1.m_x, atom1.m_y, atom1.m_z); );
  DEBUG( printf("Cooordinates of the second atom: %f %f %f\n", atom2.m_x, atom2.m_y, atom2.m_z); );
  return sqrt( ((atom1.m_x - atom2.m_x)*(atom1.m_x - atom2.m_x)) + 
	       ((atom1.m_y - atom2.m_y)*(atom1.m_y - atom2.m_y)) + 
	       ((atom1.m_z - atom2.m_z)*(atom1.m_z - atom2.m_z)) );
};

typedef map<string, Atom> AtomMap;
typedef AtomMap::iterator AtomMapIter;

// Amino-acid residues
class Residue
{
public:
  string m_name;
  unsigned long m_number;
  char m_chain;
  char m_icode;
  AtomMap m_atoms;
  long m_previous;
public:
  Residue(const string& name, char chain, char icode, unsigned long number, AtomMap atoms, long previous=-IMAX): m_name(name), m_chain(chain), m_icode(icode), m_number(number), m_previous(previous), m_atoms(atoms) {};
};

typedef vector<Residue> ResVector;
typedef ResVector::iterator ResVecIter;

/*
COLUMNS      DATA TYPE        FIELD      DEFINITION
------------------------------------------------------
 1 -  6      Record name      "ATOM    "
 7 - 11      Integer          serial     Atom serial number.
13 - 16      Atom             name       Atom name.
17           Character        altLoc     Alternate location indicator.
18 - 20      Residue name     resName    Residue name.
22           Character        chainID    Chain identifier.
23 - 26      Integer          resSeq     Residue sequence number.
27           AChar            iCode      Code for insertion of residues.
31 - 38      Real(8.3)        x          Orthogonal coordinates for X in 
                                         Angstroms
39 - 46      Real(8.3)        y          Orthogonal coordinates for Y in 
                                         Angstroms
47 - 54      Real(8.3)        z          Orthogonal coordinates for Z in 
                                         Angstroms
55 - 60      Real(6.2)        occupancy  Occupancy.
61 - 66      Real(6.2)        tempFactor Temperature factor.
77 - 78      LString(2)       element    Element symbol, right-justified.
79 - 80      LString(2)       charge     Charge on the atom.
 */

string removeSpaces(const string& s)
{
  string r("");
  for (size_t ui=0; ui<s.size(); ui++)
    if ( s[ui]!=' ' )
      r+=s[ui];
  return r;	  
};
inline AtomVector parseAtomLines(FILE* fh)
{
  char line[200];
  AtomVector result;
  while ( fgets (line, 200, fh ) != NULL)
    {
      string sline(line);
      if (sline.substr(0,4)!="ATOM")
	continue;	
      // Corrected: for >100000 atoms
      //          if len(line)==82:
      //              try:
      //                  serial_number = int(line[6:12])
      //                  line=line[0:12]+line[13:]
      //              except:
      //                  serial_number = 0
      //          else:
      //              try:
      //                  serial_number=int(line[6:11])
      //              except:
      //                  serial_number=0 
      short delta = sline.size()-81;
      if (delta<0) delta=0;	
      //printf("%d\n", sline.size());
      long atom_number = atol( sline.substr(6,5+delta).c_str() );
      string atom_name = removeSpaces(sline.substr(12+delta,4));
      char alt_loc = sline[16+delta];
      string res_name = removeSpaces(sline.substr(17+delta,3));
      char chain = sline[21+delta];
      long res_number = atol( sline.substr(22+delta,4).c_str() );
      char icode = sline[26+delta];
      double x = atof( sline.substr(30+delta, 8).c_str() );
      double y = atof( sline.substr(38+delta, 8).c_str() );
      double z = atof( sline.substr(46+delta, 8).c_str() );
      //printf("%d %d\n", atom_number, res_number );
      //printf("%s\n", sline.substr(6,4).c_str());
      result.push_back( Atom(atom_number, atom_name, alt_loc, res_name, chain, res_number, icode, x, y, z) );
    };
  return result;
};


bool atomCmp(Atom atom1, Atom atom2)
{
  if (atom1.m_chain<atom2.m_chain)
    return true;
  else if (atom1.m_chain==atom2.m_chain)
      if (atom1.m_res_no<atom2.m_res_no)
	return true;
      else if (atom1.m_res_no==atom2.m_res_no)
	  if (atom1.m_icode<atom2.m_icode)
	    return true;
  return false;
};

// Residue identifier class including information about insertion codes
class ResId
{ 
public:
  char m_icode;
  long m_resno;
public:
  ResId(const long resno, const char icode): m_resno(resno), m_icode(icode) {};
  bool operator < (const ResId& other) const 
  { 
    if (m_resno<other.m_resno)
      return true;
    else 
      if (m_resno==other.m_resno && m_icode<other.m_icode)
	return true;
    return false;
    //return (m_resno<other.m_resno && m_icode<other.m_icode); 
  };
  bool operator == (const ResId& other) const 
  { 
    return (m_resno==other.m_resno && m_icode==other.m_icode );
  };
};

typedef set<ResId> resnoset_t;
typedef map< char, resnoset_t > ChainBounds;

ChainBounds makeChainBounds(AtomVector& atomvec)
{
  ChainBounds result;
  for (AtomVector::iterator it=atomvec.begin(); it!=atomvec.end(); it++)
    {
      result[it->m_chain].insert( ResId(it->m_res_no, it->m_icode) );
    };
  return result;
};

#endif
