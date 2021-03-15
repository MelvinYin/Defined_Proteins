// Program to propose hydrogen bonds between residues in a given protein structure
/*
  Module: hbonds
  Authors: Grzegorz M. Koczyk (2008)
  ***************************************************************
  Data structures and subroutines for handling PDB data.
 */


#include <vector>
#include <string>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <utility>
#include <map>
#include <set>
#include <algorithm>
#include <math.h>
#include "pdb_data.hpp"
#define DEBUG(x) 

using namespace std;


#include <string>
#include <map>

// # Donors according to (Stickle et al. 1992)
// donors = {
//     'ALA': atomdikt( peptide_donor ),
//     'ARG': atomdikt( peptide_donor,
//                      AtomDonor(atom_name='NE', dd='CD', dd1='CZ', hybridization='Nsp2', radius=1.9, max_bonds=1),
//                      AtomDonor(atom_name='NH1', dd='CZ', dd1='NE', hybridization='Nsp2', radius=1.9, max_bonds=2),
//                      AtomDonor(atom_name='NH2', dd='CZ', dd1='NE', hybridization='Nsp2', radius=1.9, max_bonds=2)
//                      ),
//     'ASP': atomdikt( peptide_donor ),
//     'ASN': atomdikt( peptide_donor,
//                      AtomDonor(atom_name='ND2', dd='CG', dd1='CB', hybridization='Nsp2', radius=1.9, max_bonds=2),
//                      ),
//     'CYS': atomdikt( peptide_donor ),
//     'GLU': atomdikt( peptide_donor ),
//     'GLN': atomdikt( peptide_donor,
//                      AtomDonor(atom_name='NE2', dd='CD', dd1='CG', hybridization='Nsp2', radius=1.9, max_bonds=2),
//                      ),
//     'GLY': atomdikt( peptide_donor ),
//     'HIS': atomdikt( peptide_donor,
//                      AtomDonor(atom_name='ND1', dd='CG', dd1='CE1', hybridization='Nsp2', radius=1.9, max_bonds=1, excludes='NE2'),
//                      AtomDonor(atom_name='NE2', dd='CE1', dd1='CD2', hybridization='Nsp2', radius=1.9, max_bonds=1, excludes='ND1'),
//                      ),
//     'ILE': atomdikt( peptide_donor ),
//     'LEU': atomdikt( peptide_donor ),
//     'LYS': atomdikt( peptide_donor,
//                      AtomDonor(atom_name='NZ', dd='CE', dd1='CD', hybridization='Nsp3', radius=2.1, max_bonds=3),
//                      ),
//     'MET': atomdikt( peptide_donor ),
//     'PHE': atomdikt( peptide_donor ),
//     'PRO': atomdikt( ),
//     'SER': atomdikt( peptide_donor,
//                      AtomDonor(atom_name='OG', dd='CB', dd1='CA', hybridization='Osp3', radius=1.7, max_bonds=1),
//                      ),
//     'THR': atomdikt( peptide_donor,
//                      # TODO: Check radius=1.9 with Igor - whether it is on purpose (article says r=1.7)
//                      AtomDonor(atom_name='OG1', dd='CB', dd1='CA', hybridization='Osp3', radius=1.7, max_bonds=1),
//                      ),
//     'TRP': atomdikt( peptide_donor,
//                      AtomDonor(atom_name='NE1', dd='CD1', dd1='CE2', hybridization='Nsp2', radius=1.9, max_bonds=1),
//                     ) ,
//     'TYR': atomdikt( peptide_donor,
//                      AtomDonor(atom_name='OH', dd='CZ', dd1='CE1', hybridization='Osp3', radius=1.7, max_bonds=1),
//                      ),
//     'VAL': atomdikt( peptide_donor ),
//     }


class DonorData
{
public:
  // Atom data consists of the following data:
  //    atom name
  //__init__(self, atom_name, aa, aa2, hybridization, radius, max_bonds = 1, excludes=None): 
  string m_atom_name;
  //    dd    - antecedent of donor (as per Stickle et al)
  string m_dd;
  //    dd1    - second antecedent of donor (as per Stickle et al)
  string m_dd1;
  //    hybridization - hybridization of donor
  string m_hybridization;
  //    radius  - radius of influence for calculations
  double m_radius;
  //    max_hbonds - maximum number of hbonds
  int m_max_hbonds;
  // m_hybrno - number at the end of hybridization string
  short m_hybrno;
public:
  // Default constructor
  DonorData(string atom_name="", string dd="", string dd1="", string hybridization="", double radius=0., int max_hbonds=32000): m_atom_name(atom_name), m_dd(dd), m_dd1(dd1), m_hybridization(hybridization), m_radius(radius), 
																m_max_hbonds(max_hbonds)
  { if (m_hybridization.size()>0)
      this->m_hybrno = atoi(  m_hybridization.substr(m_hybridization.size()-1).c_str() ); 
  };
};

// Dictionary of HbsData objects for each atom
typedef map<string, DonorData> DonorAtomDict;
// Dictionary of VdwAtomDict objects representing residues
typedef map<string, DonorAtomDict> DonorResDict;

DonorResDict getHbsDonors()
{
  DonorResDict donors;
  DonorData peptide_donor("N", "CA", "C_prev", "Nsp2",1.9, 1);
  // ALA
  donors["ALA"]["N"] = peptide_donor;
  // ARG
  donors["ARG"]["N"] = peptide_donor;
  donors["ARG"]["NE"] = DonorData("NE", "CD", "CZ", "Nsp2", 1.9, 1);
  donors["ARG"]["NH1"] = DonorData("NH1", "CZ", "NE", "Nsp2", 1.9, 2);
  donors["ARG"]["NH2"] = DonorData("NH2", "CZ", "NE", "Nsp2", 1.9, 2);
  // ASP
  donors["ASP"]["N"] = peptide_donor;
  // ASN
  donors["ASN"]["N"] = peptide_donor;
  donors["ASN"]["ND2"] = DonorData("ND2", "CG", "CB", "Nsp2", 1.9, 2);
  // CYS
  donors["CYS"]["N"] = peptide_donor;
  // GLU
  donors["GLU"]["N"] = peptide_donor;
  // GLN
  donors["GLN"]["N"] = peptide_donor;
  donors["GLN"]["NE2"] = DonorData("NE2", "CD", "CG", "Nsp2", 1.9, 2);
  // GLY
  donors["GLY"]["N"] = peptide_donor;
  // HIS
  donors["HIS"]["N"] = peptide_donor;
  donors["HIS"]["ND1"] = DonorData("ND1", "CG", "CE1", "Nsp2", 1.9, 1);
  donors["HIS"]["NE2"] = DonorData("NE2", "CE1", "CD2", "Nsp2", 1.9, 1);
  // ILE
  donors["ILE"]["N"] = peptide_donor;
  // LEU
  donors["LEU"]["N"] = peptide_donor;
  // LYS
  donors["LYS"]["N"] = peptide_donor;
  donors["LYS"]["NZ"] = DonorData("NZ", "CE", "CD", "Nsp3", 2.1, 3);
  // MET
  donors["MET"]["N"] = peptide_donor;
  // PHE
  donors["PHE"]["N"] = peptide_donor;
  // PRO 
  // ...IS NOT A DONOR AT ALL !!!
  // SER
  donors["SER"]["N"] = peptide_donor;
  donors["SER"]["OG"] = DonorData("OG", "CB", "CA", "Osp3", 1.7, 1);
  // THR
  donors["THR"]["N"] = peptide_donor;
  donors["THR"]["OG1"] = DonorData("OG1", "CB", "CA", "Osp3", 1.7, 1);
  // TRP
  donors["TRP"]["N"] = peptide_donor;
  donors["TRP"]["NE1"] = DonorData("NE1", "CD1", "CE2", "Nsp2", 1.9, 3);
  // TYR
  donors["TYR"]["N"] = peptide_donor;
  donors["TYR"]["OH"] = DonorData("OH", "CZ", "CE1", "Osp3", 1.7, 1);
  // VAL
  donors["VAL"]["N"] = peptide_donor;
  return donors;
};

// acceptors = {
//     'ALA': atomdikt( peptide_acceptor ),
//     'ARG': atomdikt( peptide_acceptor ),
//     'ASP': atomdikt( peptide_acceptor,
//                      AtomAcceptor(atom_name='OD1', aa='CG', aa2='CB', hybridization='Osp2', radius=1.6, max_bonds=1),
//                      AtomAcceptor(atom_name='OD2', aa='CG', aa2='CB', hybridization='Osp2', radius=1.6, max_bonds=1),
//                      ),
//     'ASN': atomdikt( peptide_acceptor,
//                      AtomAcceptor(atom_name='OD1', aa='CG', aa2='CB', hybridization='Osp2', radius=1.6, max_bonds=2),
//                      ),
//     'CYS': atomdikt( peptide_acceptor,
//                      AtomAcceptor(atom_name='SG', aa='CB', aa2='CA', hybridization='Ssp3', radius=2.1, max_bonds=2),
//                      ),
//     'GLU': atomdikt( peptide_acceptor,
//                      AtomAcceptor(atom_name='OE1', aa='CD', aa2='CG', hybridization='Osp2', radius=1.6, max_bonds=2),
//                      AtomAcceptor(atom_name='OE2', aa='CD', aa2='CG', hybridization='Osp2', radius=1.6, max_bonds=2),
//                      ),
//     'GLN': atomdikt( peptide_acceptor,
//                      AtomAcceptor(atom_name='OE1', aa='CD', aa2='CG', hybridization='Nsp2', radius=1.6, max_bonds=2),
//                      ),
//     'GLY': atomdikt( peptide_acceptor ),
//     'HIS': atomdikt( peptide_acceptor,
//                      AtomAcceptor(atom_name='ND1', aa='CG', aa2='CE1', hybridization='Nsp2', radius=1.6, max_bonds=1, excludes='NE2'),
//                      AtomAcceptor(atom_name='NE2', aa='CE1', aa2='CD2', hybridization='Nsp2', radius=1.6, max_bonds=1, excludes='ND1'),
//                      ),
//     'ILE': atomdikt( peptide_acceptor ),
//     'LEU': atomdikt( peptide_acceptor ),
//     'LYS': atomdikt( peptide_acceptor ),
//     'MET': atomdikt( peptide_acceptor,
//                      AtomAcceptor(atom_name='SD', aa='CG', aa2='CB', hybridization='Ssp3', radius=1.95, max_bonds=2),
//                      ),
//     'PHE': atomdikt( peptide_acceptor ),
//     'PRO': atomdikt( peptide_acceptor ),
//     'SER': atomdikt( peptide_acceptor,
//                      AtomAcceptor(atom_name='OG', aa='CB', aa2='CA', hybridization='Osp3', radius=1.7, max_bonds=1),
//                      ),
//     'THR': atomdikt( peptide_acceptor,
//                      AtomAcceptor(atom_name='OG1', aa='CB', aa2='CA', hybridization='Osp3', radius=1.7, max_bonds=1),
//                      ),
//     'TRP': atomdikt( peptide_acceptor ),
//     'TYR': atomdikt( peptide_acceptor,
//                      AtomAcceptor(atom_name='OH', aa='CZ', aa2='CE1', hybridization='Osp3', radius=1.7, max_bonds=1),
//                      ),
//     'VAL': atomdikt( peptide_acceptor ),
//     } 
class AcceptorData
{
public:
  // Atom data consists of the following data:
  //    atom name
  //__init__(self, atom_name, aa, aa2, hybridization, radius, max_bonds = 1, excludes=None): 
  string m_atom_name;
  //    aa    - antecedent of acceptor (as per Stickle et al)
  string m_aa;
  //    aa2    - second antecedent of acceptor (as per Stickle et al)
  string m_aa2;
  //    hybridization - hybridization of donor
  string m_hybridization;
  //    radius  - radius of influence for calculations
  double m_radius;
  //    max_hbonds - maximum number of hbonds
  int m_max_hbonds;
   // m_hybrno - number at the end of hybridization string
  short m_hybrno;
public:
  // Default constructor
  AcceptorData(string atom_name="", string aa="", string aa2="", string hybridization="", double radius=0., int max_hbonds=32000): m_atom_name(atom_name), m_aa(aa), m_aa2(aa2), m_hybridization(hybridization), m_radius(radius), 
																   m_max_hbonds(max_hbonds) 
  { if (m_hybridization.size()>0)
      this->m_hybrno = atoi(  m_hybridization.substr(m_hybridization.size()-1).c_str() ); 
  };
};

// Dictionary of HbsData objects for each atom
typedef map<string, AcceptorData> AcceptorAtomDict;
// Dictionary of VdwAtomDict objects representing residues
typedef map<string, AcceptorAtomDict> AcceptorResDict;

AcceptorResDict getHbsAcceptors()
{
  AcceptorResDict acceptors;
 AcceptorData  peptide_acceptor("O", "C", "CA", "Osp2", 1.6, 1);
  // ALA
  acceptors["ALA"]["O"] = peptide_acceptor;
  // ARG
  acceptors["ARG"]["O"] = peptide_acceptor;
  // ASP
  acceptors["ASP"]["O"] = peptide_acceptor;
  acceptors["ASP"]["OD1"] = AcceptorData("OD1", "CG", "CB", "Osp2", 1.6, 2);
  acceptors["ASP"]["OD2"] = AcceptorData("OD2", "CG", "CB", "Osp2", 1.6, 2);
  // ASN
  acceptors["ASN"]["O"] = peptide_acceptor;
  acceptors["ASN"]["OD1"] = AcceptorData("OD1", "CG", "CB", "Osp2", 1.6, 2);
  // CYS
  acceptors["CYS"]["O"] = peptide_acceptor;
  acceptors["CYS"]["SG"] = AcceptorData("OD1", "CB", "CA", "Ssp3", 2.1, 2);
  // GLU
  acceptors["GLU"]["O"] = peptide_acceptor;
  acceptors["GLU"]["OE1"] = AcceptorData("OE1", "CD", "CG", "Osp2", 1.6, 2);
  acceptors["GLU"]["OE2"] = AcceptorData("OE2", "CD", "CG", "Osp2", 1.6, 2);
  // GLN
  acceptors["GLN"]["O"] = peptide_acceptor;
  acceptors["GLN"]["OE1"] = AcceptorData("OE1", "CD", "CG", "Osp2", 1.6, 2);
  // GLY
  acceptors["GLY"]["O"] = peptide_acceptor;
  // HIS
  acceptors["HIS"]["O"] = peptide_acceptor;
  acceptors["HIS"]["ND1"] = AcceptorData("ND1", "CG", "CE1", "Nsp2", 1.6, 1);
  acceptors["HIS"]["NE2"] = AcceptorData("ND2", "CE1", "CD2", "Nsp2", 1.6, 1);
  // ILE
  acceptors["ILE"]["O"] = peptide_acceptor;
  // LEU
  acceptors["LEU"]["O"] = peptide_acceptor;
  // LYS
  acceptors["LYS"]["O"] = peptide_acceptor;
  // MET
  acceptors["MET"]["O"] = peptide_acceptor;
  acceptors["MET"]["SD"] = AcceptorData("SD", "CG", "CB", "Ssp3", 1.95, 2);
  // PHE
  acceptors["PHE"]["O"] = peptide_acceptor;
  // PRO
  acceptors["PRO"]["O"] = peptide_acceptor;
  // SER
  acceptors["SER"]["O"] = peptide_acceptor;
  acceptors["SER"]["OG"] = AcceptorData("OG", "CB", "CA", "Osp3", 1.7, 1);
  // THR
  acceptors["THR"]["O"] = peptide_acceptor;
  acceptors["THR"]["OG1"] = AcceptorData("OG1", "CB", "CA", "Osp3", 1.7, 1);
  // TRP
  acceptors["TRP"]["O"] = peptide_acceptor;
  // TYR
  acceptors["TYR"]["O"] = peptide_acceptor;
  acceptors["TYR"]["OH"] = AcceptorData("OH", "CZ", "CE1", "Osp3", 1.7, 1);
  // VAL
  acceptors["VAL"]["O"] = peptide_acceptor;
  return acceptors;
};


// Define shorthand typedefs
typedef DonorResDict::iterator DonorIt;
typedef AcceptorResDict::iterator AcceptorIt;

class HBond
{
public:
  // Donor
  string m_donor_residue;
  string m_donor_atom;
  char m_donor_icode;
  char m_donor_chain;
  long m_donor_rnumber;
  // Acceptor
  string m_acceptor_residue;
  string m_acceptor_atom;
  char m_acceptor_icode;
  char m_acceptor_chain;
  long m_acceptor_rnumber;
  // Bond characteristics
  double m_distance;
  double m_donor_angle;
  double m_acceptor_angle;
  double m_planar1;
  double m_planar2;
public:
  HBond( const char donor_chain=' ', const string donor_residue="", const long donor_rnumber =0, 
	 const string donor_atom = "", const char donor_icode=' ', 
	 const char acceptor_chain=' ', const string acceptor_residue="", const long acceptor_rnumber =0,
	 const string acceptor_atom="", const char acceptor_icode=' ', 
	 const double distance=0., const double donor_angle=0., const double acceptor_angle=0., 
	 const double planar1=0., const double planar2=0.): m_donor_chain(donor_chain), m_donor_residue(donor_residue), m_donor_rnumber(donor_rnumber), m_donor_atom(donor_atom), m_donor_icode(donor_icode),
							    m_acceptor_chain(acceptor_chain),m_acceptor_residue(acceptor_residue), m_acceptor_rnumber(acceptor_rnumber), m_acceptor_atom(acceptor_atom), m_acceptor_icode(acceptor_icode),
							    m_distance(distance), m_donor_angle(donor_angle), m_acceptor_angle(acceptor_angle),
							    m_planar1(planar1), m_planar2(planar2)
  {};
};

typedef vector<HBond> HbVec;

double norm3d(double x, double y, double z)
{
  return sqrt( x*x+y*y+z*z);
};

// TODO: Unary minus for this one !
class Vec3d
{
public:
  const double m_x;
  const double m_y;
  const double m_z;
public:
  Vec3d(double x, double y, double z): m_x(x), m_y(y), m_z(z) {};
  double norm() const { return norm3d(this->m_x, this->m_y, this->m_z); };
  Vec3d operator - () const { return Vec3d(-this->m_x, -this->m_y, -this->m_z); };
};

Vec3d cross(double p1, double p2, double p3, double q1, double q2, double q3)
{
  return Vec3d( p2*q3-q2*p3, p3*q1-q3*p1, p1*q2-p2*q1);
  //return Vec3d( y1*z2 - z1*y2, z1*x2-x1*z2, x1*y2-y1*x2 );
};

Vec3d cross(const Vec3d& v1, const Vec3d& v2)
{
  return cross( v1.m_x, v1.m_y, v1.m_z, v2.m_x, v2.m_y, v2.m_z);
};

double dot(double x1, double y1, double z1, double x2, double y2, double z2 )
{
  return x1*x2+y1*y2+z1*z2;
};

double dot(const Vec3d& v1, const Vec3d& v2)
{
  return dot( v1.m_x, v1.m_y, v1.m_z, v2.m_x, v2.m_y, v2.m_z);
};

int sgn(double v)
{
  if (v>=0)
    return 1;
  return -1;
};

double todeg(double v)
{
  return (v*180)/M_PI;
};

double angle(const Vec3d& v1, const Vec3d& v2, bool absolute=true)
{
  double m = dot(v1,v2)/(v1.norm()*v2.norm());
  if (m>1.)
    m = 1.;
  if (m<-1.)
    m = -1.;
  double va = acos(m);
  if (va<0)
    va = M_PI - va;
  return todeg(va);
};

Vec3d atoms2vec(Atom& atom1, Atom& atom2)
{
  return Vec3d(atom1.m_x-atom2.m_x, atom1.m_y-atom2.m_y, atom1.m_z-atom2.m_z);
};

void checkHbsPair(HbVec& result, Residue& donor, Residue& acceptor, DonorAtomDict& donor_data, AcceptorAtomDict& acceptor_data, ResVector& resvec)
{
  // For all donor atoms in donor/all acceptor atoms in acceptor
  for (DonorAtomDict::iterator itdi=donor_data.begin(); itdi!=donor_data.end(); itdi++)
    {
      // Skipping if the atom is not represented
      AtomMap::iterator itdr = donor.m_atoms.find(itdi->first);
      if (itdr==donor.m_atoms.end())
	continue;
      
      DonorData& donor_data = itdi->second;
      Atom& d = itdr->second;
      //////////////////////////////////////////////////////////////////////////////////////////////
      // MAIN COURSE :)
      // ATOM RETRIEVAL
      // D
      itdr = donor.m_atoms.find(donor_data.m_atom_name);
      if (itdr == donor.m_atoms.end() )
	continue;
      // DD1
      if (donor_data.m_dd1=="C_prev")	  
	{
	  if (donor.m_previous==-IMAX)
	    continue;
	  else	      
	    {
	      //printf("Siamthing %d%s\n", donor.m_number, donor.m_name.c_str() );
	      itdr = resvec[donor.m_previous].m_atoms.find("C");	    
	      //printf("Afters\n");
	      if (itdr == resvec[donor.m_previous].m_atoms.end() )
		continue;	 	     
	    };
	}
      else
	{
	  itdr = donor.m_atoms.find(donor_data.m_dd1);
	  if (itdr == donor.m_atoms.end() )
	    continue;
	};
      Atom& dd1 = itdr->second;
      // DD
      itdr = donor.m_atoms.find(donor_data.m_dd);
      if (itdr == donor.m_atoms.end() )
	continue;
      Atom& dd = itdr->second;
      /*-------------------------------------*/
      // DEBUG block - output donor data
      //printf("Donor residue: %s%d\t", donor.m_name.c_str(), donor.m_number);
      //printf("Donor atom:%s\n", itdi->second.m_atom_name.c_str());

      for (AcceptorAtomDict::iterator itai=acceptor_data.begin(); itai!=acceptor_data.end(); itai++)
      {
	// Skipping if the atom is not represented
	AtomMap::iterator itar = acceptor.m_atoms.find(itai->first);
	AcceptorData& acceptor_data = itai->second;
	Atom& a = itar->second;
	if (itar==acceptor.m_atoms.end())
	  continue;
	/*-------------------------------------*/
	// DEBUG block - output donor data
	//printf("\nDonor residue: %s%d\t", donor.m_name.c_str(), donor.m_number);
	//printf("Donor atom:%s\n", itdi->second.m_atom_name.c_str());
	//printf("Acceptor residue: %s%d\t", acceptor.m_name.c_str(), acceptor.m_number);
	//printf("Acceptor atom:%s\n", itai->second.m_atom_name.c_str());
	
	// AA
	itar = acceptor.m_atoms.find(acceptor_data.m_aa);
	if (itar == acceptor.m_atoms.end() )
	  {};//continue;
	Atom& aa = itar->second;
	// AA2
	itar = acceptor.m_atoms.find(acceptor_data.m_aa2);
	if (itar == acceptor.m_atoms.end() )
	  continue;
	Atom& aa2 = itar->second;
	//----------------------------------------------------------------------------------------
	// CHECKS
	// Check distance constraint
	double dist = atomDistance(d, a);
	//printf ("Distance: %.2f\n", dist );
	if (dist>donor_data.m_radius+acceptor_data.m_radius || dist==0.)
	  continue;
	// Check angles
	// 1) D-A-AA
	double d_a_aa = angle( atoms2vec(aa,a), atoms2vec(d,a) );
	//printf ("d-a-aa: %.2f\n", d_a_aa );
	if (acceptor_data.m_hybrno==3)
	  {
	    if (d_a_aa<=60.) 
	      continue;
	  }
	else
	  if (d_a_aa<=90.)
	    continue;
	// 2) A-D-DD 
	double a_d_dd = angle( atoms2vec(a,d), atoms2vec(dd,d) );
	//printf ("a-d-dd: %.2f\n", a_d_dd );
	//printf("A=%s,D=%s,DD=%s\t", a.m_name.c_str(), d.m_name.c_str(), dd1.m_name.c_str() );
	//printf("A-D-DD: %.2f\n", a_d_dd);
	if ( donor.m_name=="LYS" && d.m_name!="N" )
	  {
	    if (a_d_dd<= 60.)
	      continue;
	  }
	else
	  if  (a_d_dd<= 90.)
	    continue;

	double planar1 = angle( cross( atoms2vec(d, dd), atoms2vec(dd1, dd) ), -cross( atoms2vec(a,d), atoms2vec(dd,d) )  );
	if ( donor_data.m_hybrno==2 )
	  {
	    if (donor.m_name=="LYS" && d.m_name!="N")
	      {
		if (planar1>90. && -planar1<=-90)
		  continue;
	      }
	    else
	      if (planar1>=60. && -planar1>=-120. )
		continue;
	  };

	double planar2 = angle( cross( atoms2vec(a, aa), atoms2vec(aa2, aa) ), cross( atoms2vec(d,a), atoms2vec(aa,a) )  );
	if ( (acceptor_data.m_hybrno==2) && (planar2>=90. && -planar2 >= -90.) )
	  continue;
	HBond hb(donor.m_chain, donor.m_name, donor.m_number, d.m_name, d.m_icode, acceptor.m_chain, acceptor.m_name, acceptor.m_number, a.m_name, a.m_icode, dist, d_a_aa, a_d_dd, planar1, planar2);
	result.push_back(hb);
      };
    };
};

//! Calculate hydrogen bonds 
HbVec doHydrogenBonds(AtomVector& atomvec)
{
  // Setup donor and acceptor constant tables
  DonorResDict donors = getHbsDonors();
  AcceptorResDict acceptors = getHbsAcceptors();
  HbVec result;
  ResVector residues;
  // Setup the atoms as residues  
  // Go through atoms sequentially
  // each time a new residue is found, or end reached - put all atoms in the stack in the residue
  // mark appropriately information about previous residue etc.
  // Sort the vector according to first the chain. then residue number then the icode
  std::sort(atomvec.begin(), atomvec.end(), atomCmp);
  // Group the atoms into residues
  AtomMap ratoms;
  int prev_resno = -IMAX;
  int resno = prev_resno;
  //Residue residue();
  for (AtomVector::iterator itat=atomvec.begin(); itat!=atomvec.end(); itat++)
    {
      Atom& atom = *itat;
      resno = atom.m_res_no;
      // If this is the first iteration or 
      if (resno!=prev_resno)
	{
	  ratoms.clear();
	  //printf("Resno %d\n", atom.m_res_no);
	  Residue residue(atom.m_res_name, atom.m_chain, atom.m_icode, atom.m_res_no, ratoms, -IMAX);
	  //printf("Resno %d\n", residue.m_number);
	  if (resno==prev_resno+1 && residues.size()>0)	    	    
	    {
	      residue.m_previous = residues.size()-1;
	    };	  
	  residues.push_back(residue);
	  prev_resno = resno;
	};
      residues.rbegin()->m_atoms[atom.m_name]=atom;
    };

  // Launch a one vs one comparison for all residue pairs (as donors/as acceptors)
  for (ResVecIter itd=residues.begin(); itd!=residues.end(); itd++)
    {
      Residue& donor = *itd;
      DonorIt itd_desc = donors.find(itd->m_name);
      if (itd_desc==donors.end())
	continue;
      for (ResVecIter ita=residues.begin(); ita!=residues.end(); ita++)
	{	
	  //if (itd==ita)
	  //  continue;
	  Residue& acceptor = *ita;
	  AcceptorIt ita_acc = acceptors.find(ita->m_name);
	  if (ita_acc==acceptors.end())
	    continue;
	  // If we have both - run checkHbsPair
	  checkHbsPair(result, donor, acceptor, itd_desc->second, ita_acc->second, residues);
	};
    };
  return result;
};

// For the openFile function
#include "basic_formatting.hpp"

bool writeHbs(FILE* fh, HbVec& hbs)
{
  for (HbVec::iterator it=hbs.begin(); it!=hbs.end(); it++)
    {
      fprintf(fh, "query\t");
      fprintf(fh, "%c\t%s\t%d\t%c\t%s\t",it->m_donor_chain, it->m_donor_residue.c_str(), it->m_donor_rnumber, it->m_donor_icode, it->m_donor_atom.c_str());
      fprintf(fh, "%c\t%s\t%d\t%c\t%s\t", it->m_acceptor_chain, it->m_acceptor_residue.c_str(), it->m_acceptor_rnumber,  it->m_acceptor_icode, it->m_acceptor_atom.c_str());
      fprintf(fh, "%.2f\t%.2f\t%.2f\t%.2f\t%.2f", it->m_distance,  it->m_acceptor_angle, it->m_donor_angle,it->m_planar1, it->m_planar2);
      char d='R';
      char a='R';
      if (it->m_donor_atom=="N")
	d='N';
      if (it->m_acceptor_atom=="O")
	a='O';
      fprintf(fh, "\t%c_%c", d, a);
      fprintf(fh, "\n");
    };
  return true;
};

bool writeHbs(const string& fname, HbVec& hbs)
{
  FILE* fh = openFile(fname.c_str(), "wt");
  if (!fh)
    return false;
  bool r = writeHbs(fh, hbs);
  fclose(fh);
  return r;
}; 

int main(int argc, char* argv[])
{

   if (argc<3)
    {
      printf("You must supply a filename to compute on and result hbonds filename.");
      return -1;
    };
   // Read the structure in
   string input_fname(argv[1]);
   string hbs_fname(argv[2]);
   FILE* fh = openFile(argv[1], "r");  
   if (! fh)
     return -1;
   AtomVector vec = parseAtomLines(fh);
   fclose(fh);
   // Compute hbonds on it
   HbVec hbs = doHydrogenBonds(vec);
   writeHbs(hbs_fname, hbs);
   return 0;
}; 

