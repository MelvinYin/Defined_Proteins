/*
  Module: vdw_data
  Authors: Grzegorz M. Koczyk (2007), Igor Berezovsky (1990-2007)
  *****************************************************************************
  Subroutines for calculating Lennard-Jones based on van der Waals-related
  constants.
 */

#ifndef __VDW_DATA_HPP__
#define __VDW_DATA_HPP__

#define DEBUG(x)

#include <vector>
#include <string>
#include <cmath>
#include "pdb_data.hpp"
#include "vdw_constants.hpp"
using namespace std;

// Min and max values for radius between atoms, for which to compute LJ
const double RMIN = 2.5;
const double RMAX = 5.0;
// Skip value for residue computations (LJ not computed for atoms of {n +/- SKIP_ADJACENT} residues)
const int SKIP_ADJACENT=2;

// Compute Lennard-Jones between two atoms based on their distance and van der Waals data object
inline double atomicLj(double r, const VdwData& vdw1, const VdwData& vdw2)
{
  // Original formula from JBSD paper
  // ckl=(362*(ha[kk])*(ha1[ll]))/((sqrt((double)((ha[kk])/(hn[kk])))+sqrt((double)((ha1[ll])/(hn1[ll])))))
  // rkl=((hr[kk])+(hr1[ll]))/2 
  // rkl6=rkl*rkl*rkl*rkl*rkl*rkl
  // ekl=((-1)*ckl)/(2*rkl6)
  // akl=(-1)*ekl*rkl6*rkl6
  // r6=r2*r2*r2
  // ev=((akl/(r6*r6))-(ckl/r6)) 
  //printf("bckl\n");
  //printf("%s\n", vdw1.m_name.c_str());
  //printf("%f %f %f %f\n", vdw1.m_ha, vdw2.m_ha, vdw1.m_hn, vdw2.m_hn);
  double ckl = (362 * (vdw1.m_ha * vdw2.m_ha) ) / ( sqrt( vdw1.m_ha/vdw1.m_hn) + sqrt(vdw2.m_ha/vdw2.m_hn) );
  //printf("ckl\n");
  double rkl = (vdw1.m_hr + vdw2.m_hr)/2.;
  double rkl6 = rkl*rkl*rkl*rkl*rkl*rkl;
  double ekl = (-ckl)/(2*rkl6);
  double akl = (-ekl*rkl6*rkl6);
  double r6 = r*r*r*r*r*r;
  double ev = (akl/(r6*r6)) - (ckl/r6);
  return ev;
};

inline double atomicLj(Atom& atom1, Atom& atom2, const VdwData& vdw1, const VdwData& vdw2)
{
  // Original formula from JBSD paper
  // ckl=(362*(ha[kk])*(ha1[ll]))/((sqrt((double)((ha[kk])/(hn[kk])))+sqrt((double)((ha1[ll])/(hn1[ll])))))
  // rkl=((hr[kk])+(hr1[ll]))/2 
  // rkl6=rkl*rkl*rkl*rkl*rkl*rkl
  // ekl=((-1)*ckl)/(2*rkl6)
  // akl=(-1)*ekl*rkl6*rkl6
  // r6=r2*r2*r2
  // ev=((akl/(r6*r6))-(ckl/r6)) 
  //printf("bckl\n");
  //printf("%s\n", vdw1.m_name.c_str());
  //printf("%f %f %f %f\n", vdw1.m_ha, vdw2.m_ha, vdw1.m_hn, vdw2.m_hn);
  double ckl = (362 * (vdw1.m_ha * vdw2.m_ha) ) / ( sqrt( vdw1.m_ha/vdw1.m_hn) + sqrt(vdw2.m_ha/vdw2.m_hn) );
  //printf("ckl\n");
  double rkl = (vdw1.m_hr + vdw2.m_hr)/2.;
  double rkl6 = rkl*rkl*rkl*rkl*rkl*rkl;
  double ekl = (-ckl)/(2*rkl6);
  double akl = (-ekl*rkl6*rkl6);
  double r2 = ( ((atom1.m_x - atom2.m_x)*(atom1.m_x - atom2.m_x)) + 
		((atom1.m_y - atom2.m_y)*(atom1.m_y - atom2.m_y)) + 
		((atom1.m_z - atom2.m_z)*(atom1.m_z - atom2.m_z)) );
  double r6 = r2*r2*r2;
  double ev = (akl/(r6*r6)) - (ckl/r6);
  return ev;
};

typedef double lj_value_t;
class LJPart
{
 public:
  unsigned long m_atom1_no;
  unsigned long m_atom2_no;
  lj_value_t m_value;
 public:
  LJPart(unsigned long atom1_no, unsigned long atom2_no, lj_value_t value): m_atom1_no(atom1_no), m_atom2_no(atom2_no), m_value(value) {};
};

typedef vector<LJPart> LJResult;

void outAtom(Atom& atom, string head="")
{
  if ( head.size()>0 )
    {
      printf("%s\t", head.c_str());
    };
  printf("%c\t%d\t%s\t%d\t%s\n", atom.m_chain, atom.m_res_no, atom.m_res_name.c_str(), atom.m_number, atom.m_name.c_str());
};

LJResult ljForAtoms(AtomVector& atomvec)
{
  LJResult result;

  // Associating a VdwData with each corresponding atom (ASSERTION)
  vector<VdwData*> vdwdata;
  vdwdata.reserve(atomvec.size());
  VdwResDict vdwconsts = getVdwConstants();
  for (AtomVector::iterator it1=atomvec.begin(); it1!=atomvec.end(); it1++)
    {
      // If records exist use it, if not use a dummy value
      if ( vdwconsts.find(it1->m_res_name)!=vdwconsts.end() && vdwconsts[it1->m_res_name].find(it1->m_name)!=vdwconsts[it1->m_res_name].end() )
	  vdwdata.push_back( &vdwconsts[it1->m_res_name][it1->m_name] );
      else
	  vdwdata.push_back(NULL);
    };
  for (size_t ui1=0; ui1<atomvec.size(); ui1++)
    {
      // Skip if the atom has no known van der Waals data
      if (vdwdata[ui1]==NULL)
	continue;
      Atom& atom1 = atomvec[ui1];
      DEBUG(outAtom(atom1,"atom1:"););
      for (size_t ui2=ui1+1; ui2<atomvec.size(); ui2++)
	{
	  if (vdwdata[ui2]==NULL)
	    continue;
	  Atom& atom2 = atomvec[ui2];
	  DEBUG(outAtom(atom2,"atom2:"););
	  // Skip condiction - atoms from adjacent/same residues in the same chain
	  if (atom1.m_chain==atom2.m_chain)
	    {
	      // Signed difference between residue positions of two atoms
	      int adj_diff = (atom1.m_res_no - atom2.m_res_no);
	      //DEBUG(printf("Adjacent difference %d with ADJ_DIFF=%d\n", adj_diff, SKIP_ADJACENT););
	      //DEBUG(printf( "-SKIP: %d\n", (int) ((-SKIP_ADJACENT)>adj_diff) ); );
	      if ( (-SKIP_ADJACENT<=adj_diff) && (adj_diff<=SKIP_ADJACENT) )
		  continue;
	    };
	  // Skip condition - too distant residues
	  const double r = atomDistance( atom1, atom2 );
	  DEBUG(printf("Atom pair: %d %d\n", atom1.m_number, atom2.m_number););	  
	  DEBUG(printf("Distance is %f with RMIN=%f and RMAX=%f\n", r, RMIN, RMAX););
	  if ( not ( (RMIN<r) && (r<=RMAX)) )
	    continue;
	  // Compute LJ between the atoms
	  result.push_back( LJPart(ui1, ui2, atomicLj(atom1, atom2, *(vdwdata[ui1]), *(vdwdata[ui2]) ) ) );	  
	};
    };
  return result;
};

#undef DEBUG
#endif
