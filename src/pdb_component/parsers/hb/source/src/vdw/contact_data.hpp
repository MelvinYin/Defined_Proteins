/*
  Module: contact_data
  Authors: Grzegorz M. Koczyk (2007), Igor Berezovsky (1990-2007)
  *****************************************************************************
  Subroutines for calculating contact data for defining loops& locks.
 */

#ifndef __CONTACT_DATA_HPP__
#define __CONTACT_DATA_HPP__

#define DEBUG(x)

#include <vector>
#include <string>
#include <cmath>
#include "pdb_data.hpp"
using namespace std;

typedef double dist_value_t;
class ContactPart
{
 public:
  unsigned long m_atom1_no;
  unsigned long m_atom2_no;
  dist_value_t m_distance;
 public:
  ContactPart(unsigned long atom1_no, unsigned long atom2_no, dist_value_t value): m_atom1_no(atom1_no), m_atom2_no(atom2_no), m_distance(value) {};
};

typedef vector<ContactPart> ContactResult;
// Constants for computation (warning: same names but a bit different from van der Waals part of the program)
// Min and max values for radius between atoms, for which to give contacts
const double RMIN = 2.5;
const double RMAX = 10.0;
// Skip value for residue computations (distance not computed for atoms of {n +/- SKIP_ADJACENT} residues)
const int SKIP_ADJACENT=4;

string ca_string("CA");
ContactResult contactForAtoms(AtomVector& atomvec)
{
  ContactResult result;
  // Filtering out atoms except CA
  vector<bool> validvec;
  validvec.reserve(atomvec.size());
  for (AtomVector::iterator it1=atomvec.begin(); it1!=atomvec.end(); it1++)
    {
      // If records exist use it, if not use a dummy value
      if ( it1->m_name==ca_string)	
	validvec.push_back(true);
      else
	validvec.push_back(false);
    };
  // Loop over the atoms, computing contacts and distance
  for (size_t ui1=0; ui1<atomvec.size(); ui1++)
    {
      Atom& atom1 = atomvec[ui1];
      if (!validvec[ui1])
	continue;
	  for (size_t ui2=ui1+1; ui2<atomvec.size(); ui2++)
	    {
	      Atom& atom2 = atomvec[ui2];
	      if (! validvec[ui2])
		continue;
	      // Skip condiction - atoms from adjacent/same residues in the same chain
	      if (atom1.m_chain==atom2.m_chain)
		{
		  // Signed difference between residue positions of two atoms
		  int adj_diff = (atom1.m_res_no - atom2.m_res_no);
		  if ( (-SKIP_ADJACENT<=adj_diff) && (adj_diff<=SKIP_ADJACENT) )
		    continue;
		};
	  // Skip condition - too distant residues
	  const double r = atomDistance( atom1, atom2 );
	  DEBUG(printf("Atom pair: %d %d\n", atom1.m_number, atom2.m_number););	  
	  DEBUG(printf("Distance is %f with RMIN=%f and RMAX=%f\n", r, RMIN, RMAX););
	  if ( not ( (RMIN<r) && (r<=RMAX)) )
	    continue;
	  // Push the resulting contact
	  result.push_back( ContactPart(ui1, ui2, r) );	  
	};
    };
  return result;
};

#undef DEBUG
#endif
