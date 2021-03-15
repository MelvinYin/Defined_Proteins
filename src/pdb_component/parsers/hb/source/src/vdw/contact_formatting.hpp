/*
  Module: vdw_formatting
  Authors: Grzegorz M. Koczyk (2007), Igor Berezovsky (1990-2007)
  *****************************************************************************
  Formatting of Lennard-Jones results for atom-atom or residue-residue contact information.
 */

#ifndef __CONTACT_FORMATTING_HPP__
#define __CONTACT_FORMATTING_HPP__
#include "contact_data.hpp"
#include "basic_formatting.hpp"


class AtomContactBasics
{
public:
  Atom m_atom1;
  Atom m_atom2;
  dist_value_t m_distance;  
public:
  AtomContactBasics(const Atom& atom1, const Atom& atom2, const dist_value_t& v): m_atom1(atom1), m_atom2(atom2), m_distance(v) {};
};
typedef vector<AtomContactBasics> AtomBasicsVector;

// Create a list of annotated contact data between atoms
AtomBasicsVector formatAtomBasics(AtomVector& atomvec, ContactResult& cresult)
{
  AtomBasicsVector result;
  for (ContactResult::iterator it=cresult.begin(); it!=cresult.end(); it++)
    result.push_back( AtomContactBasics( atomvec[ it->m_atom1_no ], atomvec[ it->m_atom2_no ], it->m_distance ) );
  return result;
};

bool writeAtomBasics(FILE* fh, AtomBasicsVector& atomcs)
{
  for (AtomBasicsVector::iterator it=atomcs.begin(); it!=atomcs.end(); it++)
    {
      Atom& atom1 = it->m_atom1;
      Atom& atom2 = it->m_atom2;
      dist_value_t& d = it->m_distance;
      fprintf(fh, "%lu\t%ld\t%c\t%s\t%c\t%s\t", atom1.m_number, atom1.m_res_no, atom1.m_icode, atom1.m_res_name.c_str(), atom1.m_chain, atom1.m_name.c_str() );
      fprintf(fh, "%lu\t%ld\t%c\t%s\t%c\t%s\t", atom2.m_number, atom2.m_res_no, atom2.m_icode, atom2.m_res_name.c_str(), atom2.m_chain, atom2.m_name.c_str() );
      fprintf(fh, "%f\n", d);
    };
  return true;
};

bool writeAtomBasics(const string& fname, AtomBasicsVector& atomcs)
{
  FILE* fh = openFile(fname.c_str(), "wt");
  if (!fh)
    return false;
  bool r = writeAtomBasics(fh, atomcs);
  fclose(fh);
  return r;
};

#endif
