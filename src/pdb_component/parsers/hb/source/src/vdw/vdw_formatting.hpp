/*
  Module: vdw_formatting
  Authors: Grzegorz M. Koczyk (2007), Igor Berezovsky (1990-2007)
  *****************************************************************************
  Formatting of Lennard-Jones results for atom-atom or residue-residue contact information.
 */

#ifndef __VDW_FORMATTING_HPP__
#define __VDW_FORMATTING_HPP__
#include "vdw_data.hpp"
#include "basic_formatting.hpp"

class AtomContactInfo
{
public:
  Atom m_atom1;
  Atom m_atom2;
  lj_value_t m_lj;  
public:
  AtomContactInfo(const Atom& atom1, const Atom& atom2, const lj_value_t& lj): m_atom1(atom1), m_atom2(atom2), m_lj(lj) {};
};
typedef vector<AtomContactInfo> AtomContacts;

// Create a list of annotated contact data between atoms
AtomContacts formatAtomContacts(AtomVector& atomvec, LJResult& ljresult)
{
  AtomContacts result;
  for (LJResult::iterator it=ljresult.begin(); it!=ljresult.end(); it++)
    result.push_back( AtomContactInfo( atomvec[ it->m_atom1_no ], atomvec[ it->m_atom2_no ], it->m_value ) );
  return result;
};


bool writeAtomContacts(FILE* fh, AtomContacts& atomcs)
{
  fprintf(fh, "%u\n", (unsigned long) atomcs.size() );
  for (AtomContacts::iterator it=atomcs.begin(); it!=atomcs.end(); it++)
    {
      Atom& atom1 = it->m_atom1;
      Atom& atom2 = it->m_atom2;
      lj_value_t& lj = it->m_lj;
      double d = atomDistance(atom1, atom2);
      fprintf(fh, "%lu\t%ld\t%c\t%s\t%c\t%s\t", atom1.m_number, atom1.m_res_no, atom1.m_icode, atom1.m_res_name.c_str(), atom1.m_chain, atom1.m_name.c_str() );
      fprintf(fh, "%lu\t%ld\t%c\t%s\t%c\t%s\t", atom2.m_number, atom2.m_res_no, atom2.m_icode, atom2.m_res_name.c_str(), atom2.m_chain, atom2.m_name.c_str() );
      fprintf(fh, "%f\t%f\n", lj, d);
    };
  return true;
};

bool writeAtomContacts(const string& fname, AtomContacts& atomcs)
{
  FILE* fh = openFile(fname.c_str(), "wt");
  if (!fh)
    return false;
  bool r = writeAtomContacts(fh, atomcs);
  fclose(fh);
  return r;
};


class ResidueContactInfo
{
public:
  string m_resname1;
  string m_resname2;
  unsigned long m_contacts;
  double m_lj;  
public:
  ResidueContactInfo(string resname1="", string resname2="", unsigned long contacts=0, double lj=0.): m_resname1(resname1), m_resname2(resname2), m_contacts(contacts), m_lj(lj) {};
};

typedef map<ResId, ResidueContactInfo> map_sno_t;
typedef map<char, map_sno_t> map_sch_t;
typedef map<ResId, map_sch_t> map_fno_t;
typedef map<char, map_fno_t> map_fch_t;


// Create a list of annotated contact data between residues
map_fch_t asResidueContacts(AtomContacts& atomcs) 
{
  map_fch_t r;
  for (AtomContacts::iterator it=atomcs.begin(); it!=atomcs.end(); it++)
    {
      Atom& atom1 = it->m_atom1;
      Atom& atom2 = it->m_atom2;
      lj_value_t& lj = it->m_lj;
      ResidueContactInfo& tmp = r[atom1.m_chain][ResId(atom1.m_res_no, atom1.m_icode)][atom2.m_chain][ResId(atom2.m_res_no, atom2.m_icode)];
      if (tmp.m_contacts==0)
	{
	  tmp.m_resname1 = atom1.m_res_name;
	  tmp.m_resname2 = atom2.m_res_name;
	};
      tmp.m_contacts+=1;
      tmp.m_lj += lj;
    };
  return r;
};

bool writeAsResidueContacts(FILE* fh, AtomContacts& atomcs)
{
  map_fch_t r = asResidueContacts(atomcs);
  // Iterate over maps of [chain]->[residue no]->[chain]->[residue no]
  for (map_fch_t::iterator it1=r.begin(); it1!=r.end(); it1++)
    {
      char chain1 = it1->first;
      map_fno_t& st2= it1->second;
      for (map_fno_t::iterator it2=st2.begin(); it2!=st2.end();it2++)
	{
	  ResId resid1 = it2->first;
	  unsigned long resno1 = resid1.m_resno;
	  char icode1 = resid1.m_icode;
	  map_sch_t& st3= it2->second;
	  for (map_sch_t::iterator it3=st3.begin(); it3!=st3.end();it3++)
	    {
	      char chain2 = it3->first;
	      map_sno_t& st4= it3->second;
	      for (map_sno_t::iterator it4=st4.begin(); it4!=st4.end();it4++)
		{
		  ResId resid2 = it4->first;
		  unsigned long resno2 = resid2.m_resno;
		  char icode2 = resid2.m_icode;
		  ResidueContactInfo& rc = it4->second;
		  fprintf(fh, "%c\t%ld\t%c\t%s\t%c\t%ld\t%c\t%s\t%u\t%f\n", chain1, resno1, icode1, rc.m_resname1.c_str(), chain2, resno2, icode2, rc.m_resname2.c_str(), rc.m_contacts, rc.m_lj);
		};	      
	    };
	};
    };
  return true;
};

bool writeAsResidueContacts(const string& fname, AtomContacts& atomcs)
{
  FILE* fh = openFile(fname.c_str(), "wt");
  if (!fh)
    return false;
  bool r = writeAsResidueContacts(fh, atomcs);
  fclose(fh);
  return r;
};

// bool asSplitEnergies(ResidueContacts& residuecs)
// {
//   // For each chain in the structure, analyze only residues of the same chain
//   for (map_fch_t::iterator it1=residuecs.begin(); it1!=residuecs.end(); it++)
//     {
//       char chain1 = it1->first;
//       map_fno_t& st2= it1->second;
//       for (map_fno_t::iterator it2=st2.begin(); it2!=st2.end();it2++)
// 	{
// 	  int resno1 = it2->first;
// 	  map_sch_t& st3= it2->second;
// 	  for (map_sch_t::iterator it3=st3.begin(); it3!=st3.end();it3++)
// 	    {
// 	      char chain2 = it3->first;
// 	      if (chain2!=chain1)
// 		continue;
// 	      map_sno_t& st4= it3->second;
// 	      for (map_sno_t::iterator it4=st4.begin(); it4!=st4.end();it4++)
// 		{		  
// 		  int resno2 = it4->first;
// 		  if (resno2<resno1)
// 		  ResidueContactInfo& rc = it4->second;
// 		  fprintf(fh, "%c\t%u\t%s\t%c\t%u\t%s\t%u\t%f\n", chain1, resno1,  rc.m_resname1.c_str(), chain2, resno2, rc.m_resname2.c_str(), rc.m_contacts, rc.m_lj);
// 		};	      
// 	    };
// 	};
//     };
// };

// bool asSplitEnergies(AtomContacts& atomcs)
// {
//   return asSplitEnergies(asResidueContacts(atomcs));
// };

// bool writeSplitEnergies(FILE* fh, AtomContacts& atomcs)
// {
// };

// bool writeSplitEnergies(FILE* fh, ResidueContacts& atomcs)
// {
// };

#endif
