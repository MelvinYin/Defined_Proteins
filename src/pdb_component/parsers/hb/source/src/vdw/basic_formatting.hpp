/*
  Module: basic_formatting
  Authors: Grzegorz M. Koczyk (2007), Igor Berezovsky (1990-2007)
  *****************************************************************************
  Helper functions for formatting output
 */

#ifndef __BASIC_FORMATTING_HPP__
#define __BASIC_FORMATTING_HPP__
#include "pdb_data.hpp"
// Utility function to open a file for reading or writing, outputting a message on error
FILE* openFile(const string& fname, const string& mode)
{
  FILE* fh = fopen(fname.c_str(), mode.c_str());
  if (!fh)
    {
      if (mode[0]=='r')	
	printf( "Input file %s non-existent or corrupted.\n", fname.c_str() );
      else
	printf( "Output file %s could not be opened.\n", fname.c_str() );
    };
  return fh;  
};


bool writeChainBounds(FILE* fh, ChainBounds& bounds)
{
  for (ChainBounds::iterator it=bounds.begin(); it!=bounds.end(); it++)
    {
      char chain = it->first;
      resnoset_t r = it->second;
      for (resnoset_t::iterator it1=r.begin(); it1!=r.end(); it1++)
	{
	  ResId resid = *it1;
	  long resno = resid.m_resno;
	  char icode = resid.m_icode;
	  fprintf(fh, "%c\t%ld\t%c\n", chain, resno, icode);
	};
    };
  return true;
};

bool writeChainBounds(const string& fname, ChainBounds& bounds)
{
  FILE* fh = openFile(fname.c_str(), "wt");
  if (!fh)
    return false;
  bool r = writeChainBounds(fh, bounds);
  fclose(fh);
  return r;
};

#endif
