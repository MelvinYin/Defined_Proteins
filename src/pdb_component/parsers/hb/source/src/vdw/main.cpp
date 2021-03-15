#include "pdb_data.hpp"
#include "vdw_data.hpp"
#include "vdw_constants.hpp"
#include "vdw_formatting.hpp"

int main(int argc, char *argv[])
{
   if (argc<5)
    {
      printf("You must supply a filename to compute on, atom contacts filename, residue contacts filename as well as the filename for residue bounds data.");
      return -1;
    };
   // Read the structure in
   string input_fname(argv[1]);
   string atomcs_fname(argv[2]);
   string res_fname(argv[3]);
   string bounds_fname(argv[4]);
   FILE* fh = openFile(argv[1], "r");  
   if (! fh)
     return -1;
   AtomVector vec = parseAtomLines(fh);
   fclose(fh);
   // Compute LJ on it
   LJResult r = ljForAtoms(vec);
   // Format LJ for individual atom contacts and output it
   AtomContacts atomcs = formatAtomContacts(vec, r);
   writeAtomContacts(atomcs_fname, atomcs);
   // Format LJ for residue-residue contacts and output it
   writeAsResidueContacts(res_fname, atomcs);
   ChainBounds bounds = makeChainBounds(vec);
   writeChainBounds(bounds_fname, bounds);
   return 0;
};
