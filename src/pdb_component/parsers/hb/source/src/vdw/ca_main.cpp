#include "pdb_data.hpp"
#include "contact_data.hpp"
#include "contact_formatting.hpp"

int main(int argc, char *argv[])
{
   if (argc<4)
    {
      printf("You must supply a filename to compute on, filename for contacts data, as well as the filename for residue bounds data.");
      return -1;
    };
   // Read the structure in
   string input_fname(argv[1]);
   string atomcs_fname(argv[2]);
   string bounds_fname(argv[3]);
   FILE* fh = openFile(argv[1], "r");  
   if (! fh)
     return -1;
   AtomVector vec = parseAtomLines(fh);
   fclose(fh);
   // Compute contacts on it
   ContactResult r = contactForAtoms(vec);
   // Format distances for individual atom contacts and output it
   AtomBasicsVector atomcs = formatAtomBasics(vec, r);
   writeAtomBasics(atomcs_fname, atomcs);
   // Write chain residues (bounds)
   ChainBounds bounds = makeChainBounds(vec);
   writeChainBounds(bounds_fname, bounds);
   return 0;
};
