/* Module prep.c for preparing set of PDB-fails or particular
   PDB-fail for using in calculations.
   ( Module find_ent.c has been included into protein.c )
   Version2.0
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Proteins.h"
PROTEIN p1;
extern char *inp_file[] ;
extern int I ;
int main(int argc, char* argv[])
{
  if (argc!=3)
    {
      printf("Incorrect number of args: %d",argc) ;
      exit(1) ;
    }
  char inpfile[15] ;
  char newfn[25] ;
  char *ptrfn ;
  int n,d,ii, no_chain ;
  char m[81];
  FILE *ll, *fp;
  int kk,t;
  t=0;
  strcpy(inpfile,argv[1]);
  t=strlen(inpfile);

  if ((fopen(inpfile,"r"))==NULL)
    {
      printf("Can't open file %s",inpfile) ;
      exit(1) ;
    }
  ptrfn=strcpy(newfn,argv[2]) ;
  //strcat(ptrfn,".hhh") ;
  n=getcoor(p1,inpfile);
  // Open for writing 
  // PDB-file with hydrogens {argv[1]}.hhh (newfn)
  if ((fp=fopen(newfn,"w"))==NULL)
    {
      printf("Can't open file %s",newfn) ;
      exit(1) ;
    }
  fclose(fp);
  for (no_chain=0; no_chain<n; no_chain++)
    {
      install_hydrogen(&(p1[no_chain]));
      writecoor(p1[no_chain],newfn);
    };
  printf("%s",m);
  return 0;
}
