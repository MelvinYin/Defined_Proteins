#define MAXCHAIN 200
typedef struct
{
  char name[5];
  char nameres[5];
  int num;
  int numres;
  float charge;
  float x,y,z;
}ATOM_STR;
typedef struct
{
    int num;
    unsigned firstatom;
    char type;     // ' '-none, 'H'-Helix, 'S'-Sheet, 'T'-Turn
    char type_num;
    char name[5];
}ST_RES;

typedef struct
{
  ATOM_STR *atom;
  int numatoms;
  int from,to;
  ST_RES *residue;
  char num_chain;
  char chain;
}PROTEIN_STR;
typedef PROTEIN_STR PROTEIN [MAXCHAIN];

extern char ion; // Ionization. ion = 0 (default) on neutral model, ion = 1 on ionic model

int getcoor(PROTEIN_STR *p,char *filename); // Read PDB file. Return the number of chains
void closeprotein(PROTEIN_STR *p); // Close protein
int getatom1(PROTEIN_STR p,int res,int n); // Get number of atom in protein from residue number and number of atom in residue
int checkatom1(PROTEIN_STR p,int res,int n); //Check atom and return -1 on error or number of atom in protein from residue number and number of atom in residue
int getatom2(PROTEIN_STR p,int res,char* a); // Get number of atom in protein from residue number and name of atom in residue
int checkatom2(PROTEIN_STR p,int res,char* a); // Check atom and return -1 on error or number of atom in protein from residue number and name of atom in residue
int getres(PROTEIN_STR p, int num,char *a);
double distatom(PROTEIN_STR p,int n1,int n2); //Get distance between atoms with number in protein n1, n2
void install_hydrogen(PROTEIN_STR *p); // Find coordinates of hydrogens for protein 'p'. Use pointer of protein for one chain
void writecoor(PROTEIN_STR p,char *filename); // Write protein 'p' to file filename
void getcharge(PROTEIN_STR *p,char *filename);

