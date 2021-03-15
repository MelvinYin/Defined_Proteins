#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "Proteins.h"

void atomcpy(ATOM_STR *a1,ATOM_STR *a2);
void subatom(ATOM_STR a1,ATOM_STR a2,double *v);
double length(double *v);
void sum(double *v1,double *v2,double *v);
void mult(double *v1,double *v2,double *v);
void create_rotation_matrix(double pivot[3],float angle0,float rm[9]);
void rotate_vector(double v[3],float rm[9],double out[3]);


char *allresname[29]={"ALA ","ARG ","ASN ","ASP ","CYS ","CYSH","GLN ","GLU ","GLY ","HISD","HISE","HIS+","ILE ","LEU ","LYS ","MET ","PHE ","SER ","THR ","TRP ","TRPF","TYR ","VAL ","PRO ","HIS ","  A ","  T ","  G ","  C "};
//char *allresname[29]={"ALA ","ARG ","ASN ","ASP ","CYS ","CYSH","GLN ","GLU ","GLY ","HISD","HISE","HIS+","ILE ","LEU ","LYS ","MET ","PHE ","SER ","THR ","TRP ","TRPF","TYR ","VAL ","PRO ","HIS ","
char *hydr_file[2]={"hi.dat","hi_ion.dat"};
float addn[2]={0.204,1.28},addc[2]={0,-1.};
char ion=0;

struct hydrogen
{
  char res_name[5];
  int n_atom;
  char name_atom[23][5];
  int hydr_link[23];
  int invoke_key[23];
  int delta_1[23];
  int delta_2[23];
  int delta_3[23];
}residue[20];

char t[150];
char t1[150];
char t2[150];
char t3[150];
char t4[150];
char t5[150];



int getatom1(PROTEIN_STR p,int res,int n)
{
   int i;
  if(p.residue[res-p.from].num!=res){printf("Error in getatom1");exit(0);}
  i=p.residue[res-p.from].firstatom;
   if(i==p.numatoms){printf("No such atom");exit(-1);}
   return i+n;
}

int checkatom1(PROTEIN_STR p,int res,int n)
{
   int i;
  if(p.residue[res-p.from].num!=res)return -1;
  i=p.residue[res-p.from].firstatom;
//   if(i>=p.numatoms)return -1;
  if(i>=p.numatoms||p.atom[i+n].numres!=res)return -1;
   return i+n;
}

int getatom2(PROTEIN_STR p,int res,char* a)
{
  int i,j;
  char t[80];
  if(p.residue[res-p.from].num!=res){printf("Error in getatom2");exit(0);}
  i=p.residue[res-p.from].firstatom;
  for(j=0;p.atom[i+j].numres==res;j++)if(strcmp(p.atom[i+j].name,a)==0)break;
  if(i==p.numatoms||p.atom[i+j].numres!=res)
  {
    printf("In %d residue atom %s NOT found",res,a);
    exit(0);
  }
  return i+j;
}

int checkatom2(PROTEIN_STR p,int res,char* a)
{
  int i,j;
  if(p.residue[res-p.from].num!=res)return -1;
  i=p.residue[res-p.from].firstatom;
  for(j=0;p.atom[i+j].numres==res;j++)if(strcmp(p.atom[i+j].name,a)==0)break;
  if(i==p.numatoms||p.atom[i+j].numres!=res)return -1;
  return i+j;
}

double distatom(PROTEIN_STR p,int n1,int n2)
{
 return sqrt((p.atom[n1].x-p.atom[n2].x)*(p.atom[n1].x-p.atom[n2].x)+(p.atom[n1].y-p.atom[n2].y)*(p.atom[n1].y-p.atom[n2].y)+(p.atom[n1].z-p.atom[n2].z)*(p.atom[n1].z-p.atom[n2].z));
}

int getcoor(PROTEIN_STR *p,char *filename)
{
   char t[200];
   FILE *in;
   int natom,nres,prevres=-1;
   char nameatom[5],nameres[5];
   char chain,chain1;
   float x,y,z;
   int i,i1,i2,n,ii,type;
   char num_type,num_helix;
   n=0;
   chain=' ';
   if((in=fopen(filename,"r"))==NULL)
   {

     printf("NO file %s\n",filename);
     exit(-3);
   }
   p[n].atom=NULL;
   p[n].residue=NULL;

   for(ii=0;!feof(in);)
   {
     fgets(t,85,in);
     if(strncmp(t,"ATOM",4)!=0)continue;
     if(t[26]!=' ')continue;
     chain1=t[21];
     t[11]=0;
     t[16]=0;
     t[21]=0;
     t[26]=0;
     t[55]=0;
     if(sscanf(t+6,"%d",&natom)!=1)printf("Bad N_ATOM in file %s\n",filename);
     strcpy(nameatom,t+12);
     strcpy(nameres,t+17);
     if(sscanf(t+22,"%d",&nres)!=1)printf("Bad N_RES in atom %d file %s\n",natom,filename);

     if(prevres<0)
     {
        if(!(p[n].residue=(ST_RES *)realloc(p[n].residue,sizeof(ST_RES)))){printf("No memory\n");exit(0);}
        p[n].from=prevres=nres;
        p[n].residue[0].num=nres;
        p[n].residue[0].type=0;
        p[n].residue[0].firstatom=0;
        chain=p[n].chain=chain1;
        strcpy(p[n].residue[0].name,nameres);
     }
     else if(nres!=prevres||chain!=chain1)
     {
      if(nres>prevres&&chain==chain1)
      {
          if(!(p[n].residue=(ST_RES *)realloc(p[n].residue,(nres-p[n].from+1)*sizeof(ST_RES)))){printf("No memory\n");exit(0);}
          if(nres>(prevres+1))
           for(i=prevres+1;i<nres;i++)p[n].residue[i-p[n].from].num=-1;
          p[n].residue[nres-p[n].from].num=nres;
          p[n].residue[nres-p[n].from].firstatom=ii;
          p[n].residue[nres-p[n].from].type=' ';
          strcpy(p[n].residue[nres-p[n].from].name,nameres);
          prevres=nres;
      }
      else
      {
       p[n].to=prevres;
       p[n].numatoms=ii;
       n++;
       if(n>=MAXCHAIN){printf("Number of chains > 20");break;}
       ii=0;
       prevres=-1;
       p[n].atom=NULL;
       p[n].residue=NULL;

      }
     }
     if(ii)
     {
      for(i=p[n].residue[nres-p[n].from].firstatom;i<ii;i++)if(strcmp(nameatom,p[n].atom[i].name)==0)break;
      if(i!=ii)continue;
     }
     if(!(p[n].atom=(ATOM_STR *)realloc(p[n].atom,(ii+1)*sizeof(ATOM_STR)))){printf("No memory\n");exit(0);}
     sscanf(t+30,"%f%f%f",&x,&y,&z);
     p[n].atom[ii].x=x;
     p[n].atom[ii].y=y;
     p[n].atom[ii].z=z;
     p[n].atom[ii].num=natom;
     p[n].atom[ii].numres=nres;
     strcpy(p[n].atom[ii].name,nameatom);
     strcpy(p[n].atom[ii].nameres,nameres);
     ii++;
   }
   p[n].to=prevres;
   p[n].chain=chain;
   p[n].numatoms=ii;
   n++;
   for(i=0;i<n;i++)p[i].num_chain=n;
   rewind(in);
   num_helix=0;
   for(;;)
   {
     if(!fgets(t,90,in))break;
     type=' ';
     if(strncmp(t,"HELIX",5)==0)
     {
        t[25]=t[37]=0;
        chain=t[19];
        if(sscanf(t+21,"%d",&i1)!=1){printf("Bad format HELIX in PDB file %s\n",filename);exit(0);}
        if(sscanf(t+33,"%d",&i2)!=1){printf("Bad format HELIX in PDB file %s\n",filename);exit(0);}
        type='H';
        num_helix++;
        num_type=num_helix;
     }
     if(strncmp(t,"SHEET",5)==0)
     {
        t[26]=t[37]=0;
        chain=t[21];
        if(sscanf(t+22,"%d",&i1)!=1){printf("Bad format SHEET in PDB file %s\n",filename);exit(0);}
        if(sscanf(t+33,"%d",&i2)!=1){printf("Bad format SHEET in PDB file %s\n",filename);exit(0);}
        type='S';
     }
     if(strncmp(t,"TURN ",5)==0)
     {
        t[24]=t[35]=0;
        chain=t[19];
        if(sscanf(t+20,"%d",&i1)!=1){printf("Bad format TURN  in PDB file %s\n",filename);exit(0);}
        if(sscanf(t+31,"%d",&i2)!=1){printf("Bad format TURN  in PDB file %s\n",filename);exit(0);}
        type='T';
     }
     if(type==' ')continue;
     for(ii=0;ii<n;ii++)
      if(p[ii].chain==chain)
      {
        for(i=i1;i<=i2;i++)
         if(checkatom1(p[ii],i,0)>=0)
         {
          p[ii].residue[i-p[ii].from].type=type;
          p[ii].residue[i-p[ii].from].type_num=num_type;
         }
      }
   }
   fclose(in);
   return n;
}
void closeprotein(PROTEIN_STR *p)
{
   int i;
   for(i=0;i<p[0].num_chain;i++)
   {
     free(p[i].residue);
     p[i].residue=NULL;
     free(p[i].atom);
     p[i].atom=NULL;
   }
}

int getdat()
{
    int i,loop;
    FILE *in;
    if((in=fopen("Hi.dat","r"))==NULL){printf("Hi.dat NOT found\n");exit(0);}
    for(loop=0;loop<20;loop++)
    {
      fgets(t,140,in);
      if(t[3]=='\n')t[3]=' ';
      t[4]=0;
      strcpy(residue[loop].res_name,t);
      fgets(t,140,in);
      sscanf(t,"%d",&residue[loop].n_atom);
      fgets(t,140,in);
      fgets(t1,140,in);
      fgets(t2,140,in);
      fgets(t3,140,in);
      fgets(t4,140,in);
      fgets(t5,140,in);
      t[strlen(t)-1]=' ';
      for(i=0;i<residue[loop].n_atom;i++)
      {
        strcpy((char *)residue[loop].name_atom[i],"    ");
        strncpy((char *)(residue[loop].name_atom[i]+1),(char *)&t[i*6],3);
        t[i*6+5]=t1[i*6+5]=t2[i*6+5]=t3[i*6+5]=t4[i*6+5]=t5[i*6+5]=0;
        sscanf(&t1[i*6],"%d",&residue[loop].hydr_link[i]);
        sscanf(&t2[i*6],"%d",&residue[loop].invoke_key[i]);
        sscanf(&t3[i*6],"%d",&residue[loop].delta_1[i]);
        sscanf(&t4[i*6],"%d",&residue[loop].delta_2[i]);
        sscanf(&t5[i*6],"%d",&residue[loop].delta_3[i]);
      }
    }
    fclose(in);
    return(0);
}

void install_hydrogen(PROTEIN_STR *p)
{
   int imain,atom_23,one_of_20,ap1,ap2,ap3,ap_delta,delta=0;
   char cur[5],an_atom[5];
   int i,j,nh;
   double h[3][3];
   char tt[180];
   double v[3];
   double f,f1,f2;
   PROTEIN_STR ta;
   cur[0]=0;
   ta.atom=NULL;
   getdat();
   for(imain=0;imain<p->numatoms;imain++)
   {
       if((ta.atom = (ATOM_STR *) realloc(ta.atom, (imain+delta+4)*sizeof(ATOM_STR)))==NULL)
       {
               printf("No Memory!!!");
                   exit(-1);
       }
       atomcpy(&(p->atom[imain]),&(ta.atom[imain+delta]));
         ta.atom[imain+delta].num=imain+delta+1;
       if(strcmp(p->atom[imain].nameres,cur)!=0)
       {
          for(i=0;i<20;i++)
          {
             if(strncmp(residue[i].res_name,p->atom[imain].nameres,3)==0)
             {
               one_of_20=i;
               strcpy(cur,p->atom[imain].nameres);
               break;
             }
          }
          if(i==20)
          {
            if(strcmp(p->atom[imain].nameres,"  A ")==0 || strcmp(p->atom[imain].nameres,"  T ")==0 || strcmp(p->atom[imain].nameres,"  G ")==0 || strcmp(p->atom[imain].nameres,"  C ")==0)continue;
        printf("Unknown residue name %s\n",p->atom[imain].nameres);
            continue;
          }
       }
       for(atom_23=0;atom_23<23;atom_23++)
         if(strcmp(residue[one_of_20].name_atom[atom_23],p->atom[imain].name)==0)break;
       if(atom_23==23)
       {
      if(strcmp(p->atom[imain].name," OXT")==0)continue;
      if(strcmp(p->atom[imain].name," OT ")==0)continue;
        printf("Unknown atom name %s\n",p->atom[imain].name);
          continue;
       }
       if(residue[one_of_20].hydr_link[atom_23]==0)continue;
       if(strcmp(p->atom[imain].name," N  ")==0)
       {
         if(imain==0||p->atom[imain].numres==1)continue;
         for(i=-1;strcmp(p->atom[imain+i].name," C  ")!=0;i--);
         ap1=imain+i;
         ap2=imain+1;
         {
           double v1[3],v2[3];
       subatom(p->atom[imain],p->atom[ap1],v1);
       subatom(p->atom[imain],p->atom[ap2],v2);
           f1=length(v1);
           f2=length(v2);
           for(i=0;i<3;i++)
           {
             v1[i]/=f1;
             v2[i]/=f2;
           }
           sum(v1,v2,h[0]);
           f=length(h[0]);
       h[0][0]=h[0][0]/f*1.08+p->atom[imain].x;
       h[0][1]=h[0][1]/f*1.08+p->atom[imain].y;
       h[0][2]=h[0][2]/f*1.08+p->atom[imain].z;
         }
       }
       else
       {
         ap_delta=residue[one_of_20].delta_1[atom_23];
         ap1=imain+ap_delta;
         if(ap1>=p->numatoms)break;
         strcpy(an_atom,residue[one_of_20].name_atom[atom_23+ap_delta]);
         if(strcmp(p->atom[ap1].name,an_atom)!=0)
         {
             j=p->atom[imain].numres;
         for(ap1=imain;p->atom[ap1].numres==j;ap1--)if(strcmp(an_atom,p->atom[ap1].name)==0)break;
             if(p->atom[ap1].numres!=j)for(ap1=imain;p->atom[ap1].numres==j;ap1++)if(strcmp(an_atom,p->atom[ap1].name)==0)break;
                if(p->atom[ap1].numres!=j)
             {
        printf("In %d %s '%s' out of order, '%s' expected. Atom structure cannot be corrected.",p->atom[ap1].numres,p->atom[ap1].nameres,p->atom[ap1].name,an_atom);

                continue;
             }
         }
         ap_delta=residue[one_of_20].delta_2[atom_23];
         ap2=imain+ap_delta;
     if(ap2>=p->numatoms)break;
         strcpy(an_atom,residue[one_of_20].name_atom[atom_23+ap_delta]);
         if(strcmp(p->atom[ap2].name,an_atom)!=0)
         {
             j=p->atom[imain].numres;
         for(ap2=imain;p->atom[ap2].numres==j;ap2--)if(strcmp(an_atom,p->atom[ap2].name)==0)break;
             if(p->atom[ap2].numres!=j)for(ap2=imain;p->atom[ap2].numres==j;ap2++)if(strcmp(an_atom,p->atom[ap2].name)==0)break;
             if(p->atom[ap2].numres!=j)
             {
        sprintf(t,"In %d %s '%s' out of order, '%s' expected.Atom structure cannot be corrected. ",p->atom[ap2].numres,p->atom[ap2].nameres,p->atom[ap2].name,an_atom);
                continue;
             }
         }
         ap3=imain+residue[one_of_20].delta_3[atom_23];
     if(ap3>=p->numatoms)break;
         switch(residue[one_of_20].invoke_key[atom_23])
         {
           case 1:
             {
               double v1[3],v2[3],v3[3];
           subatom(p->atom[imain],p->atom[ap1],v1);
               subatom(p->atom[imain],p->atom[ap2],v2);
               f=length(v1);
               for(i=0;i<3;i++)v1[i]=v1[i]/f*0.5*1.08;
               mult(v2,v1,v3);
               mult(v3,v1,v2);
               f=length(v2);
               for(i=0;i<3;i++)v2[i]=v2[i]/f*0.866*1.08;
           h[0][0]=v2[0]+v1[0]+p->atom[imain].x;
           h[0][1]=v2[1]+v1[1]+p->atom[imain].y;
           h[0][2]=v2[2]+v1[2]+p->atom[imain].z;
               break;
             }
           case 2:
             {
               double v1[3],v2[3];
           subatom(p->atom[imain],p->atom[ap1],v1);
               subatom(p->atom[imain],p->atom[ap2],v2);
               f1=length(v1);
               f2=length(v2);
               for(i=0;i<3;i++)
               {
                v1[i]/=f1;
                v2[i]/=f2;
                h[0][i]=v1[i]+v2[i];
               }
               f=length(h[0]);
               h[0][0]=h[0][0]/f*1.08+p->atom[imain].x;
               h[0][1]=h[0][1]/f*1.08+p->atom[imain].y;
               h[0][2]=h[0][2]/f*1.08+p->atom[imain].z;
               break;
             }
           case 3:
             {
               double v1[3],v2[3],v3[3];
           subatom(p->atom[imain],p->atom[ap1],v1);
           subatom(p->atom[imain],p->atom[ap2],v2);
               f=length(v1);
               for(i=0;i<3;i++)v1[i]=v1[i]/f*0.5*1.08;
               mult(v2,v1,v3);
               mult(v3,v1,v2);
               f=length(v2);
               for(i=0;i<3;i++)v2[i]=v2[i]/f*0.866*1.08;
               h[0][0]=v2[0]+v1[0]+p->atom[imain].x;
               h[0][1]=v2[1]+v1[1]+p->atom[imain].y;
               h[0][2]=v2[2]+v1[2]+p->atom[imain].z;
               h[1][0]=v1[0]-v2[0]+p->atom[imain].x;
               h[1][1]=v1[1]-v2[1]+p->atom[imain].y;
               h[1][2]=v1[2]-v2[2]+p->atom[imain].z;
               break;
             }
           case 4:
             {
                double v1[3],v2[3],v3[3];
        subatom(p->atom[imain],p->atom[ap1],v1);
        subatom(p->atom[imain],p->atom[ap2],v2);
        subatom(p->atom[imain],p->atom[ap3],v3);
                f=length(v1);
                f1=length(v2);
                f2=length(v3);
                for(i=0;i<3;i++)
                {
                  v1[i]=v1[i]/f;
                  v2[i]=v2[i]/f1;
                  v3[i]=v3[i]/f2;
                  h[0][i]=v1[i]+v2[i]+v3[i];
                }
                f=length(h[0]);
                h[0][0]=h[0][0]/f*1.08+p->atom[imain].x;
                h[0][1]=h[0][1]/f*1.08+p->atom[imain].y;
                h[0][2]=h[0][2]/f*1.08+p->atom[imain].z;
                break;
             }
           case 5:
             {
               double v1[3],v2[3],v3[3],v4[3];
           subatom(p->atom[imain],p->atom[ap1],v1);
           subatom(p->atom[imain],p->atom[ap2],v2);
               mult(v1,v2,v4);
               for(i=0;i<3;i++)v3[i]=v1[i]+v2[i];
               f1=length(v3);
               f2=length(v4);
               for(i=0;i<3;i++)
               {
                 v3[i]=v3[i]/f1*0.5*1.08;
                 v4[i]=v4[i]/f2*0.866*1.08;
               }
           h[0][0]=v3[0]+v4[0]+p->atom[imain].x;
           h[0][1]=v3[1]+v4[1]+p->atom[imain].y;
           h[0][2]=v3[2]+v4[2]+p->atom[imain].z;
           h[1][0]=v3[0]-v4[0]+p->atom[imain].x;
           h[1][1]=v3[1]-v4[1]+p->atom[imain].y;
           h[1][2]=v3[2]-v4[2]+p->atom[imain].z;
               break;
             }
           case 6:
             {
                double v1[3],v2[3],v3[3];
                float rm[9];
        subatom(p->atom[imain],p->atom[ap1],v1);
        subatom(p->atom[imain],p->atom[ap2],v2);
                f=length(v1);
                for(i=0;i<3;i++)v1[i]=v1[i]/f*0.5*1.08;
                mult(v2,v1,v3);
                mult(v1,v3,v2);
                f1=length(v2);
                f2=length(v1);
                for(i=0;i<3;i++)
                {
                  h[0][i]=v2[i]/f1*0.866*1.08+v1[i];
                  v1[i]=v1[i]/f2;
                }
                create_rotation_matrix(v1,120.,rm);
                rotate_vector(h[0],rm,h[1]);
                rotate_vector(h[1],rm,h[2]);
                for(i=0;i<3;i++)
                {
                  h[i][0]+=p->atom[imain].x;
                  h[i][1]+=p->atom[imain].y;
                  h[i][2]+=p->atom[imain].z;
                }
                break;
             }
         }
       }
       for(nh=0;nh<residue[one_of_20].hydr_link[atom_23];nh++)
       {
         delta++;
         strcpy(t,p->atom[imain].name);
         if(residue[one_of_20].hydr_link[atom_23]>1)t[0]=nh+49;
         t[1]='H';
     strcpy(ta.atom[imain+delta].name,t);
     strcpy(ta.atom[imain+delta].nameres,p->atom[imain].nameres);
         ta.atom[imain+delta].numres=p->atom[imain].numres;
         ta.atom[imain+delta].num=imain+delta+1;
     ta.atom[imain+delta].x=h[nh][0];
     ta.atom[imain+delta].y=h[nh][1];
     ta.atom[imain+delta].z=h[nh][2];
     ta.atom[imain+delta].charge=0;

       }
   }
   p->numatoms=imain+delta;
     free(p->atom);
     p->atom=ta.atom;
     j=-1;
     for(i=0;i<p->numatoms;i++)
      if(p->atom[i].numres!=j)
      {
       p->residue[p->atom[i].numres-p->from].firstatom=i;
       j=p->atom[i].numres;

      }
}

void atomcpy(ATOM_STR *a1,ATOM_STR *a2)
{
                 a2->x=a1->x;
                 a2->y=a1->y;
                 a2->z=a1->z;
         a2->num=a1->num;
         a2->numres=a1->numres;
                 a2->charge=a1->charge;
                 strcpy(a2->name,a1->name);
         strcpy(a2->nameres,a1->nameres);
}

void subatom(ATOM_STR a1,ATOM_STR a2,double *v)
{
   v[0]=a1.x-a2.x;
   v[1]=a1.y-a2.y;
   v[2]=a1.z-a2.z;
}
double length(double *v)
{
  return sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
}
void sum(double *v1,double *v2,double *v)
{
  int i;
  for(i=0;i<3;i++)v[i]=v1[i]+v2[i];
}
void mult(double *v1,double *v2,double *v)
{
  v[0]=v1[1]*v2[2]-v1[2]*v2[1];
  v[1]=v1[2]*v2[0]-v1[0]*v2[2];
  v[2]=v1[0]*v2[1]-v1[1]*v2[0];
}

void create_rotation_matrix(double pivot[3],float angle0,float rm[9])
{
  float angle,lambda,mu,nu,ro;
  angle=angle0*3.14159/180.*0.5;
  lambda=pivot[0]*sin(angle);
  mu    =pivot[1]*sin(angle);
  nu    =pivot[2]*sin(angle);
  ro    =cos(angle);
  rm[0]=lambda*lambda-mu*mu-nu*nu+ro*ro;
  rm[1]=2.*(lambda*mu+nu*ro);
  rm[2]=2.*(lambda*nu-mu*ro);
  rm[3]=2.*(lambda*mu-nu*ro);
  rm[4]=mu*mu-nu*nu-lambda*lambda+ro*ro;
  rm[5]=2.*(mu*nu+lambda*ro);
  rm[6]=2.*(lambda*nu+mu*ro);
  rm[7]=2.*(mu*nu-lambda*ro);
  rm[8]=nu*nu-mu*mu-lambda*lambda+ro*ro;
}

void rotate_vector(double v[3],float rm[9],double out[3])
{
  out[0]=v[0]*rm[0]+v[1]*rm[3]+v[2]*rm[6];
  out[1]=v[0]*rm[1]+v[1]*rm[4]+v[2]*rm[7];
  out[2]=v[0]*rm[2]+v[1]*rm[5]+v[2]*rm[8];
}

void writecoor(PROTEIN_STR p,char *filename)
{
       int i;
       FILE *in;
       if((in=fopen(filename,"a"))==NULL){printf("Cannot write coordinates to file %s\n",filename);exit(0);}
       for(i=0;i<p.numatoms;i++)
			 fprintf(in,"ATOM  %5d %s %s%c%4d    %8.3f%8.3f%8.3f\n",
           p.atom[i].num,p.atom[i].name,p.atom[i].nameres, p.chain, p.atom[i].numres,
           p.atom[i].x,p.atom[i].y,p.atom[i].z);
       fclose(in);
}

void getcharge(PROTEIN_STR *p,char *filename)
{
        char t[90],t1[102],t2[10];
        int i,j,n,k;
        float d,f;
        FILE *in;
        if((in=fopen(filename,"r"))==NULL)
        {
         printf("Cannot open charge file in Get Charge");
         exit(0);
        }
        for(i=0;i<p->numatoms;i++)p->atom[i].charge=0;
        for(;;)
        {
                if(fgets(t1,100,in)==NULL)break;
                t1[4]=0;
                for(i=0;i<29;i++)if(strcmp(t1,allresname[i])==0)break;
                if(i==29)continue;
                for(;;)
                {
                        if(fgets(t,80,in)==NULL)break;
                        if(t[1]<55)break;
                        for(j=0;j<7;j++)t2[j]=t[4+j];
                        t2[7]=0;
                        t[4]=0;
                        d=(float)atof(t2);
                        for(j=0;j<p->numatoms;j++)
                          if((strcmp(t1,p->atom[j].nameres)==0)&&(strcmp(t,p->atom[j].name)==0))
                             p->atom[j].charge=d;
                }
        }
        n=0;
        f=0;
        fclose(in);
//        for(i=0;i<p->numatoms;i++)if(p->atom[i].charge==0)printf("No charge in %d %s %s\n",p->atom[i].numres,p->atom[i].nameres,p->atom[i].name);
        for(i=0;i<50;i++)if(strcmp(p->atom[i].name," N  ")==0 && p->atom[i].numres==p->from)break;
        if(i<50 && addn[ion]!=0)
        {
          if(strcmp(p->atom[i+1].name," H  ")!=0)
             p->atom[i].charge+=addn[ion];
        }
        if(addc[ion]!=0)
        {
         n=p->numatoms-1;
         for(j=n;p->atom[j].numres==p->atom[n].numres;j--)
            if(strcmp(p->atom[j].name," O  ")==0)break;
         if(p->atom[j].numres==p->atom[n].numres)
         {
           if(strcmp(p->atom[n].name," OXT")!=0)
             p->atom[j].charge+=addc[ion];
            else
            {
              p->atom[j].charge+=addc[ion]/2;
              p->atom[n].charge+=addc[ion]/2;
            }
         }
        }

}

/*char *inp_file[1000] ;
int I ;
void find_ent()
{
	struct ffblk ffblk ;
	int done,i ;
	for(i=0;i<1000;i++)
	if((inp_file[i]=(char*)calloc(14,sizeof(char)))==NULL)
	{ printf("Not enough memory!!!\n"); exit(1); }
	done = findfirst("*.ent",&ffblk,0) ;
	for (i=0;!done;i++)
	{
	 if(i>1000)
	  {
		printf("\nSorry, Sir! You have more than 1000 files of structures.") ;
		printf("\nPlease, change value of array inp_file in module find_cds.c") ;
		exit(1) ;
	  }
	 strcpy(inp_file[i],ffblk.ff_name) ;
	 done = findnext(&ffblk) ;
	}
	I=i ;
	return ;
}*/


/*
PROTEIN p1;
void main()
{
  int i,j,n;
  n=getcoor(p1,"201l.pdb");
  for(i=p1[0].from;i<=p1[0].to;i++)
  if(p1[0].residue[i-p1[0].from].num>=0)
   printf("%d %s\n",p1[0].residue[i-p1[0].from].num,p1[0].residue[i-p1[0].from].name);
  install_hydrogen(&(p1[0]));
  writecoor(p1[0],"201l.pdb");

  printf("%d %f\n",i,(float)distatom(p1[0],getatom2(p1[0],3," CA "),getatom2(p1[0],10," CA ")));
  for(j=0;j<n;j++)
  for(i=p1[j].from;i<=p1[j].to;i++)
   if(p1[j].residue[i-p1[j].from].num>=0)
    printf("%d %s %c\n",i,p1[0].atom[getatom1(p1[j],i,0)].nameres,p1[j].residue[i-p1[j].from].type);
}
*/
