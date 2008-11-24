// hhcov.C: 
// Format hhsearch database for transitive scoring
// Compile with g++ hhcov.C -o hhcov -O3 -static

#define MAIN
#include <iostream>   // cin, cout, cerr
#include <fstream>    // ofstream, ifstream 
#include <stdio.h>    // printf
#include <stdlib.h>   // exit
#include <string>     // strcmp, strstr
#include <math.h>     // sqrt, pow
#include <limits.h>   // INT_MIN
#include <float.h>    // FLT_MIN
#include <ctype.h>    // islower, isdigit etc
#include <time.h>     // clock
#include <dirent.h>   // DIR, struct dirent, opendir, readdir, closedir, ...
#include <errno.h>    // perror()


using std::cout;
using std::cerr;
using std::endl;
using std::ios;
using std::ifstream;
using std::ofstream;

#include "util.C"        // imax, fmax, iround, iceil, ifloor, strint, strscn, strcut, substr, uprstr, uprchr, Basename etc.
#include "list.C"        // list data structure
#include "hash.C"        // hash data structure




/////////////////////////////////////////////////////////////////////////////////////
// Global variables 
/////////////////////////////////////////////////////////////////////////////////////

const char VERSION_AND_DATE[]="version 1.0 November 2005";
const char REFERENCE[]="Soding, J. To be published\n";
const char COPYRIGHT[]="(C) Johannes Soeding (see LICENSE file)\n";
const int MAXFILES=32768; // maximum number of files to process
const int LINELEN=32768;  // max length of line read in from input files; must be >= MAXCOL 
const int NAMELEN=128;    // max length of file names etc.
const int IDLEN=32;       // max length of scop hierarchy id and pdb-id
 
// Input parameters
class Parameters       // Parameters for gap penalties and pseudocounts
{
public:
  char scoresdir[NAMELEN]; // directory name with *.scores files
  char outfile[NAMELEN];   // output filename
  float Pmin;             // P-values P(k,l)<Pmin do not contribute to the negative Z score distribution
  float Pmax;             // probabilities P(k,l)>Pmax% do not contribute to the covariance calculation
  char fast;
};

int v=2;                   // verbose mode
Parameters par;



/////////////////////////////////////////////////////////////////////////////////////
// Help function
/////////////////////////////////////////////////////////////////////////////////////
void help()
{
  printf("\n");
  printf("HHcov %s\n",VERSION_AND_DATE);
  printf("Calculate a covariance matrix of an hhsearch database for transitive scoring\n");
  printf("When hhsearching a database in transitive scoring mode, the score between query\n");
  printf("HMM Q and a database HMM T will not only depend on the comparison of Q and T,\n");
  printf("but also on the transitive comparisons Q-P, P-T via intermediate database HMMs P.\n");
  printf("This transitive scoring method can increase the sensitivity considerably.\n");
  printf("HHcov needs to read the *.scores file of every database HMM.\n");
//  printf("%s",REFERENCE);
  printf("%s",COPYRIGHT);
  printf("Usage: hhcov -o <file.hhc> [options]\n");
  printf("\n");
  printf("Options:\n");
  printf(" -o <file>   output .hhc file containing weights and pairwise probabilities\n");
  printf(" -d <dir>    directory with *.scores files of all database HMMs (def='.')\n");
  printf(" -v <int>    verbose mode (def=%i)\n",v);
  printf("\n");
  printf("\n");
  return;
}


void ProcessArguments(int argc, char** argv); // Processing input options from command
inline void  ScopID(char cl[], char fold[], char sfam[], const char fam[]); // a.118.8.2 => 'a' 'a.118' 'a.118.8'
void PrintMatrix(float** V, int N, const char* str);
void ReadMatrix(float** V, int N, const char* str);
float Score2Z(float S);                      // sqrt(2)*dierfc(2*2^(-S)) transforms a -log2(P-value) bit score into a Z-score
float dierfc(float y);        



/////////////////////////////////////////////////////////////////////////////////////
//// MAIN PROGRAM
/////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{  
  // Variable declarations
  char line[LINELEN]="";  // input line
  char tname[NAMELEN];
  char qname[NAMELEN];
  char qfam[IDLEN];
  char scoresfile[NAMELEN];
  FILE* scoresf=NULL;
  int qlen;

  int N=0;             // number of scores files in directory
  int M=0;             // number of HMMs in each scores file
  int k,l,m=0;         // HMM indices: 0 < = k,l <N and 0 <= m <M
  float** Z;           // P-values transformed into normally distributed Z-values: Z = sqrt(2)*erfc^(-1)(2*P-value)
  float* Zbar;         // Zbar[k] = mean value of Z(k,m) (m=1,..,M)
  float* sigma;        // sigma[k] = sqrt(C[k][k]);
  float* w;            // weights of HMMs
  float** C;           // C[k][l] = cov( Z_k,Z_l )
  Hash<int> index(10000);  // index[name] = unique index between 0 and N-1 for each HMM
  index.Null(-1);      // value to be returned when the key is not defined
  Hash<int>* exclfold; // exclfold[k] is a hash that records the names of folds that are likely homologs of HMM k
  char** incl;         // incl[k][l]=0 iff there is member of same fold as l that has P(k,l)>Pmax; 1 otherwise
  char** name=NULL;    // name[k] = name (and description) of HMM k
  char** fold=NULL;;   // fold[k] = fold id (e.g. a.118) of HMM k
  char** sfam=NULL;;   // sfam[k] = superfamily id (e.g. a.118.8) of HMM k
  char** fam=NULL;;    // fam[k ] = family id (e.g. a.118.2) of HMM k
  Hash<int> famsz(10000);  // Hash counts number of HMMs in family
  Hash<int> sfamsz(10000); // Hash counts number of families in superfamily
  Hash<int> foldsz(10000); // Hash counts number of superfamilies in fold
  famsz.Null(0);        // Set int value to return when no data can be retrieved
  sfamsz.Null(0);       // Set int value to return when no data can be retrieved
  foldsz.Null(0);       // Set int value to return when no data can be retrieved

  // Parameter defaults
  strcpy(par.outfile,"");
  strcpy(par.scoresdir,"./");
  par.Pmax=50;    // folds with at least on member l with P(k,l)>Pmax are excluded from covariance calculation for k
  par.Pmin=0.01;
  par.fast=0;

  // Process command line options
  ProcessArguments(argc,argv);

  // Check input 
  if (!*par.outfile) {help(); cerr<<"Please specify an output file name, e.g. hhcov -o test.hhc\n"; exit(4);}
  if (par.scoresdir[strlen(par.scoresdir)-1]!='/') // if scoresdir does not end with '/' append '/'
    strcat(par.scoresdir,"/"); 


  // Open scores directory
  DIR *dirp;           // pointer to directory
  struct dirent *dirt; // file descriptor
  dirp = opendir(par.scoresdir);
  if (dirp==NULL) {cerr<<"Error opening "<<par.scoresdir<<": "; perror("opendir"); exit(2);}

  // Read all *.scores files and record HMM name in hash index[]
  if (v>=2) printf("Counting number of *.scores files in %s and number of HMMs in scores files ...\n",par.scoresdir);
  N=0; // count number of *.scores files
  while ((dirt = readdir(dirp))!=NULL)
    if (strstr(dirt->d_name,".scores") && dirt->d_name[strlen(dirt->d_name)-7]=='.' ) 
      {
	// Open scores file
	strcpy(scoresfile,par.scoresdir);
	strcat(scoresfile,dirt->d_name);
	if (v>=5) printf("%i %i %s\n",N,index.Size(),scoresfile);
	scoresf = fopen(scoresfile,"r");
	if (!scoresf) {cerr<<"\nError in hhcov: could not open \'"<<scoresfile<<"\'\n"; exit(2);}
	
	if (N==0) 
	  {
	    // Jump to first scores record and count number of scores records M
	    while (fgetline(line,LINELEN,scoresf))  // jump to first scores record
	      if (!strncmp(line,"TARGET ",7)) break;             
	    M=0;
	    while (fgetline(line,LINELEN,scoresf))
	      if (strscn(line)==NULL) continue;     // skip lines that contain only white space
	      else M++;
	    rewind(scoresf);
	    
	    // Reserve space 
	    name = new(char*[M]);
	    fam  = new(char*[M]);
	    sfam = new(char*[M]);
	    fold = new(char*[M]);
	  }
	
	// Read header
	while (fgetline(line,LINELEN,scoresf))
	  {
	    if (!strncmp(line,"NAME ",5)) strwrd(qname,line+4);      // copy first word after line+4 to qname
	    else if (!strncmp(line,"FAM  ",5)) strwrd(qfam,line+4);  // copy first word after line+4 to qfam
	    else if (!strncmp(line,"LENG ",5)) qlen=atoi(line+4);    // store integer after line+4 in qlen
	    else if (strscn(line)==NULL) continue;                   // skip lines that contain only white space
	    else if (!strncmp(line,"TARGET ",7)) break;              // => read probabilities
	  }
	k = index.Size();
	name[k] = new(char[NAMELEN]);
	fam[k]  = new(char[IDLEN]);
	sfam[k] = new(char[IDLEN]);
	fold[k] = new(char[IDLEN]);
	strcpy(name[k],qname);
	strcpy(fam[k],qfam);
	index.Add(name[k],k);     // record name of HMM and its index in hash index
	ScopID(qname, fold[k], sfam[k], fam[k]);
	if (famsz.Contains(fam[k])) (*famsz(fam[k]))++;
	else 
	  {
	    famsz.Add(fam[k],1); 
	    if (sfamsz.Contains(sfam[k])) (*sfamsz(sfam[k]))++;
	    else
	      {
		sfamsz.Add(sfam[k],1);
		if (foldsz.Contains(fold[k])) (*foldsz(fold[k]))++;
		else foldsz.Add(fold[k],1);
	      }
	  }
	fclose(scoresf);
	
	if (v>=4) printf("%i %-10.19s %s\n",k,name[k],fam[k]);
	if (v>=3 && N%100==0) printf("%i\n",N); 
	N++;
      }
  if (v>=3) printf("\n");
  if (v>=2) printf("Found a total of %i scores files with %i HMMs each\n",N,M);
  
  scoresf = fopen(scoresfile,"r");
  if (!scoresf) {cerr<<"\nError in hhcov: could not open \'"<<scoresfile<<"\'\n"; exit(2);}

  if (v>=5) printf("Reading HMMs and family codes in %s\n",scoresfile);
  scoresf = fopen(scoresfile,"r");
  if (!scoresf) {cerr<<"\nError in hhcov: could not open \'"<<scoresfile<<"\'\n"; exit(2);}
  
  // Go back to first scores record
  while (fgetline(line,LINELEN,scoresf))  // jump to first scores record
    if (!strncmp(line,"TARGET ",7)) break;             
  
  // Read HMM names and family from scores records
  while (fgetline(line,LINELEN,scoresf))
    {
      if (strscn(line)==NULL) continue; // skip lines that contain only white space
      
      // Record name and fold id (e.g. a.118) of HMM
      char* ptr;
      ptr = strwrd(tname,line); // read name 
      if (!index.Contains(tname))
	{
	  k = index.Size();
	  name[k] = new(char[NAMELEN]);
	  fam[k]  = new(char[IDLEN]);
	  sfam[k] = new(char[IDLEN]);
	  fold[k] = new(char[IDLEN]);
	  strcpy(name[k],tname);    // 
	  index.Add(name[k],k);     // record name of HMM and its index in hash index
	  strwrd(fam[k],ptr);       // read family
	  ScopID(tname, fold[k], sfam[k], fam[k]);
	  if (famsz.Contains(fam[k])) (*famsz(fam[k]))++;
	  else 
	    {
	      famsz.Add(fam[k],1); 
	      if (sfamsz.Contains(sfam[k])) (*sfamsz(sfam[k]))++;
	      else
		{
		  sfamsz.Add(sfam[k],1);
		  if (foldsz.Contains(fold[k])) (*foldsz(fold[k]))++;
		  else foldsz.Add(fold[k],1);
		}
	    }
//        strcat(strcat(strcpy(name[k],qname)," "),qfam); // d1zpda2 c.36.1.5 
	}
    }
  fclose(scoresf);
  rewinddir(dirp);

  // Reserve memory
  Z = new(float*[N]);
  for (k=0; k<N; k++) Z[k] = new(float[M]); // allocate memory for NxM matrix
  exclfold = new(Hash<int>[N]);
 
  // Calculate preliminary weights
  w = new(float[M]);
  int n = sfamsz.Size();
  for (k=0; k<M; k++) w[k] = 1.0/n/sfamsz[sfam[k]]/famsz[fam[k]];
  if (v>=3) // check if sum w[k] = 1
    {
      double one=0;
      for (k=0; k<M; k++) one += w[k];
      printf("Sum of w[k] = %6.4f (must be 1.0000) \n\n",one);
    }

  // For each file in directory
  if (v>=2) printf("Reading scores files ...\n");
  while ((dirt = readdir(dirp))!=NULL)
    {
      // Is file not a *.scores file?
      if (strstr(dirt->d_name,".scores")==NULL || dirt->d_name[strlen(dirt->d_name)-7]!='.') continue; 
      
      // Open scores file
      strcpy(scoresfile,par.scoresdir);
      strcat(scoresfile,dirt->d_name);
      if (v>=5) printf("%s\n",scoresfile);
      scoresf = fopen(scoresfile,"r");
      if (!scoresf) {cerr<<"\nError in hhcov: could not open \'"<<scoresfile<<"\'\n"; exit(2);}

      // Read header
      while (fgetline(line,LINELEN,scoresf))
	{
	  if (!strncmp(line,"NAME ",5)) strwrd(qname,line+4);      // copy first word after line+4 to qname
	  else if (!strncmp(line,"FAM  ",5)) strwrd(qfam,line+4);  // copy first word after line+4 to qfam
	  else if (!strncmp(line,"LENG ",5)) qlen=atoi(line+4);    // store integer after line+4 in qlen
	  else if (strscn(line)==NULL) continue;                   // skip lines that contain only white space
	  else if (!strncmp(line,"TARGET ",7)) break;              // => read probabilities
	}
      k = index[qname]; // set first HMM index
 
      // Read probabilities
      Hash<int> twice(10000);
      twice.Null(-1);      
      twice.RemoveAll();
      exclfold[k].New(10,-1);               // exclfold[k] will have 10 slots and null element -1

      while (fgetline(line,LINELEN,scoresf))
	{
	  char* ptr;
	  if (strscn(line)==NULL) continue; // skip lines that contain only white space
	  ptr = strwrd(tname,line);               // read name of database HMM
 	  if (twice[tname]==1) continue;
 	  twice.Add(tname,1);
	  m=index[tname];                   // returns -1 if tname is not yet contained in index
	  if (m<0) {
	    printf("Error: found scores file %s for HMM %s that does not seem to appear in scores records\n",scoresfile,tname);
	    exit(1);
	  } 
	  // TARGET     FAMILY   REL LEN COL LOG-PVA  S-AASS PROBAB
	  ptr = strwrd(qname,ptr);               // read FAMILY
	  ptr = strwrd(qname,ptr);               // read degree of RELatedness
	  ptr = strwrd(qname,ptr);               // read LENgth of HMM
	  ptr = strwrd(qname,ptr);               // read number of COLumns in alignment
	  float log2pval=fmin(fmax(atof(ptr),1E-10),1000.0); // read -log2( P-value) ...
	  Z[k][m] = Score2Z(log2pval);           // and transform to Z-value
	  if (isnan(Z[k][m])) 
	    {
	      printf("\nError: in %s: HMM %s Z[%i][%i]=nan!  value read in: %.5g  log2pval=%.5g\n",scoresfile,tname,k,m,atof(line+32),log2pval);
	      exit(4);
	    }
	  ptr = strwrd(qname,ptr);               // read name of database HMM
	  ptr = strwrd(qname,ptr);               // read name of database HMM
	  if (atof(ptr)>par.Pmax && fold[m]!=NULL && exclfold[k].Contains(fold[m])==0 ) 
	    {
	      // Exclude folds with Prob > Pmax (= 050%)
	      exclfold[k].Add(fold[m]);
	      if (v>=5) printf("Exclude from %s: %s\n",name[k],fold[m]);
	    }
	}
      fclose(scoresf);
      
      if (v>=5) 
	{
	  printf("%4i ",k);
	  printf("Name=%-10.10s ",name[k]);
	  printf("Fold=%-10.10s \n",fold[k]);
	}
      else if (v>=3 && k%100==0) printf("%i\n",k); 

    } // while ((dirt = readdir(dirp))!=NULL)

  if (index.Size()!=M && v>=1) 
    printf("WARNING: the total number of template HMMs is %i instead of %i\n",index.Size(),M);

  if (v>=3) 
    {
      printf("\n");
    }

  // Close scores directory
  if (closedir(dirp)==-1) { cerr<<"Error closing "<<par.scoresdir<<": "; perror("closedir"); exit(2);}
  
  // Symmetrize Z[k][m]
  for (k=0; k<N; k++)
    for (m=0; m<k; m++)
      Z[k][m] = 0.5*(Z[k][m]+Z[m][k]);

  // Set incl[k][m] and calculate Zbar and sigma
  incl = new(char*[N]);
  for (k=0; k<N; k++) incl[k] = new(char[M]); // allocate memory for NxM matrix
  Zbar = new(float[N]);
  sigma = new(float[N]); // vector of standard deviations
  for (k=0; k<N; k++)
    {
      double sumw=0.0; 
      double sumZ=0.0; 
      double sumZ2=0.0; 
      for (m=0; m<M; m++)
	{
	  if (fold[m]!=NULL && exclfold[k].Contains(fold[m])) 
	    incl[k][m]=0; 
	  else 
	    {
	      incl[k][m]=1; 
	      sumw  += w[m];
	      sumZ  += w[m]*Z[k][m];
	      sumZ2 += w[m]*Z[k][m]*Z[k][m];
	    }
	}
      Zbar[k] = sumZ/sumw;  // mean value of P(k,-)
      sigma[k] = sqrt((sumZ2/sumw-sumZ*sumZ/sumw/sumw));
      if (v>=3) printf("mu[%4i] = %6.3f  sigma[%4i] = %6.3f\n\n",k,Zbar[k],k,sigma[k]);
    }

  // Normalize Z
  for (k=0; k<N; k++)
    for (m=0; m<M; m++)
      Z[k][m] = (Z[k][m]-Zbar[k])/sigma[k];
  if (v>=4) PrintMatrix(Z,N,"Z.mat");

  // Allocate memory for NxN covariance matrix
  C = new(float*[N]);
  for (k=0; k<N; k++) C[k] = new(float[N]);

  // Calculate covariance matrix C(k,l) = sum_m [ Z_km*Z_lm ]
  if (v>=2) printf("Calculating covariance matrix ...\n");
  for (k=0; k<N; k++)
    { 
      if (v>=2 && (k+1)%20==0) fprintf(stderr,".",k+1); 
      if (v>=2 && (k+1)%1000==0) fprintf(stderr," %i\n",k+1); 
      for (l=k; l<N; l++)
	{
	  double sumw=0.0; 
	  double cov=0.0;
	  for (m=0; m<M; m++)
	    if (incl[k][m] && incl[l][m])
	      {
		sumw += w[m];
		cov  += w[m]*Z[k][m]*Z[l][m];
	      }
// 	  if (cov/psumw < par.rmin)  
// 	    C[k][l] = C[l][k] = 0.0;
// 	  else 
	  C[k][l] = C[l][k] = cov/sumw;
	}
    }
  if (v>=2) printf("\n");
  if (v>=4) PrintMatrix(C,N,"C.mat");
  
  // Write binary output file
  if (v>=2) printf("Writing binary output to %s ...\n",par.outfile);
  FILE* outf=NULL;
  outf = fopen(par.outfile,"wb");
  if (!outf) {cerr<<"\nError in hhcov: could not open \'"<<par.outfile<<"\'\n"; exit(2);}
  
  // Write number of HMMs
  fwrite(&N,sizeof(int),1,outf);

  // Write names in right order
  for (k=0; k<N; k++)
    fwrite(name[k],sizeof(char),IDLEN,outf);

  // Write Z-score matrix Z_kl
  for (k=0; k<N; k++)
    fwrite(Z[k]+k, sizeof(float), N-k, outf);

  // Write covariance matrix C_kl
  for (k=0; k<N; k++)
    fwrite(C[k]+k, sizeof(float), N-k, outf);

  // Write weights?
  if (v>=3) 
    {
      printf("<weight>\n");
      for (k=0; k<N; k++)
	printf(">%s\n%8.6f\n",name[k],w[k]);
      printf("</weight>\n");
    }

  fclose(outf);
  if (v>=2) printf("Done\n");
  exit(0);
}



/////////////////////////////////////////////////////////////////////////////////////
//// Processing input options from command line and .hhdefaults file
/////////////////////////////////////////////////////////////////////////////////////
void ProcessArguments(int argc, char** argv)
{
  //Processing command line input
  for (int i=1; i<argc; i++)
    { 
      if (v>=4) cout<<i<<"  "<<argv[i]<<endl; //PRINT
      if (!strcmp(argv[i],"-d"))
	{
	  if (++i>=argc || argv[i][0]=='-') 
	    {help() ; cerr<<"\nError in hhcov: no database file following -d\n"; exit(4);}
	  else 
	    strncpy(par.scoresdir,argv[i],NAMELEN);
	}
      else if (!strcmp(argv[i],"-o"))
	{
	  if (++i>=argc || argv[i][0]=='-') 
	    {help() ; cerr<<"\nError in hhcov: no output file following -o\n"; exit(4);}
	  else strcpy(par.outfile,argv[i]);
	}
      else if (!strcmp(argv[i],"-p") && (i<argc-1)) par.Pmin=atof(argv[++i]);
      else if (!strcmp(argv[i],"-v") && (i<argc-1) && argv[i+1][0]!='-' ) v=atoi(argv[++i]);
      else if (!strcmp(argv[i],"-v"))  v=2;
      else if (!strcmp(argv[i],"-f")) par.fast=1;
      else cerr<<"\nWARNING: Ignoring unknown option "<<argv[i]<<" ...\n";
      if (v>=4) cout<<i<<"  "<<argv[i]<<endl; //PRINT
    } // end of for-loop for command line input
  return;
}


/////////////////////////////////////////////////////////////////////////////////////
//// Takes family code (eg. a.1.2.3) and returns strings 'a', 'a.1', and 'a.1.2'
/////////////////////////////////////////////////////////////////////////////////////
inline void  ScopID(char cl[], char fold[], char sfam[], const char fam[])
{
  char* ptr;

  //get scop class ID 
  strcpy(cl,fam); 
  ptr = strchr(cl,'.');               //return adress of next '.' in name
  if(ptr) ptr[0]='\0';  

  //get scop fold ID
  strcpy(fold,fam); 
  ptr = strchr(fold,'.');             //return adress of next '.' in name
  if(ptr) ptr = strchr(ptr+1,'.');    //return adress of next '.' in name
  if(ptr) ptr[0]='\0';

  //get scop superfamily ID
  strcpy(sfam,fam); 
  ptr = strchr(sfam,'.');            //return adress of next '.' in name
  if(ptr) ptr = strchr(ptr+1,'.');   //return adress of next '.' in name
  if(ptr) ptr = strchr(ptr+1,'.');   //return adress of next '.' in name
  if(ptr) ptr[0]='\0';
  return;
}



/////////////////////////////////////////////////////////////////////////////////////
void PrintMatrix(float** V, int N, const char* str)
{
  FILE* outf = fopen(str,"w");
  int k,l;

  printf("\nPrint matrix %s ... ",str);
  for (k=0; k<N; k++)
    {
      fprintf(outf,"%4i >> ",k);
      for (l=0; l<N; l++)
	{
	  if (k==l) fprintf(outf,"   *:%7.1f ",V[k][l]);
	  else fprintf(outf,"%4i:%7.1f ",l,V[k][l]);
	}
      fprintf(outf,"\n");
    }
  fprintf(outf,"\n");
  fclose(outf);
  printf("finished printing\n");
}

/////////////////////////////////////////////////////////////////////////////////////
void ReadMatrix(float** V, int N, const char* str)
{
  FILE* inf = fopen(str,"r");
  if (inf==NULL) {perror("fopen"); exit(2);}
  int k,l,m;
  float val;

  printf("\nRead matrix %s ... \n",str);
  for (k=0; k<N; k++)
    {
      fscanf(inf,"%d >> ",&m);
      for (l=0; l<N; l++)
	{
	  if (k==l) fscanf(inf,"   *:%7f",&val);
	  else fscanf(inf,"%4d:%7f",&m,&val);
	  V[k][l]=float(val);
	}
      fscanf(inf,"\n");
    }
  fscanf(inf,"\n");
  fclose(inf);
  printf("finished reading\n");
}

/////////////////////////////////////////////////////////////////////////////////////
// Score2Z transforms the -log2(P-value) score into a Z-score for 0 < S
// Score2Z(S) = sqrt(2)*dierfc(2*2^(-S)), where dierfc is the inverse of error function
/////////////////////////////////////////////////////////////////////////////////////
float Score2Z(float S)
{
  double s, t, u, w, x, y, z;
  if (S<=0) return float(-100000);
  y = ( S>200 ? 0.0 : 2.0*pow(2.0,-S) );
  if (y > 1) 
    {
      z =  (S<1e-6? 2*0.69314718056*S : 2-y);
      z = 2 - y;
      w = 0.916461398268964 - log(z);
    }
  else 
    {
      z = y; 
      w = 0.916461398268964 - 0.69314718056*(1-S); // w = 0.916461398268964 - log(z);
    }

  u = sqrt(w);
  s = (log(u) + 0.488826640273108) / w;
  t = 1 / (u + 0.231729200323405);
  x = u * (1 - s * (s * 0.124610454613712 + 0.5)) - 
    ((((-0.0728846765585675 * t + 0.269999308670029) * t + 
       0.150689047360223) * t + 0.116065025341614) * t + 
     0.499999303439796) * t;
  t = 3.97886080735226 / (x + 3.97886080735226);
  u = t - 0.5;
  s = (((((((((0.00112648096188977922 * u + 
	       1.05739299623423047e-4) * u - 0.00351287146129100025) * u - 
	     7.71708358954120939e-4) * u + 0.00685649426074558612) * u + 
	   0.00339721910367775861) * u - 0.011274916933250487) * u - 
	 0.0118598117047771104) * u + 0.0142961988697898018) * u + 
       0.0346494207789099922) * u + 0.00220995927012179067;
  s = ((((((((((((s * u - 0.0743424357241784861) * u - 
		 0.105872177941595488) * u + 0.0147297938331485121) * u + 
	       0.316847638520135944) * u + 0.713657635868730364) * u + 
	     1.05375024970847138) * u + 1.21448730779995237) * u + 
	   1.16374581931560831) * u + 0.956464974744799006) * u + 
	 0.686265948274097816) * u + 0.434397492331430115) * u + 
       0.244044510593190935) * t - 
    z * exp(x * x - 0.120782237635245222);
  x += s * (x * s + 1);
  if (y > 1) {
    x = -x;
  }
    return float (1.41421356237*x);
}

