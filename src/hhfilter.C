// hhfilterC: filter alignment in a2m format with maximum sequence identity of match states and minimum coverage
// Compile:             g++ hhfilter.C -o hhfilter -O3 -fno-strict-aliasing 
// Compile with efence: g++ hhfilter.C -o hhfilter -g -O -lefence 
#define MAIN
#include <iostream>   // cin, cout, cerr
#include <fstream>    // ofstream, ifstream
#include <stdio.h>    // printf
#include <stdlib.h>   // exit
#include <string.h>     // strcmp, strstr
#include <math.h>     // sqrt, pow
#include <limits.h>   // INT_MIN
#include <float.h>    // FLT_MIN
#include <time.h>     // clock
#include <ctype.h>    // islower, isdigit etc
#include <cassert>

using std::cout;
using std::cerr;
using std::endl;
using std::ios;
using std::ifstream;
using std::ofstream;

#include "util.C"        // imax, fmax, iround, iceil, ifloor, strint, strscn, strcut, substr, uprstr, uprchr, Basename etc.
#include "list.C"        // list data structure
#include "hash.C"        // hash data structure
#include "hhdecl.C"      // Constants, global variables, struct Parameters
#include "hhutil.C"      // MatchChr, InsertChr, aa2i, i2aa, log2, fast_log2, ScopID, WriteToScreen,
#include "hhmatrices.C"  // BLOSUM50, GONNET, HSDM

// includes needed for context specific pseudocounts
#include "amino_acid.cpp"
#include "sequence.cpp"
#include "profile.cpp"
#include "cluster.cpp"
#include "simple_cluster.cpp"
#include "matrix.cpp"
#include "cs_counts.cpp"

#include "hhhit.h"       // class Hit
#include "hhalignment.h" // class Alignment
#include "hhhalfalignment.h" // class HalfAlignment
#include "hhfullalignment.h" // class FullAlignment
#include "hhhitlist.h"   // class Hit

#include "hhhmm.C"       // class HMM
#include "hhalignment.C" // class Alignment
#include "hhhit.C"       // class Hit
#include "hhhalfalignment.C" // class HalfAlignment
#include "hhfullalignment.C" // class FullAlignment
#include "hhhitlist.C"   // class HitList
#include "hhfunc.C"      // some functions common to hh programs

float Neff=0.0;            // target diversity

/////////////////////////////////////////////////////////////////////////////////////
// Exit function
/////////////////////////////////////////////////////////////////////////////////////
void help()
{
  printf("\n");
  printf("HHfilter %s\n",VERSION_AND_DATE);
  printf("Filter an alignment by maximum sequence identity of match states and minimum coverage\n");
  printf("%s",REFERENCE);
  printf("%s",COPYRIGHT);
  printf("\n");
  printf("Usage: %s -i infile -o outfile [options]                  \n",program_name);
  printf(" -i <file>     read input file in A3M/A2M or FASTA format                 \n");
  printf(" -o <file>     write to output file in A3M format                         \n");
  printf(" -a <file>     append to output file in A3M format                        \n");
  printf("\n");
  printf("Options:                                                                  \n");
  printf(" -v <int>      verbose mode: 0:no screen output  1:only warings  2: verbose\n");
  printf(" -id   [0,100] maximum pairwise sequence identity (%%) (def=%i)   \n",par.max_seqid);
  printf(" -diff [0,inf[ filter most diverse set of sequences, keeping at least this    \n");
  printf("               many sequences in each block of >50 columns (def=%i)\n",par.Ndiff);
  printf(" -cov  [0,100] minimum coverage with query (%%) (def=%i) \n",par.coverage);
  printf(" -qid  [0,100] minimum sequence identity with query (%%) (def=%i) \n",par.qid);
  printf(" -qsc  [0,100] minimum score per column with query  (def=%.1f)\n",par.qsc);
  printf(" -neff [1,inf] target diversity of alignment (default=off)\n");
  printf(" -def          read default options from ./.hhdefaults or <home>/.hhdefault. \n");
  printf("\n");         
  printf("Input alignment format:                                                    \n");
  printf(" -M a2m        use A2M/A3M (default): upper case = Match; lower case = Insert;\n");         
  printf("               '-' = Delete; '.' = gaps aligned to inserts (may be omitted)   \n");
  printf(" -M first      use FASTA: columns with residue in 1st sequence are match states\n");
  printf(" -M [0,100]    use FASTA: columns with fewer than X%% gaps are match states   \n");
  printf("                                                                          \n");
  printf("Example: %s -id 50 -i d1mvfd_.a2m -o d1mvfd_.fil.a2m          \n\n",program_name);
  cout<<endl;
}

/////////////////////////////////////////////////////////////////////////////////////
//// Processing input options from command line and .hhconfig file
/////////////////////////////////////////////////////////////////////////////////////
void ProcessArguments(int argc, char** argv)
{
  // Read command line options
  for (int i=1; i<=argc-1; i++)
    { 
      if (v>=4) cout<<i<<"  "<<argv[i]<<endl; //PRINT
      if (!strcmp(argv[i],"-i"))
	{
	  if (++i>argc-1 || argv[i][0]=='-') 
	    {cerr<<"Error in "<<program_name<<": no input file following -f\n"; exit(4);}
	  else strcpy(par.infile,argv[i]);
	}
      else if (!strcmp(argv[i],"-o"))
	{
	  par.append=0;
	  if (++i>argc-1) 
	    {cerr<<"Error in "<<program_name<<": no output file following -o\n"; exit(4);}
	  else strcpy(par.outfile,argv[i]);
	}
      else if (!strcmp(argv[i],"-a"))
	{
	  par.append=1;
	  if (++i>argc-1) 
	    {cerr<<"Error in "<<program_name<<": no output file following -a\n"; exit(4);}
	  else strcpy(par.outfile,argv[i]); 
	}
      else if (!strcmp(argv[i],"-v") && (i+1<argc) && argv[i+1][0]!='-' ) v=atoi(argv[++i]);
      else if (!strcmp(argv[i],"-v0")) v=0;
      else if (!strcmp(argv[i],"-v1")) v=1;
      else if (!strcmp(argv[i],"-v2")) v=2;
      else if (!strcmp(argv[i],"-v"))  v=2;
      else if (!strcmp(argv[i],"-v3")) v=3;
      else if (!strcmp(argv[i],"-v4")) v=4;
      else if (!strcmp(argv[i],"-v5")) v=5;
      else if (!strcmp(argv[i],"-id") && (i<argc-1))   par.max_seqid=atoi(argv[++i]); 
      else if (!strcmp(argv[i],"-qid") && (i<argc-1))  par.qid=atoi(argv[++i]); 
      else if (!strcmp(argv[i],"-qsc") && (i<argc-1))  par.qsc=atof(argv[++i]); 
      else if (!strcmp(argv[i],"-cov") && (i<argc-1))  par.coverage=atoi(argv[++i]); 
      else if (!strcmp(argv[i],"-diff") && (i<argc-1)) par.Ndiff=atoi(argv[++i]); 
      else if (!strcmp(argv[i],"-neff") && (i<argc-1)) Neff=atof(argv[++i]); 
      else if (!strcmp(argv[i],"-Neff") && (i<argc-1)) Neff=atof(argv[++i]); 
      else if (!strcmp(argv[i],"-M") && (i<argc-1)) 
	if (!strcmp(argv[++i],"a2m") || !strcmp(argv[i],"a3m"))  par.M=1; 
	else if(!strcmp(argv[i],"first"))  par.M=3; 
	else if (argv[i][0]>='0' && argv[i][0]<='9') {par.Mgaps=atoi(argv[i]); par.M=2;}
	else cerr<<endl<<"WARNING: Ignoring unknown argument: -M "<<argv[i]<<"\n";

      else if (!strcmp(argv[i],"-def")) par.readdefaultsfile=1; 
      else if (!strcmp(argv[i],"-h")|| !strcmp(argv[i],"--help")) {help(); exit(0);}
      else cerr<<endl<<"WARNING: Ignoring unknown option "<<argv[i]<<" ...\n";
      if (v>=4) cout<<i<<"  "<<argv[i]<<endl; //PRINT
    } // end of for-loop for command line input
}

float filter_by_qsc(float qsc, Alignment& qali, char* dummy)
{
  HMM q;
  for (int k=0; k<qali.N_in; k++) qali.keep[k]=dummy[k];
  qali.Filter2(qali.keep,par.coverage,0,qsc,par.max_seqid+1,par.max_seqid,0); 
  qali.FrequenciesAndTransitions(q);
//   printf("qsc=%4.1f  N_filtered=%-3i  Neff=%6.3f\n",qsc,n,q.Neff_HMM);
  return q.Neff_HMM;
} 

/////////////////////////////////////////////////////////////////////////////////////
//// MAIN PROGRAM
/////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{  
  Alignment qali;              //Create an alignment 

  char* argv_conf[MAXOPT];        // Input arguments from .hhconfig file (first=1: argv_conf[0] is not used)
  int argc_conf;                  // Number of arguments in argv_conf 

  strcpy(par.infile,"");
  strcpy(par.outfile,"");
  
  SetDefaults();
  par.nseqdis=MAXSEQ-1;        // maximum number of sequences to be written 
  par.Ndiff=0;                 // no filtering for maximum diversity

  // Make command line input globally available
  par.argv=argv; 
  par.argc=argc;
  RemovePathAndExtension(program_name,argv[0]);
  
  // Enable changing verbose mode before defaults file and command line are processed
  for (int i=1; i<argc; i++)
    { 
      if (!strcmp(argv[i],"-def")) par.readdefaultsfile=1;
      else if (argc>1 && !strcmp(argv[i],"-v0")) v=0;
      else if (argc>1 && !strcmp(argv[i],"-v1")) v=1;
      else if (argc>2 && !strcmp(argv[i],"-v")) v=atoi(argv[i+1]);
    }

  // Read .hhdefaults file?
  if (par.readdefaultsfile) 
    {
      // Process default otpions from .hhconfig file
      ReadDefaultsFile(argc_conf,argv_conf);
      ProcessArguments(argc_conf,argv_conf);
    }

  // Process command line options (they override defaults from .hhconfig file)
  ProcessArguments(argc,argv);

  // Check command line input and default values
  if (!*par.infile) {help(); cerr<<"Error: input file missing\n"; exit(4);}   
  if (!*par.outfile) {help(); cerr<<"Error: output file missing\n"; exit(4);}   
 
  if (v>=2) 
    {
      cout<<"Input file = "<<par.infile<<"\n";
      cout<<"Output file = "<<par.outfile<<"\n";
    }

  // Reads in an alignment from par.infile into matrix X[k][l] as ASCII
  FILE* inf = fopen(par.infile, "r");
  if (!inf) OpenFileError(par.infile);
  qali.Read(inf,par.infile);
  fclose(inf);

  // Convert ASCII to int (0-20),throw out all insert states, record their number in I[k][i] 
  // and store marked sequences in name[k] and seq[k]
  qali.Compress(par.infile);

  // Filter by minimum score per column with query sequence?
  if (0>-10) SetSubstitutionMatrix();

  // Remove sequences with seq. identity larger than seqid percent (remove the shorter of two)
  qali.N_filtered = qali.Filter(par.max_seqid,par.coverage,par.qid,par.qsc,par.Ndiff);

  // Atune alignment diversity q.Neff with qsc to value Neff_goal
  if (Neff>=0.999) 
    {
      int v1=v;
      v=1;
      const float TOLX=0.001; 
      const float TOLY=0.02; 
      char dummy[qali.N_in+1];   
      for (int k=0; k<qali.N_in; k++) dummy[k]=qali.keep[k];
      float x=0.0,y=0.0;
      float x0=-1.0;
      float x1=+2.0;
      float y0=filter_by_qsc(x0,qali,dummy);
      float y1=filter_by_qsc(x1,qali,dummy);
      int i=2;
      while (y0-Neff>0 && Neff-y1>0)
	{
	  x = x0 + (Neff-y0)*(x1-x0)/(y1-y0); // linear interpolation between (x0,y0) and (x1,y1)
	  y = filter_by_qsc(x,qali,dummy);
	  if (v>=0) printf(" %3i  x0=%6.3f -> %6.3f     x=%6.3f -> %6.3f     x1=%6.3f -> %6.3f \n",++i,x0,y0,x,y,x1,y1);
	  if (y>Neff) {x0=x; y0=y;} else {x1=x; y1=y;}
	  if (fabs(Neff-y)<TOLY || x1-x0<TOLX) break;
	}
      v=v1;
      if (y0>=Neff && y1<=Neff) 
	{
	  // Write filtered alignment WITH insert states (lower case) to alignment file
	  if (v>=2) printf("Found Neff=%6.3f at filter threshold qsc=%6.3f\n",y,x);
	  qali.WriteToFile(par.outfile);
	}
      else if (v>=1) 
	printf("Diversity of unfiltered alignment %.2f is below target diversity %.2f. No alignment written\n",y0,Neff);
    }
  else 
    // Write filtered alignment WITH insert states (lower case) to alignment file
    qali.WriteToFile(par.outfile);

  // Print 'Done!'
  FILE* outf=NULL;
  if (!strcmp(par.outfile,"stdout"))
    printf("Done!\n");
  else
    {
      if (!*par.outfile)
	{
	  outf=fopen(par.outfile,"a"); //open for append
	  fprintf(outf,"Done!\n");
	  fclose(outf);
	}
      if (v>=2) printf("Done\n");
    }
 
  exit(0);

} //end main




