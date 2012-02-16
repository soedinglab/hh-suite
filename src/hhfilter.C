// hhfilterC: filter alignment in a2m format with maximum sequence identity of match states and minimum coverage
// 
//     (C) Johannes Soeding 2012

//     This program is free software: you can redistribute it and/or modify
//     it under the terms of the GNU General Public License as published by
//     the Free Software Foundation, either version 3 of the License, or
//     (at your option) any later version.

//     This program is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     GNU General Public License for more details.

//     You should have received a copy of the GNU General Public License
//     along with this program.  If not, see <http://www.gnu.org/licenses/>.

//     We are very grateful for bug reports! Please contact us at soeding@genzentrum.lmu.de

//     Reference: 
//     Remmert M., Biegert A., Hauser A., and Soding J.
//     HHblits: Lightning-fast iterative protein sequence searching by HMM-HMM alignment.
//     Nat. Methods, epub Dec 25, doi: 10.1038/NMETH.1818 (2011).

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

#ifdef HH_SSE4
#include <tmmintrin.h>   // SSSE3
#include <smmintrin.h>   // SSE4.1
#define HH_SSE3
#endif

#ifdef HH_SSE3
#include <pmmintrin.h>   // SSE3
#define HH_SSE2
#endif

#ifdef HH_SSE2
#ifndef __SUNPRO_C
#include <emmintrin.h>   // SSE2
#else
#include <sunmedia_intrin.h>
#endif
#endif

using std::cout;
using std::cerr;
using std::endl;
using std::ios;
using std::ifstream;
using std::ofstream;

#include "cs.h"          // context-specific pseudocounts
#include "context_library.h"
#include "library_pseudocounts-inl.h"

#include "list.h"        // list data structure
#include "hash.h"        // hash data structure

#include "util.C"        // imax, fmax, iround, iceil, ifloor, strint, strscn, strcut, substr, uprstr, uprchr, Basename etc.
#include "hhdecl.C"      // Constants, global variables, struct Parameters
#include "hhutil.C"      // MatchChr, InsertChr, aa2i, i2aa, log2, fast_log2, ScopID, WriteToScreen,
#include "hhmatrices.C"  // BLOSUM50, GONNET, HSDM

#include "hhhmm.h"       // class HMM
#include "hhhit.h"       // class Hit
#include "hhalignment.h" // class Alignment

#include "hhhmm.C"       // class HMM
#include "hhalignment.C" // class Alignment
#include "hhfunc.C"      // some functions common to hh programs

/////////////////////////////////////////////////////////////////////////////////////
// Exit function
/////////////////////////////////////////////////////////////////////////////////////
void help()
{
  printf("\n");
  printf("HHfilter %s\n",VERSION_AND_DATE);
  printf("Filter an alignment by maximum sequence identity of match states and minimum coverage\n");
  printf("%s",COPYRIGHT);
  printf("%s",REFERENCE);
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
      else if (!strcmp(argv[i],"-neff") && (i<argc-1)) par.Neff=atof(argv[++i]); 
      else if (!strcmp(argv[i],"-Neff") && (i<argc-1)) par.Neff=atof(argv[++i]); 
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

  par.SetDefaultPaths(program_path);

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
  if (par.Neff>=1.0) qali.FilterNeff();

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




