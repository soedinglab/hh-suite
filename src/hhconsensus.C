// hhconsensus.C: read A3M/FASTA file and calculate consensus sequence 
// Compile with g++ hhconsensus.C -o hhconsensus -O3 -static -fno-strict-aliasing -L/home/soeding/programs/electric-fence-2.1.13/

////#define WINDOWS
#define MAIN
#include <iostream>   // cin, cout, cerr
#include <fstream>    // ofstream, ifstream
#include <cstdio>     // printf
#include <algorithm>  // min,max
#include <stdlib.h>   // exit
#include <string.h>     // strcmp, strstr
#include <math.h>     // sqrt, pow
#include <limits.h>   // INT_MIN
#include <float.h>    // FLT_MIN
#include <ctype.h>    // islower, isdigit etc
#include <time.h>     // clock_gettime etc. (in realtime library (-lrt compiler option))
#include <errno.h>    // perror()
#include <cassert>
#include <stdexcept>

//#include <new>
//#include "efence.h"
//#include "efence.c"

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

/////////////////////////////////////////////////////////////////////////////////////
// Global variables 
/////////////////////////////////////////////////////////////////////////////////////
Alignment qali;              //Create an alignment 
HMM q;                       //Create a HMM with maximum of MAXRES match states
 
/////////////////////////////////////////////////////////////////////////////////////
// Help functions
/////////////////////////////////////////////////////////////////////////////////////
void help()
{
  printf("\n");
  printf("HHconsensus %s\n",VERSION_AND_DATE);
  printf("Calculate the consensus sequence for an A3M/FASTA input file.   \n");
  printf("%s",REFERENCE);
  printf("%s",COPYRIGHT);
  printf("\n");
  printf("Usage: %s -i <file> [options]                           \n",program_name);
  printf(" -i <file>     query alignment (A2M, A3M, or FASTA), or query HMM          \n");
  printf("\n");
  printf("Output options:                                                            \n");
  printf(" -s <file>     append consensus sequence in FASTA (default=<infile.seq>)   \n");
  printf(" -o <file>     write alignment with consensus sequence in A3M              \n");
  printf(" -oa3m <file>  same                                                        \n");
  printf(" -oa2m <file>  write alignment with consensus sequence in A2M              \n");
  printf(" -ofas <file>  write alignment with consensus sequence in FASTA            \n");
  printf(" -v <int>      verbose mode: 0:no screen output  1:only warings  2: verbose\n");
  printf("\n");
  printf("Filter input alignment (options can be combined):                         \n");
  printf(" -id   [0,100] maximum pairwise sequence identity (%%) (def=%i)   \n",par.max_seqid);
  printf(" -diff [0,inf[ filter most diverse set of sequences, keeping at least this    \n");
  printf("               many sequences in each block of >50 columns (def=%i)\n",par.Ndiff);
  printf(" -cov  [0,100] minimum coverage with query (%%) (def=%i) \n",par.coverage);
  printf(" -qid  [0,100] minimum sequence identity with query (%%) (def=%i) \n",par.qid);
  printf(" -qsc  [0,100] minimum score per column with query  (def=%.1f)\n",par.qsc);
//   printf(" -csc  [0,100] minimum score per column with core alignment (def=%-.2f)\n",par.coresc);
//   printf(" -qscc [0,100] minimum score per column of core sequence with query (def=%-.2f)\n",par.qsc_core);
  printf("\n");         
  printf("Input alignment format:                                                    \n");
  printf(" -M a2m        use A2M/A3M (default): upper case = Match; lower case = Insert;\n");         
  printf("               '-' = Delete; '.' = gaps aligned to inserts (may be omitted)   \n");
  printf(" -M first      use FASTA: columns with residue in 1st sequence are match states\n");
  printf(" -M [0,100]    use FASTA: columns with fewer than X%% gaps are match states   \n");
  printf("\n");    
  printf("Other options:                                                               \n");
  printf(" -def          read default options from ./.hhdefaults or <home>/.hhdefault. \n");
  printf("\n");    
  printf("Example: %s -i stdin -s stdout\n",program_name);
  printf("\n");    
}


/////////////////////////////////////////////////////////////////////////////////////
//// Processing input options from command line and .hhdefaults file
/////////////////////////////////////////////////////////////////////////////////////
void ProcessArguments(int argc,char** argv)
{
  // Read command line options
  for (int i=1; i<=argc-1; i++)
    { 
      if (v>=4) cout<<i<<"  "<<argv[i]<<endl; //PRINT
      if (!strcmp(argv[i],"-i"))
	{
	  if (++i>=argc || argv[i][0]=='-') 
	    {help(); cerr<<endl<<"Error in "<<program_name<<": no input file following -i\n"; exit(4);}
	  else strcpy(par.infile,argv[i]);
	}
      else if (!strcmp(argv[i],"-s"))
	{
	  if (++i>=argc) 
	    {help(); cerr<<endl<<"Error in "<<program_name<<": no output file following -s\n"; exit(4);}
	  else strcpy(par.outfile,argv[i]);
	}
      else if (!strcmp(argv[i],"-o"))
	{
	  par.outformat=3;
	  if (++i>=argc) 
	    {help(); cerr<<endl<<"Error in "<<program_name<<": no output file following -o\n"; exit(4);}
	  else strcpy(par.alnfile,argv[i]);
	}
      else if (!strcmp(argv[i],"-ofas"))
	{
	  par.outformat=1;
	  if (++i>=argc || argv[i][0]=='-') 
	    {help() ; cerr<<endl<<"Error in "<<program_name<<": no output file following -o\n"; exit(4);}
	  else strcpy(par.alnfile,argv[i]);
	}
      else if (!strcmp(argv[i],"-oa2m"))
	{
	  par.outformat=2;
	  if (++i>=argc || argv[i][0]=='-') 
	    {help() ; cerr<<endl<<"Error in "<<program_name<<": no output file following -o\n"; exit(4);}
	  else strcpy(par.alnfile,argv[i]);
	}
      else if (!strcmp(argv[i],"-oa3m"))
	{
	  par.outformat=3;
	  if (++i>=argc || argv[i][0]=='-') 
	    {help() ; cerr<<endl<<"Error in "<<program_name<<": no output file following -o\n"; exit(4);}
	  else strcpy(par.alnfile,argv[i]);
	}
      else if (!strcmp(argv[i],"-h")|| !strcmp(argv[i],"--help")) {help(); exit(0);}
      else if (!strcmp(argv[i],"-v") && (i<argc-1) && argv[i+1][0]!='-' ) v=atoi(argv[++i]);
      else if (!strcmp(argv[i],"-v0")) v=0;
      else if (!strcmp(argv[i],"-seq") && (i<argc-1))  par.nseqdis=atoi(argv[++i]); 
      else if (!strcmp(argv[i],"-name") && (i<argc-1)) {
	strncpy(q.name,argv[++i],NAMELEN-1); //copy longname to name...
	strncpy(q.longname,argv[i],DESCLEN-1);   //copy full name to longname
      }
      else if (!strcmp(argv[i],"-id") && (i<argc-1))   par.max_seqid=atoi(argv[++i]); 
      else if (!strcmp(argv[i],"-qid") && (i<argc-1))  par.qid=atoi(argv[++i]); 
      else if (!strcmp(argv[i],"-qsc") && (i<argc-1))  par.qsc=atof(argv[++i]); 
      else if (!strcmp(argv[i],"-cov") && (i<argc-1))  par.coverage=atoi(argv[++i]); 
      else if (!strcmp(argv[i],"-diff") && (i<argc-1)) par.Ndiff=atoi(argv[++i]); 
      else if (!strcmp(argv[i],"-qscc") && (i<argc-1)) par.qsc_core=atof(argv[++i]); 
      else if (!strcmp(argv[i],"-csc") && (i<argc-1))  par.coresc=atof(argv[++i]); 
      else if (!strcmp(argv[i],"-M") && (i<argc-1)) 
	if(!strcmp(argv[++i],"a2m") || !strcmp(argv[i],"a3m"))  par.M=1; 
	else if(!strcmp(argv[i],"first"))  par.M=3; 
	else if (argv[i][0]>='0' && argv[i][0]<='9') {par.Mgaps=atoi(argv[i]); par.M=2;}
	else cerr<<endl<<"WARNING: Ignoring unknown argument: -M "<<argv[i]<<"\n";
      else if (!strcmp(argv[i],"-Gonnet")) par.matrix=0; 
      else if (!strncmp(argv[i],"-BLOSUM",7) || !strncmp(argv[i],"-Blosum",7))
	{
	  if (!strcmp(argv[i]+7,"30")) par.matrix=30; 
	  else if (!strcmp(argv[i]+7,"40")) par.matrix=40; 
	  else if (!strcmp(argv[i]+7,"50")) par.matrix=50; 
	  else if (!strcmp(argv[i]+7,"65")) par.matrix=65; 
	  else if (!strcmp(argv[i]+7,"80")) par.matrix=80; 
	  else cerr<<endl<<"WARNING: Ignoring unknown option "<<argv[i]<<" ...\n";
	}
       else if (!strcmp(argv[i],"-pcm") && (i<argc-1)) par.pcm=atoi(argv[++i]); 
      else if (!strcmp(argv[i],"-pca") && (i<argc-1)) par.pca=atof(argv[++i]); 
      else if (!strcmp(argv[i],"-pcb") && (i<argc-1)) par.pcb=atof(argv[++i]); 
      else if (!strcmp(argv[i],"-pcc") && (i<argc-1)) par.pcc=atof(argv[++i]);  
      else if (!strcmp(argv[i],"-gapb") && (i<argc-1)) { par.gapb=atof(argv[++i]); if (par.gapb<=0.01) par.gapb=0.01;} 
      else if (!strcmp(argv[i],"-gapd") && (i<argc-1)) par.gapd=atof(argv[++i]); 
      else if (!strcmp(argv[i],"-gape") && (i<argc-1)) par.gape=atof(argv[++i]); 
      else if (!strcmp(argv[i],"-gapf") && (i<argc-1)) par.gapf=atof(argv[++i]); 
      else if (!strcmp(argv[i],"-gapg") && (i<argc-1)) par.gapg=atof(argv[++i]); 
      else if (!strcmp(argv[i],"-gaph") && (i<argc-1)) par.gaph=atof(argv[++i]); 
      else if (!strcmp(argv[i],"-gapi") && (i<argc-1)) par.gapi=atof(argv[++i]); 
      else if (!strcmp(argv[i],"-def")) par.readdefaultsfile=1; 

      else cerr<<endl<<"WARNING: Ignoring unknown option "<<argv[i]<<" ...\n";
      if (v>=4) cout<<i<<"  "<<argv[i]<<endl; //PRINT
    } // end of for-loop for command line input
}

/////////////////////////////////////////////////////////////////////////////////////
//// MAIN PROGRAM
/////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{  
  char* argv_conf[MAXOPT];     // Input arguments from .hhdefaults file (first=1: argv_conf[0] is not used)
  int argc_conf;               // Number of arguments in argv_conf 

  strcpy(par.infile,"");
  strcpy(par.outfile,"");
  strcpy(par.alnfile,"");
  
  //Default parameter settings
  SetDefaults();
  par.nseqdis=MAXSEQ-1;        // maximum number of sequences to be written 
  par.showcons=0;
  par.cons=1;
  par.Ndiff=0;
  par.max_seqid=100;
  par.coverage=0;
  par.pcm=0;    // no amino acid pseudocounts
  par.pca=0.0;  // no amino acid pseudocounts
  par.gapb=0.0; // no transition pseudocounts

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

  // Process command line options (they override defaults from .hhdefaults file)
  ProcessArguments(argc,argv);

  // Check command line input and default values
  if (!*par.infile) {help(); cerr<<endl<<"Error in "<<program_name<<": input file missing\n"; exit(4);}    

  // Get basename
  RemoveExtension(q.file,par.infile);  //Get basename of infile (w/o extension):

  // Outfile not given? Name it basename.hhm
  if (!*par.outfile && !*par.alnfile)
    {
      RemoveExtension(par.outfile,par.infile); 
      strcat(par.outfile,".seq");
    }

  // Set substitution matrix; adjust to query aa distribution if par.pcm==3
  SetSubstitutionMatrix();

  // Read input file (HMM, HHM, or alignment format), and add pseudocounts etc.
  ReadAndPrepare(par.infile, q);

  // Write consensus sequence to sequence file
  // Consensus sequence is calculated in hhalignment.C, Alignment::FrequenciesAndTransitions()
  if (*par.outfile) 
    {
      FILE* outf=NULL;
      if (strcmp(par.outfile,"stdout"))
	{
	  outf=fopen(par.outfile,"a");
	  if (!outf) OpenFileError(par.outfile);
	} 
      else
	outf = stdout;
      // ">name_consensus" -> ">name consensus"
      strsubst(q.sname[q.nfirst],"_consensus"," consensus");
      fprintf(outf,">%s\n%s\n",q.sname[q.nfirst],q.seq[q.nfirst]+1); 
      fclose(outf);
    }

  // Print A3M/A2M/FASTA output alignment
  if (*par.alnfile) 
    {
      HalfAlignment qa;
      int n = imin(q.n_display,par.nseqdis+(q.nss_dssp>=0)+(q.nss_pred>=0)+(q.nss_conf>=0)+(q.ncons>=0));
      qa.Set(q.name,q.seq,q.sname,n,q.L,q.nss_dssp,q.nss_pred,q.nss_conf,q.nsa_dssp,q.ncons);

      if      (par.outformat==1) qa.BuildFASTA();
      else if (par.outformat==2) qa.BuildA2M();
      else if (par.outformat==3) qa.BuildA3M();
      qa.Print(par.alnfile);   // print alignment to outfile
    }

  // Print 'Done!'
  printf("Done!\n");

  exit(0);
} //end main




