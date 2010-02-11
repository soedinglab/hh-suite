// hhcorr.C:
// Measure the correlation of p-values for two input HMMs with respect to a representative HMM database
// compile with g++ hhcorr.C -o hhcorr -O3 -g -I/usr/include/ -L/usr/lib -lpng -lz -static -fno-strict-aliasing


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

#include "pngwriter.h"   //PNGWriter (http://pngwriter.sourceforge.net/)
#include "pngwriter.cc"  //PNGWriter (http://pngwriter.sourceforge.net/)


/////////////////////////////////////////////////////////////////////////////////////
// Global variables
/////////////////////////////////////////////////////////////////////////////////////
HMM q1,q2;                   // Create two query  HMMs
HMM t;                       // Create template HMM
Hit hit;                     // Ceate new hit object pointed at by hit
HitList hitlist1,hitlist2;   // list of hits with one Hit object for each pairwise comparison done
char corr_q2_vs_HMMs_like_PSIBLAST=0;

/////////////////////////////////////////////////////////////////////////////////////
// Help functions
/////////////////////////////////////////////////////////////////////////////////////
void help()
{
  printf("\n");
  printf("HHcorr %s\n",VERSION_AND_DATE);
  printf("Measure the correlation of p-values for two input HMMs with respect to\n");
  printf("a representative HMM database\n");
  printf("%s",REFERENCE);
  printf("%s",COPYRIGHT);
  printf("\n");
  printf("Usage %s -i query1 query2 -d database  [options]            \n",program_name);
  printf(" -i file file  input query alignments/HMMs              \n");
  printf(" -d file       input HMM database                                  \n");
  printf("\n");
  printf("Example: %s -i a.1.1.1.1_10.ext a.1.1.1.1_10.a3m -d cal.hhm \n",program_name);
  printf("\n");
  printf("More help:                                                         \n");
  printf(" -h out        options for formatting ouput                        \n");
  printf(" -h hmm        options for building HMM from multiple alignment    \n");
  printf(" -h gap        options for setting gap penalties                   \n");
  printf(" -h ali        options for HMM-HMM alignment                       \n");
  printf(" -h all        all options \n");
  cout<<endl;
}

void help_out()
{
  printf("\n");
  printf("Output options:                                                           \n");
  printf(" -v            verbose mode (default: show only warnings)                 \n");
  printf(" -v 0          suppress all screen output                                 \n");
  printf(" -o file       write output to file (def=%s)\n",par.outfile);
  printf(" -png file     write dotplot into PNG-file                         \n");
}

void help_hmm()
{
  printf("\n");
  printf("Filter input alignment (options can be combined):                         \n");
  printf(" -id   [0,100] maximum pairwise sequence identity (%%) (def=%i)   \n",par.max_seqid);
  printf(" -diff [0,inf[ filter most diverse set of sequences, keeping at least this    \n");
  printf("               many sequences in each block of >50 columns (def=%i)\n",par.Ndiff);
  printf(" -cov  [0,100] minimum coverage with query (%%) (def=%i) \n",par.coverage);
  printf(" -qid  [0,100] minimum sequence identity with query (%%) (def=%i) \n",par.qid);
  printf(" -qsc  [0,100] minimum score per column with query  (def=%.1f)\n",par.qsc);
  printf("                                                                          \n");
  printf("HMM-building options:                                                     \n");
  printf(" -M a2m        use A2M/A3M (default): upper case = Match; lower case = Insert;\n");
  printf("               '-' = Delete; '.' = gaps aligned to inserts (may be omitted)   \n");
  printf(" -M first      use FASTA: columns with residue in 1st sequence are match states\n");
  printf(" -M [0,100]    use FASTA: columns with fewer than X%% gaps are match states   \n");
  printf("                                                                          \n");
  printf("Pseudocount options:                                                      \n");
  printf(" -Gonnet       use the Gonnet substitution matrix (default)               \n");
  printf(" -Blosum50     use the Blosum50 substitution matrix                       \n");
  printf(" -Blosum62     use the Blosum62 substitution matrix                       \n");
  printf(" -HSDM         use the structure-derived HSDM substitution matrix         \n");
  printf(" -pcm  0-3     Pseudocount mode (default=%-i)                             \n",par.pcm);
  printf("               tau = substitution matrix pseudocount admixture            \n");
  printf("               0: no pseudo counts:     tau = 0                           \n");
  printf("               1: constant              tau = a                           \n");
  printf("               2: divergence-dependent: tau = a/(1 + ((Neff-1)/b)^c)       \n");
  printf("                  Neff=( (Neff_q^d+Neff_t^d)/2 )^(1/d)                       \n");
  printf("                  Neff_q = av number of different AAs per column in query  \n");
  printf("               3: column-specific:      tau = \'2\' * (Neff(i)/Neff)^(w*Neff/20)\n");
  printf(" -pca  [0,1]   set a (overall admixture) (def=%-.1f)                      \n",par.pca);
  printf(" -pcb  [1,inf[ set b (threshold for Neff) (def=%-.1f)                      \n",par.pcb);
  printf(" -pcc  [0,3]   set c (extinction exponent for tau(Neff))  (def=%-.1f)      \n",par.pcc);
  printf(" -pcw  [0,3]   set w (weight of pos-specificity for pcs) (def=%-.1f)      \n",par.pcw);
}

void help_gap()
{
  printf("\n");
  printf("Gap cost options:                                                         \n");
  printf(" -gapb [0,inf[ transition pseudocount admixture (def=%-.1f)               \n",par.gapb);
  printf(" -gapd [0,inf[ Transition pseudocount admixture for opening gap (default=%-.1f)\n",par.gapd);
  printf(" -gape [0,1.5] Transition pseudocount admixture for extending gap (def=%-.1f)\n",par.gape);
  printf(" -gapf ]0,inf] factor for increasing/reducing the gap open penalty for deletes (def=%-.1f)\n",par.gapf);
  printf(" -gapg ]0,inf] factor for increasing/reducing the gap open penalty for deletes (def=%-.1f)\n",par.gapg);
  printf(" -gaph ]0,inf] factor for increasing/reducing the gap extension penalty for deletes(def=%-.1f)\n",par.gaph);
  printf(" -gapi ]0,inf] factor for increasing/reducing the gap extension penalty for inserts(def=%-.1f)\n",par.gapi);
  printf(" -egq  [0,inf[ penalty (bits) for end gaps aligned to query residues (def=%-.1f)\n",par.egq);
  printf(" -egt  [0,inf[ penalty (bits) for end gaps aligned to template residues (def=%-.1f)\n",par.egt);
  //  printf(" -eg   [0,inf[ penalty (bits) for end gaps aligned to shorter of two HMMs (def=off)\n");
}

void help_ali()
{
  printf("\n");
  printf("Alignment options:  \n");
  printf(" -glob/-loc    global or local alignment mode (def=global)         \n");
  printf(" -sc   int     column score S(i,j)    (tja=template HMM at column j) (def=%i)\n",par.columnscore);
  printf("        0      = log2 Sum(tja*qia/pa)   (pa: aa background frequencies)   \n");
  printf("        1      = log2 Sum(tja*qia/ra)   (ra: optimized frequencies)       \n");
  printf("        2      = log2 Sum(tja*qia/ta)   (ta: av. aa freqs in template)      \n");
  printf("        3      = log2 Sum(tja*qia/qa)   (qa: av. aa freqs in query)       \n");
  printf("        4      = log2 Sum(tja*qia/pa) - log2 Sum(tja*qa/pa)               \n");
  printf("        6      = log2 Sum(tja*qia/ta) - log2 Sum(tja*qa/ta)               \n");
  printf("        7      = log2 Sum(tja*qia/qa) - log2 Sum(ta*qia/qa)               \n");
  printf("        8      = log2 Sum(tja*qia/ra) - log2 Sum(tja*qa/ra)               \n");
  printf("               - log2 Sum(ta *qia/ra) + log2 Sum(ta *qa/ra)               \n");
  printf("        9      Sum of pairs score (substitution matrix admixture a)       \n");
  printf(" -corr [0,1]   weight of term for pair correlations                       \n");
  printf(" -shift [-1,1] score offset (def=%-.3f)                                   \n",par.shift);
  printf(" -r            repeat identification: multiple hits not treated as independent\n");
  printf(" -ssm  0-2     0:no ss scoring [default=%i]               \n",par.ssm);
  printf("               1:ss scoring after alignment                               \n");
  printf("               2:ss scoring during alignment                              \n");
  printf("               3:ss scoring after alignment; use only psipred (not dssp)\n");
  printf("               4:ss scoring during alignment use onlu psipred (not dssp)\n");
  printf(" -ssw  [0,1]   weight of ss score compared to column score (def=%-.2f)    \n",par.ssw);
  printf(" -ssa  [0,1]   ss confusion matrix = (1-ssa)*I + ssa*psipred-confusion-matrix [def=%-.2f)\n",par.ssa);
}

void help_all()
{
  help();
  help_out();
  help_hmm();
  help_gap();
  help_ali();
  printf("\n");
  printf(" -def          read default options from ./.hhdefaults or HOME/.hhdefault. \n");
}

/////////////////////////////////////////////////////////////////////////////////////
//// Processing input options from command line and .hhdefaults file
/////////////////////////////////////////////////////////////////////////////////////
void ProcessArguments(int argc, char** argv)
{
  //Processing command line input
  for (int i=1; i<argc; i++)
    {
      if (v>=4) cerr<<i<<"  "<<argv[i]<<endl; //PRINT
      if (!strcmp(argv[i],"-i"))
        {
          if (++i>=argc || argv[i][0]=='-')
            {help() ; cerr<<endl<<"Error in "<<program_name<<": no input file following -i\n"; exit(4);}
          else strcpy(par.infile,argv[i]);
          if (++i>=argc || argv[i][0]=='-')
            {help() ; cerr<<endl<<"Error in "<<program_name<<": -i needs two input files as arguments\n"; exit(4);}
          else strcpy(par.tfile,argv[i]);
        }
      else if (!strcmp(argv[i],"-d"))
        {
          if (++i>=argc || argv[i][0]=='-')
            {help() ; cerr<<endl<<"Error in "<<program_name<<": no database file following -d\n"; exit(4);}
          else
	    {
              par.dbfiles = new(char[strlen(argv[i])+1]);
              strcpy(par.dbfiles,argv[i]);
            }
        }
      else if (!strcmp(argv[i],"-o"))
        {
          par.append=0;
          if (++i>=argc || argv[i][0]=='-')
            {help() ; cerr<<endl<<"Error in "<<program_name<<": no output file following -o\n"; exit(4);}
          else strcpy(par.outfile,argv[i]);
        }
      else if (!strcmp(argv[i],"-png"))
        {
          if (++i>=argc || argv[i][0]=='-')
            {help(); cerr<<endl<<"Error in "<<program_name<<": no output file following -p\n"; exit(4);}
          else {pngfile = new(char[strlen(argv[i])+1]); strcpy(pngfile,argv[i]);}
        }
       else if (!strcmp(argv[i],"-h"))
         {
           if (++i>=argc || argv[i][0]=='-') {help(); exit(0);}
           if (!strcmp(argv[i],"out")) {help_out(); exit(0);}
           if (!strcmp(argv[i],"hmm")) {help_hmm(); exit(0);}
           if (!strcmp(argv[i],"gap")) {help_gap(); exit(0);}
           if (!strcmp(argv[i],"ali")) {help_ali(); exit(0);}
           else {help(); exit(0);}
         }
      else if (!strcmp(argv[i],"-v") && (i<argc-1) && argv[i+1][0]!='-' ) v=atoi(argv[++i]);
      else if (!strcmp(argv[i],"-v0")) v=0;
      else if (!strcmp(argv[i],"-v1")) v=1;
      else if (!strcmp(argv[i],"-v2")) v=2;
      else if (!strcmp(argv[i],"-v"))  v=2;
      else if (!strcmp(argv[i],"-v3")) v=3;
      else if (!strcmp(argv[i],"-v4")) v=4;
      else if (!strcmp(argv[i],"-v5")) v=5;
      else if (!strcmp(argv[i],"-id") && (++i<=argc-1))   par.max_seqid=atoi(argv[i]);
      else if (!strcmp(argv[i],"-qid") && (++i<=argc-1))  par.qid=atoi(argv[i]);
      else if (!strcmp(argv[i],"-qsc") && (++i<=argc-1))  par.qsc=atof(argv[i]);
      else if (!strcmp(argv[i],"-cov") && (++i<=argc-1))  par.coverage=atoi(argv[i]);
      else if (!strcmp(argv[i],"-diff") && (++i<=argc-1)) par.Ndiff=atoi(argv[i]);
      else if (!strcmp(argv[i],"-Gonnet")) par.matrix=0;
      else if (!strcmp(argv[i],"-HSDM")) par.matrix=1;
      else if (!strcmp(argv[i],"-BLOSUM50")) par.matrix=2;
      else if (!strcmp(argv[i],"-Blosum50")) par.matrix=2;
      else if (!strcmp(argv[i],"-B50")) par.matrix=2;
      else if (!strcmp(argv[i],"-BLOSUM62")) par.matrix=3;
      else if (!strcmp(argv[i],"-Blosum62")) par.matrix=3;
      else if (!strcmp(argv[i],"-B62")) par.matrix=3;
      else if (!strcmp(argv[i],"-pcm") && (i<argc-1)) par.pcm=atoi(argv[++i]);
      else if (!strcmp(argv[i],"-pca") && (i<argc-1)) par.pca=atof(argv[++i]);
      else if (!strcmp(argv[i],"-pcb") && (i<argc-1)) par.pcb=atof(argv[++i]);
      else if (!strcmp(argv[i],"-pcc") && (i<argc-1)) par.pcc=atof(argv[++i]);
      else if (!strcmp(argv[i],"-pcw") && (i<argc-1)) par.pcw=atof(argv[++i]);
      else if (!strcmp(argv[i],"-gapb") && (i<argc-1)) { par.gapb=atof(argv[++i]); if (par.gapb<=0.01) par.gapb=0.01;}
      else if (!strcmp(argv[i],"-gapd") && (i<argc-1)) par.gapd=atof(argv[++i]);
      else if (!strcmp(argv[i],"-gape") && (i<argc-1)) par.gape=atof(argv[++i]);
      else if (!strcmp(argv[i],"-gapf") && (i<argc-1)) par.gapf=atof(argv[++i]);
      else if (!strcmp(argv[i],"-gapg") && (i<argc-1)) par.gapg=atof(argv[++i]);
      else if (!strcmp(argv[i],"-gaph") && (i<argc-1)) par.gaph=atof(argv[++i]);
      else if (!strcmp(argv[i],"-gapi") && (i<argc-1)) par.gapi=atof(argv[++i]);
      //      else if (!strcmp(argv[i],"-eg") && (i<argc-1))  par.eg =atof(argv[++i]);
      else if (!strcmp(argv[i],"-egq") && (i<argc-1)) par.egq=atof(argv[++i]);
      else if (!strcmp(argv[i],"-egt") && (i<argc-1)) par.egt=atof(argv[++i]);
      else if (!strcmp(argv[i],"-ssm") && (i<argc-1)) par.ssm=atoi(argv[++i]);
      else if (!strcmp(argv[i],"-ssw") && (i<argc-1)) par.ssw=atof(argv[++i]);
      else if (!strcmp(argv[i],"-ssa") && (i<argc-1)) par.ssa=atof(argv[++i]);
      else if (!strncmp(argv[i],"-glo",3)) par.loc=0;
      else if (!strncmp(argv[i],"-loc",3)) par.loc=1;
      else if (!strcmp(argv[i],"-r")) par.repmode=1;
      else if (!strcmp(argv[i],"-M") && (i<argc-1))
        if (!strcmp(argv[++i],"a2m") || !strcmp(argv[i],"a3m"))  par.M=1;
        else if(!strcmp(argv[i],"first"))  par.M=3;
        else if (argv[i][0]>='0' && argv[i][0]<='9') {par.Mgaps=atoi(argv[i]); par.M=2;}
        else cerr<<endl<<"WARNING: Ignoring unknown argument: -M "<<argv[i]<<"\n";
      else if (!strncmp(argv[i],"-cal",3)) par.calibrate=1;
      else if (!strcmp(argv[i],"-shift") && (i<argc-1)) par.shift=atof(argv[++i]);
      else if (!strcmp(argv[i],"-sc") && (i<argc-1)) par.columnscore=atoi(argv[++i]);
      else if (!strcmp(argv[i],"-corr") && (i<argc-1)) par.corr=atof(argv[++i]);
      //else if (!strcmp(argv[i],"-best") && (i<argc-1)) par.best=atof(argv[++i]);
      else if (!strcmp(argv[i],"-def")) par.readdefaultsfile=1;
      else if (!strcmp(argv[i],"-psi2")) corr_q2_vs_HMMs_like_PSIBLAST=1;
      else cerr<<endl<<"WARNING: Ignoring unknown option "<<argv[i]<<" ...\n";
      if (v>=4) cerr<<i<<"  "<<argv[i]<<endl; //PRINT
    } // end of for-loop for command line input
}



// Prepare query and template HMM for comparison
void Prepare(HMM& q, HMM& t, char pcm=par.pcm)
{
  // Generate an amino acid frequency matrix from t.f[i][a] with max pseudocounts (tau=1) -> t.g[i][a]
  t.PreparePseudocounts();

  // Add amino acid pseudocounts to template: t.p[i][a] = (1-tau)*f[i][a] + tau*g[i][a]
  char pcm_in=par.pcm;
  par.pcm = pcm;
  t.AddAminoAcidPseudocounts();
  t.CalculateAminoAcidBackground();
  par.pcm = pcm_in;

  // Factor Null model into HMM t
  // ATTENTION! t.p[i][a] is divided by pnul[a] (for reasons of efficiency) => do not reuse t.p
  t.IncludeNullModelInHMM(q,t);

  return;
}


/////////////////////////////////////////////////////////////////////////////////////
//// MAIN PROGRAM
/////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
  char readstatus = 0;         // 1:ok.continue. 2:Skip current HMM. 0:stop reading
  int N_searched;              // Number of HMMs searched
  char* argv_conf[MAXOPT];     // Input arguments from .hhdefaults file (first=1: argv_conf[0] is not used)
  int argc_conf;               // Number of arguments in argv_conf
  char in1ext[IDLEN];          // Extension of input file (hhm or a3m)
  char in2ext[IDLEN];          // Extension of input file (hhm or a3m)
  char* pngfile=NULL;

  SetDefaults();
  par.p=0.0;                   // minimum threshold for inclusion in hit list and alignment listing
  par.E=1e6;                   // maximum threshold for inclusion in hit list and alignment listing
  par.b=1;                     // min number of alignments
  par.B=1;                     // max number of alignments
  par.z=1;                     // min number of lines in hit list
  par.Z=1;                     // max number of lines in hit list
  par.altali=1;                // find up to FIVE (possibly overlapping) subalignments
  strcpy(par.outfile,"stdout");


  //make command line input globally available
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
  if (!*par.infile || !*par.tfile)
    {help(); cerr<<endl<<"Error in "<<program_name<<": no input alignment files given (-i <file> <file>)\n"; exit(4);}
  if (!*par.dbfiles)
    {help(); cerr<<endl<<"Error in "<<program_name<<": no HMM database file given (-d <file>)\n"; exit(4);}

  // Get rootname (no directory path, no extension) and extensions of infiles
  RemoveExtension(q1.file,par.infile);
  Extension(in1ext,par.infile);
  RemoveExtension(q2.file,par.tfile);
  Extension(in2ext,par.tfile);

  if (par.nseqdis>MAXSEQDIS-3) par.nseqdis=MAXSEQDIS-3; //3 reserved for secondary structure
  if (par.aliwidth<50) par.aliwidth=50;

  // Input parameters
  if (v>=3)
    {
      cerr<<"Input files :  "<<par.infile<<", "<<par.tfile<<"\n";
      cerr<<"Database file: "<<par.dbfiles<<"\n";
      cerr<<"Output file:   "<<par.outfile<<"\n";
    }

  // Set secondary structure substitution matrix
  SetSecStrucSubstitutionMatrix();

  // Set (global variable) substitution matrix with derived matrices and background frequencies
  SetSubstitutionMatrix();

  // Read input file (HMM, HHM, or alignment format), and add pseudocounts etc.
  ReadAndPrepare(par.infile,q1);

  // Read input file (HMM, HHM, or alignment format), and add pseudocounts etc.
  ReadAndPrepare(par.tfile,q2);


  // Allocate memory for dynamic programming matrix
  hit.AllocateBacktraceMatrix(imax(q1.L,q2.L)+2, MAXRES);

  // Open HMM database
  FILE* dbf = fopen(par.dbfiles,"r");
  if (!dbf) OpenFileError("<par.dbfiles");

  ////////////////////////////////////////////////////////////////////////////////
  // Read in one template HMM after the other and compare with both query HMMs
  N_searched=0; readstatus=1;
  int v1=v;
  if (v>=2) v=1;
  while (readstatus) //readstatus 0: End of db file  1: ok, continue  2: skip this HMM
    {
      readstatus = t.Read(dbf);
      if (readstatus!=1) continue; //skip current HMM of reached end of database

      // Add transition pseudocounts to template
      t.AddTransitionPseudocounts();

      // Prepare q ant t and compare

      // Do HMM-HMM comparisons
      hit.irep=1;
      Prepare(q1,t);
      hit.Viterbi(q1,t,0);
      hit.Backtrace(q1,t);
//       printf ("%-12.12s  %-12.12s   irep=%-2i  score=%6.2f\n",hit.name,hit.fam,hit.irep,hit.score);
      hitlist1.Push(hit);

      if (corr_q2_vs_HMMs_like_PSIBLAST)
        {
          // Set all transition probabilities to zero in template
          t.AddTransitionPseudocounts(0.0,0.0,1.0,1.0,1.0,1.0,10000.0);
          Prepare(q2,t,0);
        } else {
          Prepare(q2,t);
        }

      hit.irep=1;
      hit.Viterbi(q2,t,0);
      hit.Backtrace(q2,t);
//       printf ("%-12.12s  %-12.12s   irep=%-2i  score=%6.2f\n",hit.name,hit.fam,hit.irep,hit.score);
      hitlist2.Push(hit);

      N_searched++;
      if (v1>=2 & !(N_searched%20))
        {
          printf(".");
          if (!(N_searched%1000)) printf(" %-4i HMMs searched\n",N_searched);
          cout.flush();
        }
    } // read next HMM
  // End: Read in one template HMM after the other and do p-p comparison
  ////////////////////////////////////////////////////////////////////////////////
  fclose(dbf);
  if (v1>=2) printf("\n");
  v=v1;

  hitlist1.N_searched=N_searched; //hand over number of HMMs scanned to hitlist (for E-value calculation)
  hitlist2.N_searched=N_searched; //hand over number of HMMs scanned to hitlist (for E-value calculation)

  // Fit EVD (with lamda, mu) to score distribution
  if (N_searched>=100)
    {
      hitlist1.MaxLikelihoodEVD(q1,3); // first ML fit: exclude 3 best superfamilies from fit
      hitlist1.MaxLikelihoodEVD(q1,0); // second ML fit: exclude superfamilies with E-value<MINEVALEXCL
      hitlist2.MaxLikelihoodEVD(q2,3); // first ML fit: exclude 3 best superfamilies from fit
      hitlist2.MaxLikelihoodEVD(q2,0); // second ML fit: exclude superfamilies with E-value<MINEVALEXCL
      if (v>=2 && par.loc==0)
        fprintf(stderr,"WARNING: E-values for global alignment option may be unreliable.\n");
    }
  else if (v>=2)
    {
      fprintf(stderr,"Error: no E-values could be determined from only %i searched HMMs.\n",N_searched);
      exit(1);
    }

  // Write HMM to output file without pseudocounts
  if (par.calibrate && !strcmp(in1ext,"hhm")) q1.InsertCalibration(par.infile);
  if (par.calibrate && !strcmp(in2ext,"hhm")) q2.InsertCalibration(par.tfile);

  // Calculate the correlation alpha between q1 and q2:
  // Let S1 and S2 be the scores for q1 and q2 giving a P-value of pval, i.e.pval=P(s1>S1)=P(s2>S2)
  // We define f(pval) = P(s1>S1 & s2>S2)/P(s1>S1)/P(s2>S2) =  P(s1>S1|s2>S2)/P(s1>S1).
  // If s1 and s2 are perfectly uncorrelated, i.e. if P(s1>S1 & s2>S2)=P(s1>S1)*P(s2>S2), we have  f(pval) = 1 = pval^0.
  // If s1 and s2 are perfectly correlated, i.e. if P(s1>S1 & s2>S2)=P(s1>S1)=P(s2>S2), we have f(pval) = pval^(-1)
  // We define a correlation coefficient corr by f(pval) ~ pval^(-corr)

  const int NPVAL=6;
  float pval[]={1.0, 0.4, 0.2, 0.15, 0.1, 0.05};
  float pval1;     // pvalue of template with q1
  float pval2;     // pvalue of template with q2
  float c1[NPVAL], c2[NPVAL], c12[NPVAL]; // #counts with pval1<pval, pval2<pval and (pval1<pval & pval2<pval)
  float f[NPVAL];
  float corr;      // correlation coefficient
  int i;           // index for pval[i]
  const int NPIX=500;
  for (int i=0; i<NPVAL; i++) c1[i]=c2[i]=c12[i]=0;
  hitlist1.Reset();
  hitlist2.Reset();
  while (!hitlist1.End())
    {
      pval1=(hitlist1.ReadNext()).Pval;
      pval2=(hitlist2.ReadNext()).Pval;
//       Hit hit1, hit2;
//       hit1=hitlist1.ReadNext();
//       hit2=hitlist2.ReadNext();
//       pval1=hit1.Pval;
//       pval2=hit2.Pval;
//       printf ("%-12.12s  %-12.12s  %-12.12s  pval1=%6.4f  pval2=%6.4f \n",hit1.name,hit2.name,hit1.fam,pval1,pval2);
//       printf("%6.3f  %6.3f\n",pval1,pval2);
      if (pval1<1E-4 || pval2<1E-4) continue; // reject this point because of probable homology
      for (int i=0; i<NPVAL; i++)
        {
          if (pval1<=pval[i]) c1[i]++;
          if (pval2<=pval[i]) c2[i]++;
          if (pval1<=pval[i] && pval2<=pval[i]) c12[i]++;
        }
    }
  for (int i=0; i<NPVAL; i++) f[i]=c12[i]*c1[0]/c1[i]/c2[i];

  // Write png file?
  if (*pngfile)
    {
      pngwriter png(NPIX,NPIX,0.0,pngfile);
      hitlist1.Reset();
      hitlist2.Reset();
      while (!hitlist1.End())
        {
          pval1=(hitlist1.ReadNext()).Pval;
          pval2=(hitlist2.ReadNext()).Pval;
          png.plot(1+ifloor(pval1*NPIX),1+ifloor(pval2*NPIX),1.,1.,1.);
        }
      png.close();
      if (v>=2) printf("Dot plot written to %s\n",pngfile);
    }

  // Write results into par.outfile?
  if (*par.outfile)
    {
      FILE *outf;
      outf=fopen(par.outfile,"w");
      if (!outf) OpenFileError("par.outfile");
      for (int i=0; i<NPVAL; i++)
        {
          fprintf(outf,"%6.4f %3i %3i %3i  %6.3f  %6.3f  %6.3f  %6.3f\n",pval[i],(int) c12[i],(int)c1[i],(int)c2[i],-log2(f[i]*(1.+sqrt(1./c12[i])))/log2(pval[i]),-log2(pval[i]),log2(f[i]),log2(1+sqrt(1./c12[i])));
        }
      fclose(outf);
      if (v>=2) printf("Statistical results written to %s\n",par.outfile);
    }

  // Calculate corr from point 0 to point 3 (pval=1 and pval=0.1) and add the standard error sqrt(1./c12[3])
  for (i=NPVAL-1; i>=1; i--) if(c12[i]>=15) break;
  corr = -log2(f[i]*(1.+sqrt(1./c12[i])))/log2(pval[i]);
  if (v>=2) printf("Correlation coefficient <~ %4.2f\n",corr);

  // Delete memory for dynamic programming matrix
  hit.DeleteBacktraceMatrix(imax(q1.L,q2.L)+2);

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

//////////////////////////////////////////////////////////////////////////////////////////////////////
// END OF MAIN
//////////////////////////////////////////////////////////////////////////////////////////////////////



