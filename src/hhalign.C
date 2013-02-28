// hhalign.C: 
// Align a multiple alignment to an alignment or HMM 
// Print out aligned input sequences in a3m format
// Error codes: 0: ok  1: file format error  2: file access error  3: memory error  4: internal numeric error  5: command line error

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

////#define WINDOWS
#define MAIN

#include <iostream>   // cin, cout, cerr
#include <fstream>    // ofstream, ifstream
#include <cstdio>     // printf
#include <algorithm>  // min,max
#include <stdlib.h>   // exit
#include <string.h>   // strcmp, strstr
#include <math.h>     // sqrt, pow
#include <limits.h>   // INT_MIN
#include <float.h>    // FLT_MIN
#include <ctype.h>    // islower, isdigit etc
#include <time.h>     // clock_gettime etc. (in realtime library (-lrt compiler option))
#include <errno.h>    // perror()
#include <cassert>
#include <stdexcept>

#include <sys/time.h>
//#include <new>
//#include "efence.h"
//#include "efence.c"

#ifdef HH_SSE41
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
#include "crf_pseudocounts-inl.h"

#include "util.C"        // imax, fmax, iround, iceil, ifloor, strint, strscn, strcut, substr, uprstr, uprchr, Basename etc.
#include "list.C"        // list data structure
#include "hash.C"        // hash data structure
#include "hhdecl.C"      // Constants, global variables, struct Parameters
#include "hhutil.C"      // MatchChr, InsertChr, aa2i, i2aa, log2, fast_log2, ScopID, WriteToScreen,
#include "hhmatrices.C"  // BLOSUM50, GONNET, HSDM

#include "hhhmm.h"       // class HMM
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

#ifdef HH_PNG
#include "pngwriter.h"   //PNGWriter (http://pngwriter.sourceforge.net/)
#include "pngwriter.cc"  //PNGWriter (http://pngwriter.sourceforge.net/)
#endif	    

/////////////////////////////////////////////////////////////////////////////////////
// Global variables 
/////////////////////////////////////////////////////////////////////////////////////
HMM* q = new HMM;            // Create query    HMM with maximum of par.maxres match states
HMM* t = new HMM;            // Create template HMM with maximum of par.maxres match states
Alignment qali;              // (query alignment might be needed outside of hhfunc.C for -a option)
Hit hit;                     // Ceate new hit object pointed at by hit
HitList hitlist;             // list of hits with one Hit object for each pairwise comparison done
char aliindices[256];        // hash containing indices of all alignments which to show in dot plot
char* dmapfile=NULL;         // where to write the coordinates for the HTML map file (to click the alignments)
char* pngfile=NULL;          // pointer to pngfile
char* tcfile=NULL;           // TCoffee output file name
float probmin_tc=0.05;       // 5% minimum posterior probability for printing pairs of residues for TCoffee

int dotW=10;                 // average score of dot plot over window [i-W..i+W]
float dotthr=0.5;            // probability/score threshold for dot plot
int dotscale=600;            // size scale of dotplot
char dotali=0;               // show no alignments in dotplot
float dotsat=0.3;            // saturation of grid and alignments in dot plot
float pself=0.001;           // maximum p-value of 2nd and following self-alignments

/////////////////////////////////////////////////////////////////////////////////////
// Help functions
/////////////////////////////////////////////////////////////////////////////////////
void help()
{
  printf("\n");
  printf("HHalign %s\n",VERSION_AND_DATE);
  printf("Align a query alignment/HMM to a template alignment/HMM by HMM-HMM alignment\n");
  printf("If only one alignment/HMM is given it is compared to itself and the best\n");
  printf("off-diagonal alignment plus all further non-overlapping alignments above \n");
  printf("significance threshold are shown.\n");
  printf("%s",REFERENCE);
  printf("%s",COPYRIGHT);
  printf("\n");
  printf("Usage: %s -i query [-t template] [options]  \n",program_name);
  printf(" -i <file>      input query alignment  (fasta/a2m/a3m) or HMM file (.hhm)\n");
  printf(" -t <file>      input template alignment (fasta/a2m/a3m) or HMM file (.hhm)\n");
#ifdef HH_PNG
  printf(" -png <file>    write dotplot into PNG-file (default=none)           \n");
#endif
  printf("\n");         
  printf("Output options:                                                           \n");
  printf(" -o <file>      write output alignment to file\n"); 
  printf(" -ofas <file>   write alignments in FASTA, A2M (-oa2m) or A3M (-oa3m) format   \n"); 
  printf(" -Oa3m <file>   write query alignment in a3m format to file (default=none)\n");
  printf(" -Aa3m <file>   append query alignment in a3m format to file (default=none)\n");
  printf(" -atab <file>   write alignment as a table (with posteriors) to file (default=none)\n");
  printf(" -index <file>  use given alignment to calculate Viterbi score (default=none)\n");
  printf(" -v <int>       verbose mode: 0:no screen output  1:only warings  2: verbose\n");
  printf(" -seq  [1,inf[  max. number of query/template sequences displayed  (def=%i)  \n",par.nseqdis);
  printf(" -nocons        don't show consensus sequence in alignments (default=show) \n");
  printf(" -nopred        don't show predicted 2ndary structure in alignments (default=show) \n");
  printf(" -nodssp        don't show DSSP 2ndary structure in alignments (default=show) \n");
  printf(" -ssconf        show confidences for predicted 2ndary structure in alignments\n");
  printf(" -aliw int      number of columns per line in alignment list (def=%i)\n",par.aliwidth);
  printf(" -P <float>     for self-comparison: max p-value of alignments (def=%.2g\n",pself);
  printf(" -p <float>     minimum probability in summary and alignment list (def=%G) \n",par.p);
  printf(" -E <float>     maximum E-value in summary and alignment list (def=%G)     \n",par.E);
  printf(" -Z <int>       maximum number of lines in summary hit list (def=%i)       \n",par.Z);
  printf(" -z <int>       minimum number of lines in summary hit list (def=%i)       \n",par.z);
  printf(" -B <int>       maximum number of alignments in alignment list (def=%i)    \n",par.B);
  printf(" -b <int>       minimum number of alignments in alignment list (def=%i)    \n",par.b);
  printf(" -rank int      specify rank of alignment to write with -Oa3m or -Aa3m option (default=1)\n");
  printf("\n");         
#ifdef HH_PNG
  printf("Dotplot options:\n");
  printf(" -dthr <float>  probability/score threshold for dotplot (default=%.2f)        \n",dotthr);
  printf(" -dsca <int>    if value <= 20: size of dot plot unit box in pixels           \n");
  printf("                if value > 20: maximum dot plot size in pixels (default=%i)   \n",dotscale);
  printf(" -dwin <int>    average score over window [i-W..i+W] (for -norealign) (def=%i)\n",dotW);
  printf(" -dali <list>   show alignments with indices in <list> in dot plot            \n");
  printf("                <list> = <index1> ... <indexN>  or  <list> = all              \n");
  printf("\n");         
#endif
  printf("Filter input alignment (options can be combined):                         \n");
  printf(" -id   [0,100]  maximum pairwise sequence identity (%%) (def=%i)   \n",par.max_seqid);
  printf(" -diff [0,inf[  filter most diverse set of sequences, keeping at least this    \n");
  printf("                many sequences in each block of >50 columns (def=%i)\n",par.Ndiff);
  printf(" -cov  [0,100]  minimum coverage with query (%%) (def=%i) \n",par.coverage);
  printf(" -qid  [0,100]  minimum sequence identity with query (%%) (def=%i) \n",par.qid);
  printf(" -qsc  [0,100]  minimum score per column with query  (def=%.1f)\n",par.qsc);
  printf("\n");         
  printf("Input alignment format:                                                     \n");
  printf(" -M a2m         use A2M/A3M (default): upper case = Match; lower case = Insert;\n");         
  printf("                '-' = Delete; '.' = gaps aligned to inserts (may be omitted)   \n");
  printf(" -M first       use FASTA: columns with residue in 1st sequence are match states\n");
  printf(" -M [0,100]     use FASTA: columns with fewer than X%% gaps are match states   \n");
  printf("\n");    
  printf("HMM-HMM alignment options:                                                  \n");
  printf(" -glob/-loc     global or local alignment mode (def=local)         \n");
  printf(" -alt <int>     show up to this number of alternative alignments (def=%i)    \n",par.altali);
  //  printf(" -vit          use Viterbi algorithm for alignment instead of MAC algorithm \n");
  //  printf(" -mac          use Maximum Accuracy (MAC) alignment (default)  \n");
  printf(" -realign       realign displayed hits with max. accuracy (MAC) algorithm \n");
  printf(" -norealign     do NOT realign displayed hits with MAC algorithm (def=realign)\n");
  printf(" -mact [0,1[    posterior prob threshold for MAC realignment controlling greedi- \n");
  printf("                ness at alignment ends: 0:global >0.1:local (default=%.2f)       \n",par.mact);
  printf(" -macins [0,1[  posterior prob threshold for MAC realignment controlling greedi- \n");
  printf("                ness for aligning nonhomologous inserts to each other (def=%.2f)\n",par.macins);
  printf(" -excl <range>  exclude query positions from the alignment, e.g. '1-33,97-168'\n");
  printf(" -shift [-1,1]  score offset (def=%-.3f)                                      \n",par.shift);
  printf(" -corr [0,1]    weight of term for pair correlations (def=%.2f)               \n",par.corr);
  printf(" -ssm  0-4      0:no ss scoring [default=%i]               \n",par.ssm);
  printf("                1:ss scoring after alignment                                  \n");
  printf("                2:ss scoring during alignment                                 \n");
  printf(" -ssw  [0,1]    weight of ss score  (def=%-.2f)                               \n",par.ssw);
  printf("\n");
  printf(" -def           read default options from ./.hhdefaults or <home>/.hhdefault. \n");
  printf("\n");
  printf("Example: %s -i T0187.a3m -t d1hz4a_.hhm -png T0187pdb.png \n",program_name);
  cout<<endl;
//   printf("More help:                                                         \n");
//   printf(" -h out        output options                                      \n");
//   printf(" -h hmm        options for building HMM from multiple alignment    \n");
//   printf(" -h gap        options for setting gap penalties                   \n");
//   printf(" -h ali        options for HMM-HMM alignment                       \n");
//   printf(" -h all        all options \n");
}

void help_out()
{
  printf("\n");
  printf("Output options:                                                           \n");
  printf(" -o <file>      write output alignment to file\n"); 
  printf(" -ofas <file>   write alignments in FASTA, A2M (-oa2m) or A3M (-oa3m) format   \n"); 
  printf(" -Oa3m <file>   write query alignment in a3m format to file (default=none)\n");
  printf(" -Aa3m <file>   append query alignment in a3m format to file (default=none)\n");
  printf(" -atab <file>   write alignment as a table (with posteriors) to file (default=none)\n");
  printf(" -v <int>       verbose mode: 0:no screen output  1:only warings  2: verbose\n");
  printf(" -seq  [1,inf[  max. number of query/template sequences displayed  (def=%i)  \n",par.nseqdis);
  printf(" -nocons        don't show consensus sequence in alignments (default=show) \n");
  printf(" -nopred        don't show predicted 2ndary structure in alignments (default=show) \n");
  printf(" -nodssp        don't show DSSP 2ndary structure in alignments (default=show) \n");
  printf(" -ssconf        show confidences for predicted 2ndary structure in alignments\n");
  printf(" -aliw int      number of columns per line in alignment list (def=%i)\n",par.aliwidth);
  printf(" -P <float>     for self-comparison: max p-value of alignments (def=%.2g\n",pself);
  printf(" -p <float>     minimum probability in summary and alignment list (def=%G) \n",par.p);
  printf(" -E <float>     maximum E-value in summary and alignment list (def=%G)     \n",par.E);
  printf(" -Z <int>       maximum number of lines in summary hit list (def=%i)       \n",par.Z);
  printf(" -z <int>       minimum number of lines in summary hit list (def=%i)       \n",par.z);
  printf(" -B <int>       maximum number of alignments in alignment list (def=%i)    \n",par.B);
  printf(" -b <int>       minimum number of alignments in alignment list (def=%i)    \n",par.b);
  printf(" -rank int      specify rank of alignment to write with -Oa3m or -Aa3m option (default=1)\n");
  printf(" -tc <file>     write a TCoffee library file for the pairwise comparison   \n");         
  printf(" -tct [0,100]   min. probobability of residue pairs for TCoffee (def=%i%%)\n",iround(100*probmin_tc));         
  printf("\n");         
#ifdef HH_PNG
  printf("Dotplot options:\n");
  printf(" -dwin int      average score in dotplot over window [i-W..i+W] (def=%i)   \n",dotW);
  printf(" -dthr float    score threshold for dotplot (default=%.2f)                 \n",dotthr);
  printf(" -dsca int      size of dot plot box in pixels  (default=%i)               \n",dotscale);
  printf(" -dali <list>   show alignments with indices in <list> in dot plot\n");
  printf("                <list> = <index1> ... <indexN>  or  <list> = all              \n");
  printf(" -dmap <file>   print list of coordinates in png plot  \n");
#endif
}

void help_hmm()
{
  printf("\n");
  printf("Options to filter input alignment (options can be combined):              \n");
  printf(" -id   [0,100]  maximum pairwise sequence identity (%%) (def=%i)   \n",par.max_seqid);
  printf(" -diff [0,inf[  filter most diverse set of sequences, keeping at least this    \n");
  printf("                many sequences in each block of >50 columns (def=%i)\n",par.Ndiff);
  printf(" -cov  [0,100]  minimum coverage with query (%%) (def=%i) \n",par.coverage);
  printf(" -qid  [0,100]  minimum sequence identity with query (%%) (def=%i) \n",par.qid);
  printf(" -qsc  [0,100]  minimum score per column with query  (def=%.1f)\n",par.qsc);
  printf("                                                                          \n");
  printf("HMM-building options:                                                     \n");
  printf(" -M a2m         use A2M/A3M (default): upper case = Match; lower case = Insert;\n");         
  printf("                '-' = Delete; '.' = gaps aligned to inserts (may be omitted)   \n");
  printf(" -M first       use FASTA: columns with residue in 1st sequence are match states\n");
  printf(" -M [0,100]     use FASTA: columns with fewer than X%% gaps are match states   \n");
  printf(" -tags          do NOT neutralize His-, C-myc-, FLAG-tags, and \n");
  printf("                trypsin recognition sequence to background distribution    \n");
  printf("                                                                          \n");
  printf("Pseudocount (pc) options:                                                        \n");
  printf(" -pcm {0,..,3}      position dependence of pc admixture 'tau' (pc mode, default=%-i) \n",par.pc.admix);
  printf("                    0: no pseudo counts:    tau = 0                                  \n");
  printf("                    1: constant             tau = a                                  \n");
  printf("                    2: diversity-dependent: tau = a/(1+((Neff[i]-1)/b)^c)            \n");
  printf("                    3: CSBlast admixture:   tau = a(1+b)/(Neff[i]+b)                 \n");
  printf("                    (Neff[i]: number of effective seqs in local MSA around column i) \n");
  printf(" -pca  [0,1]        overall pseudocount admixture (def=%-.1f)                        \n",par.pc.pca);
  printf(" -pcb  [1,inf[      Neff threshold value for -pcm 2 (def=%-.1f)                      \n",par.pc.pcb);
  printf(" -pcc  [0,3]        extinction exponent c for -pcm 2 (def=%-.1f)                     \n",par.pc.pcc);
  // HHsearch option should be the same as HHblits option!!
  printf("Context-specific pseudo-counts:                                                  \n");
  printf(" -nocontxt      use substitution-matrix instead of context-specific pseudocounts \n");
  printf(" -contxt <file> context file for computing context-specific pseudocounts (default=%s)\n",par.clusterfile);
  printf(" -cslib  <file> column state file for fast database prefiltering (default=%s)\n",par.cs_library);
}

void help_gap()
{
  printf("\n");
  printf("Gap cost options:                                                                      \n");
  printf(" -gapb [0,inf[  Transition pseudocount admixture (def=%-.2f)                           \n",par.gapb);
  printf(" -gapd [0,inf[  Transition pseudocount admixture for open gap (default=%-.2f)          \n",par.gapd);
  printf(" -gape [0,1.5]  Transition pseudocount admixture for extend gap (def=%-.2f)            \n",par.gape);
  printf(" -gapf ]0,inf]  factor to increase/reduce the gap open penalty for deletes (def=%-.2f) \n",par.gapf);
  printf(" -gapg ]0,inf]  factor to increase/reduce the gap open penalty for inserts (def=%-.2f) \n",par.gapg);
  printf(" -gaph ]0,inf]  factor to increase/reduce the gap extend penalty for deletes(def=%-.2f)\n",par.gaph);
  printf(" -gapi ]0,inf]  factor to increase/reduce the gap extend penalty for inserts(def=%-.2f)\n",par.gapi);
  printf(" -egq  [0,inf[  penalty (bits) for end gaps aligned to query residues (def=%-.2f)      \n",par.egq);
  printf(" -egt  [0,inf[  penalty (bits) for end gaps aligned to template residues (def=%-.2f)   \n",par.egt);
  printf("\n");
}

void help_ali()
{
  printf("\n");
  printf("Alignment options:  \n");
  printf(" -glob/-loc     global or local alignment mode (def=global)                \n");
  printf(" -mac           use Maximum Accuracy (MAC) alignment instead of Viterbi\n");
  printf(" -mact [0,1[    posterior prob threshold for MAC realignment controlling greedi- \n");
  printf("                ness at alignment ends: 0:global >0.1:local (default=%.2f)       \n",par.mact);
  printf(" -macins [0,1[  posterior prob threshold for MAC realignment controlling greedi- \n");
  printf("                ness for aligning nonhomologous inserts to each other (def=%.2f)\n",par.macins);
  printf(" -sto <int>     use global stochastic sampling algorithm to sample this many alignments\n");
  printf(" -sc   <int>    amino acid score         (tja: template HMM at column j) (def=%i)\n",par.columnscore);
  printf("        0       = log2 Sum(tja*qia/pa)   (pa: aa background frequencies)    \n");
  printf("        1       = log2 Sum(tja*qia/pqa)  (pqa = 1/2*(pa+ta) )               \n");
  printf("        2       = log2 Sum(tja*qia/ta)   (ta: av. aa freqs in template)     \n");
  printf("        3       = log2 Sum(tja*qia/qa)   (qa: av. aa freqs in query)        \n");
  printf(" -corr [0,1]    weight of term for pair correlations (def=%.2f)            \n",par.corr);
  printf(" -shift [-1,1]  score offset (def=%-.3f)                                   \n",par.shift);
  printf(" -r             repeat identification: multiple hits not treated as independent\n");
  printf(" -ssm  0-2      0:no ss scoring [default=%i]               \n",par.ssm);
  printf("                1:ss scoring after alignment                               \n");
  printf("                2:ss scoring during alignment                              \n");
  printf(" -ssw  [0,1]    weight of ss score compared to column score (def=%-.2f)    \n",par.ssw);
  printf(" -ssa  [0,1]    ss confusion matrix = (1-ssa)*I + ssa*psipred-confusion-matrix [def=%-.2f)\n",par.ssa);
  printf(" -maxmem [1,inf[ limit memory for realignment (in GB) (def=%.1f)          \n",par.maxmem);
}

void help_all()
{
  help();
  help_out();
  help_hmm();
  help_gap();
  help_ali();
  printf(" -calm 0-3      empirical score calibration of 0:query 1:template 2:both (def=off)\n");
  printf("\n");
  printf("Default options can be specified in './.hhdefaults' or '~/.hhdefaults'\n");
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
      if (!strcmp(argv[i],"-i"))
	{
	  if (++i>=argc || argv[i][0]=='-') 
	    {help(); cerr<<endl<<"Error in "<<program_name<<": no query file following -i\n"; exit(4);}
	  else strcpy(par.infile,argv[i]);
	}
      else if (!strcmp(argv[i],"-t"))
	{
	  if (++i>=argc || argv[i][0]=='-') 
	    {help(); cerr<<endl<<"Error in "<<program_name<<": no template file following -d\n"; exit(4);}
	  else strcpy(par.tfile,argv[i]);
	}
      else if (!strcmp(argv[i],"-o"))
	{
	  if (++i>=argc) 
	    {help(); cerr<<endl<<"Error in "<<program_name<<": no filename following -o\n"; exit(4);}
	  else strcpy(par.outfile,argv[i]);
	}
      else if (!strcmp(argv[i],"-ofas"))
	{
	  par.outformat=1;
	  if (++i>=argc || argv[i][0]=='-') 
	    {help() ; cerr<<endl<<"Error in "<<program_name<<": no output file following -o\n"; exit(4);}
	  else strcpy(par.pairwisealisfile,argv[i]);
	}
      else if (!strcmp(argv[i],"-oa2m"))
	{
	  par.outformat=2;
	  if (++i>=argc || argv[i][0]=='-') 
	    {help() ; cerr<<endl<<"Error in "<<program_name<<": no output file following -o\n"; exit(4);}
	  else strcpy(par.pairwisealisfile,argv[i]);
	}
      else if (!strcmp(argv[i],"-oa3m"))
	{
	  par.outformat=3;
	  if (++i>=argc || argv[i][0]=='-') 
	    {help() ; cerr<<endl<<"Error in "<<program_name<<": no output file following -o\n"; exit(4);}
	  else strcpy(par.pairwisealisfile,argv[i]);
	}
      else if (!strcmp(argv[i],"-rank") && (i<argc-1)) par.hitrank=atoi(argv[++i]); 
      else if (!strcmp(argv[i],"-Oa3m"))
	{
	  par.append=0;
	  if (++i>=argc || argv[i][0]=='-') 
	    {help() ; cerr<<endl<<"Error in "<<program_name<<": no output file following -Oa3m\n"; exit(4);}
	  else strcpy(par.alnfile,argv[i]);
	}
      else if (!strcmp(argv[i],"-Aa3m"))
	{
	  par.append=1;
	  if (++i>=argc || argv[i][0]=='-') 
	    {help() ; cerr<<endl<<"Error in "<<program_name<<": no output file following -Aa3m\n"; exit(4);}
	  else strcpy(par.alnfile,argv[i]);
	}
      else if (!strcmp(argv[i],"-Opsi"))
	{
	  par.append=0;
	  if (++i>=argc || argv[i][0]=='-') 
	    {help() ; cerr<<endl<<"Error in "<<program_name<<": no output file following -Opsi\n"; exit(4);}
	  else strcpy(par.psifile,argv[i]);
	}
      else if (!strcmp(argv[i],"-Apsi"))
	{
	  par.append=1;
	  if (++i>=argc || argv[i][0]=='-') 
	    {help() ; cerr<<endl<<"Error in "<<program_name<<": no output file following -Apsi\n"; exit(4);}
	  else strcpy(par.psifile,argv[i]);
	}
      else if (!strcmp(argv[i],"-png"))
	{
	  if (++i>=argc) 
	    {help(); cerr<<endl<<"Error in "<<program_name<<": no filename following -png\n"; exit(4);}
	  else 
	    {
	      pngfile = new(char[strlen(argv[i])+1]);
	      strcpy(pngfile,argv[i]);
	    }
	}
      else if (!strcmp(argv[i],"-atab") || !strcmp(argv[i],"-Aliout"))
	{
	  if (++i>=argc || argv[i][0]=='-') 
	    {help(); cerr<<endl<<"Error in "<<program_name<<": no query file following -atab\n"; exit(4);}
	  else strmcpy(par.alitabfile,argv[i],NAMELEN-1);
	}
      else if (!strcmp(argv[i],"-index"))
	{
	  if (++i>=argc || argv[i][0]=='-') 
	    {help(); cerr<<endl<<"Error in "<<program_name<<": no index file following -index\n"; exit(4);}
	  else strcpy(par.indexfile,argv[i]);
	}
      else if (!strcmp(argv[i],"-tc"))
	{
	  if (++i>=argc || argv[i][0]=='-') 
	    {help() ; cerr<<endl<<"Error in "<<program_name<<": no output file following -Opsi\n"; exit(4);}
	  else 
	    {
	      tcfile = new(char[strlen(argv[i])+1]);
	      strcpy(tcfile,argv[i]);
	    }
	}
      else if (!strcmp(argv[i],"-h")|| !strcmp(argv[i],"--help"))
	{
	  if (++i>=argc) {help(); exit(0);} 
	  if (!strcmp(argv[i],"out")) {help_out(); exit(0);} 
	  if (!strcmp(argv[i],"hmm")) {help_hmm(); exit(0);} 
	  if (!strcmp(argv[i],"gap")) {help_gap(); exit(0);} 
	  if (!strcmp(argv[i],"ali")) {help_ali(); exit(0);} 
	  if (!strcmp(argv[i],"all")) {help_all(); exit(0);} 
	  else {help(); exit(0);}
	}
      else if (!strcmp(argv[i],"-excl"))
	{
	  if (++i>=argc) {help(); exit(4);} 
	  par.exclstr = new(char[strlen(argv[i])+1]);
	  strcpy(par.exclstr,argv[i]);
	}
      else if (!strcmp(argv[i],"-v") && (i<argc-1) && argv[i+1][0]!='-' ) v=atoi(argv[++i]);
      else if (!strcmp(argv[i],"-v"))  v=2;
      else if (!strcmp(argv[i],"-v0")) v=0;
      else if (!strcmp(argv[i],"-v1")) v=1;
      else if (!strcmp(argv[i],"-v2")) v=2;
      else if (!strcmp(argv[i],"-v3")) v=3;
      else if (!strcmp(argv[i],"-v4")) v=4;
      else if (!strcmp(argv[i],"-v5")) v=5;
      else if (!strcmp(argv[i],"-P") && (i<argc-1)) pself=atof(argv[++i]);
      else if (!strcmp(argv[i],"-p") && (i<argc-1)) par.p = atof(argv[++i]);
      else if (!strcmp(argv[i],"-e") && (i<argc-1)) par.E = atof(argv[++i]);
      else if (!strcmp(argv[i],"-E") && (i<argc-1)) par.E = atof(argv[++i]);
      else if (!strcmp(argv[i],"-b") && (i<argc-1)) par.b = atoi(argv[++i]);
      else if (!strcmp(argv[i],"-B") && (i<argc-1)) par.B = atoi(argv[++i]);
      else if (!strcmp(argv[i],"-z") && (i<argc-1)) par.z = atoi(argv[++i]);
      else if (!strcmp(argv[i],"-Z") && (i<argc-1)) par.Z = atoi(argv[++i]);
      else if (!strncmp(argv[i],"-nocons",7)) par.showcons=0;
      else if (!strncmp(argv[i],"-nopred",7)) par.showpred=0;
      else if (!strncmp(argv[i],"-nodssp",7)) par.showdssp=0;
      else if (!strncmp(argv[i],"-ssconf",7)) par.showconf=1;
      else if (!strncmp(argv[i],"-mark",7)) par.mark=1;
      else if (!strcmp(argv[i],"-seq") && (i<argc-1))  par.nseqdis=atoi(argv[++i]); 
      else if (!strcmp(argv[i],"-aliw") && (i<argc-1)) par.aliwidth=atoi(argv[++i]); 
      else if (!strcmp(argv[i],"-id") && (i<argc-1))   par.max_seqid=atoi(argv[++i]); 
      else if (!strcmp(argv[i],"-tct") && (i<argc-1))  probmin_tc=atoi(argv[++i]); 
      else if (!strcmp(argv[i],"-dwin") && (i<argc-1)) dotW=atoi(argv[++i]); 
      else if (!strcmp(argv[i],"-dsca") && (i<argc-1)) dotscale=atoi(argv[++i]); 
      else if (!strcmp(argv[i],"-dthr") && (i<argc-1)) dotthr=atof(argv[++i]); 
      else if (!strcmp(argv[i],"-dali") && (i<argc-1))  
	{
	  dotali=1; 
	  for (int index=0; index<256; index++) aliindices[index]=0;
	  while (i+1<argc && argv[i+1][0]!='-') // adds index to hash aliindices
	    {
	      i++;
	      if (strcmp(argv[i],"all")) aliindices[atoi(argv[i])]=1;	      
	      else dotali=2;
	    }
	}
      else if (!strcmp(argv[i],"-dmap"))  
	{
	  if (++i>=argc) 
	    {help(); cerr<<endl<<"Error in "<<program_name<<": no filename following -o\n"; exit(4);}
	  else 
	    {
	      dmapfile = new(char[strlen(argv[i])+1]);
	      strcpy(dmapfile,argv[i]);
	    }
	}
      else if (!strcmp(argv[i],"-dsat") && (i<argc-1)) dotsat=atof(argv[++i]); 
      else if (!strcmp(argv[i],"-qid") && (i<argc-1))  par.qid=atoi(argv[++i]); 
      else if (!strcmp(argv[i],"-qsc") && (i<argc-1))  par.qsc=atof(argv[++i]); 
      else if (!strcmp(argv[i],"-cov") && (i<argc-1))  par.coverage=atoi(argv[++i]); 
      else if (!strcmp(argv[i],"-diff") && (i<argc-1)) par.Ndiff=atoi(argv[++i]); 
      else if (!strcmp(argv[i],"-Gonnet")) par.matrix=0; 
      else if (!strcmp(argv[i],"-HSDM")) par.matrix=1; 
      else if (!strcmp(argv[i],"-BLOSUM50")) par.matrix=2; 
      else if (!strcmp(argv[i],"-Blosum50")) par.matrix=2; 
      else if (!strcmp(argv[i],"-B50")) par.matrix=2; 
      else if (!strcmp(argv[i],"-BLOSUM62")) par.matrix=3; 
      else if (!strcmp(argv[i],"-Blosum62")) par.matrix=3; 
      else if (!strcmp(argv[i],"-B62")) par.matrix=3; 
      else if (!strcmp(argv[i],"-pcm") && (i<argc-1)) par.pc.admix=(Pseudocounts::Admix)atoi(argv[++i]);
      else if (!strcmp(argv[i],"-pca") && (i<argc-1)) par.pc.pca=atof(argv[++i]);
      else if (!strcmp(argv[i],"-pcb") && (i<argc-1)) par.pc.pcb=atof(argv[++i]);
      else if (!strcmp(argv[i],"-pcc") && (i<argc-1)) par.pc.pcc=atof(argv[++i]);
      else if (!strcmp(argv[i],"-gapb") && (i<argc-1)) { par.gapb=atof(argv[++i]); if (par.gapb<=0.01) par.gapb=0.01;} 
      else if (!strcmp(argv[i],"-gapd") && (i<argc-1)) par.gapd=atof(argv[++i]); 
      else if (!strcmp(argv[i],"-gape") && (i<argc-1)) par.gape=atof(argv[++i]); 
      else if (!strcmp(argv[i],"-gapf") && (i<argc-1)) par.gapf=atof(argv[++i]); 
      else if (!strcmp(argv[i],"-gapg") && (i<argc-1)) par.gapg=atof(argv[++i]); 
      else if (!strcmp(argv[i],"-gaph") && (i<argc-1)) par.gaph=atof(argv[++i]); 
      else if (!strcmp(argv[i],"-gapi") && (i<argc-1)) par.gapi=atof(argv[++i]); 
      else if (!strcmp(argv[i],"-egq") && (i<argc-1)) par.egq=atof(argv[++i]); 
      else if (!strcmp(argv[i],"-egt") && (i<argc-1)) par.egt=atof(argv[++i]); 
      else if (!strcmp(argv[i],"-ssgap")) par.ssgap=1;
      else if (!strcmp(argv[i],"-ssgapd") && (i<argc-1)) par.ssgapd=atof(argv[++i]); 
      else if (!strcmp(argv[i],"-ssgape") && (i<argc-1)) par.ssgape=atof(argv[++i]); 
      else if (!strcmp(argv[i],"-ssgapi") && (i<argc-1)) par.ssgapi=atoi(argv[++i]); 
      else if (!strcmp(argv[i],"-ssm") && (i<argc-1)) par.ssm=atoi(argv[++i]); 
      else if (!strcmp(argv[i],"-ssw") && (i<argc-1)) par.ssw=atof(argv[++i]); 
      else if (!strcmp(argv[i],"-ssa") && (i<argc-1)) par.ssa=atof(argv[++i]); 
      else if (!strncmp(argv[i],"-glo",3)) {par.loc=0; if (par.mact>0.35 && par.mact<0.3502) {par.mact=0;} }
      else if (!strncmp(argv[i],"-loc",3)) par.loc=1;
      else if (!strncmp(argv[i],"-alt",4) && (i<argc-1)) par.altali=atoi(argv[++i]); 
      else if (!strcmp(argv[i],"-map") || !strcmp(argv[i],"-MAP") || !strcmp(argv[i],"-mac") || !strcmp(argv[i],"-MAC")) 
	SyntaxError("Please note that this option has been replaced by the '-realign' option."); 
      else if (!strcmp(argv[i],"-vit")) 
	SyntaxError("Please note that this option has been replaced by the '-norealign' option."); 
      else if (!strcmp(argv[i],"-realign")) par.realign=1;
      else if (!strcmp(argv[i],"-norealign")) par.realign=0;
      else if (!strcmp(argv[i],"-M") && (i<argc-1)) 
	if (!strcmp(argv[++i],"a2m") || !strcmp(argv[i],"a3m"))  par.M=1; 
	else if(!strcmp(argv[i],"first"))  par.M=3; 
	else if (argv[i][0]>='0' && argv[i][0]<='9') {par.Mgaps=atoi(argv[i]); par.M=2;}
	else cerr<<endl<<"WARNING: Ignoring unknown argument: -M "<<argv[i]<<"\n";
      else if (!strcmp(argv[i],"-calm") && (i<argc-1)) par.calm=atoi(argv[++i]);
      else if (!strcmp(argv[i],"-shift") && (i<argc-1)) par.shift=atof(argv[++i]); 
      else if (!strcmp(argv[i],"-mact") && (i<argc-1)) {par.mact=atof(argv[++i]);}
      else if (!strcmp(argv[i],"-macins") && (i<argc-1)) par.macins=atof(argv[++i]);
      else if (!strcmp(argv[i],"-opt") && (i<argc-1)) par.opt=atoi(argv[++i]); 
      else if (!strcmp(argv[i],"-scwin") && (i<argc-1)) {par.columnscore=5; par.half_window_size_local_aa_bg_freqs = imax(1,atoi(argv[++i]));}
      else if (!strcmp(argv[i],"-sc") && (i<argc-1)) par.columnscore=atoi(argv[++i]); 
      else if (!strcmp(argv[i],"-def")) par.readdefaultsfile=1; 
      else if (!strcmp(argv[i],"-maxres") && (i<argc-1)) {
	par.maxres=atoi(argv[++i]);
	par.maxcol=2*par.maxres;
      }
      else if (!strcmp(argv[i],"-maxmem") && (i<argc-1)) {par.maxmem=atof(argv[++i]);}
      else if (!strcmp(argv[i],"-corr") && (i<argc-1)) par.corr=atof(argv[++i]); 
      else if (!strcmp(argv[i],"-ovlp") && (i<argc-1)) par.min_overlap=atoi(argv[++i]);
      else if (!strcmp(argv[i],"-tags")) par.notags=0;
      else if (!strcmp(argv[i],"-notags")) par.notags=1;
      else if (!strcmp(argv[i],"-nocontxt")) par.nocontxt=1;
      else if (!strcmp(argv[i],"-csb") && (i<argc-1)) par.csb=atof(argv[++i]);
      else if (!strcmp(argv[i],"-csw") && (i<argc-1)) par.csw=atof(argv[++i]);
      else if (!strcmp(argv[i],"-cs"))
        {
          if (++i>=argc || argv[i][0]=='-')
            {help() ; cerr<<endl<<"Error in "<<program_name<<": no query file following -cs\n"; exit(4);}
          else strcpy(par.clusterfile,argv[i]);
        }
      else cerr<<endl<<"WARNING: Ignoring unknown option "<<argv[i]<<" ...\n";
      if (v>=4) cout<<i<<"  "<<argv[i]<<endl; //PRINT
    } // end of for-loop for command line input
}


///////////////////////////////////////////////////////////////////////////////////////
//// Realign q and *(t[bin]) with the MAC algorithm (after the database search)
//////////////////////////////////////////////////////////////////////////////////////
void RealignByWorker(Hit& hit)
{
  Hit hit_cur;
  int nhits=0;
  hit.irep=1;

  // Prepare MAC comparison(s)
  q->Log2LinTransitionProbs(1.0);
  t->Log2LinTransitionProbs(1.0);

  // Allocate space
  if (par.forward==0)
    hit.AllocateForwardMatrix(q->L+2,t->L+2);
  
  // Search positions in hitlist with correct index of template
  hitlist.Reset();
  while (!hitlist.End())
    {
      hit_cur = hitlist.ReadNext();
      
      // Realign only around previous Viterbi hit
      hit.i1 = hit_cur.i1;
      hit.i2 = hit_cur.i2;
      hit.j1 = hit_cur.j1;
      hit.j2 = hit_cur.j2;
      hit.nsteps = hit_cur.nsteps;
      hit.i = hit_cur.i;
      hit.j = hit_cur.j;
      hit.realign_around_viterbi=true;

      //fprintf(stderr,"  t->name=%s   hit_cur.irep=%i  hit.irep=%i  nhits=%i\n",t->name,hit_cur.irep,hit.irep,nhits);
      // Align q to template in *hit[bin]
      hit.Forward(q,t);
      hit.Backward(q,t);
      hit.MACAlignment(q,t);
      hit.BacktraceMAC(q,t);
      
      // Overwrite *hit[bin] with Viterbi scores, Probabilities etc. of hit_cur
      hit.score      = hit_cur.score;
      hit.score_aass = hit_cur.score_aass;
      hit.score_ss   = hit_cur.score_ss; // comment out?? => Andrea
      hit.Pval       = hit_cur.Pval;
      hit.Pvalt      = hit_cur.Pvalt;
      hit.logPval    = hit_cur.logPval;
      hit.logPvalt   = hit_cur.logPvalt;
      hit.Eval       = hit_cur.Eval;
      hit.logEval    = hit_cur.logEval;
      hit.Probab     = hit_cur.Probab;

      // Replace original hit in hitlist with realigned hit
      //hitlist.ReadCurrent().Delete();
      //hitlist.Delete().Delete();                // delete list record and hit object
      hit_cur = hitlist.Delete();                // delete list record and hit object
      if (hit_cur.irep == 1) {
	delete[] hit_cur.seq;
	delete[] hit_cur.sname;
	hit_cur.seq = NULL;      // Don't delete sname and seq if flat copy from template HMM
	hit_cur.sname = NULL;
      }
      hit_cur.Delete();
      
      hitlist.Insert(hit);
      hit.irep++;
      nhits++;
    }

  // Delete all hitlist entries with too short alignments
  hitlist.Reset();
  while (!hitlist.End())
    {
      hit_cur = hitlist.ReadNext();
      if (hit_cur.matched_cols < MINCOLS_REALIGN && nhits > 1 && nhits > par.hitrank)
	{
	  if (v>=3) printf("Deleting alignment of %s with length %i\n",hit_cur.name,hit_cur.matched_cols);
	  hitlist.Delete().Delete();               // delete the list record and hit object
	  nhits--;
	}
    }

  if (hit.irep==1)
    {
      fprintf(stderr,"*************************************************\n");
      fprintf(stderr,"\nError in %s: could not find template %s in hit list \n\n",par.argv[0],hit.name);
      fprintf(stderr,"*************************************************\n");
    }

  return;
}



/////////////////////////////////////////////////////////////////////////////////////
//// MAIN PROGRAM
/////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{  
  char* argv_conf[MAXOPT];     // Input arguments from .hhdefaults file (first=1: argv_conf[0] is not used)
  int argc_conf;               // Number of arguments in argv_conf 
  char inext[IDLEN]="";        // Extension of query input file (hhm or a3m) 
  char text[IDLEN]="";         // Extension of template input file (hhm or a3m) 
#ifdef HH_PNG
  int** ali=NULL;              // ali[i][j]=1 if (i,j) is part of an alignment
#endif
  int Nali;                    // number of normally backtraced alignments in dot plot

  strcpy(par.tfile,"");
  strcpy(par.alnfile,"");
  par.p=0.0 ;                  // minimum threshold for inclusion in hit list and alignment listing
  par.E=1e6;                   // maximum threshold for inclusion in hit list and alignment listing
  par.b=1;                     // min number of alignments
  par.B=100;                   // max number of alignments
  par.z=1;                     // min number of lines in hit list
  par.Z=100;                   // max number of lines in hit list
  par.append=0;                // append alignment to output file with -a option
  par.altali=1;                // find only ONE (possibly overlapping) subalignment 
  par.hitrank=0;               // rank of hit to be printed as a3m alignment (default=0)
  par.outformat=3;             // default output format for alignment is a3m
  hit.self=0;                  // no self-alignment
  par.forward=0;               // 0: Viterbi algorithm; 1: Viterbi+stochastic sampling; 2:Maximum Accuracy (MAC) algorithm
  par.realign=1;               // default: realign

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

  // Process command line options (they override defaults from .hhdefaults file)
  ProcessArguments(argc,argv);

  // Check command line input and default values
  if (!*par.infile) 
    {help(); cerr<<endl<<"Error in "<<program_name<<": no query alignment file given (-i file)\n"; exit(4);}

  // Get rootname (no directory path, no extension) and extension of infile
  RemoveExtension(q->file,par.infile);
  RemoveExtension(t->file,par.tfile);
  Extension(inext,par.infile); 
  Extension(text,par.tfile); 

  // Check option compatibilities
  if (par.nseqdis>MAXSEQDIS-3-par.showcons) par.nseqdis=MAXSEQDIS-3-par.showcons; //3 reserved for secondary structure
  if (par.aliwidth<20) par.aliwidth=20; 
  if (par.pc.pca<0.001) par.pc.pca=0.001; // to avoid log(0)
  if (par.b>par.B) par.B=par.b;
  if (par.z>par.Z) par.Z=par.z;
  if (par.hitrank>0) par.altali=0;
  if (par.mact>=1.0) par.mact=0.999; else if (par.mact<0) par.mact=0.0;
  if (par.macins>=1.0) par.macins=0.999; else if (par.macins<0) par.macins=0.0;

  // Input parameters
  if (v>=3) 
    {
      cout<<"query file : "<<par.infile<<"\n";
      cout<<"template file: "<<par.tfile<<"\n";
      cout<<"Output file:  "<<par.outfile<<"\n";
      cout<<"Alignment file:  "<<par.alnfile<<"\n";
   }

  // Prepare CS pseudocounts lib
  if (!par.nocontxt && *par.clusterfile) {
    InitializePseudocountsEngine();
  }

  // Set (global variable) substitution matrix and derived matrices
  SetSubstitutionMatrix();

  // Set secondary structure substitution matrix
  SetSecStrucSubstitutionMatrix();

  // Read input file (HMM, HHM, or alignment format), and add pseudocounts etc.
  char input_format=0;
  ReadQueryFile(par.infile,input_format,q,&qali); 
  PrepareQueryHMM(input_format,q);

  // Set query columns in His-tags etc to Null model distribution
  if (par.notags) q->NeutralizeTags();

  // Do self-comparison?
  if (!*par.tfile) 
    {
      if (par.loc==0) {
	cerr<<"WARNING: global alignment not allowed in self-comparison mode. Setting alignment mode to local\n";
	par.loc=1;
      }

      // Deep-copy q into t
      *t = *q;
      
      // Find overlapping alternative alignments
      hit.self=1;

      // Factor Null model into HMM t
      t->IncludeNullModelInHMM(q,t); 
    } 
  // Read template alignment/HMM t and add pseudocounts
  else 
    {
      // Read input file (HMM, HHM, or alignment format), and add pseudocounts etc.
      char input_format=0;
      ReadQueryFile(par.tfile,input_format,t); 
      PrepareTemplateHMM(q,t,input_format);
    }
  
  

  //////////////////////////////////////////////////////////////
  // Calculate Score for given alignment?
  if (*par.indexfile) {

    char line[LINELEN]="";    // input line
    char* ptr;                // pointer for string manipulation
    Hit hit;    
    int step = 0;
    int length = 0;

    // read in indices from indexfile
    FILE* indexf=NULL;
    indexf = fopen(par.indexfile, "r");
    fgetline(line,LINELEN-1,indexf);
    if (!strncmp("#LEN",line,4)) 
      {
	ptr=strscn(line+4);              //advance to first non-white-space character
	length = strint(ptr);
      }
    if (length == 0)
      {
	cerr<<endl<<"Error in "<<program_name<<": first line of index file must contain length of alignment (#LEN ...)\n"; 
	exit(4);
      }

    hit.AllocateIndices(length);

    while (fgetline(line,LINELEN-1,indexf)) 
      {
	if (strscn(line)==NULL) continue;
	if (!strncmp("#QNAME",line,6)) {
	  ptr=strscn(line+6);             // advance to first non-white-space character
	  strmcpy(q->name,ptr,NAMELEN-1);    // copy full name to name
	  strcut(q->name);
	  continue;
	} 
	else if (!strncmp("#TNAME",line,6)) {
	  ptr=strscn(line+6);             // advance to first non-white-space character
	  strmcpy(t->name,ptr,NAMELEN-1);    // copy full name to name
	  strcut(t->name); 
	  continue;
	} 
	else if (line[0] == '#') continue;
	ptr = line;
	hit.i[step] = strint(ptr);
	hit.j[step] = strint(ptr);
	step++;
      }
    
    fclose(indexf);

    // calculate score for each pair of aligned residues
    hit.ScoreAlignment(q,t,step);

    printf("\nAligned %s with %s: Score = %-7.2f \n",q->name,t->name,hit.score);

    if (par.outfile && v>=1) fprintf(stderr,"\nWARNING: no output file is written when -index option is used.\n");
    hit.DeleteIndices();

    exit(0);
  } 
  ////////////////////////////////////////////////////////////////


  // Allocate memory for dynamic programming matrix
  // Longest allowable length of database HMM (backtrace: 5 chars, fwd, bwd: 1 double
  long int Lmaxmem=(par.maxmem*1024*1024*1024)/sizeof(double)/q->L;
  if (par.forward==2 && t->L+2>=Lmaxmem) 
    {
      if (v>=1) {
	cerr<<"WARNING: Not sufficient memory to realign with MAC algorithm. Using Viterbi algorithm."<<endl;
	cerr<<"This is genarally unproboblematic but may lead to slightly sub-optimal alignments."<<endl;
	cerr<<"You can increase available memory for realignment using the -maxmem <GB> option (currently "<<par.maxmem<<" GB)."<<endl; // still to be implemented
	cerr<<"The maximum length realignable is approximately maxmem/query_length/(cpus+1)/8B."<<endl;
      }
      par.forward=0;
    }
  hit.AllocateBacktraceMatrix(q->L+2,t->L+2); // ...with a separate dynamic programming matrix (memory!!)
  if (par.forward>=1) 
    hit.AllocateForwardMatrix(q->L+2,t->L+2);
  
  // Do (self-)comparison, store results if score>SMIN, and try next best alignment
  if (v>=2) 
    {
      if (par.forward==2)       printf("Using maximum accuracy (MAC) alignment algorithm ...\n");
      else if (par.forward==0) printf("Using Viterbi algorithm ...\n");
      else printf("\nWhat alignment algorithm are we using??\n");
    }
  hit.irep=1; 
  while (1)
    {
      if (par.forward==0)        // generate Viterbi alignment
	{
	  hit.Viterbi(q,t);
	  if (hit.irep>1 && hit.irep>par.hitrank && hit.score<=SMIN && !(hit.Pvalt<pself && hit.score>0 )) {
	    hit = hitlist.ReadLast(); // last alignment was not significant => read last (significant) hit from list
	    break;
	  }
	  hit.Backtrace(q,t);
	} 
      else if (par.forward==2)   // generate forward alignment
	{
	  hit.Forward(q,t); 
	  hit.Backward(q,t); 
	  hit.MACAlignment(q,t);
	  if (hit.irep>1 && hit.irep>par.hitrank && hit.score<=SMIN && !(hit.Pvalt<pself && hit.score>0 )) {
	    hit = hitlist.ReadLast(); // last alignment was not significant => read last (significant) hit from list
	    break;
	  }
	  hit.BacktraceMAC(q,t);
	} 
      //fprintf (stderr,"%-12.12s  %-12.12s   irep=%-2i  score=%6.2f hit.Pvalt=%.2g\n",hit.name,hit.fam,hit.irep,hit.score,hit.Pvalt);

      hitlist.Push(hit);      // insert hit at beginning of list (last repeats first!) and do next alignment
      if (hit.irep>=par.hitrank && hit.score<=SMIN && !(hit.Pvalt<pself && hit.score>0 )) break; // last score too bad
      if (hit.irep>=imax(par.hitrank,par.altali)) break; // max number of alignments reached
      hit.irep++;
    } 
  Nali = hit.irep;

  if (par.realign) {
    printf("Realigning using HMM-HMM Maximum Accuracy algorithm\n");
    RealignByWorker(hit);
  }
  
  // Write posterior probability matrix as TCoffee library file
  if (tcfile) 
    {
      if (v>=2) printf("Writing TCoffee library file to %s\n",tcfile);
      int i,j; 
      FILE* tcf=NULL;
      if (strcmp(tcfile,"stdout")) tcf = fopen(tcfile, "w"); else tcf = stdout;
      if (!tcf) OpenFileError(tcfile);
      fprintf(tcf,"! TC_LIB_FORMAT_01\n");
      fprintf(tcf,"%i\n",2); // two sequences in library file
      fprintf(tcf,"%s %i %s\n",q->name,q->L,q->seq[q->nfirst]+1);
      fprintf(tcf,"%s %i %s\n",hit.name,hit.L,hit.seq[hit.nfirst]+1);
      fprintf(tcf,"#1 2\n");
      for (i=1; i<=q->L; i++)  // print all pairs (i,j) with probability above PROBTCMIN
	for (j=1; j<=t->L; j++)
	  if (hit.P_MM[i][j]>probmin_tc) 
	    fprintf(tcf,"%5i %5i %5i\n",i,j,iround(100.0*hit.P_MM[i][j]));
      for (int step=hit.nsteps; step>=1; step--)  // print all pairs on MAC alignment which were not yet printed
	{
	  i=hit.i[step]; j=hit.j[step];
// 	  printf("%5i %5i %5i  %i\n",i,j,iround(100.0*hit.P_MM[i][j]),hit.states[step]);
	  if (hit.states[step]>=MM && hit.P_MM[i][j]<=probmin_tc) 
	    fprintf(tcf,"%5i %5i %5i\n",i,j,iround(100.0*hit.P_MM[i][j]));
	}


      fprintf(tcf,"! SEQ_1_TO_N\n");
      fclose(tcf);
//       for (i=1; i<=q->L; i++)
//        	{
//        	  double sum=0.0;
//        	  for (j=1; j<=t->L; j++) sum+=hit.P_MM[i][j];
// 	  printf("i=%-3i sum=%7.4f\n",i,sum);
//        	}
//        printf("\n");
    }

  // Append last alignment to alitabfile
  if (*par.alitabfile) 
    {
      FILE* alitabf=NULL;
      if (strcmp(par.alitabfile,"stdout")) alitabf = fopen(par.alitabfile, "w"); else alitabf = stdout;
      if (!alitabf) OpenFileError(par.alitabfile);
      WriteToAlifile(alitabf,&hit);
      fclose(alitabf);
    }

  // Fit EVD (with lamda, mu) to score distribution?
  if (par.forward==0)
    {
      if (par.calm==3)  
	hitlist.CalculatePvalues(q);  // Use NN prediction of lamda and mu
      else 
	hitlist.GetPvalsFromCalibration(q);
    }
  else 
    printf("WARNING: E-values and Probabilities can only be calculated when the default Viterbi algorithm is used (with or without -norealign option)\n");

  // Print FASTA or A2M alignments?
  if (*par.pairwisealisfile) {
    if (v>=2) cout<<"Printing alignments in "<<(par.outformat==1? "FASTA" : par.outformat==2?"A2M" :"A3M")<<" format to "<<par.pairwisealisfile<<"\n"; 
    hitlist.PrintAlignments(q,par.pairwisealisfile,par.outformat);
  }

  // Print hit list and alignments
  if (*par.outfile) 
    {
      hitlist.PrintHitList(q,par.outfile);
      hitlist.PrintAlignments(q,par.outfile);
      if (v==2 && strcmp(par.outfile,"stdout")) WriteToScreen(par.outfile,1009); // write only hit list to screen
      // Write whole output file to screen? (max 10000 lines)
      if (v>=3 && strcmp(par.outfile,"stdout")) WriteToScreen(par.outfile,10009);
    }



  //////////////////////////////////////////////////////////////////////////////////////

  // Show results for hit with rank par.hitrank
  if (par.hitrank==0) hit=hitlist.Read(1); else hit=hitlist.Read(par.hitrank);
  
  // Generate output alignment or HMM file?
  if (*par.alnfile || *par.psifile) 
    {
      if (par.append==0) 
	{
	  if (v>=2 && *par.alnfile) printf("Merging template to query alignment and writing resulting alignment in A3M format to %s...\n",par.alnfile);
	  if (v>=2 && *par.psifile) printf("Merging template to query alignment and writing resulting alignment in PSI format to %s...\n",par.psifile);
	}
      else 
	{
	  if (v>=2 && *par.alnfile) printf("Merging template to query alignment and appending template alignment in A3M format to %s...\n",par.alnfile);
	  if (v>=2 && *par.psifile) printf("Merging template to query alignment and appending template alignment in PSI format to %s...\n",par.psifile);
	}

      // Read query alignment into Qali
      Alignment Qali;  // output A3M generated by merging A3M alignments for significant hits to the query alignment
      char qa3mfile[NAMELEN];
      RemoveExtension(qa3mfile,par.infile); // directory??
      strcat(qa3mfile,".a3m");
      FILE* qa3mf=fopen(qa3mfile,"r");
      if (!qa3mf) OpenFileError(qa3mfile);
      Qali.Read(qa3mf,qa3mfile);
      fclose(qa3mf);
      
      // If par.append==1 do not print query alignment
      if (par.append) Qali.MarkSeqsAsNonPrintable();

      // Align query with template in master-slave mode 
      Alignment Tali;
      FILE* ta3mf=fopen(par.tfile,"r");
      if (!ta3mf) OpenFileError(par.tfile);
      Tali.Read(ta3mf,par.tfile); // Read template alignment into Tali
      fclose(ta3mf);
      Tali.Compress(par.tfile); // Filter database alignment
      Qali.MergeMasterSlave(hit,Tali,par.tfile);
      
      // Write output A3M alignment?
      if (*par.alnfile) Qali.WriteToFile(par.alnfile,"a3m");
      
      if (*par.psifile) 
	{
	  // Convert ASCII to int (0-20),throw out all insert states, record their number in I[k][i] 
	  Qali.Compress("merged A3M file");
	  
	  // Write output PSI-BLAST-formatted alignment?
	  Qali.WriteToFile(par.psifile,"psi");
	}
    }

  //////////////////////////////////////////////////////////////////////////////////////


#ifdef HH_PNG
  // Write dot plot into a png file
  if (pngfile)
    {
      // Calculate score[i][j]
      float sum;
      float r,g,b;
      int i,j,l;
      float** s=new(float*[q->L+2]);
      for (i=0; i<q->L+2; i++) 
	if(!(s[i]=new(float[t->L+2]))) MemoryError("image map");
      for(i=1; i<=q->L; i++)
	for (j=1; j<=t->L; j++) // Loop through template positions j
	  {
	    s[i][j]=hit.Score(q->p[i],t->p[j]) + hit.ScoreSS(q,t,i,j) + par.shift;	  }
 	    //printf("%-3i %-3i %7.3f %7.3f\n",i,j,s[i][j],hit.Score(q->p[i],t->p[j]));
      
      // Choose scale automatically
      if (dotscale>20) 
	dotscale=imin(5,imax(1,dotscale/imax(q->L,t->L)));
      
      // Set alignment matrix
      if (dotali) 
	{
	  if (dotali) 
	    {
	      ali = new(int*[q->L+2]);
	      for(i=0; i<q->L+2; i++)
		{
		  ali[i] = new(int[t->L+2]);
		  for(j=1; j<=t->L; j++) ali[i][j]=0;
		}
	    } 
	  int nhits=1;
	  hitlist.Reset();
	  while (!hitlist.End() && nhits<256) 
	    {
	      hit = hitlist.ReadNext();

	      if (nhits>par.z) 
		{
		  if (nhits>=par.Z) continue;       //max number of lines reached?
		  if (hit.Probab < par.p && nhits>=par.z) continue;
		  if (hit.Eval > par.E && nhits>=par.z) continue;
		}
	      if (nhits<Nali && (dotali==2 || aliindices[nhits])) 
		{
		  for (int step=hit.nsteps; step>=1; step--)
		    ali[ hit.i[step] ][ hit.j[step] ]++;
		}
	      nhits++;
	    }
	} // end if (dotali) 
      
      // Write to dmapfile? (will contain regions around alignment traces (clickable in web browser))
      if (dmapfile) 
	{
	  if (v>=2) printf("Printing self-alignment coordinates in png plot to %s\n",dmapfile);
	  const int W=3;
	  int nhits=1;
	  
	  FILE* dmapf = fopen(dmapfile, "w");
	  if (!dmapf) OpenFileError(dmapfile);
	  
	  hitlist.Reset();
	  while (!hitlist.End() && nhits<256) 
	    {
	      hit = hitlist.ReadNext(); // Delete content of hit object
	      int d=dotscale;
	      int i0 = hit.i[hit.nsteps];
	      int j0 = hit.j[hit.nsteps];
	      int i1,j1;
	      if (i0-W<1) {j0+=-i0+W+1; i0+=-i0+W+1; } // avoid underflow 
	      for (int step=hit.nsteps; step>=1; step--)
		{
		  while (step>=1 && hit.states[step]==MM) step--;
		  i1=hit.i[step+1];
		  j1=hit.j[step+1];
		  if (i1+W>q->L) {j1+=-i1+W+q->L; i1+=-i1-W+q->L;} // avoid overflow
		  fprintf(dmapf,"%i COORDS=\"%i,%i,%i,%i,%i,%i,%i,%i\"\n", nhits,
			  d*(j0-1)+1, d*(i0-1-W)+1, d*(j0-1)+1, d*(i0-1+W)+1,
			  d*j1,       d*(i1+W),     d*j1,       d*(i1-W) );
		  while (step>=1 && hit.states[step]>MM) step--;
		  i0=hit.i[step];
		  j0=hit.j[step];
		}
	      nhits++;
	    }
	  fclose(dmapf);
	}

      // Print out dot plot for scores averaged over window of length W
//      printf("x=%i   y=%i,  %s\n",dotscale * t->L,dotscale * q->L,par.pngfile);

      pngwriter png(dotscale * t->L, dotscale * q->L , 1 ,pngfile);  // pngwriter: open png plot
      if (v>=2) cout<<"Writing dot plot to "<<pngfile<<"\n";
      for(i=1; i<=q->L; i++)
	for (j=1; j<=t->L; j++) // Loop through template positions j
	  {
	    float dotval=0.0;
	    sum=0; l=0;
	    if (par.forward<=1 && !par.realign) 
	      {
		for (int w=-dotW; w<=dotW; w++) 
		  if (i+w>=1 && i+w<=q->L && j+w>=1 && j+w<=t->L) 
		    {	
		      sum+=s[i+w][j+w];
		      l++;
		    }		 
		dotval=0.0;
	      } 
	    else 
	      {
		sum = hit.P_MM[i][j];
		dotval = fmax(0.0, 1.0 - 1.0*sum/dotthr); 
		l=1; 

	      }

	    if (i==j && hit.self) {b=r=g=0.0;} 
	    else if ((sum<=0.05 && par.realign) || (sum<=dotthr*l && !par.realign)) 
	      {
	    	if (dotali && ali[i][j]) {r=g=1-dotsat; b=1.0;}
	     	else 
	     	  {
	     	    // Score below threshold
	     	    r=g=b=1.0;
	     	    g -= dotsat/3*(0.7*(!(i%10) || !(j%10)) + (!(i%50) || !(j%50)) + (!(i%100) || !(j%100)));
	     	    b -= dotsat/3*(0.7*(!(i%10) || !(j%10)) + (!(i%50) || !(j%50)) + (!(i%100) || !(j%100)));
	     	  }
	      }	    
	    else
	      {
		// Score above threshold
		if (dotali && ali[i][j]) {r=g=0.0; b=1.0;} 
		else r=g=b=dotval;
	      }

// 	    sum = sum/float(l)*dotthr;
	    for (int ii=dotscale*(q->L-i)+1; ii<=dotscale*(q->L-i+1); ii++)
	      for (int jj=dotscale*(j-1)+1; jj<=dotscale*j; jj++)
		{
		  png.plot(jj,ii,r,g,b); // pngwriter: write to png plot
		}
	  }

      png.close();  // pngwriter: close png plot
      for (i=0; i<q->L+2; i++) delete[] s[i];
      delete[] s;

      // Delete alignment matrix?
      if (dotali) 
	{
	  for(i=0; i<q->L+2; i++) delete[] ali[i];
	  delete[] ali;
	}

    } // if (*par.pngfile)
#endif

//   double log2Pvalue;
//   if (par.ssm && (par.ssm1 || par.ssm2))
//     {
//       log2Pvalue=hit.logPval/0.693147181+0.45*(4.0*par.ssw/0.15-hit.score_ss);
//       if (v>=2) 
// 	printf("Aligned %s with %s:\nApproximate P-value INCLUDING SS SCORE = %7.2g\n",q->name,t->name,pow(2.0,log2Pvalue));
//     } else {
//       if (v>=2) 
// 	printf("Aligned %s with %s:\nApproximate P-value (without SS score) = %7.2g\n",q->name,t->name,hit.Pval);
//    }

  if (v>=2) 
    {
      if (par.hitrank==0) printf("Aligned %s with %s: Score = %-7.2f  P-value = %-7.2g\n",q->name,t->name,hit.score,hit.Pval);
      else printf("Aligned %s with %s (rank %i): Score = %-7.2f  P-value = %-7.2g\n",q->name,t->name,par.hitrank,hit.score,hit.Pval);
    }


  // Delete memory for dynamic programming matrix
  hit.DeleteBacktraceMatrix(q->L+2);
  if (par.forward>=1 || par.realign) 
    hit.DeleteForwardMatrix(q->L+2);
  
  DeletePseudocountsEngine();

  // Delete content of hits in hitlist
  hitlist.Reset();
  while (!hitlist.End()) 
    hitlist.ReadNext().Delete(); // Delete content of hit object

  delete q;
  delete t;


  if (pngfile) delete[] pngfile;
  if (tcfile) delete[] tcfile;
  if (par.exclstr) delete[] par.exclstr;

  // Print 'Done!'
  FILE* outf=NULL;
  if (!strcmp(par.outfile,"stdout")) printf("Done!\n");
  else
    {
      if (*par.outfile)
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



