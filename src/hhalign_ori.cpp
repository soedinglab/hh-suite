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

#include "hhalign_ori.h"

HHalign::HHalign(Parameters& par, std::vector<HHblitsDatabase*>& databases) :
		HHblits(par, databases) {
}

HHalign::~HHalign() {

}

/////////////////////////////////////////////////////////////////////////////////////
// Help functions
/////////////////////////////////////////////////////////////////////////////////////
void HHalign::help(Parameters& par, char all) {
  printf("\n");
  printf("HHalign %s\n", VERSION_AND_DATE);
  printf(
      "Align a query alignment/HMM to a template alignment/HMM by HMM-HMM alignment\n");
  printf(
      "If only one alignment/HMM is given it is compared to itself and the best\n");
  printf(
      "off-diagonal alignment plus all further non-overlapping alignments above \n");
  printf("significance threshold are shown.\n");
  printf("%s", REFERENCE);
  printf("%s", COPYRIGHT);
  printf("\n");
  printf("Usage: hhalign -i query [-t template] [options]  \n");
  printf(
      " -i <file>      input query alignment  (fasta/a2m/a3m) or HMM file (.hhm)\n");
  printf(
      " -t <file>      input template alignment (fasta/a2m/a3m) or HMM file (.hhm)\n");
  printf("\n");
  printf(
      "Output options:                                                           \n");
  printf(" -o <file>      write output alignment to file\n");
  printf(
      " -ofas <file>   write alignments in FASTA, A2M (-oa2m) or A3M (-oa3m) format   \n");
  printf(
      " -Oa3m <file>   write query alignment in a3m format to file (default=none)\n");
  printf(
      " -Aa3m <file>   append query alignment in a3m format to file (default=none)\n");
  printf(
      " -atab <file>   write alignment as a table (with posteriors) to file (default=none)\n");
  printf(
      " -index <file>  use given alignment to calculate Viterbi score (default=none)\n");
  printf(
      " -v <int>       verbose mode: 0:no screen output  1:only warings  2: verbose\n");
  printf(
      " -seq  [1,inf[  max. number of query/template sequences displayed  (def=%i)  \n",
      par.nseqdis);
  printf(
      " -nocons        don't show consensus sequence in alignments (default=show) \n");
  printf(
      " -nopred        don't show predicted 2ndary structure in alignments (default=show) \n");
  printf(
      " -nodssp        don't show DSSP 2ndary structure in alignments (default=show) \n");
  printf(
      " -ssconf        show confidences for predicted 2ndary structure in alignments\n");
  printf(
      " -aliw int      number of columns per line in alignment list (def=%i)\n",
      par.aliwidth);
  printf(
      " -p <float>     minimum probability in summary and alignment list (def=%G) \n",
      par.p);
  printf(
      " -E <float>     maximum E-value in summary and alignment list (def=%G)     \n",
      par.E);
  printf(
      " -Z <int>       maximum number of lines in summary hit list (def=%i)       \n",
      par.Z);
  printf(
      " -z <int>       minimum number of lines in summary hit list (def=%i)       \n",
      par.z);
  printf(
      " -B <int>       maximum number of alignments in alignment list (def=%i)    \n",
      par.B);
  printf(
      " -b <int>       minimum number of alignments in alignment list (def=%i)    \n",
      par.b);
  printf(
      " -rank int      specify rank of alignment to write with -Oa3m or -Aa3m option (default=1)\n");
  printf("\n");
  printf(
      "Filter input alignment (options can be combined):                         \n");
  printf(
      " -id   [0,100]  maximum pairwise sequence identity (%%) (def=%i)   \n",
      par.max_seqid);
  printf(
      " -diff [0,inf[  filter most diverse set of sequences, keeping at least this    \n");
  printf(
      "                many sequences in each block of >50 columns (def=%i)\n",
      par.Ndiff);
  printf(" -cov  [0,100]  minimum coverage with query (%%) (def=%i) \n",
      par.coverage);
  printf(
      " -qid  [0,100]  minimum sequence identity with query (%%) (def=%i) \n",
      par.qid);
  printf(" -qsc  [0,100]  minimum score per column with query  (def=%.1f)\n",
      par.qsc);
  printf("\n");
  printf(
      "Input alignment format:                                                     \n");
  printf(
      " -M a2m         use A2M/A3M (default): upper case = Match; lower case = Insert;\n");
  printf(
      "                '-' = Delete; '.' = gaps aligned to inserts (may be omitted)   \n");
  printf(
      " -M first       use FASTA: columns with residue in 1st sequence are match states\n");
  printf(
      " -M [0,100]     use FASTA: columns with fewer than X%% gaps are match states   \n");
  printf("\n");
  printf(
      "HMM-HMM alignment options:                                                  \n");
  printf(
      " -cs_t <file>  use given column states of the template for scoring   \n");
  printf(
      " -glob/-loc     global or local alignment mode (def=local)         \n");
  printf(
      " -alt <int>     show up to this number of alternative alignments (def=%i)    \n",
      par.altali);
  printf(
      " -realign       realign displayed hits with max. accuracy (MAC) algorithm \n");
  printf(
      " -norealign     do NOT realign displayed hits with MAC algorithm (def=realign)\n");
  printf(
      " -mact [0,1[    posterior prob threshold for MAC realignment controlling greedi- \n");
  printf(
      "                ness at alignment ends: 0:global >0.1:local (default=%.2f)       \n",
      par.mact);
  printf(
      " -macins [0,1[  controls the cost of internal gap positions in the MAC algorithm.\n");
  printf(
      "                0:dense alignments  1:gappy alignments (default=%.2f)\n",
      par.macins);
  printf(
      " -shift [-1,1]  score offset (def=%-.3f)                                      \n",
      par.shift);
  printf(
      " -corr [0,1]    weight of term for pair correlations (def=%.2f)               \n",
      par.corr);
  printf(" -ssm  0-4      0:no ss scoring [default=%i]               \n",
      par.ssm);
  printf(
      "                1:ss scoring after alignment                                  \n");
  printf(
      "                2:ss scoring during alignment                                 \n");
  printf(
      " -ssw  [0,1]    weight of ss score  (def=%-.2f)                               \n",
      par.ssw);
  printf("\n");
  printf(
      " -def           read default options from ./.hhdefaults or <home>/.hhdefault. \n");
  printf("\n");
  printf("Example: hhalign -i T0187.a3m -t d1hz4a_.hhm -png T0187pdb.png \n");

  printf("\n");
  printf(
      "Output options:                                                           \n");
  printf(" -o <file>      write output alignment to file\n");
  printf(
      " -ofas <file>   write alignments in FASTA, A2M (-oa2m) or A3M (-oa3m) format   \n");
  printf(
      " -Oa3m <file>   write query alignment in a3m format to file (default=none)\n");
  printf(
      " -Aa3m <file>   append query alignment in a3m format to file (default=none)\n");
  printf(
      " -atab <file>   write alignment as a table (with posteriors) to file (default=none)\n");
  printf(
      " -v <int>       verbose mode: 0:no screen output  1:only warings  2: verbose\n");
  printf(
      " -seq  [1,inf[  max. number of query/template sequences displayed  (def=%i)  \n",
      par.nseqdis);
  printf(
      " -nocons        don't show consensus sequence in alignments (default=show) \n");
  printf(
      " -nopred        don't show predicted 2ndary structure in alignments (default=show) \n");
  printf(
      " -nodssp        don't show DSSP 2ndary structure in alignments (default=show) \n");
  printf(
      " -ssconf        show confidences for predicted 2ndary structure in alignments\n");
  printf(
      " -aliw int      number of columns per line in alignment list (def=%i)\n",
      par.aliwidth);
  printf(
      " -p <float>     minimum probability in summary and alignment list (def=%G) \n",
      par.p);
  printf(
      " -E <float>     maximum E-value in summary and alignment list (def=%G)     \n",
      par.E);
  printf(
      " -Z <int>       maximum number of lines in summary hit list (def=%i)       \n",
      par.Z);
  printf(
      " -z <int>       minimum number of lines in summary hit list (def=%i)       \n",
      par.z);
  printf(
      " -B <int>       maximum number of alignments in alignment list (def=%i)    \n",
      par.B);
  printf(
      " -b <int>       minimum number of alignments in alignment list (def=%i)    \n",
      par.b);
  printf(
      " -rank int      specify rank of alignment to write with -Oa3m or -Aa3m option (default=1)\n");
  printf(
      " -tc <file>     write a TCoffee library file for the pairwise comparison   \n");
  printf("\n");

  printf("\n");
  printf(
      "Options to filter input alignment (options can be combined):              \n");
  printf(
      " -id   [0,100]  maximum pairwise sequence identity (%%) (def=%i)   \n",
      par.max_seqid);
  printf(
      " -diff [0,inf[  filter most diverse set of sequences, keeping at least this    \n");
  printf(
      "                many sequences in each block of >50 columns (def=%i)\n",
      par.Ndiff);
  printf(" -cov  [0,100]  minimum coverage with query (%%) (def=%i) \n",
      par.coverage);
  printf(
      " -qid  [0,100]  minimum sequence identity with query (%%) (def=%i) \n",
      par.qid);
  printf(" -qsc  [0,100]  minimum score per column with query  (def=%.1f)\n",
      par.qsc);
  printf(
      "                                                                          \n");
  printf(
      "HMM-building options:                                                     \n");
  printf(
      " -M a2m         use A2M/A3M (default): upper case = Match; lower case = Insert;\n");
  printf(
      "                '-' = Delete; '.' = gaps aligned to inserts (may be omitted)   \n");
  printf(
      " -M first       use FASTA: columns with residue in 1st sequence are match states\n");
  printf(
      " -M [0,100]     use FASTA: columns with fewer than X%% gaps are match states   \n");
  printf(" -tags          do NOT neutralize His-, C-myc-, FLAG-tags, and \n");
  printf(
      "                trypsin recognition sequence to background distribution    \n");
  printf(
      "                                                                          \n");

  printf(
      "Pseudocount (pc) options:                                                        \n");
  printf(" Context specific hhm pseudocounts:\n");
  printf(
      "  -pc_hhm_contxt_mode {0,..,3}      position dependence of pc admixture 'tau' (pc mode, default=%-i) \n",
      par.pc_hhm_context_engine.admix);
  printf(
      "               0: no pseudo counts:    tau = 0                                  \n");
  printf(
      "               1: constant             tau = a                                  \n");
  printf(
      "               2: diversity-dependent: tau = a/(1+((Neff[i]-1)/b)^c)            \n");
  printf(
      "               3: CSBlast admixture:   tau = a(1+b)/(Neff[i]+b)                 \n");
  printf(
      "               (Neff[i]: number of effective seqs in local MSA around column i) \n");
  printf(
      "  -pc_hhm_contxt_a  [0,1]        overall pseudocount admixture (def=%-.1f)                        \n",
      par.pc_hhm_context_engine.pca);
  printf(
      "  -pc_hhm_contxt_b  [1,inf[      Neff threshold value for mode 2 (def=%-.1f)                      \n",
      par.pc_hhm_context_engine.pcb);
  printf(
      "  -pc_hhm_contxt_c  [0,3]        extinction exponent c for mode 2 (def=%-.1f)                     \n\n",
      par.pc_hhm_context_engine.pcc);

  printf(
      " Context independent hhm pseudocounts (used for templates; used for query if contxt file is not available):\n");
  printf(
      "  -pc_hhm_nocontxt_mode {0,..,3}      position dependence of pc admixture 'tau' (pc mode, default=%-i) \n",
      par.pc_hhm_nocontext_mode);
  printf(
      "               0: no pseudo counts:    tau = 0                                  \n");
  printf(
      "               1: constant             tau = a                                  \n");
  printf(
      "               2: diversity-dependent: tau = a/(1+((Neff[i]-1)/b)^c)            \n");
  printf(
      "               (Neff[i]: number of effective seqs in local MSA around column i) \n");
  printf(
      "  -pc_hhm_nocontxt_a  [0,1]        overall pseudocount admixture (def=%-.1f)                        \n",
      par.pc_hhm_nocontext_a);
  printf(
      "  -pc_hhm_nocontxt_b  [1,inf[      Neff threshold value for mode 2 (def=%-.1f)                      \n",
      par.pc_hhm_nocontext_b);
  printf(
      "  -pc_hhm_nocontxt_c  [0,3]        extinction exponent c for mode 2 (def=%-.1f)                     \n\n",
      par.pc_hhm_nocontext_c);

  printf(
      " Context-specific pseudo-counts:                                                  \n");
  printf(
      "  -nocontxt      use substitution-matrix instead of context-specific pseudocounts \n");
  printf(
      "  -contxt <file> context file for computing context-specific pseudocounts (default=%s)\n",
      par.clusterfile);

  printf("\n");
  printf(
      "Gap cost options:                                                                      \n");
  printf(
      " -gapb [0,inf[  Transition pseudocount admixture (def=%-.2f)                           \n",
      par.gapb);
  printf(
      " -gapd [0,inf[  Transition pseudocount admixture for open gap (default=%-.2f)          \n",
      par.gapd);
  printf(
      " -gape [0,1.5]  Transition pseudocount admixture for extend gap (def=%-.2f)            \n",
      par.gape);
  printf(
      " -gapf ]0,inf]  factor to increase/reduce the gap open penalty for deletes (def=%-.2f) \n",
      par.gapf);
  printf(
      " -gapg ]0,inf]  factor to increase/reduce the gap open penalty for inserts (def=%-.2f) \n",
      par.gapg);
  printf(
      " -gaph ]0,inf]  factor to increase/reduce the gap extend penalty for deletes(def=%-.2f)\n",
      par.gaph);
  printf(
      " -gapi ]0,inf]  factor to increase/reduce the gap extend penalty for inserts(def=%-.2f)\n",
      par.gapi);
  printf(
      " -egq  [0,inf[  penalty (bits) for end gaps aligned to query residues (def=%-.2f)      \n",
      par.egq);
  printf(
      " -egt  [0,inf[  penalty (bits) for end gaps aligned to template residues (def=%-.2f)   \n",
      par.egt);
  printf("\n");

  printf("\n");
  printf("Alignment options:  \n");
  printf(
      " -glob/-loc     global or local alignment mode (def=global)                \n");
  printf(
      " -mac           use Maximum Accuracy (MAC) alignment instead of Viterbi\n");
  printf(
      " -mact [0,1[    posterior prob threshold for MAC realignment controlling greedi- \n");
  printf(
      "                ness at alignment ends: 0:global >0.1:local (default=%.2f)       \n",
      par.mact);
  printf(
      " -macins [0,1[  posterior prob threshold for MAC realignment controlling greedi- \n");
  printf(
      "                ness for aligning nonhomologous inserts to each other (def=%.2f)\n",
      par.macins);
  printf(
      " -sto <int>     use global stochastic sampling algorithm to sample this many alignments\n");
  printf(
      " -sc   <int>    amino acid score         (tja: template HMM at column j) (def=%i)\n",
      par.columnscore);
  printf(
      "        0       = log2 Sum(tja*qia/pa)   (pa: aa background frequencies)    \n");
  printf(
      "        1       = log2 Sum(tja*qia/pqa)  (pqa = 1/2*(pa+ta) )               \n");
  printf(
      "        2       = log2 Sum(tja*qia/ta)   (ta: av. aa freqs in template)     \n");
  printf(
      "        3       = log2 Sum(tja*qia/qa)   (qa: av. aa freqs in query)        \n");
  printf(
      " -corr [0,1]    weight of term for pair correlations (def=%.2f)            \n",
      par.corr);
  printf(
      " -shift [-1,1]  score offset (def=%-.3f)                                   \n",
      par.shift);
  printf(
      " -r             repeat identification: multiple hits not treated as independent\n");
  printf(" -ssm  0-2      0:no ss scoring [default=%i]               \n",
      par.ssm);
  printf(
      "                1:ss scoring after alignment                               \n");
  printf(
      "                2:ss scoring during alignment                              \n");
  printf(
      " -ssw  [0,1]    weight of ss score compared to column score (def=%-.2f)    \n",
      par.ssw);
  printf(
      " -ssa  [0,1]    ss confusion matrix = (1-ssa)*I + ssa*psipred-confusion-matrix [def=%-.2f)\n",
      par.ssa);
  printf(
      " -maxmem [1,inf[ limit memory for realignment (in GB) (def=%.1f)          \n",
      par.maxmem);

  printf(
      " -calm 0-3      empirical score calibration of 0:query 1:template 2:both (def=off)\n");
  printf("\n");
  printf(
      "Default options can be specified in './.hhdefaults' or '~/.hhdefaults'\n");
}

/////////////////////////////////////////////////////////////////////////////////////
//// Processing input options from command line and .hhdefaults file
/////////////////////////////////////////////////////////////////////////////////////
void HHalign::ProcessAllArguments(int argc, char** argv, Parameters& par) {
  par.argv = argv;
  par.argc = argc;

  strcpy(par.tfile, "");
  strcpy(par.alnfile, "");
  par.p = 0.0; // minimum threshold for inclusion in hit list and alignment listing
  par.E = 1e6; // maximum threshold for inclusion in hit list and alignment listing
  par.b = 1;                     // min number of alignments
  par.B = 100;                   // max number of alignments
  par.z = 1;                     // min number of lines in hit list
  par.Z = 100;                   // max number of lines in hit list
  par.append = 0;              // append alignment to output file with -a option
  par.altali = 1;           // find only ONE (possibly overlapping) subalignment
  par.hitrank = 0;     // rank of hit to be printed as a3m alignment (default=0)
  par.outformat = 3;             // default output format for alignment is a3m
  par.forward = 0; // 0: Viterbi algorithm; 1: Viterbi+stochastic sampling; 2:Maximum Accuracy (MAC) algorithm
  par.realign = 1;               // default: realign

  par.num_rounds = 1;

  // Enable changing verbose mode before command line are processed
  int v = 2;
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-v") == 0) {
      v = atoi(argv[i + 1]);
      break;
    }
  }
  par.v = Log::from_int(v);
  Log::reporting_level() = par.v;

  par.SetDefaultPaths();

  // Process default otpions from .hhdefaults file
  char* argv_conf[MAXOPT];
  int argc_conf = 0;

  ReadDefaultsFile(argc_conf, argv_conf, argv[0]);
  ProcessArguments(argc_conf, argv_conf, par);

  for (int n = 1; n < argc_conf; n++)
    delete[] argv_conf[n];

  // Process command line options (they override defaults from .hhdefaults file)
  ProcessArguments(argc, argv, par);

  // Check needed files
  // Check command line input and default values
  if (!*par.infile) {
    help(par);
    std::cerr << "Error in " << __FILE__ << ":" << __LINE__ << ": " << __func__
        << ":" << std::endl;
    std::cerr << "\tno query alignment file given (-i file)\n";
    exit(4);
  }

  // Check option compatibilities
  if (par.nseqdis > MAXSEQDIS - 3 - par.showcons)
    par.nseqdis = MAXSEQDIS - 3 - par.showcons; //3 reserved for secondary structure
  if (par.aliwidth < 20)
    par.aliwidth = 20;
  if (par.pc_hhm_context_engine.pca < 0.001)
    par.pc_hhm_context_engine.pca = 0.001; // to avoid log(0)
  if (par.b > par.B)
    par.B = par.b;
  if (par.z > par.Z)
    par.Z = par.z;
  if (par.hitrank > 0)
    par.altali = 0;
  if (par.mact >= 1.0)
    par.mact = 0.999;
  else if (par.mact < 0)
    par.mact = 0.0;
  if (par.macins >= 1.0)
    par.macins = 0.999;
  else if (par.macins < 0)
    par.macins = 0.0;
}

void HHalign::ProcessArguments(int argc, char** argv, Parameters& par) {
  for (int i = 1; i < argc; i++) {
    HH_LOG(LogLevel::DEBUG1) << i << "  " << argv[i] << endl;
    if (!strcmp(argv[i], "-i")) {
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        std::cerr << "Error in " << __FILE__ << ":" << __LINE__ << ": "
            << __func__ << ":" << std::endl;
        std::cerr << "\tno query file following -i\n";
        exit(4);
      }
      else
        strcpy(par.infile, argv[i]);
    }
    else if (!strcmp(argv[i], "-t")) {
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        std::cerr << "Error in " << __FILE__ << ":" << __LINE__ << ": "
            << __func__ << ":" << std::endl;
        std::cerr << "\tno template file following -d\n";
        exit(4);
      }
      else
        strcpy(par.tfile, argv[i]);
    }
    else if (!strcmp(argv[i], "-q2t")) {
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        std::cerr << "Error in " << __FILE__ << ":" << __LINE__ << ": "
            << __func__ << ":" << std::endl;
        std::cerr << "\tno query2template file following -d\n";
        exit(4);
      }
      else
        strcpy(par.queries_to_template_file, argv[i]);
    }
    else if (!strcmp(argv[i], "-d")) {
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        std::cerr << "Error in " << __FILE__ << ":" << __LINE__ << ": "
            << __func__ << ":" << std::endl;
        std::cerr << "\tno database file following -d\n";
        exit(4);
      }
      else
    	par.db_bases.push_back(std::string(argv[i]));
    }
    else if (!strcmp(argv[i], "-o")) {
      if (++i >= argc) {
        help(par);
        std::cerr << "Error in " << __FILE__ << ":" << __LINE__ << ": "
            << __func__ << ":" << std::endl;
        std::cerr << "\tno filename following -o\n";
        exit(4);
      }
      else
        strcpy(par.outfile, argv[i]);
    }
    else if (!strcmp(argv[i], "-ored")) {
      if (++i >= argc) {
        help(par);
        std::cerr << "Error in " << __FILE__ << ":" << __LINE__ << ": "
            << __func__ << ":" << std::endl;
        std::cerr << "\tno filename following -o\n";
        exit(4);
      }
      else
        strcpy(par.reduced_outfile, argv[i]);
    }
    else if (!strcmp(argv[i], "-ofas")) {
      par.outformat = 1;
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        std::cerr << "Error in " << __FILE__ << ":" << __LINE__ << ": "
            << __func__ << ":" << std::endl;
        std::cerr << "\tno output file following -o\n";
        exit(4);
      }
      else
        strcpy(par.pairwisealisfile, argv[i]);
    }
    else if (!strcmp(argv[i], "-oa2m")) {
      par.outformat = 2;
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        std::cerr << "Error in " << __FILE__ << ":" << __LINE__ << ": "
            << __func__ << ":" << std::endl;
        std::cerr << "\tno output file following -o\n";
        exit(4);
      }
      else
        strcpy(par.pairwisealisfile, argv[i]);
    }
    else if (!strcmp(argv[i], "-oa3m")) {
      par.outformat = 3;
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        std::cerr << "Error in " << __FILE__ << ":" << __LINE__ << ": "
            << __func__ << ":" << std::endl;
        std::cerr << "\tno output file following -o\n";
        exit(4);
      }
      else
        strcpy(par.pairwisealisfile, argv[i]);
    }
    else if (!strcmp(argv[i], "-rank") && (i < argc - 1))
      par.hitrank = atoi(argv[++i]);
    else if (!strcmp(argv[i], "-Oa3m")) {
      par.append = 0;
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        std::cerr << "Error in " << __FILE__ << ":" << __LINE__ << ": "
            << __func__ << ":" << std::endl;
        std::cerr << "\tno output file following -Oa3m\n";
        exit(4);
      }
      else
        strcpy(par.alnfile, argv[i]);
    }
    else if (!strcmp(argv[i], "-Aa3m")) {
      par.append = 1;
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        std::cerr << "Error in " << __FILE__ << ":" << __LINE__ << ": "
            << __func__ << ":" << std::endl;
        std::cerr << "\tno output file following -Aa3m\n";
        exit(4);
      }
      else
        strcpy(par.alnfile, argv[i]);
    }
    else if (!strcmp(argv[i], "-Opsi")) {
      par.append = 0;
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        std::cerr << "Error in " << __FILE__ << ":" << __LINE__ << ": "
            << __func__ << ":" << std::endl;
        std::cerr << "\tno output file following -Opsi\n";
        exit(4);
      }
      else
        strcpy(par.psifile, argv[i]);
    }
    else if (!strcmp(argv[i], "-Apsi")) {
      par.append = 1;
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        std::cerr << "Error in " << __FILE__ << ":" << __LINE__ << ": "
            << __func__ << ":" << std::endl;
        std::cerr << "\tno output file following -Apsi\n";
        exit(4);
      }
      else
        strcpy(par.psifile, argv[i]);
    }
    else if (!strcmp(argv[i], "-atab") || !strcmp(argv[i], "-Aliout")) {
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        std::cerr << "Error in " << __FILE__ << ":" << __LINE__ << ": "
            << __func__ << ":" << std::endl;
        std::cerr << "\tno query file following -atab\n";
        exit(4);
      }
      else
        strmcpy(par.alitabfile, argv[i], NAMELEN - 1);
    }
    else if (!strcmp(argv[i], "-index")) {
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        std::cerr << "Error in " << __FILE__ << ":" << __LINE__ << ": "
            << __func__ << ":" << std::endl;
        std::cerr << "\tno index file following -index\n";
        exit(4);
      }
      else
        strcpy(par.indexfile, argv[i]);
    }
    else if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) {
      if (++i >= argc) {
        help(par);
        exit(0);
      }
      if (!strcmp(argv[i], "all")) {
        help(par, 1);
        exit(0);
      }
      else {
        help(par);
        exit(0);
      }
    }
    else if (!strcmp(argv[i], "-v") && (i < argc - 1)
        && argv[i + 1][0] != '-') {
      int v = atoi(argv[++i]);
      par.v = Log::from_int(v);
      Log::reporting_level() = par.v;
    }
    else if (!strcmp(argv[i], "-p") && (i < argc - 1))
      par.p = atof(argv[++i]);
    else if (!strcmp(argv[i], "-e") && (i < argc - 1))
      par.E = atof(argv[++i]);
    else if (!strcmp(argv[i], "-E") && (i < argc - 1))
      par.E = atof(argv[++i]);
    else if (!strcmp(argv[i], "-b") && (i < argc - 1))
      par.b = atoi(argv[++i]);
    else if (!strcmp(argv[i], "-B") && (i < argc - 1))
      par.B = atoi(argv[++i]);
    else if (!strcmp(argv[i], "-z") && (i < argc - 1))
      par.z = atoi(argv[++i]);
    else if (!strcmp(argv[i], "-Z") && (i < argc - 1))
      par.Z = atoi(argv[++i]);
    else if (!strncmp(argv[i], "-nocons", 7))
      par.showcons = 0;
    else if (!strncmp(argv[i], "-nopred", 7))
      par.showpred = 0;
    else if (!strncmp(argv[i], "-nodssp", 7))
      par.showdssp = 0;
    else if (!strncmp(argv[i], "-ssconf", 7))
      par.showconf = 1;
    else if (!strncmp(argv[i], "-mark", 7))
      par.mark = 1;
    else if (!strcmp(argv[i], "-seq") && (i < argc - 1))
      par.nseqdis = atoi(argv[++i]);
    else if (!strcmp(argv[i], "-aliw") && (i < argc - 1))
      par.aliwidth = atoi(argv[++i]);
    else if (!strcmp(argv[i], "-id") && (i < argc - 1))
      par.max_seqid = atoi(argv[++i]);
    else if (!strcmp(argv[i], "-qid") && (i < argc - 1))
      par.qid = atoi(argv[++i]);
    else if (!strcmp(argv[i], "-qsc") && (i < argc - 1))
      par.qsc = atof(argv[++i]);
    else if (!strcmp(argv[i], "-cov") && (i < argc - 1))
      par.coverage = atoi(argv[++i]);
    else if (!strcmp(argv[i], "-diff") && (i < argc - 1))
      par.Ndiff = atoi(argv[++i]);
    else if (!strcmp(argv[i], "-Gonnet"))
      par.matrix = 0;
    else if (!strcmp(argv[i], "-HSDM"))
      par.matrix = 1;
    else if (!strcmp(argv[i], "-BLOSUM50"))
      par.matrix = 2;
    else if (!strcmp(argv[i], "-Blosum50"))
      par.matrix = 2;
    else if (!strcmp(argv[i], "-B50"))
      par.matrix = 2;
    else if (!strcmp(argv[i], "-BLOSUM62"))
      par.matrix = 3;
    else if (!strcmp(argv[i], "-Blosum62"))
      par.matrix = 3;
    else if (!strcmp(argv[i], "-B62"))
      par.matrix = 3;
    else if (!strcmp(argv[i], "-pcm") && (i < argc - 1))
      par.pc_hhm_context_engine.admix = (Pseudocounts::Admix) atoi(argv[++i]);
    else if (!strcmp(argv[i], "-pca") && (i < argc - 1))
      par.pc_hhm_context_engine.pca = atof(argv[++i]);
    else if (!strcmp(argv[i], "-pcb") && (i < argc - 1))
      par.pc_hhm_context_engine.pcb = atof(argv[++i]);
    else if (!strcmp(argv[i], "-pcc") && (i < argc - 1))
      par.pc_hhm_context_engine.pcc = atof(argv[++i]);
    else if (!strcmp(argv[i], "-gapb") && (i < argc - 1)) {
      par.gapb = atof(argv[++i]);
      if (par.gapb <= 0.01)
        par.gapb = 0.01;
    }
    else if (!strcmp(argv[i], "-gapd") && (i < argc - 1))
      par.gapd = atof(argv[++i]);
    else if (!strcmp(argv[i], "-gape") && (i < argc - 1))
      par.gape = atof(argv[++i]);
    else if (!strcmp(argv[i], "-gapf") && (i < argc - 1))
      par.gapf = atof(argv[++i]);
    else if (!strcmp(argv[i], "-gapg") && (i < argc - 1))
      par.gapg = atof(argv[++i]);
    else if (!strcmp(argv[i], "-gaph") && (i < argc - 1))
      par.gaph = atof(argv[++i]);
    else if (!strcmp(argv[i], "-gapi") && (i < argc - 1))
      par.gapi = atof(argv[++i]);
    else if (!strcmp(argv[i], "-egq") && (i < argc - 1))
      par.egq = atof(argv[++i]);
    else if (!strcmp(argv[i], "-egt") && (i < argc - 1))
      par.egt = atof(argv[++i]);
    else if (!strcmp(argv[i], "-ssgap"))
      par.ssgap = 1;
    else if (!strcmp(argv[i], "-ssgapd") && (i < argc - 1))
      par.ssgapd = atof(argv[++i]);
    else if (!strcmp(argv[i], "-ssgape") && (i < argc - 1))
      par.ssgape = atof(argv[++i]);
    else if (!strcmp(argv[i], "-ssgapi") && (i < argc - 1))
      par.ssgapi = atoi(argv[++i]);
    else if (!strcmp(argv[i], "-ssm") && (i < argc - 1))
      par.ssm = atoi(argv[++i]);
    else if (!strcmp(argv[i], "-ssw") && (i < argc - 1))
      par.ssw = atof(argv[++i]);
    else if (!strcmp(argv[i], "-ssa") && (i < argc - 1))
      par.ssa = atof(argv[++i]);
    else if (!strncmp(argv[i], "-glo", 3)) {
      par.loc = 0;
      if (par.mact > 0.35 && par.mact < 0.3502) {
        par.mact = 0;
      }
    }
    else if (!strncmp(argv[i], "-loc", 3))
      par.loc = 1;
    else if (!strncmp(argv[i], "-alt", 4) && (i < argc - 1))
      par.altali = atoi(argv[++i]);
    else if (!strcmp(argv[i], "-map") || !strcmp(argv[i], "-MAP")
        || !strcmp(argv[i], "-mac") || !strcmp(argv[i], "-MAC"))
      SyntaxError(__FILE__, __LINE__, __func__,
          "Please note that this option has been replaced by the '-realign' option.");
    else if (!strcmp(argv[i], "-vit"))
      SyntaxError(__FILE__, __LINE__, __func__,
          "Please note that this option has been replaced by the '-norealign' option.");
    else if (!strcmp(argv[i], "-realign"))
      par.realign = 1;
    else if (!strcmp(argv[i], "-norealign"))
      par.realign = 0;
    else if (!strcmp(argv[i], "-M") && (i < argc - 1))
      if (!strcmp(argv[++i], "a2m") || !strcmp(argv[i], "a3m"))
        par.M = 1;
      else if (!strcmp(argv[i], "first"))
        par.M = 3;
      else if (argv[i][0] >= '0' && argv[i][0] <= '9') {
        par.Mgaps = atoi(argv[i]);
        par.M = 2;
      }
      else
        cerr << endl << "WARNING: Ignoring unknown argument: -M " << argv[i]
            << "\n";
    else if (!strcmp(argv[i], "-shift") && (i < argc - 1))
      par.shift = atof(argv[++i]);
    else if (!strcmp(argv[i], "-mact") && (i < argc - 1)) {
      par.mact = atof(argv[++i]);
    }
    else if (!strcmp(argv[i], "-macins") && (i < argc - 1))
      par.macins = atof(argv[++i]);
    else if (!strcmp(argv[i], "-opt") && (i < argc - 1))
      par.opt = atoi(argv[++i]);
    else if (!strcmp(argv[i], "-scwin") && (i < argc - 1)) {
      par.columnscore = 5;
      par.half_window_size_local_aa_bg_freqs = imax(1, atoi(argv[++i]));
    }
    else if (!strcmp(argv[i], "-sc") && (i < argc - 1))
      par.columnscore = atoi(argv[++i]);
    else if (!strcmp(argv[i], "-def"))
      par.readdefaultsfile = 1;
    else if (!strcmp(argv[i], "-maxres") && (i < argc - 1)) {
      par.maxres = atoi(argv[++i]);
      par.maxcol = 2 * par.maxres;
    }
    else if (!strcmp(argv[i], "-maxmem") && (i < argc - 1)) {
      par.maxmem = atof(argv[++i]);
    }
    else if (!strcmp(argv[i], "-corr") && (i < argc - 1))
      par.corr = atof(argv[++i]);
    else if (!strcmp(argv[i], "-ovlp") && (i < argc - 1))
      par.min_overlap = atoi(argv[++i]);
    else if (!strcmp(argv[i], "-tags"))
      par.notags = 0;
    else if (!strcmp(argv[i], "-notags"))
      par.notags = 1;
    else if (!strcmp(argv[i], "-nocontxt"))
      par.nocontxt = 1;
    else if (!strcmp(argv[i], "-csb") && (i < argc - 1))
      par.csb = atof(argv[++i]);
    else if (!strcmp(argv[i], "-csw") && (i < argc - 1))
      par.csw = atof(argv[++i]);
    else if (!strcmp(argv[i], "-cs")) {
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        std::cerr << "Error in " << __FILE__ << ":" << __LINE__ << ": "
            << __func__ << ":" << std::endl;
        std::cerr << "\tno query file following -cs\n";
        exit(4);
      }
      else
        strcpy(par.clusterfile, argv[i]);
    }
    else if (!strncmp(argv[i], "-cpu", 4) && (i < argc - 1)) {
      par.threads = atoi(argv[++i]);
    }
    else {
      HH_LOG(LogLevel::WARNING) << endl << "WARNING: Ignoring unknown option "
          << argv[i] << " ...\n";
    }
    HH_LOG(LogLevel::DEBUG1) << i << "  " << argv[i] << endl;
  } // end of for-loop for command line input
}

//void wiggleQSC(Alignment& orig_qali, char query_input_format,
//    Alignment& orig_tali, char template_input_format, size_t nqsc, float* qsc,
//    HitList& recalculated_hitlist) {
//  const int COV_ABS = 25;
//  int cov_tot = std::max(
//      std::min((int) (COV_ABS / orig_qali.L * 100 + 0.5), 70), par.coverage);
//
//  Alignment qali;
//  qali = orig_qali;
//
//  Alignment tali;
//  tali = orig_tali;
//
//  HMM* q = new HMM();
//  HMM* t = new HMM();
//
//  //TODO:
//  recalculated_hitlist.N_searched = 1;
//  HitList realigned_viterbi_hitlist;
//  realigned_viterbi_hitlist.N_searched = 1;
//
//  for (size_t qsc_index = 0; qsc_index < nqsc; qsc_index++) {
//    float actual_qsc = qsc[qsc_index];
//
//    qali.Compress("filtered A3M file", par.cons, par.maxres, par.maxcol, par.M,
//        par.Mgaps);
//    qali.N_filtered = qali.Filter(par.max_seqid, S, cov_tot, par.qid,
//        actual_qsc, par.Ndiff);
//    qali.FrequenciesAndTransitions(q, par.wg, par.mark, par.cons, par.showcons,
//        par.maxres, pb, Sim, NULL, false);
//    PrepareQueryHMM(par, query_input_format, q, pc_hhm_context_engine,
//        pc_hhm_context_mode, pb, R);
//
//    tali.Compress("filtered A3M file", par.cons, par.maxres, par.maxcol, par.M,
//        par.Mgaps);
//    tali.N_filtered = qali.Filter(par.max_seqid, S, cov_tot, par.qid,
//        actual_qsc, par.Ndiff);
//    tali.FrequenciesAndTransitions(t, par.wg, par.mark, par.cons, par.showcons,
//        par.maxres, pb, Sim, NULL, false);
//    PrepareTemplateHMM(par, q, t, template_input_format, pb, R);
//
//    //run viterbi
//    Hit hit;
//    hit.AllocateBacktraceMatrix(q->L + 2, par.maxres + 1);
//    hit.self = 0;
//    hit.realign_around_viterbi = false;
//
//    for (int irep = 1; irep <= par.altali; irep++) {
//      hit.irep = irep;
//      hit.Viterbi(q, t, par.loc, par.ssm, par.maxres, par.min_overlap,
//          par.shift, par.egt, par.egq, par.ssw, par.exclstr, S73, S33);
//
//      if (hit.irep > 1 && hit.score <= SMIN)
//        break;
//
//      hit.Backtrace(q, t, par.corr, par.ssw, S73, S33);
//      realigned_viterbi_hitlist.Push(hit);
//    }
//
//    hit.DeleteBacktraceMatrix(q->L + 2);
//
//    realigned_viterbi_hitlist.CalculatePvalues(q, par.loc, par.ssm, par.ssw);
//    realigned_viterbi_hitlist.CalculateHHblitsEvalues(q, 1, par.alphaa,
//        par.alphab, par.alphac, par.prefilter_evalue_thresh);
//
//    //run mac alignment
//    q->Log2LinTransitionProbs(1.0);
//    t->Log2LinTransitionProbs(1.0);
//
//    realigned_viterbi_hitlist.Reset();
//    while (!realigned_viterbi_hitlist.End()) {
//      Hit hit_ref = realigned_viterbi_hitlist.ReadNext();
//
//      Hit hit;
//      hit = hit_ref;
//
//      hit.AllocateForwardMatrix(q->L + 2, par.maxres + 1);
//      hit.AllocateBacktraceMatrix(q->L + 2, par.maxres + 1);
//      hit.irep = 1;
//      hit.self = 0;
//      hit.i1 = hit_ref.i1;
//      hit.i2 = hit_ref.i2;
//      hit.j1 = hit_ref.j1;
//      hit.j2 = hit_ref.j2;
//      hit.nsteps = hit_ref.nsteps;
//      hit.i = hit_ref.i;
//      hit.j = hit_ref.j;
//      hit.realign_around_viterbi = false;
//
//      // Align q to template in *hit[bin]
//      hit.Forward(q, t, par.ssm, par.min_overlap, par.loc, par.shift, par.ssw,
//          par.exclstr, S73, S33);
//      hit.Backward(q, t, par.loc, par.shift, par.ssw, S73, S33);
//      hit.MACAlignment(q, t, par.loc, par.mact, par.macins);
//      hit.BacktraceMAC(q, t, par.corr, par.ssw, S73, S33);
//
//      // Overwrite *hit[bin] with Viterbi scores, Probabilities etc. of hit_cur
//      hit.score = hit_ref.score;
//      hit.score_ss = hit_ref.score_ss;
//      hit.score_aass = hit_ref.score_aass;
//      hit.score_sort = hit_ref.score_sort;
//      hit.Pval = hit_ref.Pval;
//      hit.Pvalt = hit_ref.Pvalt;
//      hit.logPval = hit_ref.logPval;
//      hit.logPvalt = hit_ref.logPvalt;
//      hit.Eval = hit_ref.Eval;
//      hit.logEval = hit_ref.logEval;
//      hit.Probab = hit_ref.Probab;
//
//      hit.DeleteForwardMatrix(q->L + 2);
//      hit.DeleteBacktraceMatrix(q->L + 2);
//
//      std::cout << hit.sum_of_probs << std::endl;
//
//      if (hit.matched_cols >= MINCOLS_REALIGN) {
//        recalculated_hitlist.Insert(hit);
//      }
//    }
//
//    realigned_viterbi_hitlist.Reset();
//    while (!realigned_viterbi_hitlist.End()) {
//      realigned_viterbi_hitlist.Delete();
//    }
//  }
//
//  recalculated_hitlist.SortList(&Hit::compare_sum_of_probs);
//}

void HHalign::optimizeQSC(HitList& input_list, HMMSimd& q_vec, char query_input_format, HitList& output_list) {
	const int COV_ABS = 25;
	const int cov_tot = std::max(std::min((int) (COV_ABS / q->L * 100 + 0.5), 70), par.coverage);

	const int nqsc = 6;
	float qscs[nqsc] = {-20, 0, 0.1, 0.2, 0.3, 0.4};


	input_list.Reset();
	while (!input_list.End()) {
		Hit hit_cur = input_list.ReadNext();

		std::vector<HHDatabaseEntry*> selected_entries;
		selected_entries.push_back(hit_cur.entry);

		float best_alignment_quality = -FLT_MAX;
		HitList best_alignments;

		for(int i = 0; i < nqsc; i++) {
			float actual_qsc = qscs[i];

			HitList tmp_list;
			tmp_list.N_searched = input_list.N_searched;

			Qali.Compress("filtered A3M file", par.cons, par.maxres, par.maxcol, par.M, par.Mgaps);
			Qali.N_filtered = Qali.Filter(par.max_seqid, S, cov_tot, par.qid, actual_qsc, par.Ndiff);
			Qali.FrequenciesAndTransitions(q, par.wg, par.mark, par.cons, par.showcons, par.maxres, pb, Sim, NULL, false);
			PrepareQueryHMM(par, query_input_format, q, pc_hhm_context_engine, pc_hhm_context_mode, pb, R);

			q_vec.MapOneHMM(q);

			ViterbiRunner viterbirunner(viterbiMatrices, dbs, par.threads);
			std::vector<Hit> hits_to_add = viterbirunner.alignment(par, &q_vec, selected_entries, actual_qsc, pb, S, Sim, R);

			add_hits_to_hitlist(hits_to_add, tmp_list);

			std::vector<Hit *> hit_vector;
			std::vector<HHDatabaseEntry*> hits_to_realign;

			int t_maxres = 0;
			int n_realignments = 0;
			tmp_list.Reset();
			while (!tmp_list.End()) {
				Hit hit_cur = tmp_list.ReadNext();
				t_maxres = std::max(t_maxres, hit_cur.L + 2);
				hits_to_realign.push_back(hit_cur.entry);
				hit_vector.push_back(tmp_list.ReadCurrentAddress());
				n_realignments++;
			}

			std::qsort(&hit_vector[0], hit_vector.size(), sizeof(Hit*), compareHitLengths);

			std::map<short int, std::vector<Hit *> > alignments;
			for (int elem = 0; elem < (int) hit_vector.size(); elem++) {
				alignments[hit_vector.at(elem)->irep].push_back(hit_vector.at(elem));
			}

			PosteriorDecoderRunnerInputData input_data(dbs, hits_to_realign, alignments, n_realignments, t_maxres);

			// Initialize a Null-value as a return value if not items are available anymore
			PosteriorDecoderRunner runner(input_data, q_vec, posteriorMatrices, viterbiMatrices, par.threads);

			runner.executeComputation(par, actual_qsc, pb, S, Sim, R);

			//check if actual qsc produced the best alignment so far
			bool new_best_hit = false;

			tmp_list.Reset();
			while (!tmp_list.End()) {
				Hit hit_cur = tmp_list.ReadNext();
				float cur_ali_quality = hit_cur.estimateAlignmentQuality(q);
				std::cout << hit_cur.name << "_" << hit_cur.irep << "\t" << actual_qsc << "\t" << cur_ali_quality << std::endl;
				if(cur_ali_quality > best_alignment_quality) {
					best_alignment_quality = cur_ali_quality;
					new_best_hit = true;
					break;
				}
			}

			if(new_best_hit) {
				//delete old hits
				best_alignments.Reset();
				while (!best_alignments.End()) {
					best_alignments.Delete().Delete();
				}

				//copy new hits
				tmp_list.Reset();
				while (!tmp_list.End()) {
					Hit hit_cur = tmp_list.ReadNext();
					best_alignments.Push(hit_cur);
				}
			}
			else {
				//delete hits for actual qsc
				tmp_list.Reset();
				while (!tmp_list.End()) {
					tmp_list.Delete().Delete();
				}
			}
		}

		//copy best qsc hit to output_list
		best_alignments.Reset();
		while (!best_alignments.End()) {
			Hit hit_cur = best_alignments.ReadNext();
			output_list.Push(hit_cur);
			best_alignments.Delete();
		}
	}
}

void HHalign::run(FILE* query_fh, char* query_path, std::vector<std::string>& templates) {
  int cluster_found = 0;
  int seqs_found = 0;
  int premerge = par.premerge;

  Hit hit_cur;
  Hash<Hit>* previous_hits = new Hash<Hit>(1631, hit_cur);
  Hash<char>* premerged_hits = new Hash<char>(1631);

  q = new HMM;
  HMMSimd q_vec(par.maxres);
  q_tmp = new HMM;

  // Read input file (HMM, HHM, or alignment format), and add pseudocounts etc.
  Qali.N_in = 0;
  char input_format = 0;
  ReadQueryFile(par, query_fh, input_format, par.wg, q, Qali, query_path, pb, S, Sim);
  PrepareQueryHMM(par, input_format, q, pc_hhm_context_engine, pc_hhm_context_mode, pb, R);
  q_vec.MapOneHMM(q);

  // Set query columns in His-tags etc to Null model distribution
  if (par.notags)
    q->NeutralizeTags(pb);

  std::vector<HHDatabaseEntry*> selected_entries;
  for (size_t i = 0; i < dbs.size(); i++) {
    dbs[i]->initSelected(templates, selected_entries);
  }

  ViterbiRunner viterbirunner(viterbiMatrices, dbs, par.threads);
  std::vector<Hit> hits_to_add = viterbirunner.alignment(par, &q_vec, selected_entries, par.qsc_db, pb, S, Sim, R);

  hitlist.N_searched = selected_entries.size();
  add_hits_to_hitlist(hits_to_add, hitlist);

  // Set new ss weight for realign
  par.ssw = par.ssw_realign;

  // Realign hits with MAC algorithm
  if (par.realign && par.forward != 2) {
    perform_realign(q_vec, selected_entries, premerge, premerged_hits);
  }

  mergeHitsToQuery(previous_hits, premerged_hits, seqs_found, cluster_found);

  // Calculate pos-specific weights, AA frequencies and transitions -> f[i][a], tr[i][a]
  Qali.FrequenciesAndTransitions(q, par.wg, par.mark, par.cons, par.showcons, par.maxres, pb, Sim, NULL, true);

  if (par.notags)
      q->NeutralizeTags(pb);

  //DEBUG
  HitList opt_hitlist;
  optimizeQSC(hitlist, q_vec, input_format, opt_hitlist);

  for(size_t i = 0; i < selected_entries.size(); i++) {
    delete selected_entries[i];
  }
  selected_entries.clear();
}

void HHalign::writeHHRFile(HHalign& hhalign, std::stringstream& out) {
	hhalign.hitlist.PrintHHR(hhalign.q, out, hhalign.par.maxdbstrlen,
			hhalign.par.showconf, hhalign.par.showcons, hhalign.par.showdssp,
			hhalign.par.showpred, hhalign.par.b, hhalign.par.B, hhalign.par.z,
			hhalign.par.Z, hhalign.par.aliwidth, hhalign.par.nseqdis, hhalign.par.p,
			hhalign.par.E, hhalign.par.argc, hhalign.par.argv, hhalign.S);
}
