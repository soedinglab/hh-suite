// hhalign.C:
// Align a multiple alignment to an alignment or HMM
// Print out aligned input sequences in a3m format
// Error codes: 0: ok  1: file format error  2: file access error  3: memory error  4: internal numeric error  5: command line error

/*
     (C) Johannes Soeding 2012

     This program is free software: you can redistribute it and/or modify
     it under the terms of the GNU General Public License as published by
     the Free Software Foundation, either version 3 of the License, or
     (at your option) any later version.

     This program is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU General Public License for more details.

     You should have received a copy of the GNU General Public License
     along with this program.  If not, see <http://www.gnu.org/licenses/>.

     We are very grateful for bug reports! Please contact us at soeding@mpibpc.mpg.de

     Reference:
     Remmert M., Biegert A., Hauser A., and Soding J.
     HHblits: Lightning-fast iterative protein sequence searching by HMM-HMM alignment.
     Nat. Methods, epub Dec 25, doi: 10.1038/NMETH.1818 (2011).
*/

#include "hhalign.h"
#include "hhsuite_config.h"

std::vector<HHblitsDatabase*> empty;
HHalign::HHalign(Parameters& par) : HHblits(par, empty), tfiles(par.tfiles) {

}

HHalign::~HHalign() {

}


void HHalign::help(Parameters& par, char all) {
  printf("HHalign %i.%i.%i\n", HHSUITE_VERSION_MAJOR, HHSUITE_VERSION_MINOR, HHSUITE_VERSION_PATCH);
  printf("Align a query alignment/HMM to a template alignment/HMM by HMM-HMM alignment\n");
  printf("If only one alignment/HMM is given it is compared to itself and the best\n");
  printf("off-diagonal alignment plus all further non-overlapping alignments above \n");
  printf("significance threshold are shown.\n");
  printf("%s", REFERENCE);
  printf("%s", COPYRIGHT);
  printf("\n");

  printf("Usage: hhalign -i query -t template [options]  \n");
  printf(" -i <file>      input/query: single sequence or multiple sequence alignment (MSA)\n");
  printf("                in a3m, a2m, or FASTA format, or HMM in hhm format\n");
  printf(" -t <file>      input/template: single sequence or multiple sequence alignment (MSA)\n");
  printf("                in a3m, a2m, or FASTA format, or HMM in hhm format\n");

  if (all) {
    printf("\n");
    printf("<file> may be 'stdin' or 'stdout' throughout.\n");
  }

  printf("\n");

  printf("Input alignment format:                                                       \n");
  printf(" -M a2m         use A2M/A3M (default): upper case = Match; lower case = Insert;\n");
  printf("               '-' = Delete; '.' = gaps aligned to inserts (may be omitted)   \n");
  printf(" -M first       use FASTA: columns with residue in 1st sequence are match states\n");
  printf(" -M [0,100]     use FASTA: columns with fewer than X%% gaps are match states   \n");
  printf(" -tags/-notags  do NOT / do neutralize His-, C-myc-, FLAG-tags, and trypsin \n");
  printf("                recognition sequence to background distribution (def=-notags)  \n");
  printf("\n");

  printf("Output options: \n");
  printf(" -o <file>      write results in standard format to file (default=<infile.hhr>)\n");
  printf(" -oa3m <file>   write query alignment in a3m or PSI-BLAST format (-opsi) to file (default=none)\n");
  printf(" -aa3m <file>   append query alignment in a3m (-aa3m) or PSI-BLAST format (-apsi )to file (default=none)\n");
  printf(" -Ofas <file>   write pairwise alignments in FASTA xor A2M (-Oa2m) xor A3M (-Oa3m) format   \n");

  printf(" -add_cons      generate consensus sequence as master sequence of query MSA (default=don't)\n");
  printf(" -hide_cons     don't show consensus sequence in alignments (default=show)     \n");
  printf(" -hide_pred     don't show predicted 2ndary structure in alignments (default=show)\n");
  printf(" -hide_dssp     don't show DSSP 2ndary structure in alignments (default=show)  \n");
  printf(" -show_ssconf   show confidences for predicted 2ndary structure in alignments\n");

  if (all) {
    printf(" -seq <int>     max. number of query/template sequences displayed (default=%i)  \n", par.nseqdis);
    printf(" -aliw <int>    number of columns per line in alignment list (default=%i)       \n", par.aliwidth);
    printf(" -p [0,100]     minimum probability in summary and alignment list (default=%G)  \n", par.p);
    printf(" -E [0,inf[     maximum E-value in summary and alignment list (default=%G)      \n", par.E);
    printf(" -Z <int>       maximum number of lines in summary hit list (default=%i)        \n", par.Z);
    printf(" -z <int>       minimum number of lines in summary hit list (default=%i)        \n", par.z);
    printf(" -B <int>       maximum number of alignments in alignment list (default=%i)     \n", par.B);
    printf(" -b <int>       minimum number of alignments in alignment list (default=%i)     \n", par.b);
  }
  printf("\n");

  printf("Filter options applied to query MSA, template MSA, and result MSA              \n");
  printf(" -id   [0,100]  maximum pairwise sequence identity (def=%i)\n", par.max_seqid);
  printf(" -diff [0,inf[  filter MSAs by selecting most diverse set of sequences, keeping \n");
  printf("                at least this many seqs in each MSA block of length 50 \n");
  printf("                Zero and non-numerical values turn off the filtering. (def=%i) \n", par.Ndiff);
  printf(" -cov  [0,100]  minimum coverage with master sequence (%%) (def=%i)             \n", par.coverage);
  printf(" -qid  [0,100]  minimum sequence identity with master sequence (%%) (def=%i)    \n", par.qid);
  printf(" -qsc  [0,100]  minimum score per column with master sequence (default=%.1f)    \n", par.qsc);
  printf(" -mark          do not filter out sequences marked by \">@\"in their name line  \n");
  printf("\n");

  printf("HMM-HMM alignment options:                                                       \n");
  printf(" -norealign     do NOT realign displayed hits with MAC algorithm (def=realign)   \n");
  printf(" -mact [0,1[    posterior prob threshold for MAC realignment controlling greedi- \n");
  printf("                ness at alignment ends: 0:global >0.1:local (default=%.2f)       \n", par.mact);
  printf(" -glob/-loc     use global/local alignment mode for searching/ranking (def=local)\n");

  if (all) {
    printf(" -realign       realign displayed hits with max. accuracy (MAC) algorithm \n");
    printf(" -excl <range>  exclude query positions from the alignment, e.g. '1-33,97-168' \n");
    printf(" -template_excl <range>  exclude template positions from the alignment, e.g. '1-33,97-168' \n");
    printf(" -ovlp <int>    banded alignment: forbid <ovlp> largest diagonals |i-j| of DP matrix (def=%i)\n", par.min_overlap);
    printf(" -alt <int>     show up to this many alternative alignments with raw score > smin(def=%i)  \n", par.altali);
    printf(" -smin <float>  minimum raw score for alternative alignments (def=%.1f)  \n", par.smin);
    printf(" -shift [-1,1]  profile-profile score offset (def=%-.2f)                         \n", par.shift);
    printf(" -corr [0,1]    weight of term for pair correlations (def=%.2f)                \n", par.corr);
    printf(" -sc   <int>    amino acid score         (tja: template HMM at column j) (def=%i)\n", par.columnscore);
    printf("        0       = log2 Sum(tja*qia/pa)   (pa: aa background frequencies)    \n");
    printf("        1       = log2 Sum(tja*qia/pqa)  (pqa = 1/2*(pa+ta) )               \n");
    printf("        2       = log2 Sum(tja*qia/ta)   (ta: av. aa freqs in template)     \n");
    printf("        3       = log2 Sum(tja*qia/qa)   (qa: av. aa freqs in query)        \n");
    printf("        5       local amino acid composition correction                     \n");
    printf(" -ssm {0,..,4}  secondary structure scoring [default=%1i]             \n", par.ssm);
    printf("          0:    = no ss scoring           \n");
    printf("        1,2:    = ss scoring after or during alignment         \n");
    printf("        3,4:    = ss scoring after or during alignment, predicted vs. predicted\n");
    printf(" -ssw [0,1]     weight of ss score  (def=%-.2f)                                  \n", par.ssw);
    printf(" -ssa [0,1]     ss confusion matrix = (1-ssa)*I + ssa*psipred-confusion-matrix [def=%-.2f)\n", par.ssa);
    printf(" -wg            use global sequence weighting for realignment!                   \n");
    printf("\n");

    printf("Gap cost options:                                                                \n");
    printf(" -gapb [0,inf[  Transition pseudocount admixture (def=%-.2f)                     \n", par.gapb);
    printf(" -gapd [0,inf[  Transition pseudocount admixture for open gap (default=%-.2f)    \n", par.gapd);
    printf(" -gapd [0,inf[  Transition pseudocount admixture for open gap (default=%-.2f)          \n", par.gapd);
    printf(" -gape [0,1.5]  Transition pseudocount admixture for extend gap (def=%-.2f)            \n", par.gape);
    printf(" -gapf ]0,inf]  factor to increase/reduce the gap open penalty for deletes (def=%-.2f) \n", par.gapf);
    printf(" -gapg ]0,inf]  factor to increase/reduce the gap open penalty for inserts (def=%-.2f) \n", par.gapg);
    printf(" -gaph ]0,inf]  factor to increase/reduce the gap extend penalty for deletes(def=%-.2f)\n", par.gaph);
    printf(" -gapi ]0,inf]  factor to increase/reduce the gap extend penalty for inserts(def=%-.2f)\n", par.gapi);
    printf(" -egq  [0,inf[  penalty (bits) for end gaps aligned to query residues (def=%-.2f)      \n", par.egq);
    printf(" -egt  [0,inf[  penalty (bits) for end gaps aligned to template residues (def=%-.2f)   \n", par.egt);
    printf("\n");

    printf("Pseudocount (pc) options:                                                        \n");
    printf(" Context specific hhm pseudocounts:\n");
    printf("  -pc_hhm_contxt_mode {0,..,3}      position dependence of pc admixture 'tau' (pc mode, default=%-i) \n", par.pc_hhm_context_engine.admix);
    printf("               0: no pseudo counts:    tau = 0                                  \n");
    printf("               1: constant             tau = a                                  \n");
    printf("               2: diversity-dependent: tau = a/(1+((Neff[i]-1)/b)^c)            \n");
    printf("               3: CSBlast admixture:   tau = a(1+b)/(Neff[i]+b)                 \n");
    printf("               (Neff[i]: number of effective seqs in local MSA around column i) \n");
    printf("  -pc_hhm_contxt_a  [0,1]        overall pseudocount admixture (def=%-.1f)                        \n", par.pc_hhm_context_engine.pca);
    printf("  -pc_hhm_contxt_b  [1,inf[      Neff threshold value for mode 2 (def=%-.1f)                      \n", par.pc_hhm_context_engine.pcb);
    printf("  -pc_hhm_contxt_c  [0,3]        extinction exponent c for mode 2 (def=%-.1f)                     \n", par.pc_hhm_context_engine.pcc);
    printf("\n");

    printf(" Context independent hhm pseudocounts (used for templates; used for query if contxt file is not available):\n");
    printf("  -pc_hhm_nocontxt_mode {0,..,3}      position dependence of pc admixture 'tau' (pc mode, default=%-i) \n", par.pc_hhm_nocontext_mode);
    printf("               0: no pseudo counts:    tau = 0                                  \n");
    printf("               1: constant             tau = a                                  \n");
    printf("               2: diversity-dependent: tau = a/(1+((Neff[i]-1)/b)^c)            \n");
    printf("               (Neff[i]: number of effective seqs in local MSA around column i) \n");
    printf("  -pc_hhm_nocontxt_a  [0,1]        overall pseudocount admixture (def=%-.1f)                        \n", par.pc_hhm_nocontext_a);
    printf("  -pc_hhm_nocontxt_b  [1,inf[      Neff threshold value for mode 2 (def=%-.1f)                      \n", par.pc_hhm_nocontext_b);
    printf("  -pc_hhm_nocontxt_c  [0,3]        extinction exponent c for mode 2 (def=%-.1f)                     \n", par.pc_hhm_nocontext_c);
    printf("\n");

    printf(" Context-specific pseudo-counts:                                                  \n");
    printf("  -nocontxt      use substitution-matrix instead of context-specific pseudocounts \n");
    printf("  -contxt <file> context file for computing context-specific pseudocounts (default=%s)\n", par.clusterfile.c_str());
    printf("  -csw  [0,inf]  weight of central position in cs pseudocount mode (def=%.1f)\n", par.csw);
    printf("  -csb  [0,1]    weight decay parameter for positions in cs pc mode (def=%.1f)\n", par.csb);
  }
  printf("\n");

  printf("Other options:                                                                   \n");
  printf(" -v <int>       verbose mode: 0:no screen output  1:only warings  2: verbose (def=%i)\n", par.v);
  if (all) {
    printf(" -atab   <file> write all alignments in tabular layout to file                   \n");
    printf(" -maxseq <int>  max number of input rows (def=%5i)\n", par.maxseq);
    printf(" -maxres <int>  max number of HMM columns (def=%5i)\n", par.maxres);
    printf(" -maxmem [1,inf[ limit memory for realignment (in GB) (def=%.1f)          \n", par.maxmem);
  }
  printf("\n");

  if (!all) {
    printf("An extended list of options can be obtained by calling 'hhalign -h all'\n");
  }

  printf("Example: hhalign -i T0187.a3m -t d1hz4a_.hhm -o result.hhr\n");
}

/////////////////////////////////////////////////////////////////////////////////////
//// Processing input options from command line
/////////////////////////////////////////////////////////////////////////////////////
void HHalign::ProcessAllArguments(Parameters& par) {
  par.tfiles.resize(0);
  strcpy(par.alnfile, "");
  par.p = 0.0; // minimum threshold for inclusion in hit list and alignment listing
  par.E = 1e6; // maximum threshold for inclusion in hit list and alignment listing
  par.b = 1;                     // min number of alignments
  par.B = 100;                   // max number of alignments
  par.z = 1;                     // min number of lines in hit list
  par.Z = 100;                   // max number of lines in hit list
  par.append = 0;              // append alignment to output file with -a option
  par.altali = 1;           // find only ONE (possibly overlapping) subalignment
  par.outformat = 3;             // default output format for alignment is a3m
  par.realign = 1;               // default: realign

  par.num_rounds = 1;

  ProcessArguments(par);

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
  if (par.mact >= 1.0)
    par.mact = 0.999;
  else if (par.mact < 0)
    par.mact = 0.0;
}

void HHalign::ProcessArguments(Parameters& par) {
  const int argc = par.argc;
  const char** argv = par.argv;

  for (int i = 1; i < argc; i++) {
    HH_LOG(DEBUG1) << i << "  " << argv[i] << std::endl;
    if (!strcmp(argv[i], "-i")) {
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        HH_LOG(ERROR) << "No query file following -i" << std::endl;
        exit(4);
      }
      else
        strcpy(par.infile, argv[i]);
    }
    else if (!strcmp(argv[i], "-t")) {
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        HH_LOG(ERROR) << "No template file following -t" << std::endl;
        exit(4);
      }
      else {
        std::string tfile(argv[i]);
        par.tfiles.push_back(tfile);
      }
    }
    else if (!strcmp(argv[i], "-o")) {
      if (++i >= argc) {
        help(par);
        HH_LOG(ERROR) << "No filename following -o" << std::endl;
        exit(4);
      }
      else
        strcpy(par.outfile, argv[i]);
    }
    else if (!strcmp(argv[i], "-Ofas")) {
      par.outformat = 1;
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        HH_LOG(ERROR) << "No output file following -Ofas" << std::endl;
        exit(4);
      }
      else
        strcpy(par.pairwisealisfile, argv[i]);
    }
    else if (!strcmp(argv[i], "-Oa2m")) {
      par.outformat = 2;
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        HH_LOG(ERROR) << "No output file following -Oa2m" << std::endl;
        exit(4);
      }
      else
        strcpy(par.pairwisealisfile, argv[i]);
    }
    else if (!strcmp(argv[i], "-Oa3m")) {
      par.outformat = 3;
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        HH_LOG(ERROR) << "No output file following -Oa3m" << std::endl;
        exit(4);
      }
      else
        strcpy(par.pairwisealisfile, argv[i]);
    }
    else if (!strcmp(argv[i], "-oa3m")) {
      par.append = 0;
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        HH_LOG(ERROR) << "No output file following -oa3m" << std::endl;
        exit(4);
      }
      else
        strcpy(par.alnfile, argv[i]);
    }
    else if (!strcmp(argv[i], "-opsi")) {
      par.append = 0;
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        HH_LOG(ERROR) << "No output file following -opsi" << std::endl;
        exit(4);
      }
      else
        strcpy(par.psifile, argv[i]);
    }
    else if (!strcmp(argv[i], "-aa3m")) {
      par.append = 1;
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        HH_LOG(ERROR) << "No output file following -aa3m" << std::endl;
        exit(4);
      }
      else
        strcpy(par.alnfile, argv[i]);
    }
    else if (!strcmp(argv[i], "-apsi")) {
      par.append = 1;
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        HH_LOG(ERROR) << "No output file following -apsi" << std::endl;
        exit(4);
      }
      else
        strcpy(par.psifile, argv[i]);
    }
    else if (!strcmp(argv[i], "-wg")) {
      par.wg = 1;
    }
    else if (!strcmp(argv[i], "-atab")) {
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        HH_LOG(ERROR) << "No query file following -atab" << std::endl;
        exit(4);
      }
      else
        strmcpy(par.alitabfile, argv[i], NAMELEN - 1);
    }
    else if (!strcmp(argv[i], "-index")) {
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        HH_LOG(ERROR) << "No index file following -index" << std::endl;
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
    else if (!strncmp(argv[i], "-add_cons", 7))
      par.cons = 0;
    else if (!strncmp(argv[i], "-hide_cons", 7))
      par.showcons = 0;
    else if (!strncmp(argv[i], "-hide_pred", 7))
      par.showpred = 0;
    else if (!strncmp(argv[i], "-hide_dssp", 7))
      par.showdssp = 0;
    else if (!strncmp(argv[i], "-show_ssconf", 7))
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

    //no help required
    else if (!strcmp(argv[i], "-Gonnet"))
      par.matrix = 0;
    //no help required
    else if (!strcmp(argv[i], "-Blosum50"))
      par.matrix = 50;
    //no help required
    else if (!strcmp(argv[i], "-Blosum62"))
      par.matrix = 62;

    else if (!strcmp(argv[i], "-pc_hhm_contxt_mode") && (i < argc - 1))
      par.pc_hhm_context_engine.admix = (Pseudocounts::Admix) atoi(argv[++i]);
    else if (!strcmp(argv[i], "-pc_hhm_contxt_a") && (i < argc - 1))
      par.pc_hhm_context_engine.pca = atof(argv[++i]);
    else if (!strcmp(argv[i], "-pc_hhm_contxt_b") && (i < argc - 1))
      par.pc_hhm_context_engine.pcb = atof(argv[++i]);
    else if (!strcmp(argv[i], "-pc_hhm_contxt_c") && (i < argc - 1))
      par.pc_hhm_context_engine.pcc = atof(argv[++i]);
    else if (!strcmp(argv[i], "-pc_hhm_nocontxt_mode") && (i < argc - 1))
      par.pc_hhm_nocontext_mode = atoi(argv[++i]);
    else if (!strcmp(argv[i], "-pc_hhm_nocontxt_a") && (i < argc - 1))
      par.pc_hhm_nocontext_a = atof(argv[++i]);
    else if (!strcmp(argv[i], "-pc_hhm_nocontxt_b") && (i < argc - 1))
      par.pc_hhm_nocontext_b = atof(argv[++i]);
    else if (!strcmp(argv[i], "-pc_hhm_nocontxt_c") && (i < argc - 1))
      par.pc_hhm_nocontext_c = atof(argv[++i]);
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
    else if (!strcmp(argv[i], "-ssm") && (i < argc - 1))
      par.ssm = atoi(argv[++i]);
    else if (!strcmp(argv[i], "-ssw") && (i < argc - 1))
      par.ssw = atof(argv[++i]);
    else if (!strcmp(argv[i], "-ssa") && (i < argc - 1))
      par.ssa = atof(argv[++i]);
    else if (!strncmp(argv[i], "-glob", 5)) {
      par.loc = 0;
      if (par.mact > 0.35 && par.mact < 0.3502) {
        par.mact = 0;
      }
    }
    else if (!strncmp(argv[i], "-loc", 3))
      par.loc = 1;
    else if (!strncmp(argv[i], "-alt", 4) && (i < argc - 1))
      par.altali = atoi(argv[++i]);
    else if (!strncmp(argv[i], "-smin", 4) && (i < argc - 1))
      par.smin = atof(argv[++i]);
    else if (!strcmp(argv[i], "-realign"))
      par.realign = 1;
    else if (!strcmp(argv[i], "-norealign"))
      par.realign = 0;
    else if (!strcmp(argv[i], "-M") && (i < argc - 1)) {
      //TODO: M a3m not defined in the help
      if (!strcmp(argv[++i], "a2m") || !strcmp(argv[i], "a3m")){
        par.M = 1;
        par.M_template = 1;
      }
      else if (!strcmp(argv[i], "first")) {
        par.M = 3;
        par.M_template = 3;
      }
      else if (argv[i][0] >= '0' && argv[i][0] <= '9') {
        par.Mgaps = atoi(argv[i]);
        par.M = 2;
        par.M_template = 2;
      }
      else
        HH_LOG(WARNING) << "Ignoring unknown argument: -M " << argv[i] << std::endl;
    }
    else if (!strcmp(argv[i], "-shift") && (i < argc - 1))
      par.shift = atof(argv[++i]);
    else if (!strcmp(argv[i], "-mact") && (i < argc - 1)) {
      par.mact = atof(argv[++i]);
    }
    //not in the help - but intended
    else if (!strcmp(argv[i], "-scwin") && (i < argc - 1)) {
      par.columnscore = 5;
      par.half_window_size_local_aa_bg_freqs = imax(1, atoi(argv[++i]));
    }
    else if (!strcmp(argv[i], "-sc") && (i < argc - 1))
      par.columnscore = atoi(argv[++i]);
    else if (!strcmp(argv[i], "-maxseq") && (i < argc - 1))
      par.maxseq = atoi(argv[++i]);
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
    else if (!strcmp(argv[i], "-contxt")) {
      if (++i >= argc || argv[i][0] == '-') {
        help(par);
        HH_LOG(ERROR) << "No query file following -contxt" << std::endl;
        exit(4);
      }
      else
        par.clusterfile = argv[i];
    }
    else if (!strcmp(argv[i],"-excl")) {
      if (++i>=argc) {
        help(par);
        HH_LOG(ERROR) << "No expression following -excl" << std::endl;
        exit(4);
      }
      par.exclstr = new char[strlen(argv[i])+1];
      strcpy(par.exclstr,argv[i]);
    }
    else if (!strcmp(argv[i],"-template_excl")) {
      if (++i>=argc) {
        help(par);
        HH_LOG(ERROR) << "No expression following -excl" << std::endl;
        exit(4);
      }
      par.template_exclstr = new char[strlen(argv[i])+1];
      strcpy(par.template_exclstr,argv[i]);
    }
    else {
      HH_LOG(WARNING) << "Ignoring unknown option " << argv[i] << std::endl;
    }
    HH_LOG(DEBUG1) << i << "  " << argv[i] << std::endl;
  } // end of for-loop for command line input
}

void HHalign::run(FILE* query_fh, char* query_path) {
  HH_LOG(DEBUG) << "Query file : " << query_path << "\n";
  HH_LOG(DEBUG) << "Template files:";
  for (size_t i = 0; i < tfiles.size(); ++i) {
    HH_LOG(DEBUG) << " " << tfiles[i];
  }
  HH_LOG(DEBUG) << "\n";

  int cluster_found = 0;
  int seqs_found = 0;

  Hit hit_cur;
  Hash<Hit>* previous_hits = new Hash<Hit>(1631, hit_cur);

  Qali = new Alignment(par.maxseq, par.maxres);
  Qali_allseqs = new Alignment(par.maxseq, par.maxres);

  q = new HMM(MAXSEQDIS, par.maxres);
  HMMSimd q_vec(par.maxres);
  q_tmp = new HMM(MAXSEQDIS, par.maxres);

  // Read input file (HMM, HHM, or alignment format), and add pseudocounts etc.
  Qali->N_in = 0;
  char input_format = 0;
  ReadQueryFile(par, query_fh, input_format, par.wg, q, Qali, query_path, pb,
          S, Sim);
  PrepareQueryHMM(par, input_format, q, pc_hhm_context_engine,
          pc_hhm_context_mode, pb, R);
  q_vec.MapOneHMM(q);
  *q_tmp = *q;

  // Set query columns in His-tags etc to Null model distribution
  if (par.notags)
    q->NeutralizeTags(pb);

  std::vector<HHEntry*> new_entries;
  for (size_t i = 0; i < tfiles.size(); i++) {
    HHEntry* template_entry = new HHFileEntry(tfiles[i].c_str(), par.maxres);
    new_entries.push_back(template_entry);
  }

  int max_template_length = getMaxTemplateLength(new_entries);
  for(int i = 0; i < par.threads; i++) {
    viterbiMatrices[i]->AllocateBacktraceMatrix(q->L, max_template_length);
  }

  ViterbiRunner viterbirunner(viterbiMatrices, dbs, par.threads);
  std::vector<Hit> hits_to_add = viterbirunner.alignment(par, &q_vec, new_entries, par.qsc_db, pb, S, Sim, R, par.ssm, S73, S33, S37);

  hitlist.N_searched = new_entries.size();
  add_hits_to_hitlist(hits_to_add, hitlist);

  // Set new ss weight for realign
  par.ssw = par.ssw_realign;

  // Realign hits with MAC algorithm
  if (par.realign) {
      perform_realign(q_vec, input_format, new_entries);
  }

  mergeHitsToQuery(previous_hits, seqs_found, cluster_found);

  // Calculate pos-specific weights, AA frequencies and transitions -> f[i][a], tr[i][a]
  Qali->FrequenciesAndTransitions(q, par.wg, par.mark, par.cons, par.showcons, pb, Sim, NULL, true);

  if (par.notags)
      q->NeutralizeTags(pb);

  for(size_t i = 0; i < new_entries.size(); i++) {
    delete new_entries[i];
  }
  new_entries.clear();

  previous_hits->Reset();
  while (!previous_hits->End())
    previous_hits->ReadNext().Delete(); // Delete hit object
  delete previous_hits;
}
