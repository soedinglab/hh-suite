// hhsearch.cpp:

#include "hhsearch.h"

HHsearch::HHsearch(Parameters& par, std::vector<HHblitsDatabase*>& databases) :
		HHblits(par, databases) {
}

HHsearch::~HHsearch() {

}

void HHsearch::ProcessAllArguments(int argc, char** argv, Parameters& par) {
	par.argv = argv;
	par.argc = argc;

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
	if (!*par.infile || !strcmp(par.infile, "")) {
		help(par);
		std::cerr << "Error in " << __FILE__ << ":" << __LINE__ << ": "
				<< __func__ << ":" << std::endl;
		std::cerr << "\tinput file missing!" << std::endl;
		exit(4);
	}

	if (par.db_bases.size() == 0) {
		help(par);
		std::cerr << "Error in " << __FILE__ << ":" << __LINE__ << ": "
				<< __func__ << ":" << std::endl;
		std::cerr << "\tdatabase missing (see -d)\n";
		exit(4);
	}

	if (!par.nocontxt) {
		if (!strcmp(par.clusterfile, "")) {
			help(par);
			std::cerr << "Error in " << __FILE__ << ":" << __LINE__ << ": "
					<< __func__ << ":" << std::endl;
			std::cerr << "\tcontext-specific library missing (see -contxt)\n";
			exit(4);
		}
		if (!strcmp(par.cs_library, "")) {
			help(par);
			std::cerr << "Error in " << __FILE__ << ":" << __LINE__ << ": "
					<< __func__ << ":" << std::endl;
			std::cerr << "\tcolumn state library (see -cslib)\n";
			exit(4);
		}
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
	if (par.maxmem < 1.0) {
		cerr
				<< "Warning: setting -maxmem to its minimum allowed value of 1.0\n";
		par.maxmem = 1.0;
	}
	if (par.mact >= 1.0)
		par.mact = 0.999;
	else if (par.mact < 0)
		par.mact = 0.0;
	if (par.macins >= 1.0)
		par.macins = 0.999;
	else if (par.macins < 0)
		par.macins = 0.0;
}

void HHsearch::help(Parameters& par, char all) {
	printf("\n");
	printf("HHsearch %s\n", VERSION_AND_DATE);
	printf("Search a database of HMMs with a query alignment or query HMM\n");
	printf("%s", COPYRIGHT);
	printf("%s", HHSEARCH_REFERENCE);
	printf("\n");
	printf(
			"Usage: hhsearch -i query -d database [options]                       \n");
	printf(
			" -i <file>      input/query multiple sequence alignment (a2m, a3m, FASTA) or HMM\n");
	printf(
			" -d <file>      HMM database of concatenated HMMs in hhm, HMMER, or a3m format,\n");
	printf(
			"                OR, if file has extension pal, list of HMM file names, one per\n");
	printf(
			"                line. Multiple dbs, HMMs, or pal files with -d '<db1> <db2>...'\n");
	if (all) {
		printf("\n");
		printf("<file> may be 'stdin' or 'stdout' throughout.\n");
	}
	printf("\n");
	printf(
			"Output options:                                                              \n");
	printf(
			" -o <file>      write results in standard format to file (default=<infile.hhr>)\n");
	if (all) {
		printf(
				" -Ofas <file>   write pairwise alignments of significant matches in FASTA format\n");
		printf(
				"                Analogous for output in a3m, a2m, and psi format (e.g. -Oa3m)\n");
		printf(
				" -oa3m <file>   write MSA of significant matches in a3m format\n");
		printf(
				"                Analogous for output in a2m, psi, and hhm format (e.g. -ohhm)\n");
		printf(
				" -e [0,1]       E-value cutoff for inclusion in multiple alignment (def=%G)    \n",
				par.e);
		printf(
				" -seq <int>     max. number of query/template sequences displayed (def=%i) \n",
				par.nseqdis);
		printf(
				"                Beware of overflows! All these sequences are stored in memory.\n");
		printf(
				" -cons          show consensus sequence as master sequence of query MSA \n");
	}
	printf(
			" -nocons        don't show consensus sequence in alignments (default=show)     \n");
	printf(
			" -nopred        don't show predicted 2ndary structure in alignments (default=show)\n");
	printf(
			" -nodssp        don't show DSSP 2ndary structure in alignments (default=show)  \n");
	printf(
			" -ssconf        show confidences for predicted 2ndary structure in alignments\n");
	printf(
			" -p <float>     minimum probability in summary and alignment list (def=%G)   \n",
			par.p);
	printf(
			" -E <float>     maximum E-value in summary and alignment list (def=%G)       \n",
			par.E);
	printf(
			" -Z <int>       maximum number of lines in summary hit list (def=%i)         \n",
			par.Z);
	printf(
			" -z <int>       minimum number of lines in summary hit list (def=%i)         \n",
			par.z);
	printf(
			" -B <int>       maximum number of alignments in alignment list (def=%i)      \n",
			par.B);
	printf(
			" -b <int>       minimum number of alignments in alignment list (def=%i)      \n",
			par.b);
	if (all) {
		printf(
				" -aliw [40,..[  number of columns per line in alignment list (def=%i)\n",
				par.aliwidth);
		printf(
				" -dbstrlen      max length of database string to be printed in hhr file\n");
	}
	printf("\n");
	printf(
			"Filter query multiple sequence alignment                                     \n");
	printf(
			" -id   [0,100]  maximum pairwise sequence identity (%%) (def=%i)   \n",
			par.max_seqid);
	printf(
			" -diff [0,inf[  filter MSA by selecting most diverse set of sequences, keeping \n");
	printf(
			"                at least this many seqs in each MSA block of length 50 (def=%i) \n",
			par.Ndiff);
	printf(" -cov  [0,100]  minimum coverage with query (%%) (def=%i) \n",
			par.coverage);
	printf(
			" -qid  [0,100]  minimum sequence identity with query (%%) (def=%i) \n",
			par.qid);
	printf(" -qsc  [0,100]  minimum score per column with query  (def=%.1f)\n",
			par.qsc);
	printf(" -neff [1,inf]  target diversity of alignment (default=off)\n");
	printf("\n");
	printf(
			"Input alignment format:                                                       \n");
	printf(
			" -M a2m         use A2M/A3M (default): upper case = Match; lower case = Insert;\n");
	printf(
			"                '-' = Delete; '.' = gaps aligned to inserts (may be omitted)   \n");
	printf(
			" -M first       use FASTA: columns with residue in 1st sequence are match states\n");
	printf(
			" -M [0,100]     use FASTA: columns with fewer than X%% gaps are match states   \n");
	if (all) {
		printf(
				" -tags          do NOT neutralize His-, C-myc-, FLAG-tags, and trypsin \n");
		printf(
				"                recognition sequence to background distribution    \n");
	}
	printf("\n");
	printf(
			"HMM-HMM alignment options:                                                    \n");
	printf(
			" -norealign     do NOT realign displayed hits with MAC algorithm (def=realign)   \n");
	printf(
			" -mact [0,1[    posterior prob threshold for MAC realignment controlling greedi- \n");
	printf(
			"                ness at alignment ends: 0:global  >0.1:local (default=%.2f)       \n",
			par.mact);
	printf(
			" -macins [0,1[  controls the cost of internal gap positions in the MAC algorithm.\n");
	printf(
			"                0:dense alignments  1:gappy alignments (default=%.2f)\n",
			par.macins);

	printf(
			" -glob/-loc     use global/local alignment mode for searching/ranking (def=local)\n");
//   printf(" -vit          use Viterbi algorithm for searching/ranking (default)          \n");
//   printf(" -mac          use Maximum Accuracy MAC algorithm for searching/ranking\n");
//   printf(" -forward      use Forward probability for searching                       \n");
	printf(
			" -alt <int>     show up to this many significant alternative alignments(def=%i)\n",
			par.altali);
	if (all) {
		printf(
				" -vit           use Viterbi algorithm for searching/ranking (default)       \n");
		printf(
				" -mac           use Maximum Accuracy (MAC) algorithm for searching/ranking\n");
		printf(
				" -forward       use Forward probability for searching                       \n");
		printf(
				" -excl <range>  exclude query positions from the alignment, e.g. '1-33,97-168' \n");
		printf(
				" -shift [-1,1]  score offset (def=%-.2f)                                       \n",
				par.shift);
		printf(
				" -corr [0,1]    weight of term for pair correlations (def=%.2f)                \n",
				par.corr);
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
				"        5       local amino acid composition correction                     \n");
	}
	printf(
			" -ssm {0,..,4}  0:   no ss scoring                                             \n");
	printf(
			"                1,2: ss scoring after or during alignment  [default=%1i]       \n",
			par.ssm);
	printf(
			"                3,4: ss scoring after or during alignment, predicted vs. predicted \n");
	if (all) {
		printf(
				" -ssw  [0,1]    weight of ss score compared to column score (def=%-.2f)     \n",
				par.ssw);
		printf(
				" -ssa  [0,1]    SS substitution matrix = (1-ssa)*I + ssa*full-SS-substition-matrix [def=%-.2f)\n",
				par.ssa);
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
//  printf("               3: CSBlast admixture:   tau = a(1+b)/(Neff[i]+b)                 \n");
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
	}

	printf("\n");
	printf(
			"Context-specific pseudo-counts:                                                  \n");
	printf(
			" -nocontxt      use substitution-matrix instead of context-specific pseudocounts \n");
	printf(
			" -contxt <file> context file for computing context-specific pseudocounts (default=%s)\n",
			par.clusterfile);

	if (all) {
		printf(
				" -csw  [0,inf]  weight of central position in cs pseudocount mode (def=%.1f)\n",
				par.csw);
		printf(
				" -csb  [0,1]    weight decay parameter for positions in cs pc mode (def=%.1f)\n",
				par.csb);
	}
	printf("\n");
	printf("Other options: \n");
	printf(
			" -cpu <int>     number of CPUs to use (for shared memory SMPs) (default=1)\n");
	printf(
			" -v <int>       verbose mode: 0:no screen output  1:only warings  2: verbose   \n");
	if (all) {
		printf(
				" -maxres <int>  max number of HMM columns (def=%5i)             \n",
				par.maxres);
		printf(
				" -maxmem [1,inf[ limit memory for realignment (in GB) (def=%.1f)          \n",
				par.maxmem);
		printf(
				" -scores <file> write scores for all pairwise comparisions to file         \n");
		printf(
				" -calm {0,..,3} empirical score calibration of 0:query 1:template 2:both   \n");
		printf(
				"                default 3: neural network-based estimation of EVD params   \n");
		// printf(" -opt  <file>   parameter optimization mode (def=off): return sum of ranks \n");
		// printf("                of true positives (same superfamily) for minimization      \n");
		// printf("                and write result into file                                 \n");
		printf("\n");
	} else {
		printf(
				"An extended list of options can be obtained by calling 'hhblits -help'\n");
	}
	printf("\n");
	printf("Example: hhsearch -i a.1.1.1.a3m -d scop70_1.71.hhm \n");
	cout << endl;

//   printf("More help:                                                         \n");
//   printf(" -h out        options for formatting ouput                        \n");
//   printf(" -h hmm        options for building HMM from multiple alignment    \n");
//   printf(" -h gap        options for setting gap penalties                   \n");
//   printf(" -h ali        options for HMM-HMM alignment                       \n");
//   printf(" -h other      other options \n");
//   printf(" -h all        all options \n");
}

/////////////////////////////////////////////////////////////////////////////////////
//// Processing input options from command line
/////////////////////////////////////////////////////////////////////////////////////
void HHsearch::ProcessArguments(int argc, char** argv, Parameters& par) {
	char program_name[NAMELEN];
	RemovePathAndExtension(program_name, par.argv[0]);

	//Processing command line input
	for (int i = 1; i < argc; i++) {
		HH_LOG(LogLevel::DEBUG1) << i << "  " << argv[i] << endl;
		if (!strcmp(argv[i], "-i")) {
			if (++i >= argc || argv[i][0] == '-') {
				help(par);
				cerr << endl << "Error in " << program_name
						<< ": no query file following -i\n";
				exit(4);
			} else
				strcpy(par.infile, argv[i]);
		} else if (!strcmp(argv[i], "-d")) {
			if (++i >= argc || argv[i][0] == '-') {
				help(par);
				cerr << endl << "Error in " << program_name
						<< ": no database file following -d\n";
				exit(4);
			} else {
				std::string db(argv[i]);
				par.db_bases.push_back(db);
			}
		} else if (!strcmp(argv[i], "-o")) {
			par.append = 0;
			if (++i >= argc || argv[i][0] == '-') {
				help(par);
				cerr << endl << "Error in " << program_name
						<< ": no output file following -o\n";
				exit(4);
			} else
				strcpy(par.outfile, argv[i]);
		} else if (!strcmp(argv[i], "-Ofas")) {
			par.append = 0;
			par.outformat = 1;
			if (++i >= argc || argv[i][0] == '-') {
				help(par);
				cerr << endl << "Error in " << program_name
						<< ": no output file following -o\n";
				exit(4);
			} else
				strcpy(par.pairwisealisfile, argv[i]);
		} else if (!strcmp(argv[i], "-Oa2m")) {
			par.append = 0;
			par.outformat = 2;
			if (++i >= argc || argv[i][0] == '-') {
				help(par);
				cerr << endl << "Error in " << program_name
						<< ": no output file following -o\n";
				exit(4);
			} else
				strcpy(par.pairwisealisfile, argv[i]);
		} else if (!strcmp(argv[i], "-Oa3m")) {
			par.append = 0;
			par.outformat = 3;
			if (++i >= argc || argv[i][0] == '-') {
				help(par);
				cerr << endl << "Error in " << program_name
						<< ": no output file following -o\n";
				exit(4);
			} else
				strcpy(par.pairwisealisfile, argv[i]);
		} else if (!strcmp(argv[i], "-oa3m")) {
			if (++i >= argc || argv[i][0] == '-') {
				help(par);
				cerr << endl << "Error in " << program_name
						<< ": no output file following -Oa3m\n";
				exit(4);
			} else
				strcpy(par.alnfile, argv[i]);
		} else if (!strcmp(argv[i], "-ohhm")) {
			if (++i >= argc || argv[i][0] == '-') {
				help(par);
				cerr << endl << "Error in " << program_name
						<< ": no output file following -ohhm\n";
				exit(4);
			} else
				strcpy(par.hhmfile, argv[i]);
		} else if (!strcmp(argv[i], "-opsi")) {
			if (++i >= argc || argv[i][0] == '-') {
				help(par);
				cerr << endl << "Error in " << program_name
						<< ": no output file following -Opsi\n";
				exit(4);
			} else
				strcpy(par.psifile, argv[i]);
		} else if (!strcmp(argv[i], "-scores")) {
			if (++i >= argc || argv[i][0] == '-') {
				help(par);
				cerr << endl << "Error in " << program_name
						<< ": no file following -scores\n";
				exit(4);
			} else {
				strcpy(par.scorefile, argv[i]);
			}
		} else if (!strcmp(argv[i], "-atab") || !strcmp(argv[i], "-Aliout")) {
			if (++i >= argc || argv[i][0] == '-') {
				help(par);
				cerr << endl << "Error in " << program_name
						<< ": no query file following -atab\n";
				exit(4);
			} else
				strmcpy(par.alitabfile, argv[i], NAMELEN - 1);
		} else if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "-help")) {
			help(par, 1);
			exit(0);
		} else if (!strcmp(argv[i], "-v")) {
			if (++i >= argc) {
				help(par);
				exit(4);
			}
			int v = atoi(argv[++i]);
			par.v = Log::from_int(v);
			Log::reporting_level() = par.v;
		} else if (!strcmp(argv[i], "-p") && (i < argc - 1))
			par.p = atof(argv[++i]);
		else if (!strcmp(argv[i], "-P") && (i < argc - 1))
			par.p = atof(argv[++i]);
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
		else if (!strcmp(argv[i], "-e") && (i < argc - 1))
			par.e = atof(argv[++i]);
		else if (!strncmp(argv[i], "-nocons", 7))
			par.showcons = 0;
		else if (!strncmp(argv[i], "-nopred", 7))
			par.showpred = 0;
		else if (!strncmp(argv[i], "-nodssp", 7))
			par.showdssp = 0;
		else if (!strncmp(argv[i], "-ssconf", 7))
			par.showconf = 1;
		else if (!strncmp(argv[i], "-cons", 5))
			par.cons = 1;
		else if (!strncmp(argv[i], "-mark", 5))
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
		else if (!strcmp(argv[i], "-neff") && (i < argc - 1))
			par.Neff = atof(argv[++i]);
		else if (!strcmp(argv[i], "-Neff") && (i < argc - 1))
			par.Neff = atof(argv[++i]);
		else if (!strcmp(argv[i], "-Gonnet"))
			par.matrix = 0;
		else if (!strncmp(argv[i], "-BLOSUM", 7)
				|| !strncmp(argv[i], "-Blosum", 7)) {
			if (!strcmp(argv[i] + 7, "30"))
				par.matrix = 30;
			else if (!strcmp(argv[i] + 7, "40"))
				par.matrix = 40;
			else if (!strcmp(argv[i] + 7, "50"))
				par.matrix = 50;
			else if (!strcmp(argv[i] + 7, "65"))
				par.matrix = 65;
			else if (!strcmp(argv[i] + 7, "80"))
				par.matrix = 80;
			else
				cerr << endl << "WARNING: Ignoring unknown option " << argv[i]
						<< " ...\n";
		} else if (!strcmp(argv[i], "-wg")) {
			par.wg = 1;
		} else if (!strcmp(argv[i], "-pcm") && (i < argc - 1))
			par.pc_hhm_context_engine.admix = (Pseudocounts::Admix) atoi(
					argv[++i]);
		else if (!strcmp(argv[i], "-pca") && (i < argc - 1))
			par.pc_hhm_context_engine.pca = atof(argv[++i]);
		else if (!strcmp(argv[i], "-pcb") && (i < argc - 1))
			par.pc_hhm_context_engine.pcb = atof(argv[++i]);
		else if (!strcmp(argv[i], "-pcc") && (i < argc - 1))
			par.pc_hhm_context_engine.pcc = atof(argv[++i]);
		//else if (!strcmp(argv[i],"-pcw") && (i<argc-1)) par.pcw=atof(argv[++i]);
		else if (!strcmp(argv[i], "-gapb") && (i < argc - 1)) {
			par.gapb = atof(argv[++i]);
			if (par.gapb <= 0.01)
				par.gapb = 0.01;
		} else if (!strcmp(argv[i], "-gapd") && (i < argc - 1))
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
		else if (!strcmp(argv[i], "-ssw_mac") && (i < argc - 1))
			par.ssw_realign = atof(argv[++i]);
		else if (!strcmp(argv[i], "-ssa") && (i < argc - 1))
			par.ssa = atof(argv[++i]);
		else if (!strcmp(argv[i], "-realign"))
			par.realign = 1;
		else if (!strcmp(argv[i], "-norealign"))
			par.realign = 0;
		else if (!strcmp(argv[i], "-mac") || !strcmp(argv[i], "-MAC"))
			par.forward = 2;
		else if (!strcmp(argv[i], "-map") || !strcmp(argv[i], "-MAP"))
			par.forward = 2;
		else if (!strcmp(argv[i], "-vit"))
			par.forward = 0;
		else if (!strncmp(argv[i], "-glo", 3)) {
			par.loc = 0;
			if (par.mact > 0.35 && par.mact < 0.3502) {
				par.mact = 0;
			}
		} else if (!strncmp(argv[i], "-loc", 4))
			par.loc = 1;
		else if (!strncmp(argv[i], "-alt", 4) && (i < argc - 1))
			par.altali = atoi(argv[++i]);
		else if (!strcmp(argv[i], "-M") && (i < argc - 1))
			if (!strcmp(argv[++i], "a2m") || !strcmp(argv[i], "a3m"))
				par.M = 1;
			else if (!strcmp(argv[i], "first"))
				par.M = 3;
			else if (argv[i][0] >= '0' && argv[i][0] <= '9') {
				par.Mgaps = atoi(argv[i]);
				par.M = 2;
			} else
				cerr << endl << "WARNING: Ignoring unknown argument: -M "
						<< argv[i] << "\n";
		else if (!strcmp(argv[i], "-shift") && (i < argc - 1))
			par.shift = atof(argv[++i]);
		else if ((!strcmp(argv[i], "-mact") || !strcmp(argv[i], "-mapt"))
				&& (i < argc - 1))
			par.mact = atof(argv[++i]);
		else if (!strcmp(argv[i], "-macins") && (i < argc - 1))
			par.macins = atof(argv[++i]);
		else if (!strcmp(argv[i], "-sc") && (i < argc - 1))
			par.columnscore = atoi(argv[++i]);
		else if (!strcmp(argv[i], "-scwin") && (i < argc - 1)) {
			par.columnscore = 5;
			par.half_window_size_local_aa_bg_freqs = imax(1, atoi(argv[++i]));
		}
//      TODO: is this necessary?
//      else if (!strcmp(argv[i],"-def")) ;
		else if (!strcmp(argv[i], "-maxres") && (i < argc - 1)) {
			par.maxres = atoi(argv[++i]);
			par.maxcol = 2 * par.maxres;
		} else if (!strncmp(argv[i], "-cpu", 4) && (i < argc - 1))
			par.threads = atoi(argv[++i]);
		else if (!strcmp(argv[i], "-maxmem") && (i < argc - 1)) {
			par.maxmem = atof(argv[++i]);
		} else if (!strcmp(argv[i], "-corr") && (i < argc - 1))
			par.corr = atof(argv[++i]);
		else if (!strcmp(argv[i], "-ovlp") && (i < argc - 1))
			par.min_overlap = atoi(argv[++i]);
		else if (!strcmp(argv[i], "-dbstrlen") && (i < argc - 1))
			par.maxdbstrlen = atoi(argv[++i]);
		else if (!strcmp(argv[i], "-tags"))
			par.notags = 0;
		else if (!strcmp(argv[i], "-notags"))
			par.notags = 1;
		else if (!strncmp(argv[i], "-idummy", 7) && (i < argc - 1))
			par.idummy = atoi(argv[++i]);
		else if (!strncmp(argv[i], "-premerge", 9) && (i < argc - 1))
			par.premerge = atoi(argv[++i]);
		else if (!strncmp(argv[i], "-fdummy", 7) && (i < argc - 1))
			par.fdummy = atof(argv[++i]);
		else if (!strcmp(argv[i], "-nocontxt"))
			par.nocontxt = 1;
		else if (!strcmp(argv[i], "-csb") && (i < argc - 1))
			par.csb = atof(argv[++i]);
		else if (!strcmp(argv[i], "-csw") && (i < argc - 1))
			par.csw = atof(argv[++i]);
		else if (!strcmp(argv[i], "-contxt") || !strcmp(argv[i], "-cs")) {
			if (++i >= argc || argv[i][0] == '-') {
				help(par);
				cerr << endl << "Error in " << program_name
						<< ": no query file following -contxt\n";
				exit(4);
			} else
				strcpy(par.clusterfile, argv[i]);
		} else {
			HH_LOG(LogLevel::WARNING) << endl
				<< "WARNING: Ignoring unknown option " << argv[i] << " ...\n";
		}
		HH_LOG(LogLevel::DEBUG1) << i << "  " << argv[i] << endl;
	} // end of for-loop for command line input
}

void HHsearch::run(FILE* query_fh, char* query_path) {
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
	ReadQueryFile(par, query_fh, input_format, par.wg, q, Qali, query_path, pb,
			S, Sim);
	PrepareQueryHMM(par, input_format, q, pc_hhm_context_engine,
			pc_hhm_context_mode, pb, R);
    q_vec.MapOneHMM(q);
    *q_tmp = *q;

	// Reset lamda?
	//TODO
	if (par.calibrate > 0) {
		q->lamda = LAMDA;
		q->mu = 0.0;
	}

	// Set query columns in His-tags etc to Null model distribution
	if (par.notags)
		q->NeutralizeTags(pb);

	// Search databases

	std::vector<HHEntry*> new_entries;
	if (!par.prefilter) {
		for (size_t i = 0; i < dbs.size(); i++) {
			dbs[i]->initNoPrefilter(new_entries);
		}
	}

	int max_template_length = getMaxTemplateLength(new_entries);
    if(max_template_length > par.maxres){
		HH_LOG(LogLevel::WARNING) << "database contains sequnces that exceeds maximum allowed size (maxres = "
        << par.maxres << "). Maxres can be increased with parameter -maxres." <<std::endl;
    }
    max_template_length = std::min(max_template_length, par.maxres);
	for(int i = 0; i < par.threads; i++) {
	  viterbiMatrices[i]->AllocateBacktraceMatrix(q->L, max_template_length);
	}

    ViterbiRunner viterbirunner(viterbiMatrices, dbs, par.threads);
    std::vector<Hit> hits_to_add = viterbirunner.alignment(par, &q_vec, new_entries, par.qsc_db, pb, S, Sim, R);

    hitlist.N_searched = new_entries.size();
    add_hits_to_hitlist(hits_to_add, hitlist);

//TODO
//  if (v1 >= 2)
//    cout << "\n";
//  v = v1;
//
//  // Sort list according to sortscore
//  if (v >= 3)
//    printf("Sorting hit list ...\n");
//  hitlist.SortList();
//
//  // Fit EVD (with lamda, mu) to score distribution?
//  if (par.calm == 3) {
//    hitlist.CalculatePvalues(q, par.loc, par.ssm, par.ssw); // Use NN prediction of lamda and mu
//  }
//  else if ((par.calm != 1 && q->lamda == 0) || par.calibrate > 0) {
//    if (v >= 2 && par.loc)
//      printf("Fitting scores with EVD (first round) ...\n");
//    hitlist.MaxLikelihoodEVD(q, 3, par.loc, par.ssm, par.ssw); // first ML fit: exclude 3 best superfamilies from fit
//
//    if (v >= 3)
//      printf("Number of families present in database: %i\n", hitlist.fams); // DEBUG
//    if (hitlist.fams >= 100) {
//      if (par.loc) {
//        if (v >= 2)
//          printf("Fitting scores with EVD (second round) ...\n");
//        hitlist.MaxLikelihoodEVD(q, 0, par.loc, par.ssm, par.ssw); // second ML fit: exclude superfamilies with E-value<MINEVALEXCL
//        hitlist.ResortList();
//      }
//      else {
//        if (v >= 2)
//          fprintf(stderr,
//              "WARNING: E-values for global alignment option may be unreliable.\n");
//        hitlist.ResortList();
//      }
//    }
//    else {
//      if (v) {
//        fprintf(stderr, "\nWARNING: no E-values could be calculated.\n");
//        fprintf(stderr, "To calculate E-values you have two options:\n");
//        fprintf(stderr,
//            "1. Calibrate your query profile HMM with a calibration database:\n");
//        fprintf(stderr, "   > hhsearch -i yourHMM.hhm -d cal.hhm -cal\n");
//        fprintf(stderr,
//            "   This will insert a line in yourHMM.hhm with lamda and mu of the score distribution.\n");
//        fprintf(stderr,
//            "   cal.hhm contains 1220 HMMs from different SCOP superfamilies and is supplied with HHsearch.\n");
//        fprintf(stderr,
//            "   Instead of cal.hhm you may also use any SCOP database file, e.g. scop70_1.69\n");
//        fprintf(stderr,
//            "   Note that your HMM needs to be recalibrated when changing HMM-HMM alignment options.\n");
//        fprintf(stderr, "2. Append cal.hhm to your own database:\n");
//        fprintf(stderr, "   > cat cal.hhm >> yourDB.hhm\n");
//        fprintf(stderr,
//            "   But note that HMMs contained in cal.hmm will pop up among your hits.\n");
//      }
//    }
//    if (par.calm == 2)
//      hitlist.GetPvalsFromCalibration(q, par.loc, par.calm, par.ssm, par.ssw);
//  }
//  else
//    hitlist.GetPvalsFromCalibration(q, par.loc, par.calm, par.ssm, par.ssw);
//
//TODO
//  // Optimization mode?
//  if (par.opt)
//    hitlist.Optimize(q);

	// Set new ss weight for realign
	par.ssw = par.ssw_realign;

	// Realign hits with MAC algorithm
	if (par.realign && par.forward != 2) {
		perform_realign(q_vec, new_entries, premerge, premerged_hits);
	}

	// Write HMM to output file without pseudocounts
	//TODO
	if (par.calibrate)
		q->InsertCalibration(par.infile);

	//TODO: does no longer search for a3m, but takes a3m from hmm if needs be
	mergeHitsToQuery(previous_hits, premerged_hits, seqs_found, cluster_found);

	// Calculate pos-specific weights, AA frequencies and transitions -> f[i][a], tr[i][a]
	Qali.FrequenciesAndTransitions(q, par.wg, par.mark, par.cons, par.showcons,
			par.maxres, pb, Sim, NULL, true);

	if (par.notags)
		q->NeutralizeTags(pb);

	for(size_t i = 0; i < new_entries.size(); i++) {
	  delete new_entries[i];
	}
	new_entries.clear();
}
