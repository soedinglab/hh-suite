#include "hhdecl.h"

/////////////////////////////////////////////////////////////////////////////////////
//// Global variable declarations
/////////////////////////////////////////////////////////////////////////////////////

void Parameters::SetDefaultPaths() {
    char hhlib[PATH_MAX];   // lib base path e.g. /usr/lib64/hh
    char hhdata[PATH_MAX];  // data base path e.g. /usr/lib64/hh/data

	// set hhlib
	FILE* testf = NULL;
	if (getenv("HHLIB"))
		strcpy(hhlib, getenv("HHLIB"));
	else
		strcpy(hhlib, "/usr/lib");

	strcat(strcpy(hhdata, hhlib), "/data");
	strcat(strcpy(clusterfile, hhdata), "/context_data.crf");
	strcat(strcpy(cs_library, hhdata), "/cs219.lib");

	testf = fopen(cs_library, "r");
	if (testf)
		fclose(testf);
	else {
		HH_LOG(DEBUG) << "WARNING in HHsuite: Could not open " << cs_library << "\n";

		char program_path[NAMELEN];
		Pathname(program_path, argv[0]);

		/* we did not find HHLIB, if called with full path or in dist dir, we can try relative to program path */
		if (program_path != NULL) {
			strcat(strcpy(hhlib, program_path), "../lib/hh");
			strcat(strcpy(hhdata, hhlib), "/data");
			strcat(strcpy(clusterfile, hhdata), "/context_data.crf");
			strcat(strcpy(cs_library, hhdata), "/cs219.lib");
			testf = fopen(cs_library, "r");
			if (testf)
				fclose(testf);
			else {
				HH_LOG(DEBUG) << "WARNING in HHsuite: Could not open " << cs_library << "\n";

				strcat(strcpy(hhlib, program_path), "..");
				strcat(strcpy(hhdata, hhlib), "/data");
				strcat(strcpy(clusterfile, hhdata), "/context_data.crf");
				strcat(strcpy(cs_library, hhdata), "/cs219.lib");
				testf = fopen(cs_library, "r");
				if (testf)
					fclose(testf);
				else
					HH_LOG(DEBUG) << "WARNING in HHsuite: Could not open " << cs_library << "\n";
			}
		}
	}
	if (!testf) {
	  HH_LOG(ERROR) << "Could not find context_data.crf and cs219.lib in '" << hhlib << std::endl;
	  HH_LOG(ERROR) << "Please set the HHLIB environment variable to the HH-suite directory" << std::endl;
	  HH_LOG(ERROR) << "(Linux bash: export HHLIB=<hh_dir>, csh/tcsh: setenv HHLIB=<hh_dir>)." << std::endl;
	  HH_LOG(ERROR) << "The missing files should be in $HHLIB/data/." << std::endl;
		exit(2);
	}
	return;
}

void Parameters::SetDefaults() {
	v = INFO;

	maxcol = 32765; // max number of columns in sequence/MSA input files; must be <= LINELEN and >= maxres
	maxres = 20001;           // max number of states in HMM; must be <= LINELEN
	maxnumdb = 20000;          // max number of hits allowed past prefilter

	append = 0;                // overwrite output file
	outformat = 0;             // 0: hhr  1: FASTA  2:A2M   3:A3M
	p = 20.0f; // minimum threshold for inclusion in hit list and alignment listing
	E = 1e6f; // maximum threshold for inclusion in hit list and alignment listing
	b = 10;                    // min number of alignments
	B = 500;                   // max number of alignments
	z = 10;                    // min number of lines in hit list
	Z = 500;                   // max number of lines in hit list
	e = 1e-3f; // maximum E-value for inclusion in output alignment, output HMM, and PSI-BLAST checkpoint model
	realign_max = 500;        // Maximum number of HMM hits to realign
	maxmem = 3.0;            // 3GB
	showcons = 1;              // show consensus sequence
	showdssp = 1;              // show predicted secondary structure ss_dssp
	showpred = 1;              // show predicted secondary structure ss_pred
	showconf = 0;           // don't show secondary structure confidence ss_conf
	cons = 0; // chose first non-SS sequence as main representative sequence (not consensus)
	nseqdis = 1;       // maximum number of query sequences for output alignment
	mark = 0; // 1: only marked sequences (or first) get displayed; 0: most divergent ones get displayed
	aliwidth = 80; // number of characters per line in output alignments for HMM search

	max_seqid = 90;           // default for maximum sequence identity threshold
	qid = 0;                 // default for minimum sequence identity with query
	qsc = -20.0f;             // default for minimum score per column with query
	coverage = 0;              // default for minimum coverage threshold
	Ndiff = 100;           // pick Ndiff most different sequences from alignment
	allseqs = false;    // if true, do not filter result MSA; show all sequences

	Neff = 0; // Filter alignment to a diversity (Neff) with a maximum Neff of par.Neff

	M = 1;                     // match state assignment is by A2M/A3M
	Mgaps = 50; // Above this percentage of gaps, columns are assigned to insert states (for par.M=2)

	wg = 0;                // 0: use local sequence weights   1: use global ones

	matrix = 0;      // Subst.matrix 0: Gonnet, 1: HSDM, 2: BLOSUM50 3: BLOSUM62

	pc_hhm_context_engine.admix = Pseudocounts::HHsearchAdmix;
	pc_hhm_context_engine.pca = 0.90;
	pc_hhm_context_engine.pcb = 4.00;
	pc_hhm_context_engine.pcc = 1.0;
	pc_hhm_context_engine.target_neff = 0.0;

	pc_prefilter_context_engine.admix = Pseudocounts::CSBlastAdmix;
	pc_prefilter_context_engine.pca = 0.80;
	pc_prefilter_context_engine.pcb = 2.00;
	pc_prefilter_context_engine.pcc = 1.0;
	pc_prefilter_context_engine.target_neff = 0.0;

	pc_hhm_nocontext_mode = 2;
	pc_hhm_nocontext_a = 1.0f;
	pc_hhm_nocontext_b = 1.5f;
	pc_hhm_nocontext_c = 1.0f;

	pc_prefilter_nocontext_mode = 2;
	pc_prefilter_nocontext_a = 1.0f;
	pc_prefilter_nocontext_b = 1.5f;
	pc_prefilter_nocontext_c = 1.0f;

	gapb = 1.0;                // default values for transition pseudocounts
	gapd = 0.15; // gap open penalty pseudocount; 0.25 corresponds to 7.1*gapf bits
	gape = 1.0;                // gap extension penalty pseudocount
	gapf = 0.6;            // factor for increasing gap open penalty for deletes
	gapg = 0.6;            // factor for increasing gap open penalty for inserts
	gaph = 0.6;       // factor for increasing gap extension penalty for deletes
	gapi = 0.6;       // factor for increasing gap extension penalty for inserts

	ssm = 2; // ss scoring mode: 0:no ss score  1:score after alignment  2:score during alignment
	ssw = 0.11f;               // weight of ss scoring
	ssw_realign = 0.11f;       // weight of ss scoring for realign
	ssa = 1.0f;                // weight of ss evolution matrix
	shift = -0.03f;            // Shift match score up
	mact = 0.3501f; // Probability threshold for MAC alignment in local mode for alignment ends (set to 0.3501 to track user modification)
	corr = 0.1f;               // Weight of correlations of scores for |i-j|<=4

	egq = 0.0f;                // no charge for end gaps as default
	egt = 0.0f;                // no charge for end gaps as default

	loc = 1;                   // local vs. global alignment as default
	altali = 4;                // find up to four (possibly overlapping) subalignments // JS:02 Mar 13: changed from 2 to avoid loosing domain predictions of repeated modules
	smin = 20.;                //Minimum score of hit needed to search for another repeat of same profile: p=exp(-(4-mu)/lamda)=0.01
	realign = 1;               // realign with MAC algorithm

	columnscore = 1; // Default column score is 1: null model pnul = 1/2 * (q_av(a)+p_av(a))
	half_window_size_local_aa_bg_freqs = 40;
	min_overlap = 0;           // automatic minimum overlap used
	readdefaultsfile = 0; // Default = do not read a defaults file ./.hhdefaults or HOME/.hhdefaults
	maxdbstrlen = 200; // maximum length of database string to be printed in 'Command' line of hhr file

	notags = 1;                // neutralize His-tags, FLAG-tags, C-myc-tags
	hmmer_used = false;

	// HHblits parameters
	dbsize = 0;

	// HHblits Evalue calculation  (alpha = a + b(Neff(T) - 1)(1 - c(Neff(Q) - 1)) )
	alphaa = 0.4;
	alphab = 0.02;
	alphac = 0.1;

	prefilter = false;              //true in hhblits

	early_stopping_filter = false;  //true in hhblits
  filter_thresh=0;                // 0.01 in hhblits

	// For HHblits prefiltering with SSE2
	prefilter_gap_open = 20;
	prefilter_gap_extend = 4;
	prefilter_score_offset = 50;
	prefilter_bit_factor = 4;
	prefilter_evalue_thresh = 1000;
	prefilter_evalue_coarse_thresh = 100000;
	preprefilter_smax_thresh = 10;
	min_prefilter_hits = 100;

	// For filtering database alignments in HHsearch and HHblits
	//JS: What are these used for? They are set to the options without _db anyway.
	max_seqid_db = max_seqid;
	qid_db = qid;
	qsc_db = qsc;
	coverage_db = coverage;
	Ndiff_db = Ndiff;

	// Initialize strings

	strcpy(infile, ""); // was reverted back from 'strcpy(infile,"stdin");' (to show help list when no options are given)
	strcpy(outfile, "");
	strcpy(matrices_output_file, "");
	strcpy(pairwisealisfile, "");
	strcpy(scorefile, "");
	strcpy(indexfile, "");
	strcpy(alnfile, "");
	strcpy(hhmfile, "");
	strcpy(psifile, "");
	strcpy(alitabfile, "");
	strcpy(alisbasename, "");
	exclstr = NULL;

	max_number_matrices = 100;

	// parameters for context-specific pseudocounts
	csb = 0.85;
	csw = 1.6;

	//hhblits specific variables
	num_rounds = 2;
	already_seen_filter = true;   // Perform filtering of already seen HHMs
	realign_old_hits = false; // Realign old hits in last round or use previous alignments
	neffmax = 10.0;
	threads = 2;
	nocontxt = false;
}

Parameters::Parameters() {
	SetDefaults();
}


