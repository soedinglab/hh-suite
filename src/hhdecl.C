#ifndef MAIN
#define EXTERN extern
#else
#define EXTERN
#endif


/////////////////////////////////////////////////////////////////////////////////////
//// Constants
/////////////////////////////////////////////////////////////////////////////////////

EXTERN const char VERSION_AND_DATE[]="version 2.0.15 (June 2012)";
EXTERN const char REFERENCE[]="Remmert M, Biegert A, Hauser A, and Soding J.\nHHblits: Lightning-fast iterative protein sequence searching by HMM-HMM alignment.\nNat. Methods 9:173-175 (2011).\n";
EXTERN const char COPYRIGHT[]="(C) Johannes Soeding, Michael Remmert, Andreas Biegert, Andreas Hauser\n";
EXTERN const int MAXSEQ=65535; //max number of sequences in input alignment (must be <~30000 on cluster nodes??)
EXTERN const int LINELEN=524288; //max length of line read in from input files; must be >= MAXCOL
EXTERN const int MAXSEQDIS=10238;//max number of sequences stored in 'hit' objects and displayed in output alignment
EXTERN const int IDLEN=255;     //max length of scop hierarchy id and pdb-id
EXTERN const int DESCLEN=32765;//max length of sequence description (longname)
EXTERN const int NAMELEN=(PATH_MAX>512? PATH_MAX:512); //max length of file names etc., defined in limits.h 
EXTERN const int MAXOPT=127;   //Maximum number of options to be read in from .hhconfig or command line
EXTERN const int NAA=20;       //number of amino acids (0-19)
EXTERN const int NTRANS=7;     //number of transitions recorded in HMM (M2M,M2I,M2D,I2M,I2I,D2M,D2D)
EXTERN const int NCOLMIN=10;   //min number of cols in subalignment for calculating pos-specific weights w[k][i]
EXTERN const int ANY=20;       //number representing an X (any amino acid) internally
EXTERN const int GAP=21;       //number representing a gap internally
EXTERN const int ENDGAP=22;    //Important to distinguish because end gaps do not contribute to tansition counts
EXTERN const int HMMSCALE=1000;//Scaling number for log2-values in HMMs
EXTERN const int NFAMMAX=5119; //Size of hash for counting number of HMMs in each family
EXTERN const int MAXPROF=32766;//Maximum number of HMM scores for fitting EVD
EXTERN const float MAXENDGAPFRAC=0.1; //For weighting: include only columns into subalignment i that have a max fraction of seqs with endgap
EXTERN const float SMIN= 20.;  //Minimum score of hit needed to search for another repeat of same profile: p=exp(-(4-mu)/lamda)=0.01
EXTERN const float LAMDA=0.388; //lamda in score EVD used for -local mode in length correction: S = S-log(Lq*Lt)/LAMDA)
EXTERN const float LAMDA_GLOB=0.42; //lamda in score EVD used for -global mode
EXTERN const float PMAX=1E-2;  //Maximum single-repeat p-value that can contribute to whole-protein p-value
EXTERN const float MINEVALEXCL=0.5; //above this E-value from first ML fit hits are not used for final ML fit of EVD
EXTERN const int SELFEXCL=3;   // exclude self-alignments with j-i<SELFEXCL
EXTERN const float PLTY_GAPOPEN=6.0f; // for -qsc option (filter for min similarity to query): 6 bits to open gap
EXTERN const float PLTY_GAPEXTD=1.0f; // for -qsc option (filter for min similarity to query): 1 bit to extend gap
EXTERN const int MINCOLS_REALIGN=6; // hits with MAC alignments with fewer matched columns will be deleted in hhsearch hitlist; must be at least 2 to avoid nonsense MAC alignments starting from the left/upper edge
EXTERN const float LOG1000=log(1000.0);

// Secondary structure
EXTERN const int NDSSP=8;      //number of different ss states determined by dssp: 0-7 (0: no state available)
EXTERN const int NSSPRED=4;    //number of different ss states predicted by psipred: 0-3 (0: no prediction availabe)
EXTERN const int MAXCF=11;     //number of different confidence values: 0-10 (0: no prediction availabe)
EXTERN const int NSA=7;        //number of classes relative solvent accesiblity (0:no coord,  1:<2%, 2:<14%, 3:<33%, 4:<55%, 5:>55%, 6:S-S bridge)

// const char aa[]="ARNDCQEGHILKMFPSTWYVX-";
//Amino acids Sorted by alphabet     -> internal numbers a
//                0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
//                A  C  D  E  F  G  H  I  K  L  M  N  P  Q  R  S  T  V  W  Y  X
EXTERN const int s2a[]={ 0, 4, 3, 6,13, 7, 8, 9,11,10,12, 2,14, 5, 1,15,16,19,17,18,20};
//Internal numbers a for amino acids -> amino acids Sorted by alphabet:
//                0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
//                A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  X
EXTERN const int a2s[]={ 0,14,11, 2, 1,13, 3, 5, 6, 7, 9, 8,10, 4,12,15,16,18,19,17,20};

enum transitions {M2M,M2I,M2D,I2M,I2I,D2M,D2D}; // index for transitions within a HMM
enum pair_states {STOP=0,SAME=1,GD=2,IM=3,DG=4,MI=5,MS=6,ML=7,SM=8,LM=9,MM=10};



/////////////////////////////////////////////////////////////////////////////////////
//// Global variable declarations
/////////////////////////////////////////////////////////////////////////////////////

EXTERN char v;                //=2 1: show only warnings 2:verbose mode
EXTERN char program_name[NAMELEN]; //name of program executed (e.g. hhmake of hhsearch)
EXTERN char program_path[NAMELEN]; //path of program executed

// substitution matrix flavours
EXTERN float __attribute__((aligned(16))) P[20][20]; // P[a][b] = combined probability for a aligned to b
EXTERN float __attribute__((aligned(16))) R[20][20]; // R[a][b]=P[a][b]/p[b]=P(a|b); precalculated for pseudocounts
EXTERN float __attribute__((aligned(16))) Sim[20][20]; // Similarity matrix Sim[a][b]: how similar are a and b?
EXTERN float __attribute__((aligned(16))) S[20][20]; // Substitution score matrix S[a][b] = log2(Pab/pa/pb)
EXTERN float __attribute__((aligned(16))) pb[21];    // pb[a] = background amino acid probabilities for chosen substitution matrix
EXTERN float __attribute__((aligned(16))) qav[21];   // qav[a] = background amino acid probabilities for query HMM (needed for rate matrix rescaling)

// secondary structure matrices
EXTERN float S73[NDSSP][NSSPRED][MAXCF];           // P[A][B][cf]       =  log2 P(A,B,cf)/P(A)/P(B,cf)
EXTERN float S33[NSSPRED][MAXCF][NSSPRED][MAXCF];  // P[B][cf][B'][cf'] =  log2 sum_B' P(A,B',cf)/P(A)/P(B,cf) * P_b(B'|B)
// float S77[NDSSP][DSSP];                  // P[A][B]           =  log2 P(A,B)/P(A)/P(B)


// Structure to store data for HHblits early stopping filter
EXTERN struct Early_Stopping {
  int length;       // Length of array of 1/evalues
  int counter;      // counter for evalue array
  double* evals;    // array of last 1/evalues
  double thresh;    // Threshold for early stopping
  double sum;       // sum of evalues in array
} *early_stopping=NULL;



/////////////////////////////////////////////////////////////////////////////////////
// Class declarations
/////////////////////////////////////////////////////////////////////////////////////

//container for the scores used for cs scoring
struct ColumnStateScoring {
  float** substitutionScores; //[i][j]; i: query; j:column states
  int query_length;
  int number_column_states;
};


// Pseudocounts
namespace Pseudocounts {
enum Admix {
  ConstantAdmix = 1,
  HHsearchAdmix = 2,
  CSBlastAdmix  = 3
};

struct Params {
  Params(
      Admix m     = ConstantAdmix,
      double a    = 1.0,
      double b    = 1.0,
      double c    = 1.0,
      double neff = 0.0)
    : admix(m), pca(a), pcb(b), pcc(c), target_neff(neff) {}

  cs::Admix* CreateAdmix() {
    switch (admix) {
      case ConstantAdmix: 
        return new cs::ConstantAdmix(pca);
        break;
      case HHsearchAdmix: 
        return new cs::HHsearchAdmix(pca, pcb, pcc);
        break;
      case CSBlastAdmix:
        return new cs::CSBlastAdmix(pca, pcb);
        break;
      default:
        return NULL;
    }
  }

  Admix admix;        // admixture mode
  double pca;         // admixture paramter a
  double pcb;         // admixture paramter b
  double pcc;         // admixture parameter c needed for HHsearchAdmix
  double target_neff; // target diversity adjusted by optimizing a
};

};  // Pseudocounts

// Input parameters
class Parameters          // Parameters for gap penalties and pseudocounts
{
public:
  char** argv;            //command line parameters
  char argc;              //dimension of argv

  char hhlib[PATH_MAX];   // lib base path e.g. /usr/lib64/hh
  char hhdata[PATH_MAX];  // data base path e.g. /usr/lib64/hh/data

  char infile[NAMELEN];   // input filename
  char outfile[NAMELEN];  // output filename
  char pairwisealisfile[NAMELEN]; // output filename with pairwise alignments
  char alnfile[NAMELEN];  // name of output alignment file in A3M format (for iterative search)
  char hhmfile[NAMELEN];  // name of output HHM file for (iterative search)
  char psifile[NAMELEN];  // name of output alignmen file in PSI-BLAST format (iterative search)
  char scorefile[NAMELEN];// table of scores etc for all HMMs in searched database
  char indexfile[NAMELEN];// optional file containing indeices of aligned residues in given alignment
  char tfile[NAMELEN];    // template filename (in hhalign)
  char wfile[NAMELEN];    // weights file generated with hhformat
  char alitabfile[NAMELEN]; // where to write pairs of aligned residues (-atab option)
  char* dbfiles;          // database filenames, separated by colons
  char* exclstr;          // optional string containing list of excluded residues, e.g. '1-33,97-168'
  int aliwidth;           // number of characters per line in output alignments for HMM search
  char append;            // append to output file? (hhmake)
  float p;                // minimum probability for inclusion in hit list and alignments
  double E;               // maximum E-value for inclusion in hit list and alignment list
  double e;               // maximum E-value for inclusion in output alignment, output HMM, and PSI-BLAST checkpoint model
  int Z;                  // max number of lines in hit list
  int z;                  // min number of lines in hit list
  int B;                  // max number of lines in alignment list
  int b;                  // min number of lines in alignment list
  int showcons;           // in query-template alignments  0: don't show consensus sequence   1:show
  int showdssp;           // in query-template alignments  0: don't show ss_dssp lines        1:show
  int showpred;           // in query-template alignments  0: don't show ss_pred and ss_conf lines  1:show
  int showconf;           // in query-template alignments  0: don't show ss_conf lines        1:show
  char cons;              // if set to 1, include consensus as first representative sequence of HMM
  int nseqdis;            // maximum number of query or template sequences in output alignments
  char mark;              // which sequences to mark for display in output alignments? 0: auto; 1:all
  char outformat;         // 0: hhr  1: FASTA  2:A2M   3:A3M
                          //0:MAC alignment, master-slave  1:MAC blending, master-slave  2:MAC alignment, combining

  int max_seqid;          // Maximum sequence identity with all other sequences in alignment
  int qid;                // Minimum sequence identity with query sequence (sequence 0)
  float qsc;              // Minimum score per column with query sequence (sequence 0)
  int coverage;           // Minimum coverage threshold
  int Ndiff;              // Pick Ndiff most different sequences that passed the other filter thresholds
  bool allseqs;           // if true, do not filter in output alignment; show all sequences

  int Mgaps;              // Maximum percentage of gaps for match states
  int M;                  // Match state assignment by  1:upper/lower case  2:percentage rule  3:marked sequence
  char matrix;            // Subst.matrix 0: Gonnet, 1: HSDM, 2: BLOSUM50

  char wg;                // 0: use local sequence weights   1: use local ones

  Pseudocounts::Params pc;     // Pseudocounts parameters
  Pseudocounts::Params pre_pc; // Pseudocounts parameters for prefiltering

  int aa_pcm;              // Admixture method
  float aa_pca;            // Admixture parameter a
  float aa_pcb;            // Admixture parameter b
  float aa_pcc;            // Admixture parameter c

  float gapb;             // Diversity threshold for adding pseudocounts to transitions from M state
  float gapd;             // Gap open penalty factor for deletions
  float gape;             // Gap extend penalty: factor to multiply hmmer values (def=1)
  float gapf;             // factor for increasing/reducing the gap opening penalty for deletes
  float gapg;             // factor for increasing/reducing the gap opening penalty for inserts
  float gaph;             // factor for increasing/reducing the gap extension penalty for deletes
  float gapi;             // factor for increasing/reducing the gap extension penalty for inserts

  float egq;              // penalty for end gaps when query not fully covered
  float egt;              // penalty for end gaps when template not fully covered

  float Neff;

  char ssgap;             // 1: add secondary structure-dependent gap penalties  0:off
  float ssgapd;           // secondary structure-dependent gap-opening penalty (per residue)
  float ssgape;           // secondary structure-dependent gap-extension penalty (per residue)
  char ssgapi;            // max. number of inside-integer(ii); gap-open-penalty= -ii*ssgapd

  char ssm;               // SS comparison mode: 0:no ss scoring  1:ss scoring AFTER alignment  2:ss score in column score
  float ssw;              // SS weight as compared to column score
  float ssw_realign;      // SS weight as compared to column score for realign
  float ssa;              // SS state evolution matrix M1 = (1-ssa)*I + ssa*M0

  char loc;               // 0: local alignment (wrt. query), 1: global alignement
  char forward;           // 0:Viterbi algorithm  1:Forward algorithm  2: MAC
  char realign;           // realign database hits to be displayed with MAC algorithm
  int altali;             // find up to this many possibly overlapping alignments
  int columnscore;        // 0: no aa comp corr  1: 1/2(qav+tav) 2: template av freqs 3: query av freqs 4:...
  int half_window_size_local_aa_bg_freqs; // half-window size to average local aa background frequencies
  float corr;             // Weight of correlations between scores with |i-j|<=4
  float shift;            // Score offset for match-match states
  float mact;             // Score threshold (negative offset) in MAC alignment
  int realign_max;        // Realign max ... hits
  float maxmem;           // maximum available memory in GB (approximately)
 
  char calibrate;         // calibration of query HMM?  0:no, 1:yes (write lamda,mu into query profile)
  char calm;              // derive P-values from: 0:query calibration  1:template calibration  2:both  3:Neural Network prediction
  int opt;                // for optimization: compare only every opt'th negative; 0: mode off
  int readdefaultsfile ;  // read defaults file ./.hhdefaults or HOME/.hhdefaults?
  int min_overlap;        // all cells of dyn. programming matrix with L_T-j+i or L_Q-i+j < min_overlap will be ignored
  int hitrank;            // rank of hit to be printed as a3m alignment
  char notags;            // neutralize His-tags, FLAG tags, C-myc tags?
  unsigned int maxdbstrlen; // maximum length of database string to be printed in 'Command' line of hhr file

  int maxcol;             // max number of columns in sequence/MSA input files; must be <= LINELEN and >= maxres
  int maxres;             // max number of states in HMM; must be <= LINELEN
  int maxnumdb;           // max number of hits allowed past prefilter
  int maxnumdb_no_prefilter;// max number of hits without prefiltering

  bool hmmer_used;        // True, if a HMMER database is used

  // Directories for SS-prediction
  int addss;                           // 1: calculate secondary structure 0: don't (default: 0)
  char blast[NAMELEN];                 // BLAST binaries (not needed with csBLAST)
  char psipred[NAMELEN];               // PsiPred binaries
  char psipred_data[NAMELEN];          // PsiPred data
  char dummydb [NAMELEN];

  // parameters for context-specific pseudocounts
  float csb;
  float csw;
  char clusterfile[NAMELEN];
  bool nocontxt;
  
  // HHblits
  int premerge;
  int dbsize;           // number of clusters of input database

  // HHblits Evalue calculation  (alpha = a + b(Neff(T) - 1)(1 - c(Neff(Q) - 1)) )
  float alphaa;
  float alphab;
  float alphac;

  // For filtering database alignments in HHsearch and HHblits 
  // JS: What are these used for? They are set to the values without _db anyway.
  int max_seqid_db;
  int qid_db;      
  float qsc_db;    
  int coverage_db; 
  int Ndiff_db;    

  // HHblits context state prefilter
  char cs_library[NAMELEN];

  // HHblits prefilter
  bool prefilter;             // perform prefiltering in HHblits?
  bool early_stopping_filter; // Break HMM search, when the sum of the last N HMM-hit-Evalues is below threshold

  Hash<int*>* block_shading;         // Cross out cells not covered by prefiltering hit in HHblits
  Hash<int>* block_shading_counter;  // Cross out cells not covered by prefiltering hit in HHblits
  int block_shading_space;           // space added to the rands of prefilter HSP
  char block_shading_mode[NAMELEN];

  // For HHblits prefiltering with SSE2
  short prefilter_gap_open;
  short prefilter_gap_extend;
  int prefilter_score_offset;
  int prefilter_bit_factor;
  double prefilter_evalue_thresh;
  double prefilter_evalue_coarse_thresh;
  int preprefilter_smax_thresh;
  int prefilter_min_alis;

  float wstruc;          // weight of structure scores in hhalign
  int idummy;
  float fdummy;

  bool useCSScoring;
  char cs_template_file[NAMELEN];

  void SetDefaultPaths(char *program_path);
  void SetDefaults();
  Parameters();
};


void Parameters::SetDefaultPaths(char *program_path)
{
  // set hhlib
  FILE* testf = NULL;
  if(getenv("HHLIB"))
    strcpy(hhlib, getenv("HHLIB"));
  else
    strcpy(hhlib, "/usr/lib/hh");

  strcat(strcpy(hhdata, hhlib), "/data");
  strcat(strcpy(clusterfile, hhdata), "/context_data.crf");
  strcat(strcpy(cs_library, hhdata), "/cs219.lib");
  
  testf = fopen(cs_library, "r");
  if (testf) fclose(testf); 
  else 
    {
      if (v>=3) cerr<<"WARNING in HHsuite: Could not open "<<cs_library<<"\n";
      
      /* we did not find HHLIB, if called with full path or in dist dir, we can try relative to program path */
      if(program_path != NULL)
	{
	  strcat(strcpy(hhlib, program_path), "../lib/hh");
	  strcat(strcpy(hhdata, hhlib), "/data");
	  strcat(strcpy(clusterfile, hhdata), "/context_data.crf");
	  strcat(strcpy(cs_library, hhdata), "/cs219.lib");
	  testf = fopen(cs_library, "r");
	  if (testf) fclose(testf);
	  else 
	    {
	      if (v>=3) cerr<<"WARNING in HHsuite: Could not open "<<cs_library<<"\n";
	      
	      strcat(strcpy(hhlib, program_path), "..");
	      strcat(strcpy(hhdata, hhlib), "/data");
	      strcat(strcpy(clusterfile, hhdata), "/context_data.crf");
	      strcat(strcpy(cs_library, hhdata), "/cs219.lib");	  
	      testf = fopen(cs_library, "r");
	      if (testf) fclose(testf);
	      else 
		if (v>=3) cerr<<"WARNING in HHsuite: Could not open "<<cs_library<<"\n";
	    }
	}
    }
  if (!testf)
    {
      cerr<<endl<<"Error in "<<argv[0]<<": could not find context_data.crf and cs219.lib in '" << hhlib << "'.\n"
	"Please set the HHLIB environment variable to the HH-suite directory\n"
	"(Linux bash: export HHLIB=<hh_dir>, csh/tcsh: setenv HHLIB=<hh_dir>).\n"
	"The missing files should be in $HHLIB/data/.\n ";
      exit(2);
    }
  return;
}


void Parameters::SetDefaults()
{

  // Moved from hhdecl.C 
  v=2;

  // Parameter class
  maxcol=32765;            // max number of columns in sequence/MSA input files; must be <= LINELEN and >= maxres
  maxres=15002;            // max number of states in HMM; must be <= LINELEN
  maxnumdb=20000;          // max number of hits allowed past prefilter
  maxnumdb_no_prefilter=20000;// max number of hits without prefiltering

  append=0;                // overwrite output file
  outformat=0;             // 0: hhr  1: FASTA  2:A2M   3:A3M
  p=20.0f;                 // minimum threshold for inclusion in hit list and alignment listing
  E=1e6f;                  // maximum threshold for inclusion in hit list and alignment listing
  b=10;                    // min number of alignments
  B=500;                   // max number of alignments
  z=10;                    // min number of lines in hit list
  Z=500;                   // max number of lines in hit list
  e=1e-3f;                 // maximum E-value for inclusion in output alignment, output HMM, and PSI-BLAST checkpoint model
  realign_max=1000;        // Realign max ... hits
  maxmem = 3.0;            // 3GB
  showcons=1;              // show consensus sequence
  showdssp=1;              // show predicted secondary structure ss_dssp
  showpred=1;              // show predicted secondary structure ss_pred
  showconf=0;              // don't show secondary structure confidence ss_conf
  cons=0;                  // chose first non-SS sequence as main representative sequence (not consensus)
  nseqdis=1;               // maximum number of query sequences for output alignment
  mark=0;                  // 1: only marked sequences (or first) get displayed; 0: most divergent ones get displayed
  aliwidth=80;             // number of characters per line in output alignments for HMM search

  max_seqid=90;            // default for maximum sequence identity threshold
  qid=0;                   // default for minimum sequence identity with query
  qsc=-20.0f;              // default for minimum score per column with query
  coverage=0;              // default for minimum coverage threshold
  Ndiff=100;               // pick Ndiff most different sequences from alignment
  allseqs = false;         // if true, do not filter result MSA; show all sequences

  Neff=0;                  // Filter alignment to a diversity (Neff) with a maximum Neff of par.Neff

  M=1;                     // match state assignment is by A2M/A3M
  Mgaps=50;                // Above this percentage of gaps, columns are assigned to insert states (for par.M=2)
  calibrate=0;             // default: no calibration
  calm=3;                  // derive P-values from: 0:query calibration  1:template calibration  2:both  3:Neural Network prediction

  wg=0;                    // 0: use local sequence weights   1: use local ones

  matrix=0;                // Subst.matrix 0: Gonnet, 1: HSDM, 2: BLOSUM50 3: BLOSUM62

  pc.admix       = Pseudocounts::HHsearchAdmix;
  pc.pca         = 0.90;
  pc.pcb         = 4.00;
  pc.pcc         = 1.0;
  pc.target_neff = 0.0;

  pre_pc.admix       = Pseudocounts::CSBlastAdmix;
  pre_pc.pca         = 0.80;
  pre_pc.pcb         = 2.00;
  pre_pc.pcc         = 1.0;
  pre_pc.target_neff = 0.0;

  aa_pcm=2;
  aa_pca=1.0f;
  aa_pcb=1.5f;
  aa_pcc=1.0f;

  gapb=1.0;                // default values for transition pseudocounts
  gapd=0.15;               // gap open penalty pseudocount; 0.25 corresponds to 7.1*gapf bits
  gape=1.0;                // gap extension penalty pseudocount
  gapf=0.6;                // factor for increasing gap open penalty for deletes
  gapg=0.6;                // factor for increasing gap open penalty for inserts
  gaph=0.6;                // factor for increasing gap extension penalty for deletes
  gapi=0.6;                // factor for increasing gap extension penalty for inserts

  ssm=2;                   // ss scoring mode: 0:no ss score  1:score after alignment  2:score during alignment
  ssw=0.11f;               // weight of ss scoring
  ssw_realign=0.11f;       // weight of ss scoring for realign
  ssa=1.0f;                // weight of ss evolution matrix
  shift=-0.03f;            // Shift match score up
  mact=0.3501f;            // Score threshold for MAC alignment in local mode (set to 0.3501 to track user modification)
  corr=0.1f;               // Weight of correlations of scores for |i-j|<=4

  egq=0.0f;                // no charge for end gaps as default
  egt=0.0f;                // no charge for end gaps as default

  ssgap=0;                 // 1: add secondary structure-dependent gap penalties  0:off
  ssgapd=1.0f;             // secondary structure-dependent gap-opening penalty (per residue)
  ssgape=0.0f;             // secondary structure-dependent gap-extension penalty (per residue)
  ssgapi=4;                // max. number of inside-integer(ii); gap-open-penalty= -ii*ssgapd

  loc=1;                   // local vs. global alignment as default
  altali=2;                // find up to two (possibly overlapping) subalignments
  forward=0;               // 0: Viterbi algorithm; 1: Viterbi+stochastic sampling; 3:Maximum Accuracy (MAC) algorithm
  realign=1;               // realign with MAC algorithm

  columnscore=1;           // Default column score is 1: null model pnul = 1/2 * (q_av(a)+p_av(a))
  half_window_size_local_aa_bg_freqs = 40;
  min_overlap=0;           // automatic minimum overlap used
  opt=0;                   // Default = optimization mode off
  readdefaultsfile=0;      // Default = do not read a defaults file ./.hhdefaults or HOME/.hhdefaults
  maxdbstrlen=200;         // maximum length of database string to be printed in 'Command' line of hhr file
  premerge=0;

  notags=1;                // neutralize His-tags, FLAG-tags, C-myc-tags
  hmmer_used=false;

  // Directories for SS-prediction
  addss=0;
  strcpy(psipred,"");
  strcpy(psipred_data,"");

  // HHblits parameters
  dbsize = 0;

  // HHblits Evalue calculation  (alpha = a + b(Neff(T) - 1)(1 - c(Neff(Q) - 1)) )
  alphaa = 0.4;
  alphab = 0.02;
  alphac = 0.1;

  prefilter = false;              //true in hhblits
  early_stopping_filter = false;  //true in hhblits

  block_shading=NULL;
  block_shading_counter=NULL;
  block_shading_space = 200;
  strcpy(block_shading_mode,"tube");

  // For HHblits prefiltering with SSE2
  prefilter_gap_open = 20;
  prefilter_gap_extend = 4;
  prefilter_score_offset = 50;
  prefilter_bit_factor = 4;
  prefilter_evalue_thresh = 1000;
  prefilter_evalue_coarse_thresh = 100000;
  prefilter_min_alis = 1000;
  preprefilter_smax_thresh = 10;

  // For filtering database alignments in HHsearch and HHblits 
  //JS: What are these used for? They are set to the options without _db anyway.
  max_seqid_db=max_seqid;
  qid_db=qid;            
  qsc_db=qsc;            
  coverage_db=coverage;  
  Ndiff_db=Ndiff;        

  // Initialize strings
//  strcpy(infile,"stdin");
  strcpy(infile,"");
  strcpy(outfile,"");
  strcpy(pairwisealisfile,"");
  strcpy(scorefile,"");
  strcpy(indexfile,""); 
  strcpy(wfile,"");
  strcpy(alnfile,"");
  strcpy(hhmfile,"");
  strcpy(psifile,"");
  strcpy(alitabfile,"");
  exclstr=NULL;
  useCSScoring=false;
  strcpy(cs_template_file, "");

  // parameters for context-specific pseudocounts
  csb = 0.85;
  csw = 1.6;

  wstruc=1.0f;             // Weight of structure score in hhalign
  idummy=0;

  return;
}

Parameters::Parameters()
{
  SetDefaults();
}


// Class to store data about hit to realign
class Realign_hitpos
{
public:
  int index;         // index of template in dbfile (1,2,..)
  long ftellpos;     // position of template in dbfile
  int operator<(const Realign_hitpos& realign_hitpos) {return ftellpos<realign_hitpos.ftellpos;}
};

EXTERN Parameters par;

// cs object declarations
cs::ContextLibrary<cs::AA>* context_lib = NULL; // Context library for pseudocounts generation
cs::Crf<cs::AA>* crf                    = NULL; // CRF for pseudocounts generation
cs::Pseudocounts<cs::AA>* pc            = NULL; // Pseudocounts engine
cs::Admix* pc_admix                     = NULL; // Pseudocounts admixture method
cs::Pseudocounts<cs::AA>* pre_pc        = NULL; // Pseudocounts engine for prefiltering
cs::Admix* pre_pc_admix                 = NULL; // Pseudocounts admixture method for prefiltering

void InitializePseudocountsEngine() {
  // Prepare pseudocounts engine
  FILE* fin = fopen(par.clusterfile, "r");
  if (!fin) {
    cerr<<endl<<"Error in "<<par.argv[0]<<": could not open file \'"<<par.clusterfile<<"\'\n";
    exit(2);
  }
  char ext[100];
  Extension(ext, par.clusterfile);
  if (strcmp(ext, "crf") == 0)  {
    crf = new cs::Crf<cs::AA>(fin);
    pc = new cs::CrfPseudocounts<cs::AA>(*crf);
    pre_pc = new cs::CrfPseudocounts<cs::AA>(*crf);
  } else {
    context_lib = new cs::ContextLibrary<cs::AA>(fin);
    cs::TransformToLog(*context_lib);
    pc = new cs::LibraryPseudocounts<cs::AA>(*context_lib, par.csw, par.csb);
    pre_pc = new cs::LibraryPseudocounts<cs::AA>(*context_lib, par.csw, par.csb);
  }
  fclose(fin);
  pc->SetTargetNeff(par.pc.target_neff);
  pre_pc->SetTargetNeff(par.pre_pc.target_neff);

  // Prepare pseudocounts admixture method
  pc_admix = par.pc.CreateAdmix();
  pre_pc_admix = par.pre_pc.CreateAdmix();
}

void DeletePseudocountsEngine() {
  if (context_lib != NULL) delete context_lib;
  if (crf != NULL) delete crf;
  if (pc != NULL) delete pc;
  if (pc_admix != NULL) delete pc_admix;
  if (pre_pc != NULL) delete pre_pc;
  if (pre_pc_admix != NULL) delete pre_pc_admix;
}
