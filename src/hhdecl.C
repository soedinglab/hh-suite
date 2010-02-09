/////////////////////////////////////////////////////////////////////////////////////
//// Constants
/////////////////////////////////////////////////////////////////////////////////////

const char VERSION_AND_DATE[]="version 1.6.0.0 (April 2009)";
const char REFERENCE[]="Soding, J. Protein homology detection by HMM-HMM comparison. Bioinformatics 2005, 21, 951-960.\n";
const char COPYRIGHT[]="(C) Johannes Soeding (see LICENSE file)\n";
const int MAXSEQ=65535; //max number of sequences in input alignment (must be <~30000 on cluster nodes)
const int MAXCOL=32765; //max number of residues in input files; must be <= LINELEN and >= MAXRES
const int MAXRES=15002; //max number of columns in HMM; must be <= LINELEN
const int LINELEN=262144; //max length of line read in from input files; must be >= MAXCOL
const int MAXSEQDIS=10238;//max number of sequences stored in 'hit' objects and displayed in output alignment
const int IDLEN=255;     //max length of scop hierarchy id and pdb-id
const int DESCLEN=32765;//max length of sequence description (longname)
const int NAMELEN=511;  //max length of file names etc.
const int MAXOPT=127;   //Maximum number of options to be read in from .hhconfig or command line
const int NAA=20;       //number of amino acids (0-19)
const int NTRANS=7;    //number of transitions recorded in HMM (M2M,M2I,M2D,I2M,I2I,D2M,D2D)
const int NCOLMIN=10;   //min number of cols in subalignment for calculating pos-specific weights w[k][i]
const int ANY=20;       //number representing an X (any amino acid) internally
const int GAP=21;       //number representing a gap internally
const int ENDGAP=22;    //Important to distinguish because end gaps do not contribute to tansition counts
const int HMMSCALE=1000;//Scaling number for log2-values in HMMs
const int NFAMMAX=5119; //Size of hash for counting number of HMMs in each family
const int MAXPROF=32766;//Maximum number of HMM scores for fitting EVD
const float MAXENDGAPFRAC=0.1; //For weighting: include only columns into subalignment i that have a max fraction of seqs with endgap
const float SMIN= 20.;  //Minimum score of hit needed to search for another repeat of same profile: p=exp(-(4-mu)/lamda)=0.01
const float LAMDA=0.388; //lamda in score EVD used for -local mode in length correction: S = S-log(Lq*Lt)/LAMDA)
const float LAMDA_GLOB=0.42; //lamda in score EVD used for -global mode
const float PMAX=1E-2;  //Maximum single-repeat p-value that can contribute to whole-protein p-value
const float MINEVALEXCL=0.5; //above this E-value from first ML fit hits are not used for final ML fit of EVD
const int SELFEXCL=3;   // exclude self-alignments with j-i<SELFEXCL
const float PLTY_GAPOPEN=6.0f; // for -qsc option (filter for min similarity to query): 6 bits to open gap
const float PLTY_GAPEXTD=1.0f; // for -qsc option (filter for min similarity to query): 1 bit to extend gap
const int MINCOLS_REALIGN=6; // hits with MAC alignments with fewer matched columns will be deleted in hhsearch hitlist
const float LOG1000=log(1000.0);


enum transitions {M2M,M2I,M2D,I2M,I2I,D2M,D2D}; // index for transitions within a HMM
enum pair_states {STOP=0,SAME=1,GD=2,IM=3,DG=4,MI=5,MS=6,ML=7,SM=8,LM=9,MM=10};

// const char aa[]="ARNDCQEGHILKMFPSTWYVX-";
//Amino acids Sorted by alphabet     -> internal numbers a
//                0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
//                A  C  D  E  F  G  H  I  K  L  M  N  P  Q  R  S  T  V  W  Y  X
const int s2a[]={ 0, 4, 3, 6,13, 7, 8, 9,11,10,12, 2,14, 5, 1,15,16,19,17,18,20};
//Internal numbers a for amino acids -> amino acids Sorted by alphabet:
//                0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
//                A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  X
const int a2s[]={ 0,14,11, 2, 1,13, 3, 5, 6, 7, 9, 8,10, 4,12,15,16,18,19,17,20};

// Secondary structure
const int NDSSP=8;      //number of different ss states determined by dssp: 0-7 (0: no state available)
const int NSSPRED=4;    //number of different ss states predicted by psipred: 0-3 (0: no prediction availabe)
const int MAXCF=11;     //number of different confidence values: 0-10 (0: no prediction availabe)
const int NSA=7;        //number of classes relative solvent accesiblity (0:no coord,  1:<2%, 2:<14%, 3:<33%, 4:<55%, 5:>55%, 6:S-S bridge)

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

// Input parameters
class Parameters          // Parameters for gap penalties and pseudocounts
{
public:
  char** argv;            //command line parameters
  char argc;              //dimension of argv

  char infile[NAMELEN];   // input filename
  char outfile[NAMELEN];  // output filename
  char pairwisealisfile[NAMELEN]; // output filename with pairwise alignments
  char alnfile[NAMELEN];  // name of output alignment file in A3M format (for iterative search)
  char hhmfile[NAMELEN];  // name of output HHM file for (iterative search)
  char psifile[NAMELEN];  // name of output alignmen file in PSI-BLAST format (iterative search)
  char scorefile[NAMELEN];// table of scores etc for all HMMs in searched database
  char tfile[NAMELEN];    // template filename (in hhalign)
  char buffer[NAMELEN];   // buffer to write results for other programs into
  char pngfile[NAMELEN];  // png image file for dotplot
  char wfile[NAMELEN];    // weights file generated with hhformat
  char* blafile;          // output of 'blastpgp -m 8' with PSI-BLAST E-values for HHblits
  float hhblits_prefilter_logpval;  // PSI-BLAST prefilter log P-value threshold of HHblits
  char* dbfiles;          // database filenames, separated by colons
  char* exclstr;          // optional string containing list of excluded residues, e.g. '1-33,97-168'
  int aliwidth;           // number of characters per line in output alignments for HMM search
  char append;            // append to output file? (hhmake)
  float p;                // minimum probability for inclusion in hit list and alignments
  float E;                // maximum E-value for inclusion in hit list and alignment list
  float e;                // maximum E-value for inclusion in output alignment, output HMM, and PSI-BLAST checkpoint model
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
  char mode;              //
                          //0:MAC alignment, master-slave  1:MAC blending, master-slave  2:MAC alignment, combining

  int max_seqid;          // Maximum sequence identity with all other sequences in alignment
  int qid;                // Minimum sequence identity with query sequence (sequence 0)
  float qsc;              // Minimum score per column with query sequence (sequence 0)
  int coverage;           // Minimum coverage threshold
  int Ndiff;              // Pick Ndiff most different sequences that passed the other filter thresholds

  int Mgaps;              // Maximum percentage of gaps for match states
  int M;                  // Match state assignment by  1:upper/lower case  2:percentage rule  3:marked sequence
  char matrix;            // Subst.matrix 0: Gonnet, 1: HSDM, 2: BLOSUM50

  char wg;                // 0: use local sequence weights   1: use local ones

  char pcm;               // 0:no pseudocounts, 1:pos-specific pcs, 2:PSIBLAST pcs
  float pca;              // Pseudocount matrix = (1-tau(i))*I + tau(i)*S
  float pcb;              // tau(i) = pca/(1 + ((Neff-1)/pcb)^pcc
  float pcc;              //
  float pcw;              // Decrease pseudocounts for conserved columns

  float gapb;             // Diversity threshold for adding pseudocounts to transitions from M state
  float gapd;             // Gap open penalty factor for deletions
  float gape;             // Gap extend penalty: factor to multiply hmmer values (def=1)
  float gapf;             // factor for increasing/reducing the gap opening penalty for deletes
  float gapg;             // factor for increasing/reducing the gap opening penalty for inserts
  float gaph;             // factor for increasing/reducing the gap extension penalty for deletes
  float gapi;             // factor for increasing/reducing the gap extension penalty for inserts

  float egq;              // penalty for end gaps when query not fully covered
  float egt;              // penalty for end gaps when template not fully covered

  float neffa;            // Coefficients to estimate Neff-dependent weights for HMM merging procedure
  float neffb;            // Coefficients to estimate Neff-dependent weights for HMM merging procedure

  char ssgap;             // 1: add secondary structure-dependent gap penalties  0:off
  float ssgapd;           // secondary structure-dependent gap-opening penalty (per residue)
  float ssgape;           // secondary structure-dependent gap-extension penalty (per residue)
  char ssgapi;            // max. number of inside-integer(ii); gap-open-penalty= -ii*ssgapd

  char ssm;               // SS comparison mode: 0:no ss scoring  1:ss scoring AFTER alignment  2:ss score in column score
  float ssw;              // SS weight as compared to column score
  float ssa;              // SS state evolution matrix M1 = (1-ssa)*I + ssa*M0

  char loc;               // 0: local alignment (wrt. query), 1: global alignement
  char forward;           // 0:Viterbi algorithm  1:Forward algorithm  2: MAC
  char realign;           // realign database hits to be displayed with MAC algorithm
  char altali;            // find up to this many possibly overlapping alignments
  int columnscore;        // 0: no aa comp corr  1: 1/2(qav+tav) 2: template av freqs 3: query av freqs 4:...
  float corr;             // Weight of correlations between scores with |i-j|<=4
  float shift;            // Score offset for match-match states
  float mact;             // Score threshold (negative offset) in MAC alignment

  char calibrate;         // calibration of query HMM?  0:no, 1:yes (write lamda,mu into query profile)
  char calm;              // derive P-values from: 0:query calibration  1:template calibration  2:both  3:Neural Network prediction
  int opt;                // for optimization: compare only every opt'th negative; 0: mode off
  int readdefaultsfile ;  // read defaults file ./.hhdefaults or HOME/.hhdefaults?
  int min_overlap;        // all cells of dyn. programming matrix with L_T-j+i or L_Q-i+j < min_overlap will be ignored
  int hitrank;            // rank of hit to be printed as a3m alignment
  char notags;            // neutralize His-tags, FLAG tags, C-myc tags?
  unsigned int maxdbstrlen; // maximum length of database string to be printed in 'Command' line of hhr file

  char trans;             // 0: normal pairwise scoring; 1:transitive scoring
  float Emax_trans;       // max E-value for intermediate HMMs in transitive scoring (i.e. l is intermediate HMM if E_lq, E_lk <Emax_trans)
  float wtrans;           // Ztot[k] = Zq[k] + wtrans * (Zforward[k]+Zreverse[k])

  // parameters for context-specific pseudocounts
  float csb;
  float csw;
  char clusterfile[NAMELEN];

  // For filtering database alignments in HHsearch and HHblits
  int max_seqid_db;
  int qid_db;      
  float qsc_db;    
  int coverage_db; 
  int Ndiff_db;    

  double filter_thresh;    // Threshold for early stopping
  int filter_length;       // Length of array of 1/evalues
  double *filter_evals;    // array of last 1/evalues
  double filter_sum;       // sum of evalues in array
  int filter_counter;      // counter for evalue array

  Hash<int*>* block_shading;         // Cross out cells not covered by prefiltering hit in HHblits
  Hash<int>* block_shading_counter;  // Cross out cells not covered by prefiltering hit in HHblits
  int block_shading_space;           // space added to the rands of prefilter HSP
  char block_shading_mode[NAMELEN];

  // For HHblits prefiltering with SSE2
  int prefilter_lmax;              // maximum block length of each DB-Sequence
  int sse_shading_space;           // space added to the rands of prefilter HSP
  short prefilter_gap_open;
  short prefilter_gap_extend;
  int prefilter_states;               // Anzahl der States im Alphabet
  int prefilter_db_overlap;           // Overlap if DB-Sequence > par.prefilter_lmax
  int prefilter_score_offset;
  int prefilter_bit_factor;
  int prefilter_smax_thresh;
  int prefilter_rmax_thresh;
  

  // SCRAP THE FOLLOWING VARIABLES?

  float wstruc;          // weight of structure scores
  char repmode;          // 1:repeat identification: multiple hits not treated as independent 0: repeat mode off
  int idummy;
  int jdummy;
  float fdummy;

};

/////////////////////////////////////////////////////////////////////////////////////
//// Global variable declarations
/////////////////////////////////////////////////////////////////////////////////////

char v=2;             // 1: show only warnings 2:verbose mode
Parameters par;
char program_name[NAMELEN]; //name of program executed (e.g. hhmake of hhsearch)

// substitution matrix flavours
float P[21][21];      // P[a][b] = combined probability for a aligned to b
float R[21][21];      // R[a][b]=P[a][b]/p[b]=P(a|b); precalculated for pseudocounts
float Sim[21][21];    // Similarity matrix Sim[a][b]: how similar are a and b?
float S[21][21];      // Substitution score matrix S[a][b] = log2(Pab/pa/pb)
float pb[21];         // pb[a] = background amino acid probabilities for chosen substitution matrix
float qav[21];        // qav[a] = background amino acid probabilities for query HMM (needed for rate matrix rescaling)

// secondary structure matrices
float S73[NDSSP][NSSPRED][MAXCF];           // P[A][B][cf]       =  log2 P(A,B,cf)/P(A)/P(B,cf)
float S33[NSSPRED][MAXCF][NSSPRED][MAXCF];  // P[B][cf][B'][cf'] =  log2 sum_B' P(A,B',cf)/P(A)/P(B,cf) * P_b(B'|B)
// float S77[NDSSP][DSSP];                  // P[A][B]           =  log2 P(A,B)/P(A)/P(B)

