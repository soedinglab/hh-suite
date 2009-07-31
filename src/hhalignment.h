// hhalignment.h

class Alignment
{
public:
  int L;                  // number of match states of alignment
  int N_in;               // total number of sequences in alignment
  int N_filtered;         // number of sequences after sequence identity filtering
  int N_ss;               // number of >ss_ or >sa sequences

  int kss_dssp;           // index of sequence with secondary structure by dssp      -1:no >ss_dssp line found 
  int ksa_dssp;           // index of sequence with solvent accessibility by dssp    -1:no >sa_dssp line found 
  int kss_pred;           // index of sequence with predicted secondary structure    -1:no >ss_pred line found
  int kss_conf;           // index of sequence with confidence values of prediction  -1:no >ss_conf line found
  int kfirst;             // index of first real sequence

  char* longname        ; // Full name of first sequence of original alignment (NAME field)
  char name[NAMELEN];     // HMM name = first word in longname in lower case
  char fam[NAMELEN];      // family ID (derived from name) (FAM field)
  char file[NAMELEN];     // Rootname (w/o path, with extension) of alignment file that is used to construct the HMM

  int n_display;          // number of sequences to be displayed (INCLUDING >ss_pred, >ss_conf, >ss_dssp sequences) 
  char** sname;           // names of display sequences (first seq=0, first char=0)
  char** seq;             // residues of display sequences (first char=1)
  int* l;                 // l[i] = position of i'th match state in alignment 

  char* keep;             // keep[k]=1 if sequence is included in amino acid frequencies; 0 otherwise (first=0)

  Alignment(int maxseq=MAXSEQ, int maxres=MAXRES);
  ~Alignment();

  // Read alignment into X (uncompressed) in ASCII characters
  void Read(FILE* inf, char infile[NAMELEN], char* line=NULL);
  
  // Convert ASCII to numbers between 0 and 20, throw out all insert states, 
  // record their number in I[k][i] and store sequences to be displayed in sname[k] and seq[k]
  void Compress(const char infile[NAMELEN]);

  // Apply sequence identity filter
  inline int FilterForDisplay(int max_seqid, int coverage=0, int qid=0, float qsc=0, int N=0);
  inline int Filter(int max_seqid, int coverage=0, int qid=0, float qsc=0, int N=0);
  int Filter2(char keep[], int coverage, int qid, float qsc, int seqid1, int seqid2, int Ndiff);

  // Filter alignment for min score per column with core query profile, defined by min_coverage_core and min_seqid_core
  int HomologyFilter(int coverage_core, float qsc_core, float coresc);

  // Calculate AA frequencies q.p[i][a] and transition probabilities q.tr[i][a] from alignment
  void FrequenciesAndTransitions(HMM& q, char* in=NULL, bool time=false);

  // Calculate freqs q.f[i][a] and transitions q.tr[i][a] (a=MM,MI,MD) with pos-specific subalignments
  void Amino_acid_frequencies_and_transitions_from_M_state(HMM& q, char* in);

  // Calculate transitions q.tr[i][a] (a=DM,DD) with pos-specific subalignments
  void Transitions_from_D_state(HMM& q, char* in);

  // Calculate transitions q.tr[i][a] (a=DM,DD) with pos-specific subalignments
  void Transitions_from_I_state(HMM& q, char* in);
  
  // Write alignment without insert states to alignment file
  void WriteWithoutInsertsToFile(const char* alnfile);

  // Write alignment to alignment file
  void WriteToFile(const char* alnfile, const char format[]=NULL);

  // Read a3m slave alignment of hit from ta3mfile and merge into (query) master alignment
  void MergeMasterSlave(Hit& hit, char ta3mfile[]);

  // Read a3m alignment of hit from ta3mfile and merge-combine with query alignment
  void Merge(Hit& hit, char ta3mfile[]);

  // Add a sequence to Qali
  void AddSequence(char Xk[], int Ik[]=NULL);

  // Add SS prediction to Qali
  void AddSSPrediction(char seq_pred[], char seq_conf[]);

  // Determine matrix of position-specific weights w[k][i] for multiple alignment
  void GetPositionSpecificWeights(float* w[]);


private:
  char** X;               // X[k][i] contains column i of sequence k in alignment (first seq=0, first char=1) (0-3: ARND ..., 20:X, 21:GAP)
  short unsigned int** I; // I[k][i] contains the number of inserts AFTER match state i (first=0)
  char* display;          // display[k]=1 if sequence will be displayed in output alignments; 0 otherwise (first=0)
  float* wg;              // w[k] = global weight of sequence k
  int* nseqs;             // number of sequences in subalignment i (only for DEBUGGING)
  int* nres;              // number of residues in sequence k
  int* first;             // first residue in sequence k
  int* last;              // last  residue in sequence k
  int* ksort;             // index for sorting sequences: X[ksort[k]]
  int FilterWithCoreHMM(char in[], float coresc, HMM& qcore);
  
  char readCommentLine;   // Set to 1, if a comment line with '#' is read
};
