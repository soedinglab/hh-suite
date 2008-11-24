/////////////////////////////////////////////////////////////////////////////////////
// Class representing a2m/a3m-formatted alignment corresponding to one HMM
/////////////////////////////////////////////////////////////////////////////////////
class HalfAlignment
{
public:
  HalfAlignment(int maxseqdis=MAXSEQDIS);
  ~HalfAlignment();

  // Initialize HalfAlignment; create index arrays s [],l[], m[]
  void Set(char* name, char** seq_in, char** sname_in, int n_in, int L_in, int n1, int n2, int n3, int n4, int nc);

  // Free memory in HalfAlignment arrays s[][], l[][], and m[][]
  void Unset();

  // Align query (HalfAlignment) to template (i.e. hit) match state structure
  void AlignToTemplate(Hit& hit); 

  // Write the a2m/a3m query alignment into alnfile 
  void Print(char* outfile);

  // Fill in insert states following match state i
  void AddInserts(int i);

  // Fill up alignment with gaps '.' to generate flush end (all h[k] equal)
  void FillUpGaps();

  // Fill in insert states following match state i and fill up gaps with '.' 
  void AddInsertsAndFillUpGaps(int i);

  // Add gap column '.' or character column
  void AddChar(char c);

  // Add match state column i as is
  void AddColumn(int i);

  // Add match state column i as insert state
  void AddColumnAsInsert(int i);

  // Build alignment in FASTA format
  void BuildFASTA();

  // Build alignment in A2M format
  void BuildA2M();

  // Build alignment in A3M format
  void BuildA3M();

  // Transform alignment sequences from A2M to FASTA ( lowercase to uppercase and '.' to '-')
  void ToFASTA();

  // Remove all characters c from template sequences
  void RemoveChars(char c);

private:
  friend class FullAlignment;

  int n;           //number of sequences in half-alignment
  char** seq;      //sequences (in ASCII) in alignment
  char** sname;    //names of sequences in alignment
  int nss_dssp;    //index of sequence with dssp sec structure states
  int nsa_dssp;    //index of sequence with dssp solvent accessiblity states
  int nss_pred;    //index of sequence with predicted sec structure states
  int nss_conf;    //index of sequence with prediction confidence values
  int ncons;       //index of consensus sequence
    
  int pos;         //After FillUpGaps() all h[k] have value pos 
  int L;           //number of match states in corresponding profile
  int* h;          //h[k] = next position of sequence k to be written
  char** s;        //s[k][h] = column h, sequence k of output alignment
  int** l;         //counts non-gap residues: l[k][i] = index of last residue AT OR BEFORE match state i in seq k 
  int** m;         //counts positions:        m[k][i] = position of match state i in string seq[k]
/*   int h[MAXSEQ];   //h[k] = next position of sequence k to be written */
/*   char* s[MAXSEQ]; //s[k][h] = column h, sequence k of output alignment */
/*   int* l[MAXSEQ];  //counts non-gap residues: l[k][i] = index of last residue AT OR BEFORE match state i in seq k  */
/*   int* m[MAXSEQ];  //counts positions:        m[k][i] = position of match state i in string seq[k] */
};


