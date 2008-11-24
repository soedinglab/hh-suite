/////////////////////////////////////////////////////////////////////////////////////
// Class for output alignment of query against template sequences
/////////////////////////////////////////////////////////////////////////////////////
class FullAlignment
{
public:
  FullAlignment(int maxseqdis=MAXSEQDIS);
  ~FullAlignment();
  void FreeMemory();
  void Build(HMM& q, Hit& hit);
  void PrintHeader(FILE* outf, HMM& q, Hit& hit);
  void PrintHHR(FILE* outf, Hit& hit);
  void PrintA2M(FILE* outf, Hit& hit);
  void PrintFASTA(FILE* outf, Hit& hit);
  void PrintA3M(FILE* outf, Hit& hit);
  int identities;      // number of identical residues in query and template sequence
  float score_sim;     // substitution matrix similarity score between query and template

private:
  HalfAlignment* qa; //query and template parts of the alignment
  HalfAlignment* ta; //query and template parts of the alignment
  char symbol[LINELEN];         //symbol[h] = symbol (= - . + |) indicating match score for col h of alignment    
  void ClearSymbols()      {for (int h=0; h<LINELEN-1; h++) symbol[h]=' ';}
  void AddColumns(int i, int j, char prev_state, char state, float S);
  void AddGaps();
  int ScoreChr(float S) {return (S<-1.5?'=':(S<-0.5?'-':(S<0.5?'.':(S<1.5?'+':'|'))));}
};
