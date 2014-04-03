/////////////////////////////////////////////////////////////////////////////////////
// Class for output alignment of query against template sequences
/////////////////////////////////////////////////////////////////////////////////////

#include "hhhalfalignment.h"

#ifndef HHFULLALIGNMENT_H
#define HHFULLALIGNMENT_H

#include <sstream>

class FullAlignment
{
public:
  FullAlignment(int maxseqdis=MAXSEQDIS);
  ~FullAlignment();
  void FreeMemory();
  void Build(HMM* q, Hit& hit);
  void PrintHeader(std::stringstream& out, HMM* q, Hit& hit);
  void PrintHHR(std::stringstream& out, Hit& hit);
  void PrintA2M(std::stringstream& out, Hit& hit);
  void PrintFASTA(std::stringstream& out, Hit& hit);
  void PrintA3M(std::stringstream& out, Hit& hit);
  int identities;      // number of identical residues in query and template sequence
  float score_sim;     // substitution matrix similarity score between query and template

private:
  HalfAlignment* qa; // query and template parts of the alignment
  HalfAlignment* ta; // query and template parts of the alignment
  char symbol[LINELEN];         // symbol[h] = symbol (= - . + |) indicating match score for col h of alignment    
  char posterior[LINELEN];      // posterior probability for pair of aligned columns 
  void ClearSymbols() {for (int h=0; h<LINELEN-1; h++) symbol[h]=posterior[h]=' ';}
  void AddColumns(int i, int j, char prev_state, char state, float S, float PP);
  void AddGaps();
  char ScoreChr(float S) {return S<-1.5?'=':(S<-0.5?'-':(S<0.5?'.':(S<1.5?'+':'|')));}
  char PosteriorChr(float PP) {return 48+imax(0,imin(9,(int)(10.0*PP)));}
};

#endif
