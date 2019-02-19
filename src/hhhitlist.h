// hhhitlist.h

#include <iostream>   // cin, cout, cerr
#include <fstream>    // ofstream, ifstream
#include <stdio.h>    // printf
#include <stdlib.h>   // exit
#include <string>     // strcmp, strstr
#include <math.h>     // sqrt, pow
#include <limits.h>   // INT_MIN
#include <float.h>    // FLT_MIN
#include <time.h>     // clock
#include <ctype.h>    // islower, isdigit etc
#include <sstream>
#include <set>

#include "hhhitlist-inl.h"
#include "hhhit.h"
#include "list.h"
#include "hash.h"
#include "hhfullalignment.h"
#include "log.h"

#ifndef HHHITLIST_H
#define HHHITLIST_H

class SearchCounter {
  public:
    SearchCounter();
    ~SearchCounter();

    int getCounter();
    void append(std::string id);

  private:
    std::set<std::string> already_seen;
};


/////////////////////////////////////////////////////////////////////////////////////
// HitList is a list of hits of type Hit which can be operated upon by several anaylsis methods 
/////////////////////////////////////////////////////////////////////////////////////
class HitList : public List<Hit>
{
private:
  double score[MAXPROF];        // HHsearch score of each HMM for ML fit
  double weight[MAXPROF];       // weight of each HMM = 1/(size_fam[tfam]*size_sfam[hit.sfam]) for ML fit
  int Nprof;                    // Number of HMMs for ML fit

public:
  int fams;                     // number of families found found in hitlist
  int sfams;                    // number of superfamilies found in hitlist
  int N_searched;               // number of sequences searched from HMM database
  Hash<float>* blast_logPvals;  // Hash containing names and log(P-values) read from BLAST file (needed for HHblits)

  HitList() {blast_logPvals=NULL;}
  ~HitList() {if (blast_logPvals) delete blast_logPvals;}

  // Print summary listing of hits
  void PrintHitList(HMM* q, std::stringstream& out, const unsigned int maxdbstrlen, const int z, const int Z, const float p, const double E, const int argc, char** argv);
  void PrintHitList(HMM* q, char* outfile, const unsigned int maxdbstrlen, const int z, const int Z, const float p, const double E, const int argc, char** argv);

  // Print alignments of query sequences against hit sequences 
  void PrintAlignments(HMM* q, char* outfile, const char showconf, const char showcons,
			const char showdssp, const char showpred, const float p, const int aliwidth, const int nseqdis,
			const int b, const int B, const double E, const float S[20][20], const int maxseq, char outformat);
  void PrintAlignments(HMM* q, std::stringstream& out, const char showconf, const char showcons,
			const char showdssp, const char showpred, const float p, const int aliwidth, const int nseqdis,
			const int b, const int B, const double E, const float S[20][20], const int maxseq, char outformat);

  void PrintHHR(HMM* q, char* outfile, const unsigned int maxdbstrlen,
			const char showconf, const char showcons, const char showdssp, const char showpred,
			const int b, const int B, const int z, const int Z, const int aliwidth, const int nseqdis,
			const float p, const double E, const int argc, char** argv, const float S[20][20], const int maxseq);
  void PrintHHR(HMM* q, std::stringstream& out, const unsigned int maxdbstrlen,
			const char showconf, const char showcons, const char showdssp, const char showpred,
			const int b, const int B, const int z, const int Z, const int aliwidth, const int nseqdis,
			const float p, const double E, const int argc, char** argv, const float S[20][20], const int maxseq);

  // Print score distribution into file score_dist
  void PrintScoreFile(HMM* q, char* outputfile);
  void PrintScoreFile(HMM* q, std::stringstream& outputstream);
    
  void PrintM8File(HMM* q, char* outputfile);
  void PrintM8File(HMM* q, std::stringstream& outputstream);

  void PrintMatrices(HMM* q, const char* matricesOutputFileName, const bool filter_matrices, const size_t max_number_matrices, const float S[20][20]);
  void PrintMatrices(HMM* q, std::stringstream& out, const bool filter_matrices, const size_t max_number_matrices, const float S[20][20]);
  
  // Write alignments in tabular output
  void WriteToAlifile(HMM* q, char* alitabfile,
      const int b, const int B, const int z, const int Z,
      const float p, const double E);

  void WriteToAlifile(HMM* q, std::stringstream& out,
      const int b, const int B, const int z, const int Z,
      const float p, const double E);

  // Calculate HHblits composite E-values 
  void CalculateHHblitsEvalues(HMM* q, const int dbsize,
		  const float alphaa, const float alphab, const float alphac, const double prefilter_evalue_thresh);

  // Calculate Pvalues as a function of query and template lengths and diversities
  void CalculatePvalues(HMM* q, const char loc, const char ssm, const float ssw);
};

#endif
