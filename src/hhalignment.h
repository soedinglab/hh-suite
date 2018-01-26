// hhalignment.h

#ifndef HHALIGNMENT_H_
#define HHALIGNMENT_H_

#include <iostream>  // cin, cout, cerr
#include <fstream>   // ofstream, ifstream
#include <cstdio>    // printf
#include <cstdlib>   // exit
#include <string>    // strcmp, strstr
#include <cmath>     // sqrt, pow
#include <climits>   // INT_MIN
#include <cfloat>    // FLT_MIN
#include <ctime>     // clock
#include <cctype>    // islower, isdigit etc
#include <cstring>

#include "hhdecl.h"
#include "hhhmm.h"
#include "hhhit.h"
#include "util.h"
#include "hhutil.h"
#include "simd.h"

#include "log.h"

extern "C" {
#include <ffindex.h>     // fast index-based database reading
}

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
  Alignment& operator=(Alignment&);

  // Read alignment into X (uncompressed) in ASCII characters
  void Read(FILE* inf, char infile[], const char mark, const int maxcol, const int nseqdis, char* firstline=NULL);
  void ReadCompressed(ffindex_entry_t* entry, char* data,
      ffindex_index_t* ffindex_sequence_database_index, char* ffindex_sequence_database_data,
      ffindex_index_t* ffindex_header_database_index, char* ffindex_header_database_data,
      const char mark, const int maxcol);

  // Read sequences from HHM-file into X (uncompressed) in ASCII characters
  void GetSeqsFromHMM(HMM* q);
  
  // Convert ASCII to numbers between 0 and 20, throw out all insert states, 
  // record their number in I[k][i] and store sequences to be displayed in sname[k] and seq[k]
  void Compress(const char infile[NAMELEN], const char cons, const int maxres, const int maxcol, const int par_M, const int Mgaps);

  // Apply sequence identity filter
  int FilterForDisplay(int max_seqid, const char mark, const float S[20][20], int coverage=0, int qid=0, float qsc=0, int N=0);
  int Filter(int max_seqid, const float S[20][20], int coverage=0, int qid=0, float qsc=0, int N=0);
  int Filter2(char keep[], int coverage, int qid, float qsc, int seqid1, int seqid2, int Ndiff, const float S[20][20]);
  void Shrink();


  void FilterNeff(char use_global_weights, const char mark, const char cons,
			const char showcons, const int maxres, const int max_seqid, const int coverage,
			const float Neff, const float* pb, const float S[20][20], const float Sim[20][20]);
  float filter_by_qsc(float qsc, char use_global_weights, const char mark, const char cons,
			const char showcons, const int maxres, const int max_seqid, const int coverage,
			char* keep_orig, const float* pb, const float S[20][20], const float Sim[20][20]);

  // Calculate AA frequencies q.p[i][a] and transition probabilities q.tr[i][a] from alignment
  void FrequenciesAndTransitions(HMM* q, char use_global_weights,
			const char mark, const char cons, const char showcons, const int maxres,
			const float* pb, const float Sim[20][20], char* in=NULL, bool time=false);

  // Calculate freqs q.f[i][a] and transitions q.tr[i][a] (a=MM,MI,MD) with pos-specific subalignments
  void Amino_acid_frequencies_and_transitions_from_M_state(HMM* q,
			char use_global_weights, char* in, const int maxres, const float* pb);

  // Calculate transitions q.tr[i][a] (a=DM,DD) with pos-specific subalignments
  void Transitions_from_D_state(HMM* q, char* in, const int maxres);

  // Calculate transitions q.tr[i][a] (a=DM,DD) with pos-specific subalignments
  void Transitions_from_I_state(HMM* q, char* in, const int maxres);
  
  // Write alignment without insert states to alignment file
  void WriteWithoutInsertsToFile(const char* alnfile, const char append);

  // Write alignment to alignment file
  void WriteToFile(const char* alnfile, const char append, const char format[]=NULL);
  void WriteToFile(std::stringstream& out, const char format[]=NULL);

  // Read a3m slave alignment of hit from ta3mfile and merge into (query) master alignment
  void MergeMasterSlave(Hit& hit, Alignment& Tali, char* ta3mfile, const int par_maxcol);

  // Add a sequence to Qali
  void AddSequence(char Xk[], int Ik[]=NULL);

  // Add SS prediction to Qali
  void AddSSPrediction(char seq_pred[], char seq_conf[]);

  // Determine matrix of position-specific weights w[k][i] for multiple alignment
  void GetPositionSpecificWeights(float* w[], char use_global_weights);

  // Set keep[] and display[] arrays to 0 to mark seqs as non-printable
  void MarkSeqsAsNonPrintable();

  char readCommentLine;   // Set to 1, if a comment line with '#' is read

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
  int maxseq;
  int FilterWithCoreHMM(char in[], float coresc, HMM* qcore, const float* pb);
  char * initX(int len);
  
};

#endif
