// hhhitlist.h

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
  void PrintHitList(HMM* q, char* outfile);

  // Print alignments of query sequences against hit sequences 
  void PrintAlignments(HMM* q, char* outfile, char outformat=0);

  void PrintMatrices(HMM* q, std::string outfile);

  // Return a figure of merit for distinction of the score with positive from the scores with negatives
  void Optimize(HMM* q);

  // Print score distribution into file score_dist
  void PrintScoreFile(HMM* q);
  
  // Write alignments in tabular output
  void WriteToAlifile(HMM* q, bool scop_only=false);

  // Log likelihood for fitting the EVD to the score distribution
  double LogLikelihoodEVD(double* v);
  // Static wrapper-function for calling the nonstatic member function LogLikelihoodEVD()
  static double LogLikelihoodEVD_static(void* pt2hitlist, double* v);


  // Subroutine to FindMin: new point given by highest point ihigh, fac and replaces ihigh if it is lower 
  double TryPoint(int ndim, double* p, double* y, double* psum, int ihigh, double fac, double (*Func)(void* pt2hitlist, double* v));
  
  // Find minimum with simplex method of Nelder and Mead (1965)
  float FindMin(int ndim, double* p, double* y, double tol, int& nfunc, double (*Func)(void* pt2hitlist, double* v));
  
  // Do a maximum likelihood fit of the scores with an EV distribution with parameters lamda and mu 
  void MaxLikelihoodEVD(HMM* q, int nbest);
  
  // Calculate HHblits composite E-values 
  void CalculateHHblitsEvalues(HMM* q);

  // Calculate Pvalues as a function of query and template lengths and diversities
  void CalculatePvalues(HMM* q);

  // Set P-values, E-values and scores according to q->lamda and q->mu (if calibration from database scan is impossible)  
  void GetPvalsFromCalibration(HMM* q);

  // HHblits: read PSI-BLAST E-values to determine correlation
  void ReadBlastFile(HMM* q);

  // Print first 20 hits of hitlist
  void Debug() 
    {
      Hit hit;
      int i=0;
      Reset();
      printf("TARGET     FAMILY     LEN COL LOG-PVA  S-AASS PROBAB SCORE_SORT\n");
      while (++i<=20 && !End()) 
	{
	  hit = ReadNext();
	  printf("%-10.10s %-10.10s %3i %3i %7.4f %7.2f %6.2f %6.2f\n",hit.name,hit.fam,hit.L,hit.matched_cols,-1.443*hit.logPval,-hit.score_aass,hit.Probab,hit.score_sort);
	}      
    }

  // Used only for transitive scoring??
  /* // Matrix manipulation */
  /* void PrintMatrix(float** V, int N); */
  /* void PrintMatrix(double** V, int N); */
  /* float NormalizationFactor(double** Csub,float* w, int M); */
  /* void Normalize(float* Ztq, char** fold, Hash<int>& excluded); */
  /* void InvertMatrix(double** B, double** A, int N); */
  /* void TransposeMatrix(double** V, int N); */
  /* void SVD(double **A, int n, double w[], double **V); */

};

