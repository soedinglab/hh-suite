#ifndef HHVITERBIRUNNER
#define HHVITERBIRUNNER

#include "hhhmmsimd.h"
#include "hhviterbimatrix.h"
#include "hhviterbi.h"
#include "hhhitlist.h"
#include "hhfunc.h"
#include <vector>
#include <map>

class ViterbiConsumerThread
{
	int thread_id;
    Viterbi * viterbiAlgo;
    Parameters par;
	HMMSimd* q_simd;
	HMMSimd* t_hmm_simd;
	ViterbiMatrix* viterbiMatrix;
    int job_size;

    
  public:

    ViterbiConsumerThread(int pthread_id, Parameters& par, HMMSimd* q_simd, HMMSimd* t_hmm_simd, ViterbiMatrix* pviterbiMatrix) :
    	thread_id(pthread_id),
    	par(par),
    	q_simd(q_simd),
    	t_hmm_simd(t_hmm_simd),
    	viterbiMatrix(pviterbiMatrix),
    	job_size(0) {
                viterbiAlgo = new Viterbi(par.maxres, par.loc, par.egq, par.egt,
                    par.corr, par.min_overlap, par.shift);
		}

    ~ViterbiConsumerThread(){
        delete viterbiAlgo;
    }

    std::vector<Hit> hits;
    std::vector<std::pair<char *,Viterbi::BacktraceResult> > excludeAlignments;
 	
 	void clear();
    void align(int maxres, int nseqdis);
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Wrapper Viterbi call
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
class ViterbiRunner {
public:
    ViterbiRunner(ViterbiMatrix ** viterbiMatrix, std::vector<HHblitsDatabase*> &databases, int threads)
    : viterbiMatrix(viterbiMatrix), databases(databases), thread_count(threads) { }

	std::vector<Hit> alignment(Parameters& par, HMMSimd * q_simd, std::vector<HHEntry*> dbfiles, const float qsc, float* pb, const float S[20][20], const float Sim[20][20], const float R[20][20]);
    
private:
    ViterbiMatrix** viterbiMatrix;
    std::vector<HHblitsDatabase* > databases;
    int thread_count;

	void merge_thread_results(std::vector<Hit> &all_hits,
							  std::vector<HHEntry*> &dbfiles_to_align,
							  std::map<std::string ,std::vector<Viterbi::BacktraceResult > >  &excludeAlignments,
							  std::vector<ViterbiConsumerThread *> &threads,
							  int alignment);


	void exclude_alignments(int maxResElem, HMMSimd* q_simd, HMMSimd* t_hmm_simd,
                            std::map<std::string ,std::vector<Viterbi::BacktraceResult > >  &excludeAlignments,
							ViterbiMatrix* viterbiMatrix);

	void exclude_regions(char* exclstr, int maxResElem, HMMSimd* q_hmm_simd, HMMSimd* t_hmm_simd, ViterbiMatrix* viterbiMatrix);
    
  float calculateEarlyStop(Parameters& par, HMM * q, std::vector<Hit> &all_hits, unsigned int startPos);

};

#endif
