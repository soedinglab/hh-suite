#ifndef HHVITERBIRUNNER
#define HHVITERBIRUNNER

#include "hhhmmsimd.h"
#include "hhviterbimatrix.h"
#include "hhviterbi.h"
#include "hhfunc.h"
#include <vector>
#include <map>

class ViterbiConsumerThread
{
	int thread_id;
    Viterbi * viterbiAlgo;
	HMMSimd* q_simd;
	HMMSimd* t_hmm_simd;
	Hit* hit_cur;     
    
	ViterbiMatrix* viterbiMatrix;
    int job_size;

    
  public:

    ViterbiConsumerThread(int pthread_id, Viterbi* vit, HMMSimd* q_simd, HMMSimd* t_hmm_simd, ViterbiMatrix* pviterbiMatrix) :
    	thread_id(pthread_id),
    	viterbiAlgo(vit),
    	q_simd(q_simd),
    	t_hmm_simd(t_hmm_simd),
    	viterbiMatrix(pviterbiMatrix),
    	job_size(0) {
                hit_cur = new Hit();
		}

    ~ViterbiConsumerThread(){
        std::cout << "Destruct ViterbiConsumerThread" << std::endl;
        delete viterbiAlgo;
        delete hit_cur;
    }

    std::vector<Hit> hits;
    std::vector<std::pair<char *,Viterbi::BacktraceResult> > excludeAlignments;
 	
 	void clear();
    void align(int maxres);
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Wrapper Viterbi call
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
class ViterbiRunner {
public:
    ViterbiRunner(ViterbiMatrix ** viterbiMatrix, std::vector<HHblitsDatabase*> &databases, int threads)
    : viterbiMatrix(viterbiMatrix), databases(databases), thread_count(threads) { }

	std::vector<Hit> alignment(Parameters& par, HMMSimd * q_simd, std::vector<HHDatabaseEntry*>  dbfiles, float* pb, const float S[20][20], const float Sim[20][20], const float R[20][20]);
    
private:
    ViterbiMatrix** viterbiMatrix;
    std::vector<HHblitsDatabase* > databases;
    int thread_count;

	void merge_thread_results(std::vector<Hit> &all_hits,
							  std::vector<HHDatabaseEntry*> &dbfiles_to_align,
							  std::map<std::string ,std::vector<Viterbi::BacktraceResult > >  &excludeAlignments,
							  std::vector<ViterbiConsumerThread *> &threads,
							  int alignment);


	void exclude_alignments(int maxResElem, HMMSimd* q_simd, HMMSimd* t_hmm_simd,
                            std::map<std::string ,std::vector<Viterbi::BacktraceResult > >  &excludeAlignments,
							ViterbiMatrix* viterbiMatrix);
};

#endif
