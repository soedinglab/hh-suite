// hhhit.C

#include "hhhit.h"

using std::cout;
using std::cerr;
using std::endl;
using std::ios;
using std::ifstream;
using std::ofstream;

/////////////////////////////////////////////////////////////////////////////////////
//// Constructor
/////////////////////////////////////////////////////////////////////////////////////
Hit::Hit() {
  longname = name = file = NULL;
  entry = NULL;
  sname = NULL;
  seq = NULL;
  self = 0;
  i = j = NULL;
  alt_i = alt_j = NULL;
  states = NULL;
  S = S_ss = P_posterior = NULL;
  sum_of_probs = 0.0;
  Neff_HMM = 0.0;
  realign_around_viterbi = false;
  //TODO: debug
  ssm1 = 0;
  ssm2 = 0;
  forward_profile = NULL;
  backward_profile = NULL;
}

/////////////////////////////////////////////////////////////////////////////////////
//// Free all allocated memory (to delete list of hits) 
//// Do NOT delete DP matrices
/////////////////////////////////////////////////////////////////////////////////////
void Hit::Delete() {
  if (i)
    delete[] i;
  if (j)
    delete[] j;

//  if (irep == 1) {
    if (alt_i) {
      delete alt_i;
      alt_i = NULL;
    }
    if (alt_j) {
      delete alt_j;
      alt_j = NULL;
    }
//  }
  
  if (states)
    delete[] states;
  if (S)
    delete[] S;
  if (S_ss)
    delete[] S_ss;
  if (P_posterior)
    delete[] P_posterior;
  //  delete[] l;    
  i = j = NULL;
  states = NULL;
  S = S_ss = P_posterior = NULL;

  delete[] longname;
  delete[] name;
  delete[] file;

  if (sname) {
    for (int k = 0; k < n_display; ++k)
      delete[] sname[k];
    delete[] sname;
  }
  sname = NULL;

  if (seq) {
    for (int k = 0; k < n_display; ++k)
      delete[] seq[k];
    delete[] seq;
  }
  seq = NULL;

  if(backward_profile) {
	  delete[] backward_profile;
  }

  if(forward_profile) {
	  delete[] forward_profile;
  }

  for(size_t i = 0; i < posterior_probabilities.size(); i++) {
	  delete posterior_probabilities[i];
  }
  posterior_probabilities.clear();

  longname = name = file = NULL;
}

void Hit::initHitFromHit(Hit& hit, int query_length) {
  Delete();

  this->longname = new char[strlen(hit.longname)+1];
  strcpy(this->longname, hit.longname);

  this->name = new char[strlen(hit.name)+1];
  strcpy(this->name, hit.name);

  this->file = new char[strlen(hit.file)+1];
  strcpy(this->file, hit.file);

  strcpy(fam, hit.fam);
  strcpy(sfam, hit.sfam);
  strcpy(fold, hit.fold);
  strcpy(cl, hit.cl);

  entry = hit.entry;

  if(hit.forward_profile) {
    forward_profile = new float[query_length + 1];
    for(int i = 1; i <= query_length; i++) {
      forward_profile[i] = hit.forward_profile[i];
    }
  }

  if(hit.backward_profile) {
    backward_profile = new float[query_length + 1];
    for(int i = 1; i <= query_length; i++) {
      backward_profile[i] = hit.backward_profile[i];
    }
  }

  for(size_t i = 0; i < hit.posterior_probabilities.size(); i++) {
    Posterior_Triple* triple = new Posterior_Triple(hit.posterior_probabilities[i]->query_pos,
                                                    hit.posterior_probabilities[i]->template_pos,
                                                    hit.posterior_probabilities[i]->posterior_probability);

    posterior_probabilities.push_back(triple);
  }

  predicted_alignment_quality = hit.predicted_alignment_quality;

  score = hit.score;
  score_sort = hit.score;
  score_aass = hit.score_aass;
  score_ss = hit.score_ss;
  Pval = hit.Pval;
  Pvalt = hit.Pvalt;
  logPval = hit.Pvalt;
  logPvalt = hit.Pvalt;
  Eval = hit.Eval;
  logEval = hit.Eval;
  Probab = hit.Probab;
  Pforward = hit.Pforward;

  L = hit.L;
  irep = hit.irep;
  lastrep = hit.lastrep;

  n_display = hit.n_display;

  sname = new char*[hit.n_display];
  seq   = new char*[hit.n_display];

  for (int k=0; k < hit.n_display; k++) {
    sname[k] = new char[strlen(hit.sname[k])+1];
    seq[k]   = new char[strlen(hit.seq[k])+1];
    strcpy(sname[k], hit.sname[k]);
    strcpy(seq[k], hit.seq[k]);
  }

  nss_dssp = hit.nss_dssp;
  nsa_dssp = hit.nsa_dssp;
  nss_pred = hit.nss_pred;
  nss_conf = hit.nss_conf;
  nfirst = hit.nfirst;
  ncons = hit.ncons;

  nsteps = hit.nsteps;
  i = new int[hit.nsteps];
  j = new int[hit.nsteps];
  states = new char[hit.nsteps];
  S = new float[hit.nsteps];
  S_ss = new float[hit.nsteps];
  P_posterior = new float[hit.nsteps];

  for(int index = 0; index < hit.nsteps; index++) {
    i[index] = hit.i[index];
    j[index] = hit.j[index];
    states[index] = hit.states[index];
    S[index] = hit.S[index];
    S_ss[index] = hit.S_ss[index];
    if(hit.P_posterior)
      P_posterior[index] = hit.P_posterior[index];
  }

  if(hit.alt_i) {
    alt_i = new std::vector<int>();
    for(size_t i = 0; i < hit.alt_i->size(); i++) {
      alt_i->push_back(hit.alt_i->at(i));
    }
  }
  if(hit.alt_j) {
    alt_j = new std::vector<int>();
    for(size_t i = 0; i < hit.alt_j->size(); i++) {
      alt_j->push_back(hit.alt_j->at(i));
    }
  }

  i1 = hit.i1;
  i2 = hit.i2;
  j1 = hit.j1;
  j2 = hit.j2;
  matched_cols = hit.matched_cols;
  ssm1 = hit.ssm1;
  ssm2 = hit.ssm2;
  self = hit.self;
  sum_of_probs = hit.sum_of_probs;
  Neff_HMM = hit.Neff_HMM;

  realign_around_viterbi = hit.realign_around_viterbi;

  qsc = hit.qsc;

  state = hit.state;
  min_overlap = hit.min_overlap;
}

/////////////////////////////////////////////////////////////////////////////////////
//// Allocate/delete memory for indices by given alignment
/////////////////////////////////////////////////////////////////////////////////////
void Hit::AllocateIndices(int len) {
  i = new int[len];
  j = new int[len];
}

void Hit::DeleteIndices() {
  delete[] i;
  delete[] j;
}

float Hit::calculateSimilarity(HMM* q, const float S[20][20]) {
	float alignment_similarity = 0.0;

	char template_mapping[MAXRES];
	char query_mapping[MAXRES];

	char c;
	int l = 1;
	int i = 1;
	while ((c = seq[nfirst][l++])) {
		if (c != '.' && !(c >= 'a' && c <= 'z'))
			template_mapping[i++] = c;
	}

	l = 1;
	i = 1;
	while ((c = q->seq[q->nfirst][l++])) {
		if (c != '.' && !(c >= 'a' && c <= 'z'))
			query_mapping[i++] = c;
	}

	for (int step = nsteps; step >= 1; step--) {
		char state = states[step];
		if (state == MM) {
		  char qc = query_mapping[this->i[step]];
		  char tc = template_mapping[j[step]];
		  alignment_similarity += S[aa2i(qc)][aa2i(tc)];
		}
	}

	alignment_similarity /= matched_cols;

	return alignment_similarity;
}



/////////////////////////////////////////////////////////////////////////////////////
// Allocate memory for data of new alignment (sequence names, alignment, scores,...)
/////////////////////////////////////////////////////////////////////////////////////
//void Hit::InitializeBacktrace(HMM* q, HMM* t) {
//  //Copy information about template profile to hit and reset template pointers to avoid destruction
//  longname = new char[strlen(t->longname) + 1];
//  name = new char[strlen(t->name) + 1];
//  file = new char[strlen(t->file) + 1];
//  if (!file)
//    MemoryError(
//        "space for alignments with database HMMs. \nNote that all alignments have to be kept in memory", __FILE__, __LINE__, __func__);
//  strcpy(longname, t->longname);
//  strcpy(name, t->name);
//  strcpy(fam, t->fam);
//  strcpy(sfam, t->sfam);
//  strcpy(fold, t->fold);
//  strcpy(cl, t->cl);
//  strcpy(file, t->file);
//
//  // Allocate new space
//  this->i = new int[i2 + j2 + 2];
//  this->j = new int[i2 + j2 + 2];
//  states = new char[i2 + j2 + 2];
//  S = S_ss = P_posterior = NULL;
//
//  sname = new char*[t->n_display];
//  seq = new char*[t->n_display];
//  if (!sname || !seq)
//    MemoryError(
//        "space for alignments with database HMMs.\nNote that all sequences for display have to be kept in memory", __FILE__, __LINE__, __func__);
//
//  if (irep == 1) {
//    // Make flat copy for first alignment of template seqs to save speed
//    for (int k = 0; k < t->n_display; k++) {
//      sname[k] = t->sname[k];
//      seq[k] = t->seq[k];
//    }
//    t->dont_delete_seqs = 1;
//  }
//  else {
//    // Make deep copy for all further alignments
//    for (int k = 0; k < t->n_display; k++) {
//      sname[k] = new char[strlen(t->sname[k]) + 1];
//      seq[k] = new char[strlen(t->seq[k]) + 1];
//      strcpy(sname[k], t->sname[k]);
//      strcpy(seq[k], t->seq[k]);
//    }
//  }
//
//  n_display = t->n_display;
//  ncons = t->ncons;
//  nfirst = t->nfirst;
//  nss_dssp = t->nss_dssp;
//  nsa_dssp = t->nsa_dssp;
//  nss_pred = t->nss_pred;
//  nss_conf = t->nss_conf;
//  L = t->L;
//  Neff_HMM = t->Neff_HMM;
//  Eval = 1.0;
//  logEval = 0.0;
//  Pval = 1.0;
//  Pvalt = 1.0;
//  logPval = 0.0;
//  logPvalt = 0.0;
//  Probab = 1.0;
//}

void Hit::initHitFromHMM(HMM * t, const int par_nseqdis){
    this->alt_i = NULL;
    this->alt_j = NULL;
    this->P_posterior = NULL;
    //Copy information about template profile to hit and reset template pointers to avoid destruction
    this->longname=new char[strlen(t->longname)+1];
    //    printf("%s\n",t->longname);
    this->name    =new char[strlen(t->name)+1];
    this->file    =new char[strlen(t->file)+1];
    if (!this->file) MemoryError("space for alignments with database HMMs. \nNote that all alignments have to be kept in memory", __FILE__, __LINE__, __func__);

    strcpy(this->longname,t->longname);

    strcpy(this->name,t->name);
    strcpy(this->fam ,t->fam);
    strcpy(this->sfam ,t->sfam);
    strcpy(this->fold ,t->fold);
    strcpy(this->cl ,t->cl);
    strcpy(this->file,t->file);

    this->sname = new char*[t->n_display];
    this->seq   = new char*[t->n_display];
    if (!this->sname || !this->seq)
        MemoryError("space for alignments with database HMMs.\nNote that all sequences for display have to be kept in memory", __FILE__, __LINE__, __func__);

    // Make deep copy for all further alignments
    for (int k=0; k < std::min(t->n_display, par_nseqdis); k++)
    {
        this->sname[k] = new char[strlen(t->sname[k])+1];
        this->seq[k]   = new char[strlen(t->seq[k])+1];
        strcpy(this->sname[k],t->sname[k]);
        strcpy(this->seq[k],t->seq[k]);
    }

    this->n_display = std::min(t->n_display, par_nseqdis);
    this->ncons  = t->ncons;
    this->nfirst = t->nfirst;
    this->nss_dssp = t->nss_dssp;
    this->nsa_dssp = t->nsa_dssp;
    this->nss_pred = t->nss_pred;
    this->nss_conf = t->nss_conf;
    this->L = t->L;
    this->Neff_HMM = t->Neff_HMM;
    this->Eval   = 1.0;
    this->logEval= 0.0;
    this->Pval   = 1.0;
    this->Pvalt  = 1.0;
    this->logPval = 0.0;
    this->logPvalt= 0.0;
    this->Probab = 1.0;
}

float Hit::estimateAlignmentQuality(HMM* q) {
	  const float biases[5] = {1.004681, 3.927253, -36.391424, -0.3767489, 0.0197167};
	  const float out_bias = -12.7733130;

	  const float out_weights[5] = {-0.2237661, 13.7567823, -1.0593247, 5.7120146, -5.8497739};
	  const float in_weights[10] = {21.419840, -43.118904, 20.178319, -2.598612, -1.512682, 9.677802, 6.2671908, -2.5100477, 5.7690630, -2.4830519};

	  float sum = out_bias;
	  for(int n = 0; n < 5; n++) {
	     float tmp = biases[n] + in_weights[n*2] * (sum_of_probs / q->L) + in_weights[n*2+1] * (score / q->L);
	     if(tmp < -15) {
	        tmp = 0.0;
	     }
	     else if(tmp > 15) {
	        tmp = 1.0;
	     }
	     else {
	        tmp = (1.0/(1.0+exp(-tmp)));
	     }
	     sum += out_weights[n] * tmp;
	  }

	  return sum;
}



int compareHitLengths(const void * a, const void * b) {

    Hit* pa = *(Hit**)a;
    Hit* pb = *(Hit**)b;

//  return (pb->L_template - pa->L_template);
    return (pb->L - pa->L);
}

