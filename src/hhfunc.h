/*
 * hhfunc.h
 *
 *  Created on: Apr 4, 2014
 *      Author: meiermark
 */

#ifndef HHFUNC_H_
#define HHFUNC_H_

#include <cstdio>
#include "hhhmm.h"
#include "hhhit.h"
#include "hhalignment.h"
#include "hhhitlist.h"

void ReadQueryFile(Parameters& par, FILE* inf, char& input_format, char use_global_weights, HMM* q, Alignment* qali, char infile[],
		float* pb, const float S[20][20], const float Sim[20][20]);

void ReadQueryFile(Parameters& par, char* infile, char& input_format, char use_global_weights, HMM* q, Alignment* qali,
		float* pb, const float S[20][20], const float Sim[20][20]);

// Add transition and amino acid pseudocounts to query HMM, calculate aa background etc.
void PrepareQueryHMM(Parameters& par, char& input_format, HMM* q, cs::Pseudocounts<cs::AA>* pc_hhm_context_engine, cs::Admix* pc_hhm_context_mode,
		const float* pb, const float R[20][20]);

// Do precalculations for q and t to prepare comparison
void PrepareTemplateHMM(Parameters& par, HMM* q, HMM* t, int format, float linear_tranistion_probs, const float* pb, const float R[20][20]);

// Read number of sequences in annotation, after second '|'
int SequencesInCluster(char* name);

void InitializePseudocountsEngine(Parameters& par,
    cs::ContextLibrary<cs::AA>*& context_lib, cs::Crf<cs::AA>*& crf,
    cs::Pseudocounts<cs::AA>*& pc_hhm_context_engine, cs::Admix*& pc_hhm_context_mode,
    cs::Pseudocounts<cs::AA>*& pc_prefilter_context_engine, cs::Admix*& pc_prefilter_context_mode);

void DeletePseudocountsEngine(
    cs::ContextLibrary<cs::AA>* context_lib, cs::Crf<cs::AA>* crf,
    cs::Pseudocounts<cs::AA>* pc_hhm_context_engine, cs::Admix* pc_hhm_context_mode,
    cs::Pseudocounts<cs::AA>* pc_prefilter_context_engine, cs::Admix* pc_prefilter_context_mode);

#endif /* HHFUNC_H_ */
