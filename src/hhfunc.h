/*
 * hhfunc.h
 *
 *  Created on: Apr 4, 2014
 *      Author: meiermark
 */

#ifndef HHFUNC_H_
#define HHFUNC_H_

#include <cstdio>
#include <iostream>
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
void PrepareTemplateHMM(Parameters& par, HMM* q, HMM* t, int format, const float* pb, const float R[20][20]);

// Calculate secondary structure prediction with PSIPRED
void CalculateSS(char *ss_pred, char *ss_conf, char *tmpfile, const char* psipred_data, const char* psipred);

// Calculate secondary structure for given HMM and return prediction
void CalculateSS(HMM* q, char *ss_pred, char *ss_conf, const char* psipred_data, const char* psipred, const float* pb);

// Calculate secondary structure for given HMM
void CalculateSS(HMM* q, const int maxres, const char* psipred_data, const char* psipred, const float* pb);

// Write alignment in tab format (option -atab)
void WriteToAlifile(FILE* alitabf, Hit* hit, const char forward, const char realign);

// Read number of sequences in annotation, after second '|'
int SequencesInCluster(char* name);

void InitializePseudocountsEngine(Parameters& par,
    cs::ContextLibrary<cs::AA>* context_lib, cs::Crf<cs::AA>* crf,
    cs::Pseudocounts<cs::AA>* pc_hhm_context_engine, cs::Admix* pc_hhm_context_mode,
    cs::Pseudocounts<cs::AA>* pc_prefilter_context_engine, cs::Admix* pc_prefilter_context_mode);

void DeletePseudocountsEngine(
    cs::ContextLibrary<cs::AA>* context_lib, cs::Crf<cs::AA>* crf,
    cs::Pseudocounts<cs::AA>* pc_hhm_context_engine, cs::Admix* pc_hhm_context_mode,
    cs::Pseudocounts<cs::AA>* pc_prefilter_context_engine, cs::Admix* pc_prefilter_context_mode);

void AlignByWorker(Parameters& par, Hit* hit, HMM* t, HMM* q, const int format, const float* pb, const float R[20][20], const float S73[NDSSP][NSSPRED][MAXCF], const float S33[NSSPRED][MAXCF][NSSPRED][MAXCF], HitList& hitlist);

void PerformViterbiByWorker(Parameters& par, Hit* hit, HMM* t, HMM* q, const int format, const float* pb, const float R[20][20], const float S73[NDSSP][NSSPRED][MAXCF], const float S33[NSSPRED][MAXCF][NSSPRED][MAXCF], HitList& hitlist, Hash<Hit>* previous_hits);

void RealignByWorker(Parameters& par, Hit* hit, HMM* q, HMM* t, const int format, const float* pb, const float R[20][20], const float S73[NDSSP][NSSPRED][MAXCF], const float S33[NSSPRED][MAXCF][NSSPRED][MAXCF]);

#endif /* HHFUNC_H_ */
