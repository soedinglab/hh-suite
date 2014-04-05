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

void ReadQueryFile(FILE* inf, char& input_format, char use_global_weights, HMM* q, Alignment* qali, char infile[],
		float* pb, const float S[20][20], const float Sim[20][20]);

void ReadQueryFile(char* infile, char& input_format, char use_global_weights, HMM* q, Alignment* qali,
		float* pb, const float S[20][20], const float Sim[20][20]);

// Add transition and amino acid pseudocounts to query HMM, calculate aa background etc.
void PrepareQueryHMM(char& input_format, HMM* q, cs::Pseudocounts<cs::AA>* pc_hhm_context_engine, cs::Admix* pc_hhm_context_mode,
		const float* pb, const float R[20][20]);

// Do precalculations for q and t to prepare comparison
void PrepareTemplateHMM(HMM* q, HMM* t, int format, const float* pb, const float R[20][20]);

// Calculate secondary structure prediction with PSIPRED
void CalculateSS(char *ss_pred, char *ss_conf, char *tmpfile);

// Calculate secondary structure for given HMM and return prediction
void CalculateSS(HMM* q, char *ss_pred, char *ss_conf, const float* pb);

// Calculate secondary structure for given HMM
void CalculateSS(HMM* q, const float* pb);

// Write alignment in tab format (option -atab)
void WriteToAlifile(FILE* alitabf, Hit* hit);

// Read number of sequences in annotation, after second '|'
int SequencesInCluster(char* name);

#endif /* HHFUNC_H_ */
