/*
 * PosteriorDecoderRunnerInputData.h
 *
 *  Created on: Apr 6, 2014
 *      Author: stefan
 */

#ifndef POSTERIORDECODERRUNNERINPUTDATA_H_
#define POSTERIORDECODERRUNNERINPUTDATA_H_

#include <map>
#include <vector>
#include "hhdatabase.h"

// Structure used as data container for better transparency at
// 	initialization of PosteriorDecoderRunner
class PosteriorDecoderRunnerInputData {

public:
	PosteriorDecoderRunnerInputData(std::vector<HHblitsDatabase*> & databases,
			std::vector<HHDatabaseEntry*> & dbfiles,
			std::map<short int, std::vector<Hit *> > & alignments,
			const int n_t_hmms_to_align,
			const int t_maxres) :
            databases(databases),
            dbfiles(dbfiles),
            alignments(alignments),
            n_t_hmms_to_align(n_t_hmms_to_align),
            t_maxres(t_maxres) {};

	std::vector<HHblitsDatabase*> & databases;	// databases used to read HMMs from
	std::vector<HHDatabaseEntry* > & dbfiles;					// database files
	// This map includes structured hit objects:
	// 	outer map: key = irep; inner map: key = dbfile
	// 	This is to make sure that no hits from the same template HMM are processed in one vector!
	//	--> always first hits with irep == 1, ...
	std::map<short int, std::vector<Hit *> > & alignments;

	const int n_t_hmms_to_align;	// Number of alignments (including all sub-optimal alignments)
	const int t_maxres;						// maximum template resolution (t.L + 2)


//	char** dbfiles_new;	// database files
//	int ndb_new;				// number of db-files to align

};

#endif /* POSTERIORDECODERRUNNERINPUTDATA_H_ */
