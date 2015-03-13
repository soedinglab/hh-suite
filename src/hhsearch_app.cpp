/*
 * hhsearch_app.cpp
 *
 *  Created on: May 21, 2014
 *      Author: meiermark
 */

#include "hhsearch.h"

#ifdef OPENMP
	#include <omp.h>
#endif

int main(int argc, char **argv) {
  Parameters par;
  HHsearch::ProcessAllArguments(argc, argv, par);

  std::vector<HHblitsDatabase*> databases;
  HHblits::prepareDatabases(par, databases);

#ifdef OPENMP
  omp_set_num_threads(par.threads);
#endif
  HHblits hhsearch(par, databases);

  FILE* inf;
  if(strcmp(par.infile, "stdin") == 0) {
      inf = stdin;
  }
  else {
      inf = fopen(par.infile, "r");
  }

  if(!inf) {
      HH_LOG(ERROR) << "Input file (" << par.infile << ") could not be opened!" << std::endl;
      exit(1);
  }

  hhsearch.run(inf, par.infile);
  fclose(inf);

  hhsearch.writeHHRFile(par.outfile);
  hhsearch.writeScoresFile(par.scorefile);
  hhsearch.writePairwiseAlisFile(par.pairwisealisfile, par.outformat);
  hhsearch.writeAlitabFile(par.alitabfile);
  hhsearch.writePsiFile(par.psifile);
  hhsearch.writeHMMFile(par.hhmfile);
  hhsearch.writeA3MFile(par.alnfile);
}

