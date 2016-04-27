/*
 * hhalign_app.cpp
 *
 *  Created on: Jun 24, 2014
 *      Author: meiermark
 */


#include "hhalign.h"

int main(int argc, char **argv) {
  Parameters par;
  HHalign::ProcessAllArguments(argc, argv, par);

  //is empty and will stay empty in this application
  std::vector<HHblitsDatabase*> databases;
#ifdef OPENMP
  omp_set_num_threads(par.threads);
#endif
  HHalign hhalign(par, databases);

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

  hhalign.run(inf, par.infile, par.tfile);
  fclose(inf);

  if(Log::reporting_level() >= INFO) {
    hhalign.printHHRFile();
  }

  hhalign.writeHHRFile(par.outfile);
  hhalign.writeScoresFile(par.scorefile);
  hhalign.writeM8(par.m8file);
  hhalign.writePairwiseAlisFile(par.pairwisealisfile, par.outformat);
  hhalign.writeAlitabFile(par.alitabfile);
  hhalign.writePsiFile(par.psifile);
  hhalign.writeHMMFile(par.hhmfile);
  hhalign.writeA3MFile(par.alnfile);
}



