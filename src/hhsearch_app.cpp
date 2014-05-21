/*
 * hhsearch_app.cpp
 *
 *  Created on: May 21, 2014
 *      Author: meiermark
 */

#include "hhsearch.h"

int main(int argc, char **argv) {
  Parameters par;
  HHsearch::ProcessAllArguments(argc, argv, par);

  HHblits hhsearch(par);

  FILE* inf;
  if(strcmp(par.infile, "stdin") == 0) {
      inf = stdin;
  }
  else {
      inf = fopen(par.infile, "r");
  }

  if(!inf) {
      std::cerr << "Input file (" << par.infile << ") could not be opened!" << std::endl;
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

