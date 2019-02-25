#include "hhsearch.h"
#include "hhalign.h"

void checkOutput(Parameters& par) {
  if (!*par.outfile) {
    RemoveExtension(par.outfile, par.infile);
    strcat(par.outfile, ".hhr");
    HH_LOG(INFO) << "Search results will be written to " << par.outfile << "\n";
  }
}

int main(int argc, const char **argv) {
  Parameters par(argc, argv);

#ifdef HHSEARCH
  HHsearch::ProcessAllArguments(par);
#elif HHALIGN
  HHalign::ProcessAllArguments(par);
#else
  HHblits::ProcessAllArguments(par);
#endif
  checkOutput(par);

  std::vector<HHblitsDatabase*> databases;
#ifdef HHSEARCH
  HHsearch::prepareDatabases(par, databases);
#elif HHALIGN
#else
  HHblits::prepareDatabases(par, databases);
#endif

#ifdef OPENMP
  omp_set_num_threads(par.threads);
#endif

#ifdef HHSEARCH
  HHblits app(par, databases);
#elif HHALIGN
  HHalign app(par);
#else
  HHblits app(par, databases);
#endif

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

  app.run(inf, par.infile);
  fclose(inf);

  if(Log::reporting_level() >= INFO) {
    app.printHitList();
  }

  app.writeHHRFile(par.outfile);
  app.writeAlisFile(par.alisbasename);
  app.writeScoresFile(par.scorefile);
  app.writeM8(par.m8file);
  app.writePairwiseAlisFile(par.pairwisealisfile, par.outformat);
  app.writeAlitabFile(par.alitabfile);
  app.writePsiFile(par.psifile);
  app.writeHMMFile(par.hhmfile);
  app.writeA3MFile(par.alnfile);
  app.writeMatricesFile(par.matrices_output_file);

  for(size_t i = 0; i < databases.size(); i++) {
    delete databases[i];
  }
  databases.clear();
}
