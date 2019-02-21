/*
     hhblits.cpp:
	 Iterative search for a multiple alignment in a profile HMM database

	 Error codes: 0: ok  1: file format error  2: file access error  3: memory error  4: command line error  6: internal logic error  7: internal numeric error

     This program is free software: you can redistribute it and/or modify
     it under the terms of the GNU General Public License as published by
     the Free Software Foundation, either version 3 of the License, or
     (at your option) any later version.

     This program is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU General Public License for more details.

     You should have received a copy of the GNU General Public License
     along with this program.  If not, see <http://www.gnu.org/licenses/>.

     We are very grateful for bug reports! Please contact us at soeding@mpibpc.mpg.de

     Reference:
     Remmert M., Biegert A., Hauser A., and Soding J.
     HHblits: Lightning-fast iterative protein sequence searching by HMM-HMM alignment.
     Nat. Methods 9:173-175 (2011); epub Dec 25, doi: 10.1038/NMETH.1818 (2011).

     This program contains, in file hhprefilter.C, code adapted from Michael Farrar
     (http://sites.google.com/site/farrarmichael/smith-waterman). His code is marked
     in the file hhprefilter.C.
     The copy right of his code is shown below:

     Copyright 2006, by Michael Farrar.  All rights reserved. The SWSSE2
     program and documentation may not be sold or incorporated into a
     commercial product, in whole or in part, without written consent of
     Michael Farrar.

     For further information regarding permission for use or reproduction,
     please contact Michael Farrar at:

         farrar.michael@gmail.com


     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
     EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
     MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
     IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
     CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
     TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
     SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

     Note by J. Soeding: Michael Farrar died unexpectedly in December 2010.
     Many thanks posthumously for your great code!
 */

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
