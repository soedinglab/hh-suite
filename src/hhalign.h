/*
 * hhalign.h
 *
 *  Created on: Jun 24, 2014
 *      Author: meiermark
 */

#ifndef HHALIGN_H_
#define HHALIGN_H_

#include <iostream>   // cin, cout, cerr
#include <fstream>    // ofstream, ifstream
#include <cstdio>     // printf
#include <algorithm>  // min,max
#include <stdlib.h>   // exit
#include <string.h>   // strcmp, strstr
#include <math.h>     // sqrt, pow
#include <limits.h>   // INT_MIN
#include <float.h>    // FLT_MIN
#include <ctype.h>    // islower, isdigit etc
#include <time.h>     // clock_gettime etc. (in realtime library (-lrt compiler option))
#include <errno.h>    // perror()
#include <cassert>
#include <stdexcept>
#include <map>

#include <sys/time.h>

using std::cout;
using std::cerr;
using std::endl;
using std::ios;
using std::ifstream;
using std::ofstream;

// context-specific pseudocounts
#include "cs.h"
#include "context_library.h"
#include "library_pseudocounts-inl.h"
#include "crf_pseudocounts-inl.h"

#include "util.h"        // imax, fmax, iround, iceil, ifloor, strint, strscn, strcut, substr, uprstr, uprchr, Basename etc.
#include "list.h"        // list data structure
#include "hash.h"        // hash data structure
#include "hhdecl.h"      // Constants, global variables, struct Parameters

#include "hhutil.h"      // MatchChr, InsertChr, aa2i, i2aa, log2, fast_log2, ScopID, WriteToScreen,
#include "hhmatrices.h"  // BLOSUM50, GONNET, HSDM
#include "hhhmm.h"       // class HMM
#include "hhhit.h"       // class Hit
#include "hhalignment.h" // class Alignment
#include "hhhalfalignment.h" // class HalfAlignment
#include "hhfullalignment.h" // class FullAlignment
#include "hhhitlist.h"   // class Hit
#include "hhfunc.h"      // some functions common to hh programs

#include "hhblits.h"

class HHalign : public HHblits {
  public:
    HHalign(Parameters& par, std::vector<HHblitsDatabase*>& databases);
    virtual ~HHalign();
    void run(FILE* query_fh, char* query_path, std::vector<std::string>& template_paths);
    static void ProcessAllArguments(int argc, char** argv, Parameters& par);

  private:
    static void help(Parameters& par, char all=0);
    static void ProcessArguments(int argc, char** argv, Parameters& par);
};

#endif /* HHALIGN_H_ */
