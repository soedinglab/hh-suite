#ifndef CS_CS_H_
#define CS_CS_H_

// C includes
#include <assert.h>
#include <ctype.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

// STL includes
#include <algorithm>
#include <bitset>
#include <iostream>
#include <iomanip>
#include <set>
#include <string>
#include <utility>
#include <vector>

// Special header files for Open-MP and MPI
#ifdef OPENMP
#include <omp.h>
#endif
#ifdef PARALLEL
#include <mpi.h>
#endif

// Basic includes
#include "globals.h"
#include "aa.h"
#include "as.h"
#include "assert_helpers.h"
#include "bitset.h"
#include "dna.h"
#include "exception.h"
#include "io.h"
#include "log.h"
#include "matrix.h"
#include "ran.h"
#include "scoped_ptr.h"
#include "shared_ptr.h"
#include "timer.h"
#include "utils.h"
#include "vector.h"

#endif  // CS_CS_H_
