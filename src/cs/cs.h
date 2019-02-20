/*
  Copyright 2012 Andreas Biegert

  This file is part of the CS-BLAST package.

  The CS-BLAST package is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  The CS-BLAST package is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

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
#include "io.h"
#include "log.h"
#include "matrix.h"
#include "ran.h"
#include "scoped_ptr.h"
#include "utils.h"
#include "vector.h"

#endif  // CS_CS_H_
