/*
 * hhsearch.h
 *
 * Search for a multiple alignment (transformed into HMM) in a profile HMM database

 * Error codes: 0: ok  1: file format error  2: file access error  3: memory error  4: command line error  6: internal logic error  7: internal numeric error

 *     (C) Johannes Soeding and Michael Remmert 2012

 *     This program is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.

 *     This program is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.

 *     You should have received a copy of the GNU General Public License
 *     along with this program.  If not, see <http://www.gnu.org/licenses/>.

 *     We are very grateful for bug reports! Please contact us at soeding@mpibpc.mpg.de

 *     Reference:
 *     Remmert M., Biegert A., Hauser A., and Soding J.
 *     HHblits: Lightning-fast iterative protein sequence searching by HMM-HMM alignment.
 *     Nat. Methods 9:173-175 (2011); epub Dec 25, doi: 10.1038/NMETH.1818
 */

#ifndef HHSEARCH_H_
#define HHSEARCH_H_

#include "hhblits.h"

//const char HHSEARCH_REFERENCE[] =
//    "Soding, J. Protein homology detection by HMM-HMM comparison. Bioinformatics 21:951-960 (2005).\n";

class HHsearch {
public:
    static void prepareDatabases(Parameters& par, std::vector<HHblitsDatabase*>& databases);
    static void ProcessAllArguments(Parameters& par);
private:
    static void help(Parameters& par, char all = 0);
    static void ProcessArguments(Parameters& par);

};

#endif /* HHSEARCH_H_ */
