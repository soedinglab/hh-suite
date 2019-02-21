# FFindex - A database wrapped around mmap

FFindex is a very simple index/database for huge amounts of small files. The
files are stored concatenated in one big data file, seperated by '\0'. A second
file contains a plain text index, giving name, offset and length of of the
small files. The lookup is currently done with a binary search on an array made
from the index file. The attatched binaries (see Usage below) and their source
code shall give an impression of how to use the functions supported by the library in C/C++ code.
 
## Copyright

FFindex was written by Andreas Hauser <hauser@genzentrum.lmu.de>.

Please add your name here if you distribute modified versions.
* Martin Steinegger <martin.steinegger@mpibpc.mpg.de>
* Markus Meier <markus.meier@mpibpc.mpg.de>

FFindex is provided under the Create Commons license "Attribution-ShareAlike 4.0",
which basically captures the spirit of the Gnu Public License (GPL).

See:
http://creativecommons.org/licenses/by-sa/4.0/

## Thanks

Thanks to Laszlo Kajan for creating and maintaining Debian packages
and many suggestions to improve the build and user experience.



## Installation

### Compilation
With the sourcecode ready, simply run cmake with the default settings and libraries should be auto-detected:

	mkdir build
	cd build
	cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo -G "Unix Makefiles" -DCMAKE_INSTALL_PREFIX=${INSTALL_BASE_DIR} ..
	make
	make install


**Please use a sensible value for ${INSTALL_BASE_DIR}, e.g. /usr/local or /opt/ffindex or $HOME/ffindex**


### Setting environment variables

Please note that before querying or unlinking entries a ffindex must be
sorted, although you can add to it without. So either specify -s with
ffindex_build or sorted later with ffindex_modify -s.

Setup environment:

	export PATH="${INSTALL_BASE_DIR}/bin:${PATH}"
	export LD_LIBRARY_PATH="${INSTALL_BASE_DIR}/lib:${LD_LIBRARY_PATH}"
On OS X set DYLD_LIBRARY_PATH instead of LD_LIBRARY_PATH.

## Usage

Build index from files in test/data and test/data2.

	ffindex_build -s /tmp/test.data /tmp/test.ffindex test/data test/data2

Retrieve three entries:

	ffindex_get  /tmp/test.data /tmp/test.ffindex a b foo

Unlink (Remove reference from index) an entry:

	ffindex_modify -u /tmp/test.ffindex b

Retrieve three entries, "b" should now be missing:

	ffindex_get /tmp/test.data /tmp/test.ffindex a b foo

Convert a Fasta file to ffindex, entry names are incerental IDs starting from 1:

	ffindex_from_fasta -s fasta.ffdata fasta.ffindex NC_007779.ffn

Get first entry by name:

	ffindex_get fasta.ffdata fasta.ffindex 1

Get first and third entry by entry index, this a little faster:

	ffindex_get fasta.ffdata fasta.ffindex -n 1 3

Count the characters including header in each entry:

	ffindex_apply fasta.ffdata fasta.ffindex wc -c

Count the number of characters in each sequence, without the header:

	ffindex_apply fasta.ffdata fasta.ffindex perl -ne '$x += length unless(/^>/); END{print "$x\n"}'

Parallel version for counting the characters including header in each entry:

	mpirun -np 4 ffindex_apply_mpi fasta.ffdata fasta.ffindex -- wc -c

Parallel version for counting the characters including header in each entry and
saving the output to a new ffindex:

	mpirun -np 4 ffindex_apply_mpi fasta.ffdata fasta.ffindex -i out-wc.ffindex -o out-wc.ffdata -- wc -c
