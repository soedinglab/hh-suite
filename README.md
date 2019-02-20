# HH-suite3 for sensitive sequence searching

(C) Johannes Soeding, Markus Meier, Martin Steinegger, Milot, Mirdita Michael Remmert, Andreas Hauser, Andreas Biegert

[ ![Codeship Status for soedinglab/hh-suite](https://codeship.com/projects/0936c290-2248-0133-bcb4-52bb0fef976f/status?branch=master)](https://codeship.com/projects/96085)

[ ![Build Status](https://travis-ci.org/soedinglab/hh-suite.svg?branch=master)](https://travis-ci.org/soedinglab/hh-suite)

The HH-suite is an open-source software package for sensitive protein sequence searching based on the pairwise alignment of hidden Markov models (HMMs).

## Requirements

To compile from source, you will need:
 * a recent C/C++ compiler
 * [CMake](http://cmake.org/) 2.8.12 or later

## Installation
We recommend compiling HH-suite on the machine that should run the computations so that it can be optimized for the appropriate CPU architecture.

### Compilation
With the sourcecode ready, simply run cmake with the default settings and libraries should be auto-detected:

    git clone https://github.com/soedinglab/hh-suite.git
    cd hh-suite && mkdir build && cd build
	cmake -DCMAKE_INSTALL_PREFIX=${INSTALL_BASE_DIR} ..
	make && make install

### Download Databases
Download current databases from our [server](http://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/)
To build up multiple sequences alignments using HHblits Uniclust30 is sufficient.

## Usage
For performing a single search iteration of HHblits, run HHblits with the
following command:

	hhblits -i <input-file> -o <result-file> -n 1 -d <database-basename>

For generating an alignment of homologous sequences:

	hhblits -i <input-file> -o <result-file> -oa3m <result-alignment> -d <database-basename>

You can get a detailed list of options for HHblits by running HHblits with the "-h" option.

### Specify BLAST, PSIPRED, PDB, DSSP paths

Specify paths in ${INSTALL\_BASE\_DIR}/scripts/HHPaths.pm where they are read by HHsuite's perl scripts.

In your shell set the environment variable HHLIB to ${INSTALL\_BASE\_DIR}, e.g (for bash, zsh, ksh):

	export HHLIB=${INSTALL_BASE_DIR}

Some scripts will look for the column state library file cs219.lib
and the context library file context_data.lib in the ${HHLIB}/data/ folder.
The HH-suite scripts also read HHLIB to locate the perl modules Align.pm and HHPaths.pm
in ${HHLIB]/scripts/.  

Add the location of HHsuite binaries and scripts to your search PATH variable

	export PATH=${PATH}:${INSTALL_BASE_DIR}/bin:${INSTALL_BASE_DIR}/scripts


## Links

* [Soeding lab](http://www.mpibpc.mpg.de/soeding)
* [Databases for the HH-suite](http://wwwuser.gwdg.de/~compbiol/data/hhsuite/)

## Acknowledgements

The hhsuite contains in file hhprefilter.cpp code adapted from Michael
Farrar (http://sites.google.com/site/farrarmichael/smith-waterman).
His code is marked in the file hhprefilter.cpp. For the copyright of that
code, please see the LICENSE file that comes with HH-suite.
Reference: Farrar M. Striped Smith-Waterman speeds database searches six
times over other SIMD implementations. Bioinformatics. 2007, 23:156-61.
Many posthumous thanks to Michael Farrar for his great code!

