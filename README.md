# HH-suite3 for sensitive sequence searching

(C) Johannes Soeding, Markus Meier, Martin Steinegger, Milot, Mirdita, Michael Remmert, Andreas Hauser, Andreas Biegert

[![Travis Build Status](https://travis-ci.org/soedinglab/hh-suite.svg?branch=master)](https://travis-ci.org/soedinglab/hh-suite)[![Codeship Status for soedinglab/hh-suite](https://codeship.com/projects/0936c290-2248-0133-bcb4-52bb0fef976f/status?branch=master)](https://codeship.com/projects/96085)

The HH-suite is an open-source software package for sensitive protein sequence searching based on the pairwise alignment of hidden Markov models (HMMs).

## Documentation

We provide an extensive [user guide](https://github.com/milot-mirdita/hh-suite/wiki) with many usage examples, frequently asked questions and guides to build your own databases. 

## Installation
We recommend compiling HH-suite on the machine that should run the computations so that it can be optimized for the available CPU instruction sets.

### Compilation
To compile from source, you will need a recent C/C++ compiler (at least GCC 4.8 or Clang 3.6) and [CMake](http://cmake.org/) 2.8.12 or later.

To download the source code and compile the HH-suite execute the following commands:
```
git clone https://github.com/soedinglab/hh-suite.git
cd hh-suite && mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=. ..
make -j 4 && make install
export PATH="$(pwd)/bin:$(pwd)/scripts:$PATH"
```

### Download Databases
Download current databases from our [download server](http://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/).
To build up multiple sequences alignments using HHblits, the Uniclust30 database is sufficient.

## Usage
For performing a single search iteration of HHblits, run HHblits with the following command:
```
hhblits -i <input-file> -o <result-file> -n 1 -d <database-basename>
```

For generating an alignment of homologous sequences:
```
hhblits -i <input-file> -o <result-file> -oa3m <result-alignment> -d <database-basename>
```

A detailed list of options for HHblits is available by running HHblits with the `-h` parameter.

## Reference

Steinegger M, Meier M, Mirdita M, Vöhringer H, Haunsberger S J, and Söding J (2019)
HH-suite3 for fast remote homology detection and deep protein annotation, *bioRxiv*, 560029. [doi: 10.1101/560029](https://doi.org/10.1101/560029)

## Links

* [Soeding lab](http://www.mpibpc.mpg.de/soeding)
* [Databases for the HH-suite](http://wwwuser.gwdg.de/~compbiol/data/hhsuite/)
