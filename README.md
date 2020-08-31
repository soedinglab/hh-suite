# HH-suite3 for sensitive sequence searching

(C) Johannes Soeding, Markus Meier, Martin Steinegger, Milot Mirdita, Michael Remmert, Andreas Hauser, Andreas Biegert

[![BioConda Install](https://img.shields.io/conda/dn/bioconda/hhsuite.svg?style=flag&label=BioConda%20install)](https://anaconda.org/bioconda/hhsuite)
[![Github All Releases](https://img.shields.io/github/downloads/soedinglab/hh-suite/total.svg)](https://github.com/soedinglab/hh-suite/releases/latest)
[![Docker Pulls](https://img.shields.io/docker/pulls/soedinglab/hh-suite.svg)](https://hub.docker.com/r/soedinglab/hh-suite)
[![Travis Build Status](https://travis-ci.org/soedinglab/hh-suite.svg?branch=master)](https://travis-ci.org/soedinglab/hh-suite)
[![Codeship Status](https://codeship.com/projects/0936c290-2248-0133-bcb4-52bb0fef976f/status?branch=master)](https://codeship.com/projects/96085)

The HH-suite is an open-source software package for sensitive protein sequence searching based on the pairwise alignment of hidden Markov models (HMMs).

## Documentation

We provide an extensive [user guide](https://github.com/soedinglab/hh-suite/wiki) with many usage examples, frequently asked questions and guides to build your own databases. 

### Installation

HH-suite3 can also be installed by downloading a statically compiled version, [conda](https://github.com/conda/conda) or [Docker](https://github.com/moby/moby). HH-suite3 requires a 64-bit system (check with `uname -a | grep x86_64`) with at least the SSE2 instruction set (check by executing `cat /proc/cpuinfo | grep sse2` on Linux or `sysctl -a | grep machdep.cpu.features | grep SSE2` on macOS). `AVX2` is roughly 2x faster compared to SSE2.

```
# install via conda
conda install -c conda-forge -c bioconda hhsuite 
# install docker
docker pull soedinglab/hh-suite
# static SSE2 build
wget https://github.com/soedinglab/hh-suite/releases/download/v3.3.0/hhsuite-3.3.0-SSE2-Linux.tar.gz; tar xvfz hhsuite-3.3.0-SSE2-Linux.tar.gz; export PATH="$(pwd)/bin:$(pwd)/scripts:$PATH"
# static AVX2 build
wget https://github.com/soedinglab/hh-suite/releases/download/v3.3.0/hhsuite-3.3.0-AVX2-Linux.tar.gz; tar xvfz hhsuite-3.3.0-AVX2-Linux.tar.gz; export PATH="$(pwd)/bin:$(pwd)/scripts:$PATH"
```
:exclamation: Only the self-compiled HH-suite3 version includes MPI support, since MPI configuration is specific to the local environment.

### Available Databases
List of available database for HH-suite3: 
  - [Uniclust30](https://uniclust.mmseqs.com) [[pub]](https://doi.org/10.1093/nar/gkw1081)
  - [BFD](https://bfd.mmseqs.com) (consists of 2.5 billion, mostly enviromental, protein sequences) [[pub]](https://doi.org/10.1038/s41592-019-0437-4)
  - [Pfam/SCOP/PDB70/dbCAN](http://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/)


### Compilation
To compile from source, you will need a recent C/C++ compiler (at least GCC 4.8 or Clang 3.6) and [CMake](http://cmake.org/) 2.8.12 or later.

To download the source code and compile the HH-suite execute the following commands:
```
git clone https://github.com/soedinglab/hh-suite.git
mkdir -p hh-suite/build && cd hh-suite/build
cmake -DCMAKE_INSTALL_PREFIX=. ..
make -j 4 && make install
export PATH="$(pwd)/bin:$(pwd)/scripts:$PATH"
```

:exclamation: To compile HH-suite3 on macOS, first install the `gcc` compiler from [Homebrew](https://brew.sh). The default macOS `clang` compiler does not support OpenMP and HH-suite3 will only be able to use a single thread. Then replace the `cmake` call above with the following one:

```
CC="$(brew --prefix)/bin/gcc-8" CXX="$(brew --prefix)/bin/g++-8" cmake -DCMAKE_INSTALL_PREFIX=. ..
```    


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
HH-suite3 for fast remote homology detection and deep protein annotation, *BMC Bioinformatics*, 473. [doi: 10.1186/s12859-019-3019-7](https://doi.org/10.1186/s12859-019-3019-7)

## Links

* [Soeding lab](http://www.mpibpc.mpg.de/soeding)
* [Databases for the HH-suite](http://wwwuser.gwdg.de/~compbiol/data/hhsuite/)
