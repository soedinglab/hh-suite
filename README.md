# Beta-Test

# HHsuite for sensitive sequence searching version 3.0-beta.3 (14-07-2017)

 (C) Johannes Soeding, Markus Meier, Martin Steinegger, Michael Remmert, Andreas Hauser, Andreas Biegert 2015

[ ![Codeship Status for soedinglab/hh-suite](https://codeship.com/projects/0936c290-2248-0133-bcb4-52bb0fef976f/status?branch=master)](https://codeship.com/projects/96085)

[ ![Build Status](https://travis-ci.org/soedinglab/hh-suite.svg?branch=master)](https://travis-ci.org/soedinglab/hh-suite)

The HH-suite is an open-source software package for sensitive protein sequence searching based on the pairwise alignment of hidden Markov models (HMMs).

## WARNING ##
We had to rename our repository for ffindex, so it might be necessary to run the following commands to update the submodule in your clone:

	git pull
	git submodule deinit
	git submodule init


## Requirements

To compile from source, you will need:
 * a recent C/C++ compiler
 * [CMake](http://cmake.org/) 2.8.12 or later
 * install [pdbx python library](https://github.com/soedinglab/pdbx)


## Installation
We recommend compiling HHsuite on the machine that should run the computations so that it can be optimized for the appropriate CPU architecture.

### Packages
Some distributions incorporate HHsuite on their own:
* Ubuntu/Debian/etc. [DPKGs](http://packages.debian.org/source/sid/hhsuite) are provided by Laszlo Kajan.
* For Archlinux you can find a [PKGBUILD on aur](https://aur.archlinux.org/packages/hhsuite/)

### Release tarballs
The [release tarballs](http://wwwuser.gwdg.de/~compbiol/data/hhsuite/releases/) should contain all required source files. Simply download and extract

### Cloning from GIT
If you want to compile the most recent version, simply clone the git repository. Then, from the repository root, initialize the ffindex submodule:

	git submodule init
	git submodule update


### Compilation
With the sourcecode ready, simply run cmake with the default settings and libraries should be auto-detected:

	mkdir build
	cd build
	cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo -G "Unix Makefiles" -DCMAKE_INSTALL_PREFIX=${INSTALL_BASE_DIR} ..
	make
	make install


### Setting paths

#### Setting environment variables
In your shell set environment variable HHLIB to ${INSTALL\_BASE\_DIR}, e.g (for bash, zsh, ksh):

	export HHLIB=${INSTALL_BASE_DIR}

HHsearch and HHblits look for the column state library file cs219.lib
and the context library file context_data.lib in ${HHLIB}/data/. The HHsuite
scripts also read HHLIB to locate the perl modules Align.pm and HHPaths.pm
in ${HHLIB]/scripts/.  

Add the location of HHsuite binaries and scripts to your search PATH variable

	export PATH=${PATH}:${INSTALL_BASE_DIR}/bin:${INSTALL_BASE_DIR}/scripts


#### Specify BLAST, PSIPRED, PDB, DSSP paths

Specify paths in ${INSTALL\_BASE\_DIR}/scripts/HHPaths.pm where they are read by HHsuite's perl scripts.


### Download Databases
Download current databases from our [server](http://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/)
To build up multiple sequences alignments using HHblits uniprot20 is sufficient.


## Usage
For performing a single search iteration of HHblits, run HHblits with the
following command:

	hhblits -i <input-file> -o <result-file> -n 1 -d <database-basename>

For generating an alignment of homologous sequences:

	hhblits -i <input-file> -o <result-file> -oa3m <result-alignment> -d <database-basename>

You can get a detailed list of options for HHblits by running HHblits with the "-h" option.

## Building packages for a Release
It might be good to do the following steps in a fresh VM with Ubuntu.
1. Adjust `/README.md` (see TODOs)
	The following line has to be updated

  		HHsuite for sensitive sequence searching version <b>current version</b> (<b>release date</b>)

1. Adjust `/CMakeLists.txt`
	1. Update the version number

			set (HHSUITE_VERSION_MAJOR 3)
			set (HHSUITE_VERSION_MINOR 0)
			set (HHSUITE_VERSION_PATCH 3)

	1. Update the release date

			set (HHSUITE_DATE "14-07-2017")

	1. Check the package version

  			set (CPACK_PACKAGE_VERSION_MAJOR "${HHSUITE_VERSION_MAJOR}")

		The beta is required for the current beta releases

  			set (CPACK_PACKAGE_VERSION_MINOR "${HHSUITE_VERSION_MINOR}-beta")
  			set (CPACK_PACKAGE_VERSION_MINOR "${HHSUITE_VERSION_PATCH}")

1. Build the packages

        rm -rf build
        mkdir build
        cd build

        cmake -DCMAKE_INSTALL_PREFIX=/home/mmeier/opt/hh-suite -DHAVE_SSSE3=1 \
        -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=OFF \
        -DCMAKE_EXE_LINKER_FLAGS_RELEASE=-static -static-libgcc \
        -static-libstdc++ -DCMAKE_FIND_LIBRARY_SUFFIXES=.a ..

        make -j 16
        make package
        make package_source

The generated packages can be found in /build

## TODO
- Build README.md with cmake to update release version and release date (compare to src/hhsuite_config.h.in)

## License

The HHsearch/HHblits software package is distributed under Gnu Public Licence, Version 3.
This means that the HH-suite is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software Foundation, either
version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE. See the GNU General Public License for more details.

See the copy of the GNU General Public License in the LICENSE file.
If you do not have this file, see http://www.gnu.org/licenses/


## Notes
For full documentation see the user guide in hhsuite-userguide.pdf


We are very grateful for bug reports!
Please contact us at soeding@mpibpc.mpg.de

## Links

* [soeding lab](http://www.mpibpc.mpg.de/soeding)
* [databases and precompiled hh-suite versions](http://wwwuser.gwdg.de/~compbiol/data/hhsuite/)


## Acknowledgements

The hhsuite contains in file hhprefilter.cpp code adapted from Michael
Farrar (http://sites.google.com/site/farrarmichael/smith-waterman).
His code is marked in the file hhprefilter.cpp. For the copy right of that
code, please see the LICENSE file that comes with HHsuite.
Reference: Farrar M. Striped Smith-Waterman speeds database searches six
times over other SIMD implementations. Bioinformatics. 2007, 23:156-61.
Many posthumous thanks to Michael Farrar for his great code!
