Summary: Protein sequenging analysing tools using HMM-HMM comparison
Name: hhsuite
Version: 2.0.13
Release: 1
License: GPL
Group: Utilities/System
Source: ftp://toolkit.genzentrum.lmu.de/HH-suite/hhsuite-2.0.13.tar.gz
Requires: libpng, perl
BuildRequires: libpng-devel
BuildRoot: /tmp/hhsuite-build

%description
Protein homology detection tools using HMM-HMM comparison.

%prep
%setup

%build
make all

%install
mkdir -p "$RPM_BUILD_ROOT/usr"
make install INSTALL_DIR="$RPM_BUILD_ROOT/usr"

%files
%doc CHANGES LICENSE README hhsuite-userguide.pdf

/usr/bin/hhmake
/usr/bin/hhalign
/usr/bin/hhsearch
/usr/bin/hhfilter
/usr/bin/hhconsensus
/usr/bin/hhblits
/usr/bin/ffindex_build
/usr/bin/ffindex_get
/usr/bin/ffindex_modify
/usr/include/ffindex.h
/usr/include/ffutil.h
/usr/lib/hh/bin/cstranslate
/usr/lib/hh/data/context_data.lib
/usr/lib/hh/data/cs219.lib
%{_libdir}/libffindex.a
%{_libdir}/libffindex.so
%{_libdir}/libffindex.so.0.1
/usr/lib/hh/scripts
/usr/lib/hh/scripts/create_profile_from_hmmer.pl
/usr/lib/hh/scripts/hhmakemodel.pl
/usr/lib/hh/scripts/Align.pm
/usr/lib/hh/scripts/create_profile_from_hhm.pl
/usr/lib/hh/scripts/HHPaths.pm
/usr/lib/hh/scripts/reformat.pl
/usr/lib/hh/scripts/splitfasta.pl
/usr/lib/hh/scripts/addss.pl
/usr/lib/hh/data/do_not_delete
/usr/lib/hh/data/do_not_delete.phr
/usr/lib/hh/data/do_not_delete.pin
/usr/lib/hh/data/do_not_delete.psq
