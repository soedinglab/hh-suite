package MyPaths;

use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION);
use Exporter;
our $VERSION = 2.00;
our @ISA     = qw(Exporter);
our @EXPORT  = qw($hhlib $hhdata $dummydb $perl $hhblits $dsspdir $dssp $pdbdir $ncbidir $execdir $datadir);

# Please set HHLIB 
our $hhlib   = $ENV{"HHLIB"}; # hh perl scripts (addss.pl, reformat.pl etc.)
our $hhdata  = $hhlib."/data"

our $perl    = $ENV{"HHLIB"}."/scripts"; # hh perl scripts (addss.pl, reformat.pl etc.)
$ENV{"PATH"} = $perl.":".$ENV{"PATH"}; # Add Perl scripts directory to PATH

# PSIPRED etc
our $execdir = ".../psipred/bin";        # Where the PSIPRED V2 programs have been installed
our $datadir = ".../psipred/data";       # Where the PSIPRED V2 data files have been installed    
our $hmmerdir= ".../hmmer/binaries";
our $dummydb = "$hhdata/dummydb"; # Name of a dummy blast database (single sequence formatted with formatdb)

# BLAST
our $ncbidir = ".../blast/bin";          # Where the NCBI programs have been installed (for PSIPRED in addss.pl)

# Structures
our $pdbdir  =  ".../pdb/all";            # where are the pdb files? Used in hhmakemodel.pl.
our $dsspdir =  ".../dssp/data";          # where are the dssp files? Used in addss.pl.
our $dssp    =  ".../dssp/bin/dsspcmbi";  # where is the dssp binary? Used in addss.pl.

return 1;
