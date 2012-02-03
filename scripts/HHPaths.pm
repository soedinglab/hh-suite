# HHPaths.pm 
# (c) J. Soeding, A. Hauser 2012

# HHsuite version 2.0

# PLEASE INSERT CORRECT PATHS AT POSITIONS INDICATED BY ... BELOW
# THE ENVIRONMENT VARIABLE HHLIB NEEDS TO BE SET TO YOUR LOCAL HH-SUITE DIRECTORY, 
# AS DESCRIBED IN THE HH-SUITE USER GUIDE AND README FILE

package HHPaths;

# This block can stay unmodified
use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION);
use Exporter;
our $VERSION = "version 2.0.10 (Jan 2012)";
our @ISA     = qw(Exporter);
our @EXPORT  = qw($VERSION $hhlib $hhdata $hhbin $hhscripts $execdir $datadir $ncbidir $dummydb $pdbdir $dsspdir $dssp $cs_lib $context_lib);

##############################################################################################
# PLEASE COMPLETE THE PATHS ... TO PSIPRED AND OLD-STYLE BLAST (NOT BLAST+) (NEEDED FOR PSIPRED) 
#our $execdir = ".../psipred/bin";         # path to PSIPRED V2 binaries
#our $datadir = ".../psipred/data";        # path to PSIPRED V2 data files
#our $ncbidir = ".../blast/bin";           # path to NCBI binaries (for PSIPRED in addss.pl)
our $execdir = "/cluster/bioprogs/psipred/bin";         # Where the PSIPRED V2 programs have been installed
our $datadir = "/cluster/bioprogs/psipred/data";        # Where the PSIPRED V2 data files have been installed
our $ncbidir = "/cluster/bioprogs/blast/bin";           # Where the NCBI programs have been installed (for PSIPRED in addss.pl)

##############################################################################################
# PLEASE COMPLETE THE PATHS ... TO YOUR LOCAL PDB FILES, DSSP FILES ETC.
#our $pdbdir  =  ".../pdb/all";            # where are the pdb files? Used in hhmakemodel.pl.
#our $dsspdir =  ".../dssp/data";          # where are the dssp files? Used in addss.pl.
#our $dssp    =  ".../dssp/bin/dsspcmbi";  # where is the dssp binary? Used in addss.pl.
our $pdbdir  =  "/cluster/databases/pdb/all";            # where are the pdb files? Used in hhmakemodel.pl
our $dsspdir =  "/cluster/databases/dssp/data";          # where are the dssp files? Used in addss.pl
our $dssp    =  "/cluster/databases/dssp/bin/dsspcmbi";  # where is the dssp binary? Used in addss.pl
##############################################################################################

# The lines below probably do not need to be changed

# Setting paths for hh-suite perl scripts
our $hhlib    = $ENV{"HHLIB"};     # main hh-suite directory
our $hhdata   = $hhlib."/data";    # path to data directory for hhblits, example files
our $hhbin    = $hhlib."/bin";     # path to cstranslate (path to hhsearch, hhblits etc. should be in $PATH)
our $hhscripts= $hhlib."/scripts"; # path to hh perl scripts (addss.pl, reformat.pl, hhblitsdb.pl etc.)
our $dummydb  = $hhdata."/do_not_delete"; # Name of dummy blast db for PSIPRED (single sequence formatted with NCBI formatdb)

# HHblits data files
our $cs_lib = "$hhdata/cs219.lib";
our $context_lib = "$hhdata/context_data.lib";

# Add hh-suite scripts directory to search path
$ENV{"PATH"} = $hhscripts.":".$ENV{"PATH"}; # Add hh scripts directory to PATH

return 1;
