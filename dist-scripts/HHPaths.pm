#
# Please edit and adjust to your local environment.
#
package HHPaths;

# Skip this
use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION);
use Exporter;
our $VERSION = 2.00;
our @ISA     = qw(Exporter);
our @EXPORT  = qw($hhlib $hhdata $dummydb $perl $hhblits $dsspdir $dssp $pdbdir $ncbidir $execdir $datadir);


# Please set HHLIB 
our $hhlib   = $ENV{"HHLIB"};     # hh perl scripts (addss.pl, reformat.pl etc.)
our $hhdata  = $hhlib."/data";
our $hhlibbin= $hhlib."/bin";    # path to cstranslate

# PSIPRED etc
our $execdir = ".../psipred/bin";         # Where the PSIPRED V2 programs have been installed
our $datadir = ".../psipred/data";        # Where the PSIPRED V2 data files have been installed    
our $dummydb = "$hhdata/do_not_delete";   # Name of a dummy blast database (single sequence formatted with formatdb)
our $hmmerdir= ".../hmmer/binaries";      # HMMER suite by Sean Eddy
our $ncbidir = ".../blast/bin";           # Where the NCBI programs have been installed (for PSIPRED in addss.pl)

# Structures
our $pdbdir  =  ".../pdb/all";            # where are the pdb files? Used in hhmakemodel.pl.
our $dsspdir =  ".../dssp/data";          # where are the dssp files? Used in addss.pl.
our $dssp    =  ".../dssp/bin/dsspcmbi";  # where is the dssp binary? Used in addss.pl.


# Below probably does not need to be edited.

# Essential databases
our $cs_lib = "$hhdata/cs219.lib";
our $context_lib = "$ahhdata/context_data.lib";

# Perl
our $perl_dir= $ENV{"HHLIB"}."/scripts";   # hh perl scripts (addss.pl, reformat.pl etc.)
$ENV{"PATH"} = $perl_dir.":".$ENV{"PATH"}; # Add Perl scripts directory to PATH

return 1;
