# HHPaths.pm 

#     HHsuite version 2.0.14 (May 2012)
#     (C) J. Soeding, A. Hauser 2012

#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.

#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.

#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.

#     We are very grateful for bug reports! Please contact us at soeding@genzentrum.lmu.de
#     HHsuite version 2.0

# PLEASE INSERT CORRECT PATHS AT POSITIONS INDICATED BY ... BELOW
# THE ENVIRONMENT VARIABLE HHLIB NEEDS TO BE SET TO YOUR LOCAL HH-SUITE DIRECTORY, 
# AS DESCRIBED IN THE HH-SUITE USER GUIDE AND README FILE

package HHPaths;

# This block can stay unmodified
use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION);
use Exporter;
our $VERSION = "version 2.0.14 (Mar 2012)";
our @ISA     = qw(Exporter);
our @EXPORT  = qw($VERSION $hhlib $hhdata $hhbin $hhscripts $execdir $datadir $ncbidir $dummydb $pdbdir $dsspdir $dssp $cs_lib $context_lib $v);

##############################################################################################
# PLEASE COMPLETE THE PATHS ... TO PSIPRED AND OLD-STYLE BLAST (NOT BLAST+) (NEEDED FOR PSIPRED) 
#our $execdir = ".../psipred/bin";         # path to PSIPRED V2 binaries
#our $datadir = ".../psipred/data";        # path to PSIPRED V2 data files
#our $ncbidir = ".../blast/bin";           # path to NCBI binaries (for PSIPRED in addss.pl)
our $execdir = "/cluster/toolkit/production/bioprogs/psipred/bin";  # path to PSIPRED V2 binaries
our $datadir = "/cluster/toolkit/production/bioprogs/psipred/data"; # path to PSIPRED V2 data files
our $ncbidir = "/cluster/toolkit/production/bioprogs/blast/bin";    # path to NCBI binaries (for PSIPRED in addss.pl)

##############################################################################################
# PLEASE COMPLETE THE PATHS ... TO YOUR LOCAL PDB FILES, DSSP FILES ETC.
#our $pdbdir  =  ".../pdb/all";            # where are the pdb files? (pdb/divided directory will also work)
#our $dsspdir =  ".../dssp/data";          # where are the dssp files? Used in addss.pl.
#our $dssp    =  ".../dssp/bin/dsspcmbi";  # where is the dssp binary? Used in addss.pl.
our $pdbdir  =  "/cluster/databases/pdb/all";            # where are the pdb files? (pdb/divided directory will also work)
our $dsspdir =  "/cluster/databases/dssp/data";          # where are the dssp files? Used in addss.pl
our $dssp    =  "/cluster/databases/dssp/bin/dsspcmbi";  # where is the dssp binary? Used in addss.pl
##############################################################################################

# The lines below probably do not need to be changed

# Setting paths for hh-suite perl scripts
our $hhlib    = $ENV{"HHLIB"};     # main hh-suite directory
our $hhdata   = $hhlib."/data";    # path to data directory for hhblits, example files
our $hhbin    = $hhlib."/bin";     # path to cstranslate (path to hhsearch, hhblits etc. should be in environment variable PATH)
our $hhscripts= $hhlib."/scripts"; # path to hh perl scripts (addss.pl, reformat.pl, hhblitsdb.pl etc.)
our $dummydb  = $hhdata."/do_not_delete"; # Name of dummy blast db for PSIPRED (single sequence formatted with NCBI formatdb)

# HHblits data files
our $cs_lib = "$hhdata/cs219.lib";
our $context_lib = "$hhdata/context_data.lib";

# Add hh-suite scripts directory to search path
$ENV{"PATH"} = $hhscripts.":".$ENV{"PATH"}; # Add hh scripts directory to environment variable PATH

################################################################################################
### System command with return value parsed from output
################################################################################################
sub System()
{
    if ($v>=2) {printf("\$ %s\n",$_[0]);} 
    system($_[0]);
    if ($? == -1) {
	die("\nError: failed to execute '$_[0]': $!\n\n");	
    } elsif ($? != 0) {
	printf("\nError: command '$_[0]' returned error code %d\n\n", $? >> 8);
	return 1;
   }
    return $?;
}

return 1;
