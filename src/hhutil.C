/////////////////////////////////////////////////////////////////////////////////////
// Transform a character to lower case and '.' to '-' and vice versa
/////////////////////////////////////////////////////////////////////////////////////
inline char MatchChr(char c)  {return ((c>='a' && c<='z')? c-'a'+'A' : (c=='.'? '-':c) );}
inline char InsertChr(char c) {return ((c>='A' && c<='Z')? c+'a'-'A' : ((c>='0' && c<='9') || c=='-')? '.':c );}
inline int  WordChr(char c) {return (int)((c>='A' && c<='Z') || (c>='a' && c<='z'));}


/////////////////////////////////////////////////////////////////////////////////////
// Transforms the one-letter amino acid code into an integer between 0 and 22
/////////////////////////////////////////////////////////////////////////////////////
inline char aa2i(char c)
{
  //A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V
  if (c>='a' && c<='z') c+='A'-'a';
  switch (c)
    {
    case 'A': return 0;
    case 'R': return 1;
    case 'N': return 2;
    case 'D': return 3;
    case 'C': return 4;
    case 'Q': return 5;
    case 'E': return 6;
    case 'G': return 7;
    case 'H': return 8;
    case 'I': return 9;
    case 'L': return 10;
    case 'K': return 11;
    case 'M': return 12;
    case 'F': return 13;
    case 'P': return 14;
    case 'S': return 15;
    case 'T': return 16;
    case 'W': return 17;
    case 'Y': return 18;
    case 'V': return 19;
    case 'X': return ANY;
    case 'J': return ANY;
    case 'O': return ANY;
    case 'U': return 4;  //Selenocystein -> Cystein
    case 'B': return 3;  //D (or N)
    case 'Z': return 6;  //E (or Q)
    case '-': return GAP;
    case '.': return GAP;
    case '_': return GAP;
    }
  if (c>=0 && c<=32) return -1; // white space and control characters
  return -2;
}

/////////////////////////////////////////////////////////////////////////////////////
// Transforms integers between 0 and 22 into the one-letter amino acid code
/////////////////////////////////////////////////////////////////////////////////////
inline char i2aa(char c)
{
  //A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V
  switch (c)
    {
    case 0: return 'A';
    case 1: return 'R';
    case 2: return 'N';
    case 3: return 'D';
    case 4: return 'C';
    case 5: return 'Q';
    case 6: return 'E';
    case 7: return 'G';
    case 8: return 'H';
    case 9: return 'I';
    case 10: return 'L';
    case 11: return 'K';
    case 12: return 'M';
    case 13: return 'F';
    case 14: return 'P';
    case 15: return 'S';
    case 16: return 'T';
    case 17: return 'W';
    case 18: return 'Y';
    case 19: return 'V';
    case ANY: return 'X';
    case GAP: return '-';
    case ENDGAP: return '-';
    }
  return '?';
}

/////////////////////////////////////////////////////////////////////////////////////
// Transforms the dssp/psipred secondary structure code into an integer number
/////////////////////////////////////////////////////////////////////////////////////
inline char ss2i(char c)
{
  //- H E C S T G B
  if (c>='a' && c<='z') c+='A'-'a';
  switch (c)
    {
    case '.': return 0;
    case '-': return 0;
    case 'X': return 0;
    case 'H': return 1;
    case 'E': return 2;
    case 'C': return 3;
    case '~': return 3;
    case 'S': return 4;
    case 'T': return 5;
    case 'G': return 6;
    case 'B': return 7;
    case 'I': return 3;
    case ' ': return -1;
    case '\t': return -1;
    case '\n': return -1;
    }
  return -2;
}

/////////////////////////////////////////////////////////////////////////////////////
// Transforms integers between 0 and 8 into the dssp/psipred secondary structure code
/////////////////////////////////////////////////////////////////////////////////////
inline char i2ss(int c)
{
  //- H E C S T G B
  switch (c)
    {
    case 0: return '-';
    case 1: return 'H';
    case 2: return 'E';
    case 3: return 'C';
    case 4: return 'S';
    case 5: return 'T';
    case 6: return 'G';
    case 7: return 'B';
    case 8: return 'I';
    }
  return '?';
}


/////////////////////////////////////////////////////////////////////////////////////
// Transforms the solvend accessiblity code into an integer number
/////////////////////////////////////////////////////////////////////////////////////
inline char sa2i(char c)
{
  //- A B C D E
  if (c>='a' && c<='z') c+='A'-'a';
  switch (c)
    {
    case '.': return 0;
    case '-': return 0;
    case 'A': return 1;
    case 'B': return 2;
    case 'C': return 3;
    case 'D': return 4;
    case 'E': return 5;
    case 'F': return 6;
    case ' ': return -1;
    case '\t': return -1;
    case '\n': return -1;
    }
  return -2;
}

/////////////////////////////////////////////////////////////////////////////////////
// Transforms integers between 0 and 5 into the solvent accessibility code
/////////////////////////////////////////////////////////////////////////////////////
inline char i2sa(int c)
{
  //- H E C S T G B
  switch (c)
    {
    case 0: return '-';
    case 1: return 'A';
    case 2: return 'B';
    case 3: return 'C';
    case 4: return 'D';
    case 5: return 'E';
    case 6: return 'F';
    }
  return '?';
}


/////////////////////////////////////////////////////////////////////////////////////
// Transforms alternative secondary structure symbols into symbols
/////////////////////////////////////////////////////////////////////////////////////
inline char ss2ss(char c)
{
  //- H E C S T G B
  switch (c)
    {
    case '~': return 'C';
    case 'I': return 'C';
    case 'i': return 'c';
    case 'H':
    case 'E':
    case 'C':
    case 'S':
    case 'T':
    case 'G':
    case 'B':
    case 'h':
    case 'e':
    case 'c':
    case 's':
    case 't':
    case 'g':
    case 'b':
    case '.':
      return c;
    }
  return '-';
}

/////////////////////////////////////////////////////////////////////////////////////
// Transforms confidence values of psipred into internal code
/////////////////////////////////////////////////////////////////////////////////////
inline char cf2i(char c)
{
  switch (c)
    {
    case '-': return 0;
    case '.': return 0;
    case '0': return 1;
    case '1': return 2;
    case '2': return 3;
    case '3': return 4;
    case '4': return 5;
    case '5': return 6;
    case '6': return 7;
    case '7': return 8;
    case '8': return 9;
    case '9': return 10;
    }
  return 0;
}

/////////////////////////////////////////////////////////////////////////////////////
// Transforms internal representation of psipred confidence values into printable chars
/////////////////////////////////////////////////////////////////////////////////////
inline char i2cf(char c)
{
  switch (c)
    {
    case 0: return '-';
    case 1: return '0';
    case 2: return '1';
    case 3: return '2';
    case 4: return '3';
    case 5: return '4';
    case 6: return '5';
    case 7: return '6';
    case 8: return '7';
    case 9: return '8';
    case 10: return '9';
    }
  return '-';
}


/////////////////////////////////////////////////////////////////////////////////////
// Fast lookup of log2(1+2^(-x)) for x>=0 (precision < 0.35%)
/////////////////////////////////////////////////////////////////////////////////////
inline float fast_addscore(float x)
{
  static float val[2001];         // val[i]=log2(1+2^(-x))
  static char initialized=0;
  if (x>20) return 0.0;
  if (x<0)
    {
      fprintf(stderr,"Error in function fast_addscore: argument %g is negative\n",x);
      exit(7);
    }
  if (!initialized)   //First fill in the log2-vector
    {
      for (int i=0; i<=2000; i++) val[i]=log2(1.0+pow(2,-0.01*(i+0.5)));
      initialized=1;
    }
  return val[(int)(100.0*x)];
}



/////////////////////////////////////////////////////////////////////////////////////
// Little utilities for output
/////////////////////////////////////////////////////////////////////////////////////
inline void fout(FILE* outf, int d)
{
  if (d>=99999) fprintf(outf,"*\t"); else fprintf(outf,"%i\t",d);
  return;
}

/////////////////////////////////////////////////////////////////////////////////////
// Errors
/////////////////////////////////////////////////////////////////////////////////////
int FormatError(const char infile[], const char details[]="")
{
  cerr<<"Error in "<<par.argv[0]<<": wrong format while reading file \'"<<infile<<". "<<details<<"\n";
  exit(1);
}

int OpenFileError(const char outfile[])
{
  cerr<<endl<<"Error in "<<par.argv[0]<<": could not open file \'"<<outfile<<"\'\n";
  exit(2);
}

int MemoryError(const char arrayname[])
{
  cerr<<"Error in "<<par.argv[0]<<": Memory overflow while creating \'"<<arrayname<<"\'. Please report this bug to developers\n";
  exit(3);
}

int NoMemoryError(const char arrayname[])
{
  cerr<<"Error in "<<par.argv[0]<<": Could not allocate memory in \'"<<arrayname<<"\'.\n";
  exit(3);
}

int SyntaxError(const char details[]="")
{
  cerr<<"Error in "<<par.argv[0]<<" on command line: "<<details<<"\n";
  exit(4);
}

int InternalError(const char errstr[])
{
  cerr<<"Error in "<<par.argv[0]<<":  "<<errstr<<". Please report this bug to developers\n";
  exit(6);
}


/////////////////////////////////////////////////////////////////////////////////////
//// Takes family code (eg. a.1.2.3) and returns strings 'a', 'a.1', and 'a.1.2'
/////////////////////////////////////////////////////////////////////////////////////
inline void  ScopID(char cl[], char fold[], char sfam[], const char fam[])
{
  char* ptr;

  //get scop class ID
  strcpy(cl,fam);
  ptr = strchr(cl,'.');               //return adress of next '.' in name
  if(ptr) ptr[0]='\0';

  //get scop fold ID
  strcpy(fold,fam);
  ptr = strchr(fold,'.');             //return adress of next '.' in name
  if(ptr) ptr = strchr(ptr+1,'.');    //return adress of next '.' in name
  if(ptr) ptr[0]='\0';

  //get scop superfamily ID
  strcpy(sfam,fam);
  ptr = strchr(sfam,'.');            //return adress of next '.' in name
  if(ptr) ptr = strchr(ptr+1,'.');   //return adress of next '.' in name
  if(ptr) ptr = strchr(ptr+1,'.');   //return adress of next '.' in name
  if(ptr) ptr[0]='\0';
  return;
}

/////////////////////////////////////////////////////////////////////////////////////
// Read up to n lines of outfile and write to screen (STDERR)
/////////////////////////////////////////////////////////////////////////////////////
void WriteToScreen(char* outfile, int n)
{
  char line[LINELEN]="";
  ifstream outf;
  outf.open(outfile, ios::in);
  if (!outf) {OpenFileError(outfile);}
  cout<<"\n";
  for(; n>0 && outf.getline(line,LINELEN); n--) cout<<line<<"\n";
  outf.close();
  cout<<"\n";
}

inline void WriteToScreen(char* outfile) {WriteToScreen(outfile,INT_MAX);}



/////////////////////////////////////////////////////////////////////////////////////
// Read .hhdefaults file into array argv_conf (beginning at argv_conf[1])
/////////////////////////////////////////////////////////////////////////////////////
void ReadDefaultsFile(int& argc_conf, char** argv_conf, char* path=NULL)
{
  char line[LINELEN]="";
  char filename[NAMELEN];
  char* c_first;   //pointer to first character of argument string
  char* c;         //pointer to scan line read in for end of argument
  //  ifstream configf;
  FILE* configf=NULL;
  argc_conf=1;     //counts number of arguments read in

  // Open config file
  strcpy(filename,"./.hhdefaults");
  configf = fopen(filename,"r");
  if (!configf && path) 
    {
      strcpy(filename,path);
      strcat(filename,".hhdefaults");
      configf = fopen(filename,"r");
    }
  if (!configf && getenv("HOME"))
    {
      strcpy(filename,getenv("HOME"));
      strcat(filename,"/.hhdefaults");
      configf = fopen(filename,"r");
      if (!configf)
        {
          if (v>=3) cerr<<"Warning: could not find ./.hhdefaults or "<<filename<<"\n";
          return;
        }
    }
  else if (!configf) return; // only webserver has no home directory => need no warning

  // Scan file until line 'program_nameANYTHING'
  while (fgets(line,LINELEN,configf))
    if (!strncmp(line,program_name,6)) break;
  // Found line 'program_nameANYTHING'?
  if (!strncmp(line,program_name,6))
    {
      // Read in options until end-of-file or empty line
      //while (fgets(line,LINELEN,configf) && strcmp(line,"\n"))
      while (fgets(line,LINELEN,configf))
        {
          // Analyze line
          c=line;
          do
            {
              // Find next word
              while (*c==' ' || *c=='\t') c++; //Advance until next non-white space
              if ((*c=='h' && *(c+1)=='h') || *c=='\0' || *c=='\n' || *c=='#' || *c==13) break;  //Is next word empty string? (char 13 needed for Windows!)
              c_first=c;
              while (*c!=' ' && *c!='\t'  && *c!='#' && *c!='\0' && *c!='\n' && *c!=13) c++; //Advance until next white space or '#' (char 13 needed for Windows!)
              if (*c=='\0' || *c=='\n' || *c=='#' || *c==13)         //Is end of line reached? (char 13 needed for Windows!)
                {
                  *c='\0';
                  argv_conf[argc_conf]=new(char[strlen(c_first)+1]);
                  strcpy(argv_conf[argc_conf++],c_first);
                  break;
                }
              *c='\0';
              argv_conf[argc_conf]=new(char[strlen(c_first)+1]);
              strcpy(argv_conf[argc_conf++],c_first);
              if (v>2) printf("Default argument: %s\n",c_first);
              c++;
            } while (1);
	  if (*c=='h' && *(c+1)=='h') break; // Next program found
        } //end read line
      if (v>=3)
        {
          cout<<"Arguments read in from .hhdefaults ("<<filename<<"):";
          for (int argc=1; argc<argc_conf; argc++) cout<<(argv_conf[argc][0]=='-'? " ":"")<<argv_conf[argc]<<" ";
          cout<<"\n";
        }
      else if (v>=3) cout<<"Read in "<<argc_conf<<" default arguments for "<<program_name<<" from "<<filename<<"\n";
    }
  else //found no line 'program_name   anything"
    {
      if (v>=3) cerr<<endl<<"Warning: no default options for \'"<<program_name<<"\' found in "<<filename<<"\n";
      return; //no line 'program_name   anything' found
    }
  //   configf.close();
  fclose(configf);
}


/////////////////////////////////////////////////////////////////////////////////////
// Set default parameter values
/////////////////////////////////////////////////////////////////////////////////////
void SetDefaults()
{

  par.append=0;                // overwrite output file
  par.outformat=0;             // 0: hhr  1: FASTA  2:A2M   3:A3M
  par.p=20.0f;                 // minimum threshold for inclusion in hit list and alignment listing
  par.E=1e6f;                  // maximum threshold for inclusion in hit list and alignment listing
  par.b=10;                    // min number of alignments
  par.B=500;                   // max number of alignments
  par.z=10;                    // min number of lines in hit list
  par.Z=500;                   // max number of lines in hit list
  par.e=1e-3f;                 // maximum E-value for inclusion in output alignment, output HMM, and PSI-BLAST checkpoint model
  par.realign_max=1000;
  par.showcons=1;              // show consensus sequence
  par.showdssp=1;              // show predicted secondary structure ss_dssp
  par.showpred=1;              // show predicted secondary structure ss_pred
  par.showconf=0;              // don't show secondary structure confidence ss_conf
  par.cons=0;                  // chose first non-SS sequence as main representative sequence (not consensus)
  par.nseqdis=1;               // maximum number of query sequences for output alignment
  par.mark=0;                  // 1: only marked sequences (or first) get displayed; 0: most divergent ones get displayed
  par.aliwidth=80;             // number of characters per line in output alignments for HMM search

  par.max_seqid=90;            // default for maximum sequence identity threshold
  par.qid=0;                   // default for minimum sequence identity with query
  par.qsc=-20.0f;              // default for minimum score per column with query
  par.coverage=0;              // default for minimum coverage threshold
  par.Ndiff=100;               // pick Ndiff most different sequences from alignment

  par.Neff=0;                 // Filter alignment to a diversity (Neff) with a maximum Neff of par.Neff

  par.M=1;                     // match state assignment is by A2M/A3M
  par.Mgaps=50;                // Above this percentage of gaps, columns are assigned to insert states (for par.M=2)
  par.calibrate=0;             // default: no calibration
  par.calm=3;                  // derive P-values from: 0:query calibration  1:template calibration  2:both  3:Neural Network prediction
  par.mode=0;                  //

  par.wg=0;                    // 0: use local sequence weights   1: use local ones

  par.matrix=0;                // Subst.matrix 0: Gonnet, 1: HSDM, 2: BLOSUM50 3: BLOSUM62
  par.pcm=2;                   // pseudocount mode: default=divergence-dependent (but not column-specific)
  par.pca=1.0f;                // default values for substitution matrix pseudocounts
  par.pcb=1.5f;                // significant reduction of pcs by Neff_M starts around Neff_M-1=pcb
  par.pcc=1.0f;                // pcs are reduced prop. to 1/Neff^pcc
  par.pcw=0.0f;                // wc>0 weighs columns according to their intra-clomun similarity


  par.pre_pca=0.75f;            // PREFILTER - default values for substitution matrix pseudocounts 
  par.pre_pcb=1.75f;            // PREFILTER - significant reduction of pcs by Neff_M starts around Neff_M-1=pcb

  par.gapb=1.0;                // default values for transition pseudocounts
  par.gapd=0.15;               // gap open penalty pseudocount; 0.25 corresponds to 7.1*gapf bits
  par.gape=1.0;                // gap extension penalty pseudocount
  par.gapf=0.6;                // factor for increasing gap open penalty for deletes
  par.gapg=0.6;                // factor for increasing gap open penalty for inserts
  par.gaph=0.6;                // factor for increasing gap extension penalty for deletes
  par.gapi=0.6;                // factor for increasing gap extension penalty for inserts

  par.ssm=2;                   // ss scoring mode: 0:no ss score  1:score after alignment  2:score during alignment
  par.ssw=0.11f;               // weight of ss scoring
  par.ssw_realign=0.11f;       // weight of ss scoring for realign
  par.ssa=1.0f;                // weight of ss evolution matrix
  par.shift=-0.01f;            // Shift match score up
  par.mact=0.3001f;            // Score threshold for MAC alignment in local mode (set to 0.5001 to track user modification)
  par.corr=0.1f;               // Weight of correlations of scores for |i-j|<=4
  par.wstruc=1.0f;             // Weight of structure scores

  par.egq=0.0f;                // no charge for end gaps as default
  par.egt=0.0f;                // no charge for end gaps as default

  par.ssgap=0;                 // 1: add secondary structure-dependent gap penalties  0:off
  par.ssgapd=1.0f;             // secondary structure-dependent gap-opening penalty (per residue)
  par.ssgape=0.0f;             // secondary structure-dependent gap-extension penalty (per residue)
  par.ssgapi=4;                // max. number of inside-integer(ii); gap-open-penalty= -ii*ssgapd

  par.loc=1;                   // local vs. global alignment as default
  par.altali=2;                // find up to two (possibly overlapping) subalignments
  par.forward=0;               // 0: Viterbi algorithm; 1: Viterbi+stochastic sampling; 3:Maximum Accuracy (MAC) algorithm
  par.realign=1;               // realign with MAC algorithm

  par.repmode=0;               // repeats score independently of one another
  par.columnscore=1;           // Default column score is 1: null model pnul = 1/2 * (q_av(a)+p_av(a))
  par.half_window_size_local_aa_bg_freqs = 40;
  par.min_overlap=0;           // automatic minimum overlap used
  par.opt=0;                   // Default = optimization mode off
  par.readdefaultsfile=0;      // Default = do not read a defaults file ./.hhdefaults or HOME/.hhdefaults
  par.maxdbstrlen=200;         // maximum length of database string to be printed in 'Command' line of hhr file
  par.mode=0;
  par.idummy=0;
  par.premerge=0;

  par.notags=1;                // neutralize His-tags, FLAG-tags, C-myc-tags
  par.hmmer_used=false;

  // Directories for SS-prediction
  par.addss=0;
  strcpy(par.psipred,"");
  strcpy(par.psipred_data,"");

  // HHblits parameters
  par.hhblits_prefilter_logpval=0;

  par.dbsize = 0;

  // HHblits Evalue calculation  (alpha = a + b(Neff(T) - 1)(1 - c(Neff(Q) - 1)) )
  par.alphaa = 0.4;
  par.alphab = 0.02;
  par.alphac = 0.1;

  par.prefilter = false;              //true in hhblits
  par.early_stopping_filter = false;  //true in hhblits

  par.filter_thresh=0;                // 0.01 in hhblits
  par.filter_length=200;
  par.filter_evals=NULL;
  par.filter_sum=0.0;
  par.filter_counter=0;

  par.block_shading=NULL;
  par.block_shading_counter=NULL;
  par.block_shading_space = 200;
  strcpy(par.block_shading_mode,"tube");

  // For HHblits prefiltering with SSE2
  par.prefilter_gap_open = 20;
  par.prefilter_gap_extend = 4;
  par.prefilter_states = cs::AS219::kSize;
  par.prefilter_score_offset = 50;
  par.prefilter_bit_factor = 4;
  par.prefilter_evalue_thresh = 1000;
  par.preprefilter_smax_thresh = 10;

  // for filtering database alignments in HHsearch and HHblits
  par.max_seqid_db=par.max_seqid;
  par.qid_db=par.qid;            
  par.qsc_db=par.qsc;            
  par.coverage_db=par.coverage;  
  par.Ndiff_db=par.Ndiff;        

  // Initialize strings
  strcpy(par.infile,"stdin");
  strcpy(par.outfile,"");
  strcpy(par.pairwisealisfile,"");
  strcpy(par.buffer,"buffer.txt");
  strcpy(par.scorefile,"");
  strcpy(par.indexfile,""); 
  strcpy(par.wfile,"");
  strcpy(par.alnfile,"");
  strcpy(par.hhmfile,"");
  strcpy(par.psifile,"");
  strcpy(par.alitabfile,"");
  par.exclstr=NULL;

  // parameters for context-specific pseudocounts
  par.csb = 0.85;
  par.csw = 1.6;
  strcpy(par.clusterfile,""); // default in config-file: /cluster/user/michael/hh/cs/data/K4000.lib
  strcpy(par.cs_library,""); // default in config-file: /cluster/scripts/update_scripts/nr20/nr20_sampled_clusters_neff1.2_W1_N10M_n0_nopc_K62_wcenter1000_gauss_init.lib

  return;
}
