// coemit.C: 
// Generate random sequences that are similar to *two* given HMMs at the same time
// Compile with g++ coemit.C -o coemit -O3 

// hhhmm.C

#include <iostream>   // cin, cout, cerr
#include <fstream>    // ofstream, ifstream 
#include <stdio.h>    // printf
using std::cout;
using std::cerr;
using std::endl;
using std::ios;
using std::ifstream;
using std::ofstream;
#include <stdlib.h>   // exit
#include <string>     // strcmp, strstr
#include <math.h>     // sqrt, pow
#include <limits.h>   // INT_MIN
#include <float.h>    // FLT_MIN
#include <time.h>     // clock
#include <ctype.h>    // islower, isdigit etc

#define CALCULATE_MAX6(max, var1, var2, var3, var4, var5, var6, varb) \
if (var1>var2) { max=var1; varb=STOP;} \
else           { max=var2; varb=MM;}; \
if (var3>max)  { max=var3; varb=GD;}; \
if (var4>max)  { max=var4; varb=IM;}; \
if (var5>max)  { max=var5; varb=DG;}; \
if (var6>max)  { max=var6; varb=MI;}; 

// Generate random number in [0,1[
#define frand() ((float) rand()/(RAND_MAX+1.0))


const int MAXRES=15002; //max number of columns in HMM; must be <= LINELEN
const int MAXCOL=32766; //max number of residues in input files; must be <= LINELEN and >= MAXRES
const int LINELEN=32766; //max length of line read in from input files; must be >= MAXCOL 
const int MAXSEQALI=10238; //max number of sequences stored in 'hit' objects and displayed in output alignment 
const int NAA=20;       //number of amino acids (0-19)
const int IDLEN=31;     //max length of scop hierarchy id and pdb-id
const int DESCLEN=4095; //max length of sequence description (longname)
const int NAMELEN=255;  //max length of file names etc.
const int NTRANS=10;    //number of transitions recorded in HMM (M2M,M2I,M2D,I2M,I2I,D2M,D2D,M2M_GAPOPEN,GAPOPEN,GAPEXTD)
const int ANY=20;       //number representing an X (any amino acid) internally
const int GAP=21;       //number representing a gap internally 
const int ENDGAP=22;    //Important to distinguish because end gaps do not contribute to tansition counts 
const int HMMSCALE=1000;//Scaling number for log2-values in HMMs
enum transitions {M2M,M2I,M2D,I2M,I2I,D2M,D2D,M2M_GAPOPEN,GAPOPEN,GAPEXTD}; // index for transitions within a HMM


float pb[21];         // pb[a] = background probability for chosen substitution matrix
int v=3;

//Amino acids Sorted by alphabet     -> internal numbers a 
//                0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
//                A  C  D  E  F  G  H  I  K  L  M  N  P  Q  R  S  T  V  W  Y  X
const int s2a[]={ 0, 4, 3, 6,13, 7, 8, 9,11,10,12, 2,14, 5, 1,15,16,19,17,18,20};
//Internal numbers a for amino acids -> amino acids Sorted by alphabet: 
//                0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
//                A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  X
const int a2s[]={ 0,14,11, 2, 1,13, 3, 5, 6, 7, 9, 8,10, 4,12,15,16,18,19,17,20};
// Removes the newline and other control characters at the end of a string (if present)
// and returns the new length of the string (-1 if str is NULL)

inline int imin(int x, int y) { return (x<y? x : y);}

// transforms chr into an uppercase character
inline char uprchr(char chr)
{
  return (chr>='a' && chr<='z')? chr+'A'-'a' : chr;
}

// transforms chr into an lowercase character
inline char lwrchr(char chr)
{
  return (chr>='A' && chr<='Z')? chr-'A'+'a' : chr;
}

inline int chomp(char* str) 
{
  if (!str) return -1;
  int l;
  for (l=strlen(str)-1; l>=0 && str[l]<32; l--);
  str[++l]='\0';
  return l;
}

// Emulates the ifstream::getline method; similar to fgets(str,maxlen,FILE*), 
// but removes the newline at the end and returns NULL if at end of file or read error
inline char* fgetline(char* str, int maxlen, FILE* file) 
{
  if (fgets(str,maxlen,file)) 
    {
      chomp(str);
      return(str);
    }
  else return NULL;
}



// copies substring str[a,b] into substr and returns substr 
char *substr(char* substr, char* str, int a, int b)
{
  if (b<a) {int i=b; b=a; a=i;}
  if (b-a>1000) 
    {printf("Function substr: >1000 chars to copy. Exiting.\n"); exit(1);} 
  char* dest=substr;
  char* source=str+a;
  char* send=str+b;
  while (*source!='\0' && source<=send) *(dest++) = *(source++); 
  *dest='\0';
  return substr;
}


// Returns pointer to first non-white-space character in str OR to NULL if none found
inline char* strscn(char* str)
{
  if (!str) return NULL;
  char* ptr=str;
  while (*ptr!='\0' && *ptr<=32) ptr++;
  return (*ptr=='\0')? NULL: ptr;
}

// Cuts string at first white space character found by overwriting it with '\0'. 
// Returns pointer to next non-white-space char OR to NULL if no such char found 
inline char* strcut(char* str)
{
  if (!str) return NULL;
  char* ptr=str;
  while (*ptr!='\0' && *ptr>32) ptr++;
  if (*ptr=='\0') return NULL;
  *ptr='\0';
  ptr++;
  while (*ptr!='\0' && *ptr<=32) ptr++;
  return (*ptr=='\0')? NULL:ptr;
}

// Returns leftmost integer in ptr and sets the pointer to first char after 
// the integer. If no integer is found, returns INT_MIN and sets pt to NULL
int strint(char*& ptr)
{
  int i;
  char* ptr0=ptr;
  if (!ptr) return INT_MIN;
  while (*ptr!='\0' && !(*ptr>='0' && *ptr<='9')) ptr++;
  if (*ptr=='\0') 
    {
      ptr=0;
      return INT_MIN;
    }
  if (*(ptr-1)=='-' && ptr>ptr0) i=-atoi(ptr); else i=atoi(ptr);
  while (*ptr>='0' && *ptr<='9') ptr++;
  return i;
}

// Same as strint, but interpretes '*' as default
int strinta(char*& ptr, int deflt=99999)
{
  int i;
  if (!ptr) return INT_MIN;
  while (*ptr!='\0' && !(*ptr>='0' && *ptr<='9') && *ptr!='*') ptr++;
  if (*ptr=='\0') 
    {
      ptr=0;
      return INT_MIN;
    }
  if (*ptr=='*') 
    {
      ptr++;
      return deflt;
    }
  if (*(ptr-1)=='-') i=atoi(ptr-1);
  else i=atoi(ptr);
  while (*ptr>='0' && *ptr<='9') ptr++;
  return i;
}


// Swaps two integer elements in array k
inline void swapi(int k[], int i, int j)
{
  int temp;
  temp=k[i]; k[i]=k[j]; k[j]=temp;
}

// QSort sorting routine. time complexity of O(N ln(N)) on average
// Sorts the index array k between elements i='left' and i='right' in such a way that afterwards 
// v[k[i]] is sorted downwards (up=-1) or upwards (up=+1)
void QSortFloat(float v[], int k[], int left, int right, int up=+1)
{
  int i;      
  int last;   // last element to have been swapped
  
  if (left>=right) return;        // do nothing if less then 2 elements to sort
  // Put pivot element in the middle of the sort range to the side (to position 'left') ...
  swapi(k,left,(left+right)/2);  
  last=left; 
  // ... and swap all elements i SMALLER than the pivot 
  // with an element that is LARGER than the pivot (element last+1):
  if (up==1)
    {
    for (i=left+1; i<=right; i++)
      if (v[k[i]]<v[k[left]]) swapi(k,++last,i);
    }
  else
    for (i=left+1; i<=right; i++)
      if (v[k[i]]>v[k[left]]) swapi(k,++last,i);

  // Put the pivot to the right of the elements which are SMALLER, left to elements which are LARGER
  swapi(k,left,last);

  // Sort the elements left from the pivot and right from the pivot
  QSortFloat(v,k,left,last-1,up);
  QSortFloat(v,k,last+1,right,up);
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
//// Class HMM
/////////////////////////////////////////////////////////////////////////////////////

class HMM
{
 public:
  HMM();
  ~HMM();
  HMM& operator=(HMM&);
  
  int n_display;            // number of sequences stored for display of alignment (INCLUDING >ss_ and >cf_ sequences)
  char* sname[MAXSEQALI];   // names of stored sequences 
  char* seq[MAXSEQALI];     // residues of stored sequences (first at pos 1!)
  int ncons;                // index of consensus sequence
  int nfirst;               // index of first sequence (query sequence of HMM)
  int nss_dssp;             // index of seq[] with secondary structure by dssp
  int nsa_dssp;             // index of seq[] with solvent accessibility by dssp
  int nss_pred;             // index of seq[] with predicted secondary structure
  int nss_conf;             // index of seq[] with confidence values for secondary structure prediction

  int L;                    // length of HMM = number of match states; set in declaration of HMM object
  int N_in;                 // number of sequences in alignment
  int N_filtered;           // number of sequences after filtering
  float Neff_M[MAXRES];     // Neff_M[i] = diversity of subalignment of seqs that have residue in col i
  float Neff_I[MAXRES];     // Neff_I[i] = diversity of subalignment of seqs that have insert in col i
  float Neff_D[MAXRES];     // Neff_D[i] = diversity of subalignment of seqs that have delete in col i
  float Neff_HMM;           // average number of Neff over total length of HMM
  char annotchr[MAXRES];    // consensus amino acids in ASCII format, or, in HMMER format, the reference annotation character in insert line

  char longname[DESCLEN];   // Full name of first sequence of original alignment (NAME field)
  char name[NAMELEN];       // HMM name = first word in longname in lower case
  char file[NAMELEN];       // Basename (w/o path, with extension) of alignment file that was used to construct the HMM
  char fam[NAMELEN];        // family ID (derived from name) (FAM field)
  char sfam[NAMELEN];       // superfamily ID (derived from name) 
  char fold[NAMELEN];       // fold ID (derived from name)
  char cl[NAMELEN];         // class ID (derived from name)

  float lamda, mu;          // coefficients for score distribution of HMM using parameters in 'Parameters par'

  float f[MAXRES][NAA+3];   // f[i][a] = prob of finding amino acid a in column i WITHOUT pseudocounts
  float g[MAXRES][NAA];     // f[i][a] = prob of finding amino acid a in column i WITH pseudocounts
  float p[MAXRES][NAA];     // p[i][a] = prob of finding amino acid a in column i WITH OPTIMUM pseudocounts
  float tr[MAXRES][NTRANS]; // log2 of transition probabilities M2M M2I M2D I2M I2I D2M D2D M2M_GAPOPEN GAPOPEN GAPEXTD

  // Read an HMM from a HHsearch .hhm file and return 0 at end of file
  int Read(FILE* dbf, char firstline[LINELEN]=NULL);

  // Read an HMM from a HMMer .hmm file; return 0 at end of file
  int ReadHMMer(FILE* dbf);


private:
  char ss_dssp[MAXRES];     // secondary structure determined by dssp 0:-  1:H  2:E  3:C  4:S  5:T  6:G  7:B
  char sa_dssp[MAXRES];     // solvent accessibility state determined by dssp 0:-  1:A (absolutely buried) 2:B  3:C  4:D  5:E (exposed)
  char ss_pred[MAXRES];     // predicted secondary structure          0:-  1:H  2:E  3:C
  char ss_conf[MAXRES];     // confidence value of prediction         0:-  1:0 ... 10:9
  float pav[NAA];           // pav[a] = average freq of amino acids in HMM (including subst matrix pseudocounts)
  int l[MAXRES];            // l[i] = pos. of j'th match state in aligment

  // Utility for Read()
  int Warning(FILE* dbf, char line[], char name[])
    {
      if (v) cerr<<"\nWARNING: could not read line\n\'"<<line<<"\'\nin HMM "<<name<<" in "<<file<<"\n";
      while (fgetline(line,LINELEN,dbf) && !(line[0]=='/' && line[1]=='/'));
      if (line) return 2;  //return status: skip HMM
      return 0;            //return status: end of database file
    }

  friend class Hit;
};


/////////////////////////////////////////////////////////////////////////////////////
// Object constructor
/////////////////////////////////////////////////////////////////////////////////////
HMM::HMM()
{
  L=0; 
  Neff_HMM=0; 
  n_display=N_in=N_filtered=0; 
  nss_dssp=nsa_dssp=nss_pred=nss_conf=nfirst=ncons=-1;
//   lamda_hash.New(37,0.0); // Set size and NULL element for hash
//   mu_hash.New(37,0.0);    // Set size and NULL element for hash
  lamda=mu=0.0;
  name[0]=longname[0]=fam[0]='\0';
  for (int i=0; i<=L+1; i++) annotchr[i]=' ';
}


/////////////////////////////////////////////////////////////////////////////////////
// Object destructor
/////////////////////////////////////////////////////////////////////////////////////
HMM::~HMM()
{
  //Delete name and seq matrices
  for (int k=0; k<n_display; k++) {delete [] sname[k];}
  for (int k=0; k<n_display; k++) {delete [] seq[k];}
  n_display=0;
}


/////////////////////////////////////////////////////////////////////////////////////
//// Read an HMM from an HHsearch .hhm file; return 0 at end of file
/////////////////////////////////////////////////////////////////////////////////////
int HMM::Read(FILE* dbf, char firstline[LINELEN])
{
  char line[LINELEN]="";    // input line
  char str3[8],str4[8];     // first 3 and 4 letters of input line
  char* ptr;                // pointer for string manipulation
  int i=0;                  // index for match state (first=1)
  int a;                    // amino acid index
  
  L=0; 
  Neff_HMM=0; 
  n_display=N_in=N_filtered=0; 
  nss_dssp=nsa_dssp=nss_pred=nss_conf=nfirst=ncons=-1;
  lamda=mu=0.0;
  name[0]=longname[0]=fam[0]='\0';
  //If at the end of while-loop L is still 0 then we have reached end of db file

  //Do not delete name and seq vectors because their adresses are transferred to hitlist as part of a hit!!

  // Does firstline already contain first line of file?
  if (firstline!= NULL) strncpy(line,firstline,LINELEN-1);
  
  while (firstline || fgetline(line,LINELEN-1,dbf) && !(line[0]=='/' && line[1]=='/'))
    {
      
      firstline=NULL;
      if (strscn(line)==NULL) continue;    // skip lines that contain only white space
      substr(str3,line,0,2);               // copy the first three characters into str3
      substr(str4,line,0,3);               // copy the first four characters into str4

      if (!strncmp("HH",line,2)) continue;

      if (!strcmp("NAME",str4))
	{
	  ptr=strscn(line+4);              //advance to first non-white-space character
	  if (ptr && name[0]=='\0') 	  
	    {
	      strncpy(longname,ptr,DESCLEN-1); //copy full name to longname
	      longname[DESCLEN-1]='\0';
	      strncpy(name,ptr,NAMELEN-1);     //copy longname to name...
	      strcut(name);                    //...cut after first word...
	    }
	  else if (name[0]=='\0')
	    {
	      strcpy(longname,"");
	      strcpy(name,"");
	    }
	  if (v>=4) cout<<"Reading in HMM "<<name<<":\n";
	}

      else if (!strcmp("FAM",str3))
	{
	  ptr=strscn(line+3);              //advance to first non-white-space character
	  if (ptr) strncpy(fam,ptr,IDLEN-1); else strcpy(fam,""); //copy family name to basename
	  ScopID(cl,fold,sfam,fam);        //get scop classification from basename (e.g. a.1.2.3.4)
	}

      else if (!strcmp("FILE",str4))
	{
	  ptr=strscn(line+4);              //advance to first non-white-space character
	  if (ptr) 
	    strncpy(file,ptr,NAMELEN-1);     //copy family name to basename
	  else strcpy(file,"");
	}

      else if (!strcmp("LENG",str4)) 
	{ 
	  ptr=line+4; 
	  L=strint(ptr);                   //read next integer (number of match states)
	}
      else if (!strcmp("FILT",str4) || !strcmp("NSEQ",str4)) 
	{
	  ptr=line+4; 
	  N_filtered=strint(ptr);          //read next integer: number of sequences after filtering
	  N_in=strint(ptr);                //read next integer: number of sequences in alignment
	}

      else if (!strcmp("NEFF",str4) || !strcmp("NAA",str3)) sscanf(line+6,"%f",&Neff_HMM);

      else if (!strcmp("EVD",str3)) 
	{
// 	  char key[IDLEN];
	  sscanf(line+6,"%f %f",&lamda,&mu);
// 	  sscanf(line+22,"%s",key);
//  	  lamda_hash.Add(key,lamda);
//  	  mu_hash.Add(key,mu);
	}

      else if (!strcmp("DESC",str4)) continue;
      else if (!strcmp("COM",str3))  continue;
      else if (!strcmp("DATE",str4)) continue;

     /////////////////////////////////////////////////////////////////////////////////////
      // Read template sequences that should get displayed in output alignments 
      else if (!strcmp("SEQ",str3))
	{
	  char cur_seq[MAXCOL]; //Sequence currently read in
	  int k;                // sequence index; start with -1; after reading name of n'th sequence-> k=n
	  int h;                // index for character in input line
	  int l=1;              // index of character in sequence seq[k]
	  int i=1;              // index of match states in ss_dssp[i] and ss_pred[i] sequence 
	  int n_seq=0;          // number of sequences to be displayed EXCLUDING ss sequences 
	  cur_seq[0]='-';       // overwrite '\0' character at beginning to be able to do strcpy(*,cur_seq)
	  k=-1;
	  while (fgetline(line,LINELEN-1,dbf) && line[0]!='#')
	    {
	      if (v>=4) cout<<"Read from file:"<<line<<"\n"; //DEBUG
	      if (line[0]=='>') //line contains sequence name
		{
		  if (k>=MAXSEQALI-1) //maximum number of allowable sequences exceeded
		    {while (fgetline(line,LINELEN-1,dbf) && line[0]!='#'); break;}
		  k++; 
		  if      (!strncmp(line,">ss_dssp",8)) nss_dssp=k;
		  else if (!strncmp(line,">sa_dssp",8)) nsa_dssp=k;
		  else if (!strncmp(line,">ss_pred",8)) nss_pred=k;
		  else if (!strncmp(line,">ss_conf",8)) nss_conf=k;
		  else if (!strncmp(line,">Cons-",6) || !strncmp(line,">Consensus",10)) ncons=k;
		  else 
		    {
		      if (nfirst==-1) nfirst=k;
		      if (n_seq>=10)
			{while (fgetline(line,LINELEN-1,dbf) && line[0]!='#'); k--; break;}
		      n_seq++;
		    }

		  //If this is not the first sequence then store residues of previous sequence
		  if (k>0) {
		    seq[k-1]=new(char[strlen(cur_seq)+1]); 
		    if (!seq[k-1]) {cerr<<"Error: Memory overflow.\nConsider reducing the number of displayed sequences\n"; exit(0);}
		    strcpy(seq[k-1],cur_seq);
		  }

		  // store sequence name
		  strcut(line+1); //find next white-space character and overwrite it with end-of-string character
		  sname[k] = new (char[strlen(line+1)+1]); //+1 for terminating '\0'
		  if (!sname[k]) {cerr<<"Error: Memory overflow.\nConsider reducing the number of displayed sequences\n"; exit(0);}
		  strcpy(sname[k],line+1);           //store sequence name in **name
		  l=1; i=1;		
		}
	      else //line contains sequence residues
		{
		  if (k==-1) 
		    {
		      cerr<<endl<<"WARNING: Ignoring following line while reading HMM"<<name<<":\n\'"<<line<<"\'\n"; 
		      continue;
		    }

		  h=0; //counts characters in current line

		  // Check whether all characters are correct; store into cur_seq
		  if (k==nss_dssp) // lines with dssp secondary structure states (. - H E C S T G B)
		    {
		      while (h<LINELEN && line[h]>'\0' && l<MAXCOL-1)
			{
			  if (ss2i(line[h])>=0 && line[h]!='.') 
			    {
			      char c=ss2ss(line[h]);
			      cur_seq[l]=c; 
			      if (c!='.' && !(c>='a' && c<='z')) ss_dssp[i++]=ss2i(c); 
			      l++;
			    }
			  else if (v && ss2i(line[h])==-2) 
			    cerr<<endl<<"WARNING: invalid symbol \'"<<line[h]<<"\' at pos. "<<h<<" in line '"<<line<<"' of HMM "<<name<<"\n";
			  h++;
			} 
		    }
		  if (k==nsa_dssp) // lines with dssp secondary solvent accessibility (- A B C D E)
		    {
		      while (h<LINELEN && line[h]>'\0' && l<MAXCOL-1)
			{
			  if (sa2i(line[h])>=0) 
			    {
			      char c=line[h];
			      cur_seq[l]=c; 
			      if (c!='.' && !(c>='a' && c<='z')) sa_dssp[i++]=sa2i(c); 
			      l++;
			    }
			  else if (v && sa2i(line[h])==-2) 
			    cerr<<endl<<"WARNING: invalid symbol \'"<<line[h]<<"\' at pos. "<<h<<" in line '"<<line<<"' of HMM "<<name<<"\n";
			  h++;
			} 
		    }
		  else if (k==nss_pred) // lines with predicted secondary structure (. - H E C)
		    {
		      while (h<LINELEN && line[h]>'\0' && l<MAXCOL-1)
			{
			  if (ss2i(line[h])>=0 && ss2i(line[h])<=3 && line[h]!='.') 
			    {
			      char c=ss2ss(line[h]);
			      cur_seq[l]=c; 
			      if (c!='.' && !(c>='a' && c<='z')) ss_pred[i++]=ss2i(c); 
			      l++;
			    }
			  else if (v && ss2i(line[h])==-2) 
			    cerr<<endl<<"WARNING: invalid symbol \'"<<line[h]<<"\' at pos. "<<h<<" in line '"<<line<<"' of HMM "<<name<<"\n";
			  h++;
			} 
		    }
		  else if (k==nss_conf) // lines with confidence values should contain only 0-9, '-', or '.'
		    {
		      while (h<LINELEN && line[h]>'\0' && l<MAXCOL-1)
			{
			  if (line[h]=='-' || (line[h]>='0' && line[h]<='9')) 
			    {
			      cur_seq[l]=line[h]; 
			      ss_conf[l]=cf2i(line[h]); 
			      l++;
			    }
			  else if (v && cf2i(line[h])==-2) 
			    cerr<<endl<<"WARNING: invalid symbol \'"<<line[h]<<"\' at pos. "<<h<<" in line '"<<line<<"' of HMM "<<name<<"\n";
			  h++;
			} 
		    }
		  else // normal line containing residues
		    {
		      while (h<LINELEN && line[h]>'\0' && l<MAXCOL-1)
			{
			  if (aa2i(line[h])>=0 && line[h]!='.') // ignore '.' and white-space characters ' ', \t and \n (aa2i()==-1)
			    {cur_seq[l]=line[h]; l++;}
			  else if (aa2i(line[h])==-2 && v) 
			    cerr<<endl<<"WARNING: invalid symbol \'"<<line[h]<<"\' at pos. "<<h<<" in line '"<<line<<"' of HMM "<<name<<"\n";
			  h++;
			} 
		    }
		  cur_seq[l]='\0';  //Ensure that cur_seq ends with a '\0' character

		  if (v && l>MAXRES-2) 
		    cerr<<endl<<"WARNING: maximum number "<<MAXRES-2<<" of residues exceded in HMM "<<name<<"\n";
		} //end else
	    } //while(getline)
	  //If this is not the first sequence some residues have already been read in
	  if (k>=0) {
	    seq[k]=new(char[strlen(cur_seq)+1]); 
	    if (!seq[k]) {cerr<<"Error: Memory overflow.\nConsider reducing the number of displayed sequences\n"; exit(0);}
	    strcpy(seq[k],cur_seq);
	  }
	  n_display=k+1;
	  
	  // DEBUG
	  if (v>=4)
	    {
	      printf("nss_dssp=%i  nsa_dssp=%i  nss_pred=%i  nss_conf=%i  nfirst=%i\n",nss_dssp,nsa_dssp,nss_pred,nss_conf,nfirst);
	      for (k=0; k<n_display; k++)
		{
		  int j;
		  cout<<">"<<sname[k]<<"(k="<<k<<")\n";
		  if      (k==nss_dssp) {for (j=1; j<=L; j++) cout<<char(i2ss(ss_dssp[j]));}
		  else if (k==nsa_dssp) {for (j=1; j<=L; j++) cout<<char(i2sa(sa_dssp[j]));}
		  else if (k==nss_pred) {for (j=1; j<=L; j++) cout<<char(i2ss(ss_pred[j]));}
		  else if (k==nss_conf) {for (j=1; j<=L; j++) cout<<int(ss_conf[j]-1);}
		  else                  {for (j=1; j<=L; j++) cout<<seq[k][j];}
		  cout<<"\n";
		}
	    }

	} //end if("SEQ")

      /////////////////////////////////////////////////////////////////////////////////////
      // Read average amino acid frequencies for HMM
      else if (!strcmp("FREQ",str4)) 
	{
	  fprintf(stderr,"Error: hhm file has obsolete format.\n"); 
	  fprintf(stderr,"Please use hhmake version > 1.1 to generate hhm files.\n"); 
	  exit(0);
	}
      
      else if (!strcmp("NULL",str4))
	{
	  ptr=line+4;
	  for (a=0; a<20 && ptr; a++)
	    //s2[a]: transform amino acids Sorted by alphabet -> internal numbers for amino acids  
	    pb[s2a[a]] = (float) pow(2,float(-strinta(ptr))/HMMSCALE); 
	  if (!ptr) return Warning(dbf,line,name);
	  if (v>=4) 
	    {
	      printf("\nNULL  ");
	      for (a=0; a<20; a++) printf("%5.1f ",100.*pb[s2a[a]]); 
	      printf("\n");
	    }
	}

      else if (!strcmp("AVER",str4))
	{
	  ptr=line+4;
	  for (a=0; a<20 && ptr; a++)
	    //s2a[a]: transform amino acids Sorted by alphabet -> internal numbers for amino acids  
	    pav[s2a[a]] = (float) pow(2,float(-strinta(ptr))/HMMSCALE); 
	  if (!ptr) return Warning(dbf,line,name);
	  if (v>=4) 
	    {
	      printf("\nAVER  ");
	      for (a=0; a<20; a++) printf("%5.1f ",100.*pav[s2a[a]]); 
	      printf("\n");
	    }
	}

      /////////////////////////////////////////////////////////////////////////////////////
      // Read transition probabilities from start state
      else if (!strcmp("HMM",str3))
	{
	  fgetline(line,LINELEN-1,dbf); // Skip line with amino acid labels
	  fgetline(line,LINELEN-1,dbf); // Skip line with transition labels
	  ptr=line;
	  for (a=0; a<=D2D && ptr; a++)
	    tr[0][a] = float(-strinta(ptr))/HMMSCALE; //store transition probabilites as log2 values
	    // strinta returns next integer in string and puts ptr to first char 
	    // after the integer. Returns -99999 if '*' is found.
	    // ptr is set to 0 if no integer is found after ptr.
	  Neff_M[0] = float(strinta(ptr))/HMMSCALE;  // Read eff. number of sequences with M->? transition
	  Neff_I[0] = float(strinta(ptr))/HMMSCALE;  // Read eff. number of sequences with I->? transition
	  Neff_D[0] = float(strinta(ptr))/HMMSCALE;  // Read eff. number of sequences with D->? transition
	  if (!ptr) return Warning(dbf,line,name);
	  if (v>=4)
	    {
	      printf("       ");
	      for (a=0; a<=D2D; a++) 
	      printf("      ");
	      for (a=0; a<=D2D && ptr; a++) printf("%5.1f ",100*pow(2.,tr[i][a]));
	      printf("\n");
	    }
	  if (!ptr) return Warning(dbf,line,name);

	  /////////////////////////////////////////////////////////////////////////////////////
	  // Read columns of HMM
	  int next_i=0;  // index of next column
	  while (fgetline(line,LINELEN-2,dbf) &&  !(line[0]=='/' && line[1]=='/') && line[0]!='#')
	    {
	      if (strscn(line)==NULL) continue; // skip lines that contain only white space

 	      // Read in AA probabilities
	      ptr=line+1;
	      int prev_i = next_i;
	      next_i = strint(ptr); i++;
	      if (v && prev_i+1!=next_i) 
		cerr<<endl<<"WARNING: in HMM "<<name<<" state "<<prev_i<<" is followed by state "<<next_i<<"\n";
	      if (i>L)
		{
		  cerr<<endl<<"WARNING: in HMM "<<name<<" there are more columns than the stated length "<<L<<". Skipping HMM\n";
		  return 2;
		}
	      if (i>=MAXRES-2) 
		{
		  fgetline(line,LINELEN-1,dbf); // Skip line
		  continue;
		}

	      for (a=0; a<20 && ptr; a++)
		f[i][s2a[a]] = (float)pow(2.,float(-strinta(ptr))/HMMSCALE); 
	      //s2a[a]: transform amino acids Sorted by alphabet -> internal numbers for amino acids  
	      l[i]=strint(ptr);	  
	      if (!ptr) return Warning(dbf,line,name);
	      if (v>=4) 
		{
		  printf("%6i ",i);
		  for (a=0; a<20; a++) printf("%5.1f ",100*f[i][s2a[a]]); 
		  printf("%5i",l[i]);
		  printf("\n");
		}
	      
	      // Read transition probabilities
	      fgetline(line,LINELEN-1,dbf); // Skip line with amino acid labels
	      if (line[0]!=' ' && line[0]!='\t') return Warning(dbf,line,name);
	      ptr=line;
	      for (a=0; a<=D2D && ptr; a++)  
		tr[i][a] = float(-strinta(ptr))/HMMSCALE; //store transition prob's as log2-values 
	      Neff_M[i] = float(strinta(ptr))/HMMSCALE;  // Read eff. number of sequences with M->? transition
	      Neff_I[i] = float(strinta(ptr))/HMMSCALE;  // Read eff. number of sequences with I->? transition
	      Neff_D[i] = float(strinta(ptr))/HMMSCALE;  // Read eff. number of sequences with D->? transition
	      if (!ptr) return Warning(dbf,line,name);
	      if (v>=4) 
		{
		  printf("       ");
		  for (a=0; a<=D2D; a++) printf("%5.1f ",100*pow(2.,tr[i][a]));
		  printf("%5.1f %5.1f %5.1f \n",Neff_M[i],Neff_I[i],Neff_D[i]);
		}
	    }
	  if (line[0]=='/' && line[1]=='/') break;
	}
      else if (v) cerr<<endl<<"WARNING: Ignoring line\n\'"<<line<<"\'\nin HMM "<<name<<"\n";
      
    } //while(getline)

  if (L==0) return 0; //End of db file -> stop reading in

  // Set coefficients of EVD (= 0.0 if not calibrated for these parameters)
//   lamda = lamda_hash.Show(par.Key());
//   mu    = mu_hash.Show(par.Key());
  if (lamda && v>=3) printf("HMM %s is already calibrated: lamda=%-5.3f, mu=%-5.2f\n",name,lamda,mu);

  if (v && i!=L) cerr<<endl<<"Warning: in HMM "<<name<<" there are only "<<i<<" columns while the stated length is "<<L<<"\n";
  if (v && i>=MAXRES-2) {i=MAXRES-2; cerr<<endl<<"WARNING: maximum number "<<MAXRES-2<<" of residues exceeded while reading HMM "<<name<<"\n";}
  if (v && !i)  cerr<<endl<<"WARNING: HMM "<<name<<" contains no match states. Check the alignment that gave rise to this HMM.\n";
  if (v>=2) cout<<"Read in HMM "<<name<<" with "<<L<<" match states and effective number of sequences = "<<Neff_HMM<<"\n";
  L = i;

  // Set emission probabilities of zero'th (begin) state and L+1st (end) state to background probabilities
  for (a=0; a<20; a++) f[0][a]=f[L+1][a]=pb[a];

  return 1; //return status: ok
}


/////////////////////////////////////////////////////////////////////////////////////
//// Read an HMM from a HMMer .hmm file; return 0 at end of file
/////////////////////////////////////////////////////////////////////////////////////
int HMM::ReadHMMer(FILE* dbf)
{
  char line[LINELEN]="";    // input line
  char desc[DESCLEN];       // description of family
  char str4[5];             // first 4 letters of input line
  char* ptr;                // pointer for string manipulation
  int i=0;                  // index for match state (first=1)
  int a;                    // amino acid index
  char dssp=0;              // 1 if a consensus SS has been found in the transition prob lines
  char annot=0;             // 1 if at least one annotation character in insert lines is ne '-' or ' '
  int k=0;                  // index for seq[k]
  static char ignore_hmmer_cal = 0;

  L=0; 
  Neff_HMM=0; 
  n_display=N_in=N_filtered=0; 
  nss_dssp=nsa_dssp=nss_pred=nss_conf=nfirst=ncons=-1;
  lamda=mu=0.0;
  name[0]=longname[0]=desc[0]=fam[0]='\0';
  //If at the end of while-loop L is still 0 then we have reached end of db file

  // Do not delete name and seq vectors because there adresses are transferred to hitlist as part of a hit!!

  while (fgetline(line,LINELEN-1,dbf) && !(line[0]=='/' && line[1]=='/'))
    {
      
      if (strscn(line)==NULL) continue;   // skip lines that contain only white space
      if (!strncmp("HMMER",line,5)) continue;

      substr(str4,line,0,3);              // copy the first four characters into str4
 
      if (!strcmp("NAME",str4) && name[0]=='\0')
	{
	  ptr=strscn(line+4);             // advance to first non-white-space character
	  strncpy(name,ptr,NAMELEN-1);    // copy full name to name
	  strcut(name);                   // ...cut after first word...
	  if (v>=4) cout<<"Reading in HMM "<<name<<":\n";
	}

      else if (!strcmp("ACC ",str4))
	{
	  ptr=strscn(line+4);              // advance to first non-white-space character
	  strncpy(longname,ptr,DESCLEN-1); // copy Accession id to longname...
	}

      else if (!strcmp("DESC",str4))
	{
	  ptr=strscn(line+4);             // advance to first non-white-space character
	  if (ptr)
	    {
	      strncpy(desc,ptr,DESCLEN-1);   // copy description to name...
	      desc[DESCLEN-1]='\0';
	      strcut(ptr);                   // ...cut after first word...
	    }
	  if (!ptr || ptr[1]!='.' || strchr(ptr+3,'.')==NULL) strcpy(fam,""); else strcpy(fam,ptr); // could not find two '.' in name?
	}

      else if (!strcmp("LENG",str4)) 
	{ 
	  ptr=line+4; 
	  L=strint(ptr);                  //read next integer (number of match states)
	}

      else if (!strcmp("ALPH",str4)) continue;
      else if (!strcmp("RF  ",str4)) continue;
      else if (!strcmp("CS  ",str4)) continue;
      else if (!strcmp("MAP ",str4)) continue;
      else if (!strcmp("COM ",str4)) continue;
      else if (!strcmp("NSEQ",str4)) 
	{
	  ptr=line+4; 
	  N_in=N_filtered=strint(ptr);    //read next integer: number of sequences after filtering
	}

      else if (!strcmp("DATE",str4)) continue;
      else if (!strncmp("CKSUM ",line,5)) continue;
      else if (!strcmp("GA  ",str4)) continue;
      else if (!strcmp("TC  ",str4)) continue;
      else if (!strcmp("NC  ",str4)) continue;

      else if (!strncmp("SADSS",line,5)) 
	{
	  if (nsa_dssp<0) 
	    {
	      nsa_dssp=k++;
	      seq[nsa_dssp] = new(char[MAXRES+2]);
	      sname[nsa_dssp] = new(char[NAMELEN]);
	      strcpy(seq[nsa_dssp]," ");
	      strcpy(sname[nsa_dssp],"sa_dssp");
	      
	    }
	  ptr=strscn(line+5);
	  if (ptr) 
	    {
	      strcut(ptr);
	      if (strlen(seq[nsa_dssp])+strlen(ptr)>=(unsigned)(MAXRES)) 
		printf("\nWARNING: HMM %s has SADSS records with more than %i residues.\n",name,MAXRES);
	      else strcat(seq[nsa_dssp],ptr);
	    }
	}
      
      else if (!strncmp("SSPRD",line,5)) 
	{
	  if (nss_pred<0) 
	    {
	      nss_pred=k++;
	      seq[nss_pred] = new(char[MAXRES+2]);
	      sname[nss_pred] = new(char[NAMELEN]);
	      strcpy(seq[nss_pred]," ");
	      strcpy(sname[nss_pred],"ss_pred");
	      
	    }
	  ptr=strscn(line+5);
	  if (ptr) 
	    {
	      strcut(ptr);
	      if (strlen(seq[nss_pred])+strlen(ptr)>=(unsigned)(MAXRES)) 
		printf("\nWARNING: HMM %s has SSPRD records with more than %i residues.\n",name,MAXRES);
	      else strcat(seq[nss_pred],ptr);
	    }
	}
      
      else if (!strncmp("SSCON",line,5)) 
	{
	  if (nss_conf<0) 
	    {
	      nss_conf=k++;
	      seq[nss_conf] = new(char[MAXRES+2]);
	      sname[nss_conf] = new(char[NAMELEN]);
	      strcpy(seq[nss_conf]," ");
	      strcpy(sname[nss_conf],"ss_conf");
	    }
	  ptr=strscn(line+5);
	  if (ptr) 
	    {
	      strcut(ptr);
	      if (strlen(seq[nss_conf])+strlen(ptr)>=(unsigned)(MAXRES)) 
		printf("\nWARNING: HMM %s has SSPRD records with more than %i residues.\n",name,MAXRES);
	      else strcat(seq[nss_conf],ptr);
	    }
	}

      else if (!strncmp("SSCIT",line,5)) continue; 
      else if (!strcmp("XT  ",str4)) continue;
      else if (!strcmp("NULT",str4)) continue;

      else if (!strcmp("NULE",str4))
	{
	  ptr=line+4;
	  for (a=0; a<20 && ptr; a++)
	    //s2a[a]: transform amino acids Sorted by alphabet -> internal numbers for amino acids  
	    pb[s2a[a]] = (float) 0.05 * pow(2,float(strinta(ptr,-99999))/HMMSCALE); 
	  if (!ptr) return Warning(dbf,line,name);
	  if (v>=4) 
	    {
	      printf("\nNULL  ");
	      for (a=0; a<20; a++) printf("%5.1f ",100.*pb[s2a[a]]); 
	      printf("\n");
	    }
	}

      else if (!strcmp("EVD ",str4)) 
	{
	  char* ptr=line+4;
	  ptr = strscn(ptr);
	  sscanf(ptr,"%f",&lamda);
	  ptr = strscn(ptr);
	  sscanf(ptr,"%f",&mu);
	  if (lamda<0) 
	    {
	      if (ignore_hmmer_cal==0) 
		cerr<<endl<<"Warning: some HMMs have been calibrated with HMMER's 'hmmcalibrate'. These calibrations will be ignored\n"; 
	      ignore_hmmer_cal=1;
	      mu = lamda = 0.0;
	    }
	}

      /////////////////////////////////////////////////////////////////////////////////////
      // Read transition probabilities from start state
      else if (!strncmp("HMM",line,3))
	{
	  fgetline(line,LINELEN-1,dbf); // Skip line with amino acid labels
	  fgetline(line,LINELEN-1,dbf); // Skip line with transition labels
	  ptr=line;
	  for (a=0; a<=M2D && ptr; a++)
	    tr[0][a] = float(strinta(ptr,-99999))/HMMSCALE; //store transition probabilites as log2 values
	    // strinta returns next integer in string and puts ptr to first char 
	    // after the integer. Returns -99999 if '*' is found.
	    // ptr is set to 0 if no integer is found after ptr.
	  tr[0][I2M] = tr[0][D2M] = 0.0;
	  tr[0][I2I] = tr[0][D2D] = -99999.0;
	  if (!ptr) return Warning(dbf,line,name);
	  if (v>=4)
	    {
	      printf("       ");
	      for (a=0; a<=D2D && ptr; a++) printf("%5.1f ",100*pow(2.,tr[i][a]));
	      printf("\n");
	    }

	  // Prepare to store DSSP states (if there are none, delete afterwards)
	  nss_dssp=k++;
	  seq[nss_dssp] = new(char[MAXRES+2]);
	  sname[nss_dssp] = new(char[NAMELEN]);
	  strcpy(sname[nss_dssp],"ss_dssp");

	  /////////////////////////////////////////////////////////////////////////////////////
	  // Read columns of HMM
	  while (fgetline(line,LINELEN-1,dbf) &&  !(line[0]=='/' && line[1]=='/') && line[0]!='#')
	    {
	      if (strscn(line)==NULL) continue; // skip lines that contain only white space

 	      // Read in AA probabilities
	      int next_i;  // index of next column
	      ptr=line;
	      next_i = strint(ptr); i++;
	      if (v && next_i!=i) 
		cerr<<endl<<"WARNING: in HMM "<<name<<" state "<<i-1<<" is followed by state "<<next_i<<"\n";
	      if (i>L)
		{
		  cerr<<endl<<"Error: in HMM "<<name<<" there are more columns than the stated length "<<L<<"\n";
		  return 2;
		}
	      if (i>L && v)
		cerr<<endl<<"WARNING: in HMM "<<name<<" there are more columns than the stated length "<<L<<"\n";
	      if (i>=MAXRES-2) 
		{
		  fgetline(line,LINELEN-1,dbf); // Skip two lines
		  fgetline(line,LINELEN-1,dbf); 
		  continue;
		}

	      for (a=0; a<20 && ptr; a++)
		f[i][s2a[a]] = (float) pb[s2a[a]]*pow(2.,float(strinta(ptr,-99999))/HMMSCALE); 
	      //s2a[a]: transform amino acids Sorted by alphabet -> internal numbers for amino acids  
	      if (!ptr) return Warning(dbf,line,name);
	      if (v>=4) 
		{
		  printf("%6i ",i);
		  for (a=0; a<20; a++) printf("%5.1f ",100*f[i][s2a[a]]); 
		  printf("\n");
		}
	      
	      // Read insert emission line
	      fgetline(line,LINELEN-1,dbf); 
	      ptr = strscn(line);
	      if (!ptr) return Warning(dbf,line,name);
	      annotchr[i]=*ptr;
	      if (*ptr!='-' && *ptr!=' ') annot=1;
	      
	      // Read annotation character and seven transition probabilities
	      fgetline(line,LINELEN-1,dbf);
	      ptr = strscn(line);
	      switch (*ptr)
		{
		case 'H':
		  ss_dssp[i]=1;
		  seq[nss_dssp][i]=*ptr;
		  dssp=1;
		  break;
		case 'E':
		  ss_dssp[i]=2;
		  seq[nss_dssp][i]=*ptr;
		  dssp=1;
		  break;
		case 'C':
		  ss_dssp[i]=3;
		  seq[nss_dssp][i]=*ptr;
		  dssp=1;
		  break;
		case 'S':
		  ss_dssp[i]=4;
		  seq[nss_dssp][i]=*ptr;
		  dssp=1;
		  break;
		case 'T':
		  ss_dssp[i]=5;
		  seq[nss_dssp][i]=*ptr;
		  dssp=1;
		  break;
		case 'G':
		  ss_dssp[i]=6;
		  seq[nss_dssp][i]=*ptr;
		  dssp=1;
		  break;
		case 'B':
		  ss_dssp[i]=7;
		  seq[nss_dssp][i]=*ptr;
		  dssp=1;
		  break;
		case 'I':
		  dssp=1;
		case '~':
		  ss_dssp[i]=3;
		  seq[nss_dssp][i]=*ptr;
		  break;
		case '-':
		default: 
		  ss_dssp[i]=0;
		  seq[nss_dssp][i]=*ptr;
		  break;
		  
		}

	      ptr+=2;
	      for (a=0; a<=D2D && ptr; a++)  
		tr[i][a] = float(strinta(ptr,-99999))/HMMSCALE; //store transition prob's as log2-values 
	      if (!ptr) return Warning(dbf,line,name);
	      if (v>=4) 
		{
		  printf("       ");
		  for (a=0; a<=D2D; a++) printf("%5.1f ",100*pow(2.,tr[i][a]));
		  printf("\n");
		}
	    }

	  if (line[0]=='/' && line[1]=='/') break;

	}
      
    } //while(getline)
  
  if (L==0) return 0; //End of db file -> stop reading in
  
  // Set coefficients of EVD (= 0.0 if not calibrated for these parameters)
  //   lamda = lamda_hash.Show(par.Key());
  //   mu    = mu_hash.Show(par.Key());
  if (lamda && v>=2) printf("HMM %s is already calibrated: lamda=%-5.3f, mu=%-5.2f\n",name,lamda,mu);
  
  if (v && i!=L) cerr<<endl<<"Warning: in HMM "<<name<<" there are only "<<i<<" columns while the stated length is "<<L<<"\n";
  if (v && i>=MAXRES-2) {i=MAXRES-2; cerr<<endl<<"WARNING: maximum number "<<MAXRES-2<<" of residues exceeded while reading HMM "<<name<<"\n";}
  if (v && !i)  cerr<<endl<<"WARNING: HMM "<<name<<" contains no match states. Check the alignment that gave rise to this HMM.\n";
  L = i;
  
  if (strlen(longname)>0) strcat(longname," ");
  strncat(longname,name,DESCLEN-strlen(longname)-1);  // longname = ACC NAME DESC
  if (strlen(name)>0) strcat(longname," ");
  strncat(longname,desc,DESCLEN-strlen(longname)-1);
  longname[DESCLEN-1]='\0';
  ScopID(cl,fold,sfam,fam);// get scop classification from basename (e.g. a.1.2.3.4)
  strcpy(file,"");

  // Secondary structure
  if (!dssp)
    {
      // remove dssp sequence
      delete[] seq[nss_dssp];    // memory that had been allocated in case ss_dssp was given needs to be freed
      delete[] sname[nss_dssp];  // memory that had been allocated in case ss_dssp was given needs to be freed
      nss_dssp=-1;
      k--;
    }
  if (nss_pred>=0) 
    {
      for (i=1; i<=L; i++) ss_pred[i] = ss2i(seq[nss_pred][i]);
      if (nss_conf>=0) 
	for (i=1; i<=L; i++) ss_conf[i] = cf2i(seq[nss_conf][i]);
      else
	for (i=1; i<=L; i++) ss_conf[i] = 5;
    }

  sname[k]=new(char[10]);
  strcpy(sname[k],"Consensus");
  sname[k+1]=new(char[strlen(longname)+1]);
  strcpy(sname[k+1],longname);
  seq[k]=new(char[L+2]); 
  seq[k][0]=' '; 
  seq[k][L+1]='\0'; 
  seq[k+1]=new(char[L+2]); 
  seq[k+1][0]=' '; 
  seq[k+1][L+1]='\0'; 
  for (i=1; i<=L; i++)
    {  
      float pmax=0.0; 
      int amax=0;
      for (a=0; a<NAA; a++) 
	if (f[i][a]>pmax) {amax=a; pmax=f[i][a];}
      if (pmax>0.6) seq[k][i]=i2aa(amax);
      else if (pmax>0.4) seq[k][i]=lwrchr(i2aa(amax));
      else seq[k][i]='x';
      seq[k+1][i]=i2aa(amax);
    }
  nfirst=k++;

  n_display=k;

  // Calculate overall Neff_HMM
  Neff_HMM=0;
  for (i=1; i<=L; i++)
    {
      float S=0.0;
      for (a=0; a<20; a++) 
	if (f[i][a]>1E-10) S-=f[i][a]*log2(f[i][a]); 
      Neff_HMM+=(float) pow(2.0,S);
    }
  Neff_HMM/=L;
  for (i=0; i<=L; i++) Neff_M[i] = Neff_I[i] = Neff_D[i] = 10.0;
  if (v>=2) 
    cout<<"Read in HMM "<<name<<" with "<<L<<" match states and effective number of sequences = "<<Neff_HMM<<"\n";
  
  // Set emission probabilities of zero'th (begin) state and L+1st (end) state to background probabilities
  for (a=0; a<20; a++) f[0][a]=f[L+1][a]=pb[a];

  return 1; //return status: ok
}

/////////////////////////////////////////////////////////////////////////////////////
// Global variables
/////////////////////////////////////////////////////////////////////////////////////

HMM q, t;
int N=1000;      // number of sequences to generate (in random sequence or systematic mode)
int K=100;       // number of top-scoring sequences to print out
int L=1;         // length of HMMs
int n=0;         // index for (random or systematic) sequences generated
float w1=0.5;    // weight of first profile, second has 1-w1
char mode=0;     // default: systematically find highest-scoring sequences
float* sthr;     // sthr[i] is score threshold for prefix sequences up to residue i
float* colscore; // scores of all sequences found
char** seq;      // residues of generated sequences: seq[n][a]
float thr=1.0;   // score threshold for finding highest-scoring sequences is (score of highest-scoring sequence) - thr (bits)
char warn=1;


/////////////////////////////////////////////////////////////////////////////////////
//// Find all sequences with scores better than sthr[L]
/////////////////////////////////////////////////////////////////////////////////////
void FindHighScoringSequences(int i, float score_prev, char seq_prev[])
{  
  if (n<N) 
    {
      char* seq_this = new(char[i+1]);
      for (int j=1; j<i; j++) seq_this[j]=seq_prev[j];

      for (int a=0; a<20; a++) 
	{
	  float score_this = score_prev + w1*q.f[i][a] + (1-w1)*t.f[i][a];
	  if (score_this > sthr[i]) 
	    {
	      seq_this[i] = a;

	      if (v>=3) 
		{
		  printf("n=%-4i i=%-2i S=%6.2f Sthr=%6.2f  ",n,i,score_this,sthr[i]);
		  for (int j=1; j<=i; j++) printf("%c",i2aa(seq_this[j]));
		  if (i==L) printf(" => accept\n"); else printf("\n"); 
		}

	      if (i<L) 
		// Check next residue
		FindHighScoringSequences(i+1,score_this,seq_this);
	      else 
		{
		  // Found new sequence with score larger than sthr[L]
		  colscore[n] = score_this/L;
		  seq[n] = new(char[L+1]);
		  for (int j=1; j<=L; j++) seq[n][j] = seq_this[j];		      
		  n++;
		}
	    }
	}
      delete seq_this;
    } 
  else 
    {
      fprintf(stderr,"Error: more than %i sequences found with score > threshold. Lower threshold (with -t <int>) and restart\n",N);
      exit(1);
    }
  return;
}


/////////////////////////////////////////////////////////////////////////////////////
//// MAIN PROGRAM
/////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{  
  float *colscore1, *colscore2;  // scores for profile q and p, respectively, for all sequences found
  float mincolscore1=-10.0;
  float mincolscore2=-10.0;
  float minscore1, minscore2;
  float score1, score2;
  int i;           // index for HMM columns
  int a=0;         // index for amino acids    
  int k;           // index for sequences to print out
  char* program_name = argv[0];
  char *qfile, *tfile;
  char line[LINELEN]="";
  int* index;
  
   
  if (argc<3) 
    {
      printf("Generate random sequences that are similar to *two* given HMMs at the same time.\n");
      printf("The HMMs need to be in HHsearch format (generated with hhmake) or HMMer format.\n");
      printf("Usage:  coemit <model1.hhm> <model2.hhm> <options>\n");
      printf("\n");
      printf("Options:\n");
      printf(" -m {0,1}     0: systematically find n highest-scoring sequences\n");
      printf("              1: generate n random sequences and sort by score\n");
      printf(" -t <float>   threshold score in bits below best score (mode 0) (default=%.3f)\n",thr);
      printf(" -n <int>     number of (random or best-scoring) sequences to generate\n");
      printf(" -k <int>     number of top-scoring sequences to print out (default=%i)\n",K);
      printf(" -w <int>     weight of profile1 score for total score (default=%.2f)\n",w1);
      printf(" -s1 <float>  score-per-column threshold for first  profile (default=%.1f)\n",mincolscore1);
      printf(" -s2 <float>  score-per-column threshold for second profile (default=%.1f)\n",mincolscore2);
      printf("\n");
      exit(0);
    }
  qfile = argv[1];
  tfile = argv[2];

  //Processing command line input
  for (int i=3; i<argc; i++)
    { 
      if (v>=4) cout<<i<<"  "<<argv[i]<<endl; //PRINT
      if (!strcmp(argv[i],"-v") && (i<argc-1))  v=atoi(argv[++i]); 
      else if (!strcmp(argv[i],"-m") && (i<argc-1))  
	if ( (mode=atoi(argv[++i]))==0 ) N=1000000; else N=1000; // set mode (systematic or random) and default value for N
      else if (!strcmp(argv[i],"-n") && (i<argc-1))  N=atoi(argv[++i]); 
      else if (!strcmp(argv[i],"-k") && (i<argc-1))  K=atoi(argv[++i]); 
      else if (!strcmp(argv[i],"-w") && (i<argc-1))  w1=atof(argv[++i]); 
      else if (!strcmp(argv[i],"-t") && (i<argc-1))  thr=atof(argv[++i]); 
      else if (!strcmp(argv[i],"-s1") && (i<argc-1)) mincolscore1=atof(argv[++i]); 
      else if (!strcmp(argv[i],"-s2") && (i<argc-1)) mincolscore2=atof(argv[++i]); 

    }


  // Read first HMM
  FILE *qf;
  qf = fopen(qfile,"r");
  if (!qf) {cerr<<endl<<"Error in "<<program_name<<": could not open database file \'"<<qfile<<"\'\n"; exit(0);}
  fgetline(line,LINELEN,qf);  
  if (!strncmp(line,"HMMER",5))      // read HMMER format
      q.ReadHMMer(qf); 
  else if (!strncmp(line,"HH",2))     // read HHM format
      q.Read(qf);	      
  else {
    cerr<<endl<<"Error in "<<program_name<<": unrecognized HMM file format in \'"<<qfile<<"\'\n"; 
    cerr<<"line = "<<line<<"\n"; 
    exit(0);
  }
  if (v>=2)  cout<<"Reading in HMM "<<q.name<<":\n";
  
  // Read second HMM
  FILE *tf;
  tf = fopen(tfile,"r");
  if (!tf) {cerr<<endl<<"Error in "<<program_name<<": could not open database file \'"<<tfile<<"\'\n"; exit(0);}
  fgetline(line,LINELEN,tf);  
  if (!strncmp(line,"HMMER",5))      // read HMMER format
      t.ReadHMMer(tf); 
  else if (!strncmp(line,"HH",2))     // read HHM format
      t.Read(tf);	      
  else {
    cerr<<endl<<"Error in "<<program_name<<": unrecognized HMM file format in \'"<<tfile<<"\'\n"; 
    cerr<<"line = "<<line<<"\n"; 
    exit(0);
  }
  if (v>=2)  cout<<"Reading in HMM "<<t.name<<":\n";
  
  // Checking lengths of HMMs and allocating memory
  L = imin(q.L,t.L);
  if (q.L!=t.L) printf("WARNING: HMMs do not have same lengths. Taking minimum length %i.\n",L);
  seq = new(char*[N]);
  for (n=0; n<N; n++) seq[n] = new(char[L+1]);
  colscore1 = new(float[N]);
  colscore2 = new(float[N]);
  colscore  = new(float[N]);
  index = new(int[N]);
  minscore1 = mincolscore1*L;
  minscore2 = mincolscore2*L;
  
  // Systematically find best-scoring sequences OR score random sequences? 
  if (mode==0) 
    {
      ///////////////////////////////////////////////////////////////////////////////
      // Systematically find best-scoring sequences 

      // Transform profiles from lin to log
      for (i=1; i<=L; i++)
	for (a=0; a<20; a++) 
	  {
	    q.f[i][a] = log(q.f[i][a]+1E-5);   
	    t.f[i][a] = log(t.f[i][a]+1E-5);   
	  }

      // Find highest-scoring sequence and its cumulative scores
      float* scum = new(float[L+2]); 
      sthr = new(float[L+2]);
      scum[L+1]=0;
      for (i=L; i>=1; i--)
	{
	  float score_max = -10000;
	  char a_max=0;
	  for (a=0; a<20; a++) // find highest-scoring amino acid yi at pos i
	    {
	      float score = w1*q.f[i][a] + (1-w1)*t.f[i][a];
	      if (score>score_max) {a_max=a; score_max=score;}
	    }
	  // Calculate cumulative score of highest-scoring sequence yi: scum[i] = sum_j=i->L [ w1*qj*(yj) + (1-w1)*pj(yj) ]
	  scum[i] = scum[i+1] + score_max;
	  if (v>=3) printf("i=%-2i a=%-2i scum=%6.2f\n",i,a_max,scum[i]);
	}      
      for (i=1; i<=L; i++) { sthr[i] = scum[1]-thr - scum[i+1];
      if (v>=3) printf("i=%-2i scum=%6.2f sthr=%6.2f\n",i,scum[i],sthr[i]);}

      // Find sequences scoring higher than sthr
      n=0; // no sequences recorded so far
      FindHighScoringSequences(1,0.0,NULL);

      // Calculate colscore1[] and colscore2[]
      if (v>=2) printf("Calculating separate scores with two profiles for %i sequences ....\n",n);
      for (int m=0; m<n; m++) 
	{
	  score1 = score2 = 0.0;
	  for (i=1; i<=L; i++) 
	    {
	      score1 += q.f[i][(int)seq[m][i]];
	      score2 += t.f[i][(int)seq[m][i]];
	    }
	  colscore1[m] = score1/L; 
	  colscore2[m] = score2/L; 
	}

      // Sort sequences by colscore
      if (v>=2) printf("Sorting %i sequences ....\n",n);
      for (int m=0; m<N; m++) index[m] = m;
      // QSort sorting routine. time complexity of O(N ln(N)) on average
      // Sorts the index array k between elements i='left' and i='right' in such a way that afterwards 
      // v[k[i]] is sorted downwards (up=-1) or upwards (up=+1)
      QSortFloat(colscore, index, 0, n-1, -1);
      if (K>n) K=n;


      delete[] scum;
      delete[] sthr;

      // Systematically find best-scoring sequences 
      ///////////////////////////////////////////////////////////////////////////////
    }
  else 
    {
      ///////////////////////////////////////////////////////////////////////////////
      // Score random sequences

      // Calculate cumulative distribution for product of two profiles cumqt[i][a]
      float **cumqt;   // cumqt[i][a] is normalized cumulative distribution of q[i][a]*t[i][a]
      cumqt = new(float*[L+1]);
      for (i=0; i<=L; i++) cumqt[i] = new(float[21]);
      for (i=1; i<=L; i++)
	{
	  cumqt[i][0] = q.f[i][0]*t.f[i][0];
	  for (a=1; a<20; a++)
	    cumqt[i][a] = cumqt[i][a-1] + pow(q.f[i][a],w1)*pow(t.f[i][a],1-w1);
	  for (a=0; a<19; a++)
	    cumqt[i][a] /= cumqt[i][19];
	  cumqt[i][19]=1.00001;
	}	  
      
      /* initialize random generator */
      srand ( time(NULL) );
      
      // Transform profiles from lin to log
      for (i=1; i<=L; i++)
	for (a=0; a<20; a++) 
	  {
	    q.f[i][a] = log(q.f[i][a]+1E-5);   
	    t.f[i][a] = log(t.f[i][a]+1E-5);   
	  }
      
      // Generate N random sequences compatible with both profiles
      if (v>=2) printf("Generating %i sequences compatible with both profiles ....\n",N);
      for (n=0; n<N; n++)
	{
	  int trial=0;
	  for (trial=0; trial<10000; trial++) 
	    {
	      // Generate random sequence according to cumqt[i]
	      score1 = score2 = 0.0;
	      for (i=1; i<=L; i++)
		{
		  float rnd = frand();
		  for (a=0; a<20; a++) 
		    if (rnd<cumqt[i][a]) break;
		  if (a==20) {printf("Error: rnd=%.4f\n",rnd); exit(1);}
		  seq[n][i] = a;
		  
		  // Sum up score of random sequence with q ant t
		  score1 += q.f[i][a];
		  score2 += t.f[i][a];
		}
	      
	      // Accept sequence?
	      if (score1>=minscore1 && score2>=minscore2) break;
	    }
	  
	  if (trial>=10000)
	    printf("\nNo compatible sequence found in %i trials!\n",trial);
	  else 
	    {
	      // Found another compatible sequence  
	      colscore1[n] = score1/L; 
	      colscore2[n] = score2/L; 
	      colscore[n]  = ( w1*score1 + (1-w1)*score2 )/L; 
	      if (v>=3) fprintf(stderr," %i\n",trial); 
	      else if ((n+1)%1000 == 0) fprintf(stderr,"%i\n",n+1);
	    }
	}
      fprintf(stderr,"\n");
  
      // Sort sequences by colscore
      if (v>=2) printf("Sorting %i sequences ....\n",N);
      for (n=0; n<N; n++) index[n] = n;
      // QSort sorting routine. time complexity of O(N ln(N)) on average
      // Sorts the index array k between elements i='left' and i='right' in such a way that afterwards 
      // v[k[i]] is sorted downwards (up=-1) or upwards (up=+1)
      QSortFloat(colscore, index, 0, N-1, -1);
      
      for (i=0; i<=L; i++) delete[](cumqt[i]);
      delete[](cumqt); 

      // Score random sequences
      ///////////////////////////////////////////////////////////////////////////////
    }  

  


  // Print out sorted sequences 
  if (v>=2) printf("Printing out top %i sequences:\n\n",K);
  printf("No     Sequence                                S_tot      S1      S2   pos neg\n");
  for (k=0; k<K; k++) 
    {
      int pos=0, neg=0;
      n = index[k];
      printf("%-5i  ",k+1);
      for (i=1; i<=L; i++)
	{
	  printf("%c",i2aa(seq[n][i]));
	  if (seq[n][i]==3 || seq[n][i]==6)  neg++;
	  if (seq[n][i]==1 || seq[n][i]==11) pos++;
	}
      printf("   %6.3f  %6.3f  %6.3f    %2i  %2i\n",colscore[n],colscore1[n],colscore2[n],pos,neg);
    }

  // Free allocated memory
  for (n=0; n<N && seq[n]; n++) delete[](seq[n]);	
  delete[](seq);
  delete[](colscore1); 
  delete[](colscore2); 
  delete[](colscore); 
  delete[](index); 

  exit(0);
}
