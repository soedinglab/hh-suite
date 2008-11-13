// hhalignment.C

/////////////////////////////////////////////////////////////////////////////////////
//// Class Alignment
/////////////////////////////////////////////////////////////////////////////////////

// hhalignment.C

#ifndef MAIN
#define MAIN
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
#include "util.C"     // imax, fmax, iround, iceil, ifloor, strint, strscn, strcut, substr, uprstr, uprchr, Basename etc.
#include "list.h"     // list data structure
#include "hash.h"     // hash data structure
#include "hhdecl.C"
#include "hhutil.C"   // imax, fmax, iround, iceil, ifloor, strint, strscn, strcut, substr, uprstr, uprchr, Basename etc.
#include "hhhmm.h"
#endif


/////////////////////////////////////////////////////////////////////////////////////
// Class Alignment
/////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////
// Object constructor
/////////////////////////////////////////////////////////////////////////////////////
Alignment::Alignment(int maxseq, int maxres) 
{
  longname = new(char[DESCLEN]);
  sname = new(char*[maxseq]);
  seq = new(char*[maxseq]); 
  l = new(int[maxres]);
  X = new(char*[maxseq]); 
  I = new(short unsigned int*[maxseq]);
  keep = new(char[maxseq]);
  display = new(char[maxseq]);
  wg = new(float[maxseq]);   
  ncols = new(int[maxres]);
  nseqs = new(int[maxres]);
  N_in=0;
  nres=NULL;           // number of residues per sequence k
  first=NULL;          // first residue in sequence k
  last=NULL;           // last  residue in sequence k
  ksort=NULL;          // sequence indices sorted by descending nres[k]
  name[0]='\0';           // no name defined yet
  

}

/////////////////////////////////////////////////////////////////////////////////////
// Object destructor
/////////////////////////////////////////////////////////////////////////////////////
Alignment::~Alignment()
{
  delete[] longname;
  for(int k=0; k<N_in; k++) 
    {
      delete[] sname[k]; 
      delete[] seq[k];
      delete[] X[k]; 
      delete[] I[k];
    }
  delete[] sname; 
  delete[] seq;
  delete[] X;
  delete[] I;
  delete[] l;
  delete[] keep;
  delete[] display;
  delete[] wg;
  delete[] ncols;
  delete[] nseqs;
  delete[] nres;
  delete[] first;
  delete[] last;
  delete[] ksort;
}


/////////////////////////////////////////////////////////////////////////////////////
// Reads in an alignment from file into matrix seq[k][l] as ASCII
/////////////////////////////////////////////////////////////////////////////////////
void Alignment::Read(FILE* inf, char infile[], char* firstline)
{
  int l;                  // Postion in alignment incl. gaps (first=1) 
  int h;                  // Position in input line (first=0)
  int k;                  // Index of sequence being read currently (first=0)
  char line[LINELEN]="";  // input line
  char cur_seq[MAXCOL];   // Sequence currently read in
  char* cur_name;         // Sequence currently read in
  int linenr=0;           // current line number in input file
  char skip_sequence=0;   
  RemovePathAndExtension(file,infile);  //copy rootname (w/o path) of infile into file variable of class object

  kss_dssp=ksa_dssp=kss_pred=kss_conf=kfirst=-1;
  n_display=0;
  N_in=0;
  N_filtered=0;
  N_ss=0;
  cur_seq[0]=' ';         // overwrite '\0' character at beginning to be able to do strcpy(*,cur_seq)
  l=1; k=-1;

  // Does firstline already contain first line of file?
  if (firstline!= NULL) strcpy(line,firstline);

  /////////////////////////////////////////////////////////////////////////
  // Read infile line by line
  while(firstline || fgetline(line,LINELEN,inf) && k<MAXSEQ)
    {
      linenr++;
      firstline=NULL;
      if (line[0]=='>')             //line contains sequence name
	{
	  if (k>=MAXSEQ-1) 
	    {
	      if (v>=1 && k>=MAXSEQ) 
		cerr<<endl<<"WARNING: maximum number "<<MAXSEQ<<" of sequences exceded in file "<<infile<<"\n";
	      break;
	    }
	  cur_name=line+1;          //beginning of current sequence name
	  if (k>=0) //if this is at least the second name line
	    {
	      if (strlen(cur_seq)==0) 
		{
		  cerr<<endl<<"Error: sequence "<<sname[k]<<" contains no residues."<<endl;
		  exit(1);
		}

	      // Create space for residues and paste new sequence in
	      seq[k]=new(char[strlen(cur_seq)+2]);
	      if (!seq[k]) MemoryError("array for input sequences");
	      X[k]=new(char[strlen(cur_seq)+2]);
	      if (!X[k]) MemoryError("array for input sequences");
	      I[k]=new(short unsigned int[strlen(cur_seq)+2]);
	      if (!I[k]) MemoryError("array for input sequences");
	      strcpy(seq[k],cur_seq);
	    }
	  skip_sequence=0;

	  k++; 
	  l=1; //position in current sequence (first=1)

	  // display[k]= 0: do not show in Q-T alignments  1: show if not filtered out later     2: show in any case    (do not filter out)
	  // keep[k]   = 0: do not include in profile      1: include if not filtered out later  2: include in any case (do not filter out)
	  if (line[1]=='@') cur_name++;   //skip @-character in name  
	  if (!strncmp(line,">ss_dssp",8)) {
	    if (kss_dssp<0) {display[k]=2; n_display++; keep[k]=0; kss_dssp=k; N_ss++;} else {skip_sequence=1; k--; continue;}
	  }
	  else if (!strncmp(line,">sa_dssp",8)) {
	    if (ksa_dssp<0) {display[k]=2; n_display++; keep[k]=0; ksa_dssp=k; N_ss++;} else {skip_sequence=1; k--; continue;}
	  }
	  else if (!strncmp(line,">ss_pred",8)) {
	    if (kss_pred<0) {display[k]=2; n_display++; keep[k]=0; kss_pred=k; N_ss++;} else {skip_sequence=1; k--; continue;}
   	  }
	  else if (!strncmp(line,">ss_conf",8)) {
	    if (kss_conf<0) {display[k]=2; n_display++; keep[k]=0; kss_conf=k; N_ss++;} else {skip_sequence=1; k--; continue;}
	  }
	  else if (!strncmp(line,">ss_",4) || !strncmp(line,">sa_",4)) {
	    display[k]=2; n_display++; keep[k]=0; N_ss++;
	  }
	  else if (!strncmp(line,">aa_",4)) {  // ignore sequences beginning with ">aa_"
	    skip_sequence=1; k--; continue;
	  }
	  //store first real seq
	  else if (kfirst<0)                    {display[k]=keep[k]=2; n_display++; kfirst=k;} 
	  //store all sequences
	  else if (par.mark==0)                 {display[k]=keep[k]=1; n_display++;} 
	  //store sequences up to maxqnum 
	  else if (line[1]=='@'&& n_display-N_ss<par.maxqnum) {display[k]=keep[k]=2; n_display++;}
	  else                                  {display[k]=0; keep[k]=1;}

	  // store sequence name
	  if (v>=4) printf("Reading seq %-16.16s k=%3i  n_displ=%3i  display[k]=%i keep[k]=%i\n",cur_name,k,n_display,display[k],keep[k]);
	  sname[k] = new(char[strlen(cur_name)+1]);
	  if (!sname[k]) {MemoryError("array for sequence names");}
	  strcpy(sname[k],cur_name);
	} // end if(line contains sequence name)

      else if (line[0]=='#') // Commentary line?
	{
	  // #PF01367.9 5_3_exonuc: 5'-3' exonuclease, C-terminal SAM fold; PDB 1taq, 1bgx (T:271-174), 1taq (271-174)
	  if (name[0]) continue;  // if already name defined: skip commentary line
	  char *ptr1, *ptr2;
	  ptr1=strscn_(line+1); // set ptr1 to first non-whitespace character after '#' -> AC number
	  strncpy(longname,ptr1,DESCLEN-1); // copy whole commentary line after '# ' into longname
	  longname[DESCLEN-1]='\0';
	  ptr2=strcut_(ptr1);  // cut after AC number and set ptr2 to first non-whitespace character after AC number
// 	  strcpy(fam,ptr1);    // copy AC number to fam
// 	  if (!strncmp(fam,"PF",2)) strcut_(fam,'.'); // if PFAM identifier contains '.' cut it off
//	  strcut_(ptr2);       // cut after first word ...
	  strcpy(name,ptr1);   // ... and copy first word into name
	}

      //line contains sequence residues or SS information and does not belong to a >aa_ sequence
      else if (!skip_sequence)  
	{
	  if (v>=4) cout<<line<<"\n"; //DEBUG
	  if (k==-1 && v) 
	    {
	      cerr<<endl<<"WARNING: No sequence name preceding following line in "<<infile<<":\n\'"<<line<<"\'\n"; 
	      continue;
	    }

	  h=0; //counts characters in current line

	  // Check whether all characters are correct; store into cur_seq
	  if (keep[k]) // normal line containing residues
	    {
	      while (h<LINELEN && line[h]>'\0' && l<MAXCOL-1)
		{
		  if (aa2i(line[h])>=0) // ignore white-space characters ' ', \t and \n (aa2i()==-1)
		    {cur_seq[l]=line[h]; l++;}
		  else if (aa2i(line[h])==-2 && v) 
		    cerr<<endl<<"WARNING: invalid symbol \'"<<line[h]<<"\' at pos. "<<h<<" in line "<<linenr<<" of "<<infile<<"\n";
		  h++;
		} 
	    }
	  else if (k==kss_dssp) // lines with dssp secondary structure states (. - H E C S T G B)
	    {
	      while (h<LINELEN && line[h]>'\0' && l<MAXCOL-1)
		{
		  if (ss2i(line[h])>=0 && ss2i(line[h])<=7) 
		    {cur_seq[l]=ss2ss(line[h]); l++;}
		  else if (v) 
		    cerr<<endl<<"WARNING: invalid symbol \'"<<line[h]<<"\' at pos. "<<h<<" in line "<<linenr<<" of "<<infile<<"\n";
		  h++;
		} 
	    }
	  else if (k==ksa_dssp) // lines with dssp solvent accessibility states (. - ???)
	    {
	      while (h<LINELEN && line[h]>'\0' && l<MAXCOL-1)
		{
		  if (sa2i(line[h])>=0) 
		    cur_seq[l++]=line[h];
		  else if (v) 
		    cerr<<endl<<"WARNING: invalid symbol \'"<<line[h]<<"\' at pos. "<<h<<" in line "<<linenr<<" of "<<infile<<"\n";
		  h++;
		} 
	    }
	  else if (k==kss_pred) // lines with predicted secondary structure (. - H E C)
	    {
	      while (h<LINELEN && line[h]>'\0' && l<MAXCOL-1)
		{
		  if (ss2i(line[h])>=0 && ss2i(line[h])<=3) 
		    {cur_seq[l]=ss2ss(line[h]); l++;}
		  else if (v) 
		    cerr<<endl<<"WARNING: invalid symbol \'"<<line[h]<<"\' at pos. "<<h<<" in line "<<linenr<<" of "<<infile<<"\n";
		  h++;
		} 
	    }
	  else if (k==kss_conf) // lines with confidence values should contain only 0-9, '-', or '.'
	    {
	      while (h<LINELEN && line[h]>'\0' && l<MAXCOL-1)
		{
		  if (line[h]=='-' || line[h]=='.' || (line[h]>='0' && line[h]<='9')) 
		    {cur_seq[l]=line[h]; l++;}
		  else if (v) 
		    cerr<<endl<<"WARNING: invalid symbol \'"<<line[h]<<"\' at pos. "<<l<<" in line "<<linenr<<" of "<<infile<<"\n";
		  h++;
		} 
	    }
	  else if (display[k]) // other lines such as >sa_pred etc
	    {
	      while (h<LINELEN && line[h]>'\0' && l<MAXCOL-1)
		{
		  if (line[h]=='-' || line[h]=='.' || (line[h]>='0' && line[h]<='9') || (line[h]>='A' && line[h]<='B')) 
		    {cur_seq[l]=line[h]; l++;}
		  else if (v) 
		    cerr<<endl<<"WARNING: invalid symbol \'"<<line[h]<<"\' at pos. "<<l<<" in line "<<linenr<<" of "<<infile<<"\n";
		  h++;
		} 
	    }
	  if (v && l>=MAXCOL-1) 
	    {
	      cerr<<endl<<"WARNING: maximum number of residues "<<MAXCOL-2<<" exceded in sequence "<<sname[k]<<"\n";
	      skip_sequence=1;
	    }
	  cur_seq[l]='\0';  //Ensure that cur_seq ends with a '\0' character
	} //end else

    } 
  /////////////////////////////////////////////////////////////////////////


  if (k>=0) //if  at least one sequence was read in
    {
      seq[k]=new(char[strlen(cur_seq)+2]);
      if (!seq[k]) MemoryError("array for input sequences");
      X[k]=new(char[strlen(cur_seq)+2]);
      if (!X[k]) MemoryError("array for input sequences");
      I[k]=new(short unsigned int[strlen(cur_seq)+2]);
      if (!I[k]) MemoryError("array for input sequences");
      strcpy(seq[k],cur_seq);
    }
  else
    {cerr<<endl<<"Error: no sequences found in file "<<infile<<"\n"; exit(1);} 

  N_in = k+1;

  // Set name, longname, fam
  if (!*name)  // longname, name and family were not set by '#...' line yet -> extract from first sequence
    {
      char* ptr;
//    strtr(sname[kfirst],"~"," ");              // 'transpose': replaces the tilde with a blanc everywhere in sname[kfirst]
      strncpy(longname,sname[kfirst],DESCLEN-1); // longname is name of first sequence
      longname[DESCLEN-1]='\0';
      strncpy(name,sname[kfirst],NAMELEN-1);     // Shortname is first word of longname...
      name[NAMELEN-1]='\0';
      ptr = strcut(name);                        // ...until first white-space character
      if (ptr && islower(ptr[0]) && ptr[1]=='.' && isdigit(ptr[2])) //Scop family code present as second word? 
	{
	  lwrstr(name);                          // Transform upper case to lower case
	  strcut(ptr);                           // Non-white-space characters until next white-space character..
	  strcpy(fam,ptr);                       // ...are the SCOP familiy code
	}
      else if (name[0]=='P' && name[1]=='F' && isdigit(name[2]) &&  isdigit(name[3]) ) //Pfam code
	{
	  strcpy(fam,name);                      // set family name = Pfam code
	} 
      else strcpy(fam,"");
    }




  
  // Checking for warning messages
  if (v==0) return;
  if (v>=2) cout<<"Read "<<infile<<" with "<<N_in<<" sequences\n";
  if (v>=3) cout<<"Query sequence for alignment has number "<<kfirst<<" (0 is first)\n";
  return;
}


/////////////////////////////////////////////////////////////////////////////////////
// Convert ASCII in seq[k][l] to int (0-20), throw out all insert states, record their number in I[k][i] 
// and store sequences to be displayed in seq[k]
/////////////////////////////////////////////////////////////////////////////////////
void Alignment::Compress(char infile[])
{
  int i;                  //Index for match state (first=1)
  int l;                  //Postion in alignment incl. gaps (first=1) 
  int k;                  //Index for sequences (first=0)
  int a;                  //amino acid index
  char c;
  int unequal_lengths=0;  //k: seq k doesn't have same number of match states as seq 0 => WARNING
  static short unsigned int h[MAXSEQ];  //points to next character in seq[k] to be written

  // Initialize 
  for (k=0;k<N_in; k++) I[k][0]=0;

  if (v>=3) 
    {
      if (par.M==1) 
	cout<<"Using match state assignment by capital letters (a2m format)\n";
      else if (par.M==2) cout<<"Using percentage-rule match state assignment\n";
      else if (par.M==3) cout<<"Using residues of first sequence as match states\n";
    }

  // Create matrices X and I with amino acids represented by integer numbers
  switch(par.M) 
    {

    ///////////////////////////////////////////////////////////////////////////////////// 
    // a2m/a3m format: match states capital case, inserts lower case, delete states '-', inserted gaps '.'
    // The confidence values for ss prediction are interpreted as follows: 0-9:match states(!)  '-' :match state  '.':insert
    case 1:
    default:

      // Warn if alignment is ment to be -M first or -M NN instead of A2M/A3M
      if (v && strchr(seq[kfirst],'-') ) // Seed/query sequence contains a gap ...
	{
	  for (k=1; k<N_in; k++) 
	    if (strpbrk(seq[k],"abcdefghiklmnpqrstuvwxyz.")) break;
	  if (k==N_in) // ... but alignment contains no lower case residue 
	    printf("WARNING: input alignment %s looks like aligned FASTA instead of A2M/A3M format. Consider using '-M first' or '-M 50'\n",infile); 
	}

      // Remove '.' characters from seq[k]
      for(k=0; k<N_in; k++)
	{
	  char* ptrS=seq[k];         // pointer to source: character in seq[k]
	  char* ptrD=seq[k];         // pointer to destination: seq[k]
	  while(1) 	             // omit '.' symbols
	    {
	      if (*ptrS!='.') {*ptrD=*ptrS; ptrD++;} //leave out '.' symbols
	      if (!*ptrS) break;
	      ptrS++;
	    }      
	}
      L=MAXRES-2; // needed because L=imin(L,i)
      for (k=0; k<N_in; k++)
	{
	  i=1; l=1;
	  if (keep[k])  //skip >ss_dssp, >ss_pred, >ss_conf, >aa_... sequences
	    {
	      while((c=seq[k][l++]))  // assign residue to c at same time
		{
		  if (c>='a' && c<='z') I[k][i-1]++;//insert state = lower case character
		  else if (c!='.') //match state = upper case character
		    {
		      X[k][i]=aa2i(c);      
		      I[k][i]=0;      
		      i++;
		    }		    
		}
	    }
	  else if (k==kss_dssp || k==kss_pred) // does alignment contain sequence of secondary structure states?
	    {
	      while((c=seq[k][l++]))  // assign residue to c at same time
		if (c!='.' && !(c>='a' && c<='z')) X[k][i++]=ss2i(c); //match state = upper case character
	    }
	  else if (k==ksa_dssp)        // does alignment contain sequence of prediction confidence values?
	    {
	      while((c=seq[k][l++]))  // assign residue to c at same time
		if (c!='.' && !(c>='a' && c<='z')) X[k][i++]=sa2i(c); //match state = upper case character
	    }
	  else if (k==kss_conf)        // does alignment contain sequence of prediction confidence values?
	    {
	      while((c=seq[k][l++]))  // assign residue to c at same time
		if (c!='.') X[k][i++]=cf2i(c); //match state = 0-9 or '-'
	    }
	  else continue;
	  i--;
	  if (L!=i && L!=MAXRES-2 && !unequal_lengths) unequal_lengths=k;   //sequences have different lengths
	  L=imin(L,i);
	}
      if (unequal_lengths) break;
      
      for (i=1; i<=L; i++) this->l[i]=i; //assign column indices to match states
      if (L<=0) 
	{
	  cout<<"\nError: Alignment in "<<infile<<" contains no match states. Consider using -M first or -M <int> option"<<endl; 
	  exit(1);
	}
  
      if (v>=2) cout<<"Alignment in "<<infile<<" contains "<<L<<" match states\n";
      break;
      
    ///////////////////////////////////////////////////////////////////////////////////// 
    // gap-rule assignment of match states
    case 2: 
      int nl[NAA+2];              //nl[a] = number of seq's with amino acid a at position l 
      float percent_gaps[MAXCOL]; //percentage of gaps in column k (with weighted sequences)

      //determine number of columns L in alignment
      L=strlen(seq[kfirst])-1;

      // Conversion to integer representation, checking for unequal lengths and initialization
      if (nres==NULL) nres=new(int[N_in]);
      for (k=0; k<N_in; k++) 
	{
	  if (!keep[k]) continue;
	  int nr=0;
	  wg[k]=0; nres[k]=0;
	  for (l=1; l<=L; l++) 
	    {
	      X[k][l]=aa2i(seq[k][l]); 
	      if (X[k][l]<NAA) nr++;		
	    }
	  nres[k]=nr;
	  if (seq[k][L+1]!='\0' && !unequal_lengths) unequal_lengths=k; 
	}
      if (unequal_lengths) break;

      // Quick and dirty calculation of the weight per sequence wg[k]
      for (l=1; l<=L; l++) // for all positions l in alignment
	{
	  int naa=0;            //number of different amino acids
	  for (a=0; a<20; a++) nl[a]=0;
	  for (k=0; k<N_in; k++) if (keep[k]) nl[ (int)X[k][l]]++;
	  for (a=0; a<20; a++) if(nl[a]) naa++;
	  if (!naa) naa=1;      //naa=0 when column consists of only gaps and Xs (=ANY) 
	  for (k=0; k<N_in; k++) 
	    if (keep[k] && X[k][l]<20) 
	      wg[k]+=1.0/float(nl[ (int)X[k][l]]*naa*nres[k]);
	}
      // Add up percentage of gaps
      for (l=1; l<=L; l++) 
	{
	  float res=0;
	  float gap=0;
	  for (k=0; k<  N_in; k++)
	    if (keep[k] && X[k][l]<20) res+=wg[k]; else  gap+=wg[k]; // else: ANY, GAP or ENDGAP
	  percent_gaps[l]=100.*gap/(res+gap);
	  if (v>=3) cout<<"percent gaps["<<l<<"]="<<percent_gaps[l]<<"\n";
	}
      // Throw out insert states and keep only match states
      i=0;
      for (k=0; k<N_in; k++) {h[k]=1; seq[k][0]='-';}
      for (l=1; l<=L; l++)
	{
	  if (percent_gaps[l]<=float(par.Mgaps))
	    {
	      i++;
	      this->l[i]=l;
	      for (k=0; k<N_in; k++)
		{
		  if (keep[k])
		    {
		      seq[k][h[k]++]=MatchChr(seq[k][l]);
		      X[k][i]=X[k][l];
		      I[k][i]=0;
		    }
 		  else if (k==kss_dssp || k==kss_pred) 
		    {
		      seq[k][h[k]++]=MatchChr(seq[k][l]);
		      X[k][i]=ss2i(seq[k][l]);		      
		    }
 		  else if (k==ksa_dssp) 
		    {
		      seq[k][h[k]++]=MatchChr(seq[k][l]);
		      X[k][i]=sa2i(seq[k][l]);		      
		    }
		  else if (k==kss_conf)
		    {
		      seq[k][h[k]++]=seq[k][l];
		      X[k][i]=cf2i(seq[k][l]);		      
		    }
		}
	    }
	  else
	    {
	      for (k=0; k<N_in; k++)
		if (keep[k] && X[k][l]<GAP) 
		  {
		    I[k][i]++; 
		    seq[k][h[k]++]=InsertChr(seq[k][l]);
		  }
	    }
	}
      for (k=0; k<N_in; k++) seq[k][h[k]]='\0';
      
      if (v>=2) cout<<"Alignment in "<<infile<<" contains "<<L<<" columns and "<<i<<" match states\n";
      L = i;        //Number of match states
      break;
      

    ///////////////////////////////////////////////////////////////////////////////////// 
    // Using residues of first sequence as match states
    case 3:
      char match_state[MAXCOL]; //1: column assigned to match state 0: insert state

      //determine number of columns L in alignment
      L=strlen(seq[0]+1); printf("Length of first seq = %i\n",L);
      // Check for sequences with unequal lengths
      for (k=1; k<N_in; k++) 
	if (int(strlen(seq[k]+1))!=L) {unequal_lengths=k; break;} 
      if (unequal_lengths) break;

      // Determine match states: seq kfirst has residue at pos l -> match state
      for (l=1; l<=L; l++)
	if (isalpha(seq[kfirst][l])) match_state[l]=1; else match_state[l]=0;
      // Throw out insert states and keep only match states
      for (k=0; k<N_in; k++) {h[k]=1; seq[k][0]='-';}
      i=0;
      for (l=1; l<=L; l++)
	{
	  if (match_state[l]) // does sequence 0 have residue at position l?
	    {
	      i++;
	      this->l[i]=l;
	      for (k=0; k<N_in; k++)
		{
		  if (keep[k])
		    {
		      seq[k][h[k]++]=MatchChr(seq[k][l]);
		      X[k][i]=aa2i(seq[k][l]);
		      I[k][i]=0;
		    }
 		  else if (k==kss_dssp || k==kss_pred) 
		    {
		      seq[k][h[k]++]=MatchChr(seq[k][l]);
		      X[k][i]=ss2i(seq[k][l]);		      
		    }
 		  else if (k==ksa_dssp) 
		    {
		      seq[k][h[k]++]=MatchChr(seq[k][l]);
		      X[k][i]=sa2i(seq[k][l]);		      
		    }
		  else if (k==kss_conf)
		    {
		      seq[k][h[k]++]=seq[k][l];
		      X[k][i]=cf2i(seq[k][l]);		      
		    }
		}
	    }
	  else
	    {
	      for (k=0; k<N_in; k++)
		if (keep[k] && aa2i(seq[k][l])<GAP) 
		  {
		    I[k][i]++; 
		    seq[k][h[k]++]=InsertChr(seq[k][l]);
		  }
	    }
	}
      for (k=0; k<N_in; k++) seq[k][h[k]]='\0';
 
      if (v>=2) cout<<"Alignment in "<<infile<<" contains "<<L<<" columns and "<<i<<" match states\n";
      L = i;        //Number of match states
      break;
      
    } //end switch()
  ///////////////////////////////////////////////////////////////////////////////////// 
  

  // Error
  if (unequal_lengths) 
    {
      strcut(sname[unequal_lengths]);
      if(par.M==1)
	cerr<<endl<<"Error: sequences in "<<infile<<" do not all have the same number of match states,\ne.g. first sequence and sequence "<<sname[unequal_lengths]<<".\nCheck input format for '-M a2m' option and consider using '-M first' or '-M 50'\n";
      else  cerr<<endl<<"Error: sequences in "<<infile<<" do not all have the same number of columns, \ne.g. first sequence and sequence "<<sname[unequal_lengths]<<".\n";
      exit(1);
    }

  //Replace GAP with ENDGAP for all end gaps 
  for (k=0; k<N_in; k++)
    {
      if (!keep[k]) continue;
      for (i=1; i<=L && X[k][i]==GAP; i++) X[k][i]=ENDGAP;
      for (i=L; i>=1 && X[k][i]==GAP; i--) X[k][i]=ENDGAP;
    }

  // DEBUG
  if (v>=4)
    for (k=0; k<N_in; k++)
      {
	if (!display[k]) continue;
	cout<<">"<<sname[k]<<"\n";
	if (k==kss_dssp || k==kss_pred) {for (i=1; i<=L; i++) cout<<char(i2ss(X[k][i]));}
	else if (k==kss_conf)           {for (i=1; i<=L; i++) cout<<char(i2cf(X[k][i]));}
	else if (k==ksa_dssp)           {for (i=1; i<=L; i++) cout<<char(i2sa(X[k][i]));}
	else
	  {
	    for (i=1; i<=L; i++) cout<<char(i2aa(X[k][i]));
	    cout<<"\n";
	    for (i=1; i<=L; i++) 
	      if (I[k][i]==0) cout<<"-"; else if (I[k][i]>9) cout<<"X"; else cout<<I[k][i];
	  } 
	cout<<"\n";
     }
}


/////////////////////////////////////////////////////////////////////////////////////
// Remove sequences with seq. identity larger than seqid percent (remove the shorter of two) or coverage<cov_thr
/////////////////////////////////////////////////////////////////////////////////////
inline int Alignment::FilterForDisplay(int max_seqid, int coverage, int qid, float qsc, int N)
{
  if (par.mark) return n_display;
  char dummy[MAXSEQ];   
  int vtmp=v, seqid; 
  v=0;
  n_display=0;
  display[kss_dssp]=display[ksa_dssp]=display[kss_pred]=display[kss_conf]=0;
  for (seqid=10; n_display<N && seqid<=max_seqid; seqid++)
    {
      for (int k=0; k<N_in; k++) dummy[k]=display[k];
      n_display = Filter2(dummy,coverage,qid,qsc,20,seqid,0); 
    }
  if (n_display>N) 
    {
      for (int k=0; k<N_in; k++) dummy[k]=display[k];
      n_display = Filter2(dummy,coverage,qid,qsc,20,--(--seqid),0);
    }
  v=vtmp;
  for (int k=0; k<N_in; k++) display[k]=dummy[k];
  if (kss_dssp>=0) {display[kss_dssp]=1; n_display++;}
  if (ksa_dssp>=0) {display[ksa_dssp]=1; n_display++;}
  if (kss_pred>=0) {display[kss_pred]=1; n_display++;}
  if (kss_conf>=0) {display[kss_conf]=1; n_display++;}
  return n_display;
}

/////////////////////////////////////////////////////////////////////////////////////
// Remove sequences with seq. identity larger than seqid percent (remove the shorter of two) or coverage<cov_thr
/////////////////////////////////////////////////////////////////////////////////////
inline int Alignment::Filter(int max_seqid, int coverage, int qid, float qsc, int N)
{
  return Filter2(keep,coverage,qid,qsc,20,max_seqid,N); 
}


/////////////////////////////////////////////////////////////////////////////////////
// Select set of representative sequences in the multiple sequence alignment
// Filter criteria:
//   * Remove sequences with coverage of query less than "coverage" percent
//   * Remove sequences with sequence identity to query of less than "qid" percent 
//   * If Ndiff==0, remove sequences with seq. identity larger than seqid2(=max_seqid) percent
//   * If Ndiff>0, remove sequences with minimum-sequence-identity filter of between seqid1 
//     and seqid2 (%), where the minimum seqid threshold is determined such that, 
//     in all column blocks of at least WMIN=25 residues, at least Ndiff sequences are left.
//     This ensures that in multi-domain proteins sequences covering one domain are not
//     removed completely because sequences covering other domains are more diverse.
//
// Allways the shorter of two compared sequences is removed (=> sort sequences by length first).
// Please note: sequence identity of sequence x with y when filtering x is calculated as 
//   number of residues in sequence x that are identical to an aligned residue in y / number of residues in x
//   Example: two sequences x and y are 100% identical in their overlapping region but one overlaps by 10% of its 
//   length on the left and the other by 20% on the right. Then x has 10% seq.id with y and y has 20% seq.id. with x.
/////////////////////////////////////////////////////////////////////////////////////
int Alignment::Filter2(char keep[], int coverage, int qid, float qsc, int seqid1, int seqid2, int Ndiff)
{
  // In the beginnning, keep[k] is 1 for all regular amino acid sequences and 0 for all others (ss_conf, ss_pred,...)
  // In the end, keep[k] will be 1 for all regular representative sequences kept in the alignment, 0 for all others
  char in[MAXSEQ];          // in[k]=1: seq k has been accepted; in[k]=0: seq k has not yet been accepted at current seqid
  char inkk[MAXSEQ];        // inkk[k]=1 iff in[ksort[k]]=1 else 0;
  int idmax[MAXCOL];        // position-dependent maximum-sequence-identity threshold for filtering
  int idmaxwin[MAXCOL];     // minimum value of idmax[i-WFIL,i+WFIL]
  const int WFIL=25;        // see previous line
  int seqid_prev[MAXSEQ];   // maximum-sequence-identity threshold used in previous round of filtering (with lower seqid)
  int N[MAXCOL];            // N[i] number of already accepted sequences at position i
  int seqid;                // current maximum value for the position-dependent maximum-sequence-identity thresholds in idmax[]

  float diff_min_frac;  // minimum fraction of differing positions between sequence j and k needed to accept sequence k
  float qdiff_max_frac=0.9999-0.01*qid; // maximum allowable number of residues different from query sequence  
  int diff;                 // number of differing positions between sequences j and k (counted so far) 
  int diff_suff;            // number of differing positions between sequences j and k that would be sufficient 
  int qdiff_max;            // maximum number of residues required to be different from query
  int cov_kj;               // upper limit of number of positions where both sequence k and j have a residue
  int first_kj;             // first non-gap position in sequence j AND k
  int last_kj;              // last  non-gap position in sequence j AND k
  int kk, jj;               // indices for sequence from 1 to N_in
  int k, j;                 // kk=ksort[k], jj=ksort[j]
  int i;                    // counts residues
  int n;                    // number of sequences accepted so far

  
  // Initialize in[k]
  for (n=k=0; k<N_in; k++) if (keep[k]==2) {in[k]=2; n++;} else in[k]=0;

  // Determine first[k], last[k]?
  if (first==NULL) 
    {
      first=new(int[N_in]);// first non-gap position in sequence k
      last =new(int[N_in]);// last  non-gap position in sequence k
      for (k=0; k<N_in; k++)  // do this for ALL sequences, not only those with in[k]==1 (since in[k] may be display[k])
	{
	  for (i=1; i<=L; i++) if (X[k][i]<NAA) break;
	  first[k]=i; 
	  for (i=L; i>=1; i--) if (X[k][i]<NAA) break;
	  last[k]=i;
	}
   }
  
  // Determine number of residues nres[k]?
  if (nres==NULL) 
    {
      nres=new(int[N_in]);
      for (k=0; k<N_in; k++)  // do this for ALL sequences, not only those with in[k]==1 (since in[k] may be display[k])
	{
	  int nr=0;
	  for (i=first[k]; i<=last[k]; i++) 
	    if (X[k][i]<NAA) nr++; 
 	  nres[k]=nr;
	  //    printf("%20.20s nres=%3i  first=%3i  last=%3i\n",sname[k],nr,f,l);
	}
    }
  
  // Sort sequences according to length; afterwards, nres[ksort[kk]] is sorted by size
  if (ksort==NULL) 
    {
      ksort=new(int[N_in]);
      for (k=0; k<N_in; k++) ksort[k]=k;       
      QSortInt(nres,ksort,kfirst+1,N_in-1,-1); //Sort sequences after kfirst (query) in descending order 
    }
  for (kk=0; kk<N_in; kk++) inkk[kk]=in[ksort[kk]];

  // Initialize N[i], idmax[i], idprev[i]
  for (i=1; i<first[kfirst]; i++) N[i]=0;
  for (i=first[kfirst]; i<=last[kfirst]; i++) N[i]=1;
  for (i=last[kfirst]+1; i<=L; i++) N[i]=0; 
  for (i=1; i<=L; i++) {idmax[i]=seqid1; idmaxwin[i]=-1;}
  for (k=0; k<N_in; k++) seqid_prev[k]=-1;
  if (Ndiff<=0 || Ndiff>=N_in) {seqid1=seqid2; Ndiff=N_in;}

  // Check coverage and sim-to-query criteria for each sequence k
  for (k=0; k<N_in; k++) 
    {
      if (keep[k]==0 || keep[k]==2) continue;  // seq k not regular sequence OR is marked sequence
      if (100*nres[k]<coverage*L) {keep[k]=0; continue;} // coverage too low? => reject once and for all
      
      float qsc_sum=0.0;

      // Check if score-per-column with query is at least qsc
      if (qsc>-10) 
	{
	  float qsc_min = qsc*nres[k]; // minimum total score of seq k with query  

	  int gapq=0, gapk=0; // number of consecutive gaps in query or k'th sequence at position i
	  for (int i=first[k]; i<=last[k]; i++)
	    {
	      if (X[k][i]<20) 
		{
		  gapk=0;
		  if (X[kfirst][i]<20)
		    {
		      gapq=0;
		      qsc_sum += S[(int)X[kfirst][i]][(int)X[k][i]];
		    }
		  else if (gapq++) qsc_sum-=PLTY_GAPEXTD; else qsc_sum-=PLTY_GAPOPEN;
		}
	      else if (X[kfirst][i]<20)
		{
		  gapq=0;
		  if (gapk++) qsc_sum-=PLTY_GAPEXTD; else qsc_sum-=PLTY_GAPOPEN;
		}
	    }
// 	  printf("k=%3i qsc=%6.2f\n",k,qsc_sum);
 	  if (qsc_sum<qsc_min) {keep[k]=0; continue;} // too different from query? => reject once and for all
	}
         
      //Check if sequence similarity with query at least qid?
      if (qdiff_max_frac<0.999) 
	{
	  qdiff_max=int(qdiff_max_frac*nres[k]+0.9999); 
	  //	  printf("k=%-4i  nres=%-4i  qdiff_max=%-4i first=%-4i last=%-4i",k,nres[k],qdiff_max,first[k],last[k]);
	  diff=0;
	  for (int i=first[k]; i<=last[k]; i++)
	    // enough different residues to reject based on minimum qid with query? => break
	    if (X[k][i]<20 && X[k][i]!=X[kfirst][i] && ++diff>=qdiff_max) break;
	  //	  printf("  diff=%4i\n",diff);
	  if (diff>=qdiff_max) {keep[k]=0; continue;} // too different from query? => reject once and for all
	}
//      printf("  qsc=%6.2f     qid=%6.2f  \n",qsc_sum/nres[k],100.0*(1.0-(float)(diff)/nres[k]));
    }  


  // Successively increment idmax[i] at positons where N[i]<Ndiff
  for (seqid=seqid1; seqid<=seqid2; seqid+=1+(seqid>=50))
    {
      // Update idmax[i]
      for (i=1; i<=L; i++) if (N[i]<Ndiff) idmax[i]=seqid;
      
      // Update idmaxwin[i] as minimum of idmax[i-WFIL,i+WFIL]. If idmaxwin[] has not changed then stop
      char stop=1;
      for (i=1; i<=L; i++)
	{ 
	  int idmax_min=seqid2;
	  for (j=imax(1,imin(L-2*WFIL+1,i-WFIL)); j<=imin(L,imax(2*WFIL,i+WFIL)); j++)
	    if (idmax[j]<idmax_min) idmax_min=idmax[j];
	  if (idmax_min>idmaxwin[i]) stop=0; // idmaxwin[i] has changed => do not stop 
	  idmaxwin[i]=idmax_min;
	}
      if (stop) break;
    
//       // DEBUG
//       printf("idmax    ");
//       for (i=1; i<=L; i++) printf("%2i ",idmax[i]);
//       printf("\n");
//       printf("idmaxwin ");
//       for (i=1; i<=L; i++) printf("%2i ",idmaxwin[i]);
//       printf("\n");
//       printf("N[i]     ");
//       for (i=1; i<=L; i++) printf("%2i ",N[i]);
//       printf("\n");

      // Loop over all candidate sequences kk (-> k)
      for (kk=0; kk<N_in; kk++) 
	{
	  if (inkk[kk]) continue;   // seq k already accepted
	  k=ksort[kk];
	  if (!keep[k]) continue;  // seq k is not regular aa sequence or already suppressed by coverage or qid criterion
	  if (keep[k]==2) {in[k]=inkk[kk]=2; n++; continue;} // accept all marked sequences 

	  // Calculate max-seq-id threshold seqidk for sequence k (as maximum over idmaxwin[i])
	  if (seqid>=100) {in[k]=inkk[kk]=1; n++; continue;}
	  float seqidk=seqid1;
	  for (i=first[k]; i<=last[k]; i++)
	    if (idmaxwin[i]>seqidk) seqidk=idmaxwin[i];
	  if (seqid==seqid_prev[k]) continue; // sequence has already been rejected at this seqid threshold => reject this time
	  seqid_prev[k]=seqid;
	  diff_min_frac =0.9999-0.01*seqidk;  // min fraction of differing positions between sequence j and k needed to accept sequence k

	  // Loop over already accepted sequences 
	  for (jj=0; jj<kk; jj++)  
	    {
	      if (!inkk[jj]) continue;
	      j=ksort[jj];
	      first_kj=imax(first[k],first[j]);
	      last_kj =imin(last[k],last[j]);
	      cov_kj = last_kj-first_kj+1;
	      diff_suff=int(diff_min_frac*imin(nres[k],cov_kj));  // nres[j]>nres[k] anyway because of sorting 
	      diff=0; 
	      for (int i=first_kj; i<=last_kj; i++)
		{
		  // enough different residues to accept? => break
		  if (X[k][i]>=NAA || X[j][i]>=NAA) 
		    cov_kj--;
		  else 
		    if (X[k][i]!=X[j][i] && ++diff>=diff_suff) break; // accept (k,j)
		}
// 	      // DEBUG
//  	      printf("%20.20s with %20.20s:  diff=%i  diff_min_frac*cov_kj=%f  diff_suff=%i  nres=%i  cov_kj=%i\n",sname[k],sname[j],diff,diff_min_frac*cov_kj,diff_suff,nres[k],cov_kj);
// 	      printf("%s\n%s\n\n",seq[k],seq[j]);
	      if (float(diff)<fmin(diff_min_frac*cov_kj,diff_suff)) break; //similarity > acceptace threshold? Reject!

	    }
	  if (jj>=kk)      // did loop reach end? => accept k. Otherwise reject k (the shorter of the two)
	    {
	      in[k]=inkk[kk]=1; 
	      n++;
	      for (i=first[k]; i<=last[k]; i++) N[i]++; // update number of sequences at position i
// 	      printf("%i %20.20s accepted\n",k,sname[k]);
	    }
// 	  else 
// 	    {
// 	      printf("%20.20s rejected: too similar with seq %20.20s  diff=%i  diff_min_frac*cov_kj=%f  diff_suff=%i  nres=%i  cov_kj=%i\n",sname[k],sname[j],diff,diff_min_frac*cov_kj,diff_suff,nres[k],cov_kj);
// 	      printf("%s\n%s\n\n",seq[k],seq[j]);
// 	    } 

	} // End Loop over all candidate sequences kk

//       // DEBUG
//       printf("\n");
//       printf("seqid_prev[k]= \n");
//       for (k=0; k<N_in; k++) printf("%2i ",seqid_prev[k]);
//       printf("\n");

    } // End Loop over seqid

  if (v>=2)
    {
      printf("%i out of %i sequences passed filter (",n,N_in-N_ss);
      if (par.coverage) 
	printf("%i%% min coverage, ",coverage);
      if (qid) 
	printf("%i%% min sequence identity to query, ",qid);
      if (qsc>-10) 
	printf("%.2f bits min score per column to query, ",qsc);
      if (Ndiff<N_in && Ndiff>0) 
	printf("up to %i%% position-dependent max pairwise sequence identity)\n",seqid);
      else 
	printf("%i%% max pairwise sequence identity)\n",seqid1);
    }

  for (k=0; k<N_in; k++) keep[k]=in[k];
  return n;
}




/////////////////////////////////////////////////////////////////////////////////////
// Filter for min score per column coresc with core query profile, defined by coverage_core and qsc_core
/////////////////////////////////////////////////////////////////////////////////////
int Alignment::HomologyFilter(int coverage_core, float qsc_core, float coresc)
{
  const int seqid_core=90; //maximum sequence identity in core alignment
  const int qid_core=0;
  const int Ndiff_core=0;
  int n;
  HMM qcore;
  char* coreseq=new(char[N_in]);   // coreseq[k]=1 if sequence belongs to core of alignment (i.e. it is very similar to query)
  for (int k=0; k<N_in; k++) coreseq[k]=keep[k];   // Copy keep[] into coreseq[]

  // Remove sequences with seq. identity larger than seqid percent (remove the shorter of two)
  int v1=v; v=1; 
  n = Filter2(coreseq,coverage_core,qid_core,qsc_core,seqid_core,seqid_core,Ndiff_core); 
  v=v1;
  if (v>=2) 
    {
      printf("%i out of %i core alignment sequences passed filter (",n,N_in-N_ss);
      if (par.coverage_core) 
	printf("%i%% min coverage, ",coverage_core);
      if (qid_core) 
	printf("%i%% min sequence identity to query, ",qid_core);
      if (qsc_core>-10) 
	printf("%.2f bits min score per column to query, ",qsc_core);
      printf("%i%% max pairwise sequence identity)\n",seqid_core);
    }

  // Calculate bare AA frequencies and transition probabilities -> qcore.f[i][a], qcore.tr[i][a]
  FrequenciesAndTransitions(qcore,coreseq);

  // Add transition pseudocounts to query -> q.p[i][a] (gapd=1.0, gape=0.333, gapf=gapg=1.0, gaph=gapi=1.0, gapb=1.0
  qcore.AddTransitionPseudocounts(1.0,0.333,1.0,1.0,1.0,1.0,1.0);
  
  // Generate an amino acid frequency matrix from f[i][a] with full pseudocount admixture (tau=1) -> g[i][a] 
  qcore.PreparePseudocounts();

  // Add amino acid pseudocounts to query:  qcore.p[i][a] = (1-tau)*f[i][a] + tau*g[i][a]
  qcore.AddAminoAcidPseudocounts(2,1.5,2.0,1.0); // pcm=2, pca=1.0, pcb=2.5, pcc=1.0

  // Filter out all sequences below min score per column with qcore
  n=FilterWithCoreHMM(keep, coresc, qcore);

  if (v>=2) cout<<n<<" out of "<<N_in-N_ss<<" sequences filtered by minimum score-per-column threshold of "<<qsc_core<<"\n";
  delete[] coreseq;
  return n;
}


/////////////////////////////////////////////////////////////////////////////////////
// Filter out all sequences below a minimum score per column with profile qcore
/////////////////////////////////////////////////////////////////////////////////////
int Alignment::FilterWithCoreHMM(char in[], float coresc, HMM& qcore)
{
  int k;     // count sequences in alignment
  int i;     // column in query alignment
  int a;     // amino acid (0..19)
  int n=1;   // number of sequences that passed filter
  float** logodds=new(float*[L+1]); // log-odds ratios for HMM qcore
  char gap;  // 1: previous state in seq k was a gap  0: previous state in seq k was an amino acid
  float score; // score of sequence k aligned with qcore

  for (i=1; i<=L; i++) logodds[i]=new(float[21]);

  // Determine first[k], last[k]?
  if (first==NULL) 
    {
      first=new(int[N_in]);// first non-gap position in sequence k
      last =new(int[N_in]);// last  non-gap position in sequence k
      for (k=0; k<N_in; k++)  // do this for ALL sequences, not only those with in[k]==1 (since in[k] may be display[k])
	{
	  for (i=1; i<=L; i++) if (X[k][i]<NAA) break;
	  first[k]=i; 
	  for (i=L; i>=1; i--) if (X[k][i]<NAA) break;
	  last[k]=i;
	}
   }

  // Determine number of residues nres[k]?
  if (nres==NULL) 
    {
      nres=new(int[N_in]);
      for (k=0; k<N_in; k++)  // do this for ALL sequences, not only those with in[k]==1 (since in[k] may be display[k])
	{
	  int nr=0;
	  for (i=first[k]; i<=last[k]; i++) 
	    if (X[k][i]<NAA) nr++; 
 	  nres[k]=nr;
	  //    printf("%20.20s nres=%3i  first=%3i  last=%3i\n",sname[k],nr,f,l);
	}
    }
  
  // Precalculate the log-odds for qcore 
  for (i=1; i<=L; i++)
    {
      for (a=0; a<NAA; a++)
	logodds[i][a]=fast_log2(qcore.p[i][a]/pb[a]);
      logodds[i][ANY]=-0.5;  // half a bit penalty for X

//       printf("         A     R     N     D     C     Q     E     G     H     I     L     K     M     F     P     S     T     W     Y     V\n");
//       printf("%6i ",i);
//       for (a=0; a<20; ++a) fprintf(stdout,"%5.1f ",100*qcore.f[i][a]); 
//       printf("\n");
//       printf("       ");
//       for (a=0; a<20; ++a) fprintf(stdout,"%5.1f ",100*qcore.g[i][a]); 
//       printf("\n");
//       printf("       ");
//       for (a=0; a<20; ++a) fprintf(stdout,"%5.1f ",100*qcore.p[i][a]); 
//       printf("\n");
//       printf("       ");
//       for (a=0; a<20; ++a) fprintf(stdout,"%5.1f ",100*pb[a]); 
//       printf("\n");
//       printf("       ");
//       for (a=0; a<20; ++a) fprintf(stdout,"%5.2f ",fast_log2(qcore.p[i][a]/pb[a])); 
//       printf("\n");
   }

  // Main loop: test all sequences k
  for (k=kfirst+1; k<N_in; k++)   
    {
      if (!in[k]) continue;  // if in[k]==0 sequence k will be suppressed directly

      float score_M=0.0;
      float score_prev=0.0;

      // Calculate score of sequence k with core HMM 
      score=0; gap=0;
      for (i=first[k]; i<=last[k]; i++)
	{
	  score_M=0.0;
	  if (X[k][i]<=ANY)       // current state is Match
	    {
	      score_M=logodds[i][ (int)X[k][i]];
	      score+=logodds[i][ (int)X[k][i]];
	      if (gap) score+=qcore.tr[i][D2M]; else score+=qcore.tr[i][M2M];
	      gap=0;
	    }
	  else if (X[k][i]==GAP) // current state is Delete (ignore ENDGAPs)
	    {
	      if (gap) score+=qcore.tr[i][D2D]; else score+=qcore.tr[i][M2D];
	      gap=1;
	    }
	  if (I[k][i]) score+=qcore.tr[i][M2I]+(I[k][i]-1)*qcore.tr[i][I2I]+qcore.tr[i][I2M];
// 	  if (k==2) printf("i=%3i %c:%c   score_M=%6.2f   score=%6.2f  score_sum=%6.2f \n",i,i2aa(X[kfirst][i]),i2aa(X[k][i]),score_M,score-score_prev,score);
	  score_prev=score;
	}

      printf("k=%3i score=%6.2f\n",k,score);
      if (score<nres[k]*coresc) in[k]=0; else n++;// reject sequence k?
    }
  for (i=1; i<=L; i++) delete[] logodds[i];
  delete[] logodds;
  return n;
}



/////////////////////////////////////////////////////////////////////////////////////
// Calculate AA frequencies q.p[i][a] and transition probabilities q.tr[i][a] from alignment
/////////////////////////////////////////////////////////////////////////////////////
void Alignment::FrequenciesAndTransitions(HMM& q, char* in)
{
  int k;                // index of sequence
  int i;                // position in alignment
  int a;                // amino acid (0..19)
  int ni[NAA+3];        // number of times amino acid a occurs at position i
  int naa;              // number of different amino acids

  if (v>=3) 
     cout<<"Calculating position-dependent weights on subalignments\n";

  if (in==NULL) in=keep;  // what's this good for?

  for (k=0; k<N_in; k++) wg[k]=0.0; // initialized wg[k]

  for (i=1; i<=L; i++) // for all positions i in alignment
    {
      for (a=0; a<20; a++) ni[a]=0;
      for (k=0; k<N_in; k++) if (in[k]) ni[ (int)X[k][i]]++;
      naa=0; for (a=0; a<20; a++) if(ni[a]) naa++;
      if (!naa) naa=1; //naa=0 when column consists of only gaps and Xs (=ANY) 
      for (k=0; k<N_in; k++) 
	if (in[k] && X[k][i]<20) 
 	  wg[k]+=1.0/float(ni[ (int)X[k][i]]*naa*nres[k]);
      // ensure that each residue of a short sequence contributes as much as a residue of a long sequence:
      // contribution is proportional to one over sequence length nres[k].
    }
  NormalizeTo1(wg,N_in);
  
  // Do pos-specific sequence weighting and calculate amino acid frequencies and transitions
  for (k=0; k<N_in; k++) X[k][0]=ENDGAP;    // make sure that sequences ENTER subalignment j for j=1
  for (k=0; k<N_in; k++) X[k][L+1]=ENDGAP;  // does it have an influence?
  Amino_acid_frequencies_and_transitions_from_M_state(q,in); // use subalignments of seqs with residue in i
  if (v>=3) 
    {
      printf("\nMatches:\n");
      printf("col     Neff ncols  nseqs\n");
      for (i=1; i<=L; i++)
	printf("%3i    %5.2f   %3i    %3i\n",i,q.Neff_M[i],ncols[i],nseqs[i]);
    }
  Transitions_from_I_state(q,in);                            // use subalignments of seqs with insert in i
  if (v>=3) 
    {
      printf("\nInserts:\n");
      printf("col     Neff ncols  nseqs\n");
      for (i=1; i<=L; i++)
	printf("%3i    %5.2f   %3i    %3i\n",i,q.Neff_I[i],ncols[i],nseqs[i]);
    }
  Transitions_from_D_state(q,in);                            // use subalignments of seqs with delete in i
  if (v>=3) 
    {
      printf("\nDeletes:\n");
      printf("col     Neff ncols  nseqs\n");
      for (i=1; i<=L; i++)
	printf("%3i    %5.2f   %3i    %3i\n",i,q.Neff_D[i],ncols[i],nseqs[i]);
    }

  // Copy column information into HMM q
  q.L=L;
  q.N_in=N_in;
  q.N_filtered=N_filtered;
  for (i=1; i<=L; i++) q.l[i]=l[i];

  // Set names in HMM q
  if (strlen(q.name)==0) strcpy(q.name,name);
  if (strlen(q.longname)==0) strcpy(q.longname,longname);
  if (strlen(q.fam)==0) strcpy(q.fam,fam);
  strcpy(q.fam,fam);
  ScopID(q.cl,q.fold,q.sfam,q.fam); // derive superfamily, fold and class code from family name
  strcpy(q.file,file);              // Store basename of alignment file name in q.file

  // Copy sequences to be displayed into HMM
  q.nss_dssp=q.nsa_dssp=q.nss_pred=q.nss_conf=q.nfirst=-1;
  int n=0;
  if (kss_dssp>=0) q.nss_dssp=n++; // copy dssp sequence?
  if (ksa_dssp>=0) q.nsa_dssp=n++; // copy dssp sequence?
  if (kss_pred>=0) q.nss_pred=n++; // copy psipred sequence?
  if (kss_conf>=0) q.nss_conf=n++; // copy confidence value sequence?

  // Calculate consensus sequence?
  if (par.showcons || par.cons)    
    {
      float maxw;
      int maxa;
      float wa[NAA];
      if (par.showcons) 
	{
	  // Calculate consensus amino acids
	  q.ncons=n++;
	  q.sname[q.ncons]=new(char[10]);
	  if (!q.sname[q.ncons]) {MemoryError("array of names for displayed sequences");}
	  strcpy(q.sname[q.ncons],"Consensus");
	  q.seq[q.ncons]=new(char[L+2]);
	  if (!q.seq[q.ncons]) {MemoryError("array of names for displayed sequences");} 
	}
      if (par.cons) 
	{
	  // Calculate consensus amino acids
	  q.nfirst=n++; kfirst=-1;
	  q.sname[q.nfirst]=new(char[strlen(name)+11]);
	  if (!q.sname[q.nfirst]) {MemoryError("array of names for displayed sequences");}
	  strcpy(q.sname[q.nfirst],name);
	  strcat(q.sname[q.nfirst],"_consensus");
	  q.seq[q.nfirst]=new(char[L+2]);
	  if (!q.seq[q.nfirst]) {MemoryError("array of names for displayed sequences");} 
	}
      for (i=1; i<=L; i++) 
	{
	  maxw=0.0; maxa=0;
	  for (a=0; a<20; a++) wa[a]=0;
	  for (a=0; a<20; a++) 
	    for (int b=0; b<20; b++) wa[a]+=q.f[i][b]*Sim[a][b]*Sim[a][b];
	  
	  for (a=0; a<20; a++)
	    if (wa[a]>maxw) {maxw = wa[a]; maxa = a;}
	  
	  if (par.showcons) 
	    {
	      maxw*=q.Neff_M[i]/(q.Neff_HMM+1);  // columns with many gaps don't get consensus symbol
	      if (maxw>0.6) q.seq[q.ncons][i] = uprchr(i2aa(maxa));
	      else if (maxw>0.4) q.seq[q.ncons][i] = lwrchr(i2aa(maxa));
	      else q.seq[q.ncons][i] = 'x';
	    }
	  if (par.cons) q.seq[q.nfirst][i] = uprchr(i2aa(maxa));
	}
      if (par.showcons) 
	{
	  q.seq[q.ncons][0]='-';
	  q.seq[q.ncons][L+1]='\0';
	}
      if (par.cons) 
	{
	  q.seq[q.nfirst][0]='-';
	  q.seq[q.nfirst][L+1]='\0';
	}
    }
  
  // Copy sequences to be displayed from alignment to HMM
  for (k=0; k<N_in; k++) 
    {
      int nn;
      if (display[k])
	{
	  if (n>=MAXSEQALI) {
	    if (par.mark) cerr<<"WARNING: maximum number "<<MAXSEQALI<<" of sequences for display of alignment exceeded\n"; 
	    break;
	  }
	  if (k==kss_dssp)      nn=q.nss_dssp; // copy dssp sequence to nss_dssp
	  else if (k==ksa_dssp) nn=q.nsa_dssp; 
	  else if (k==kss_pred) nn=q.nss_pred; 
	  else if (k==kss_conf) nn=q.nss_conf; 
	  else if (k==kfirst)   nn=q.nfirst=n++; 
	  else nn=n++;
//	  strcut(sname[k],"  "); // delete rest of name line beginning with two spaces "  " // Why this?? Problem for pdb seqs without chain 
	  q.sname[nn]=new(char[strlen(sname[k])+1]);
	  if (!q.sname[nn]) {MemoryError("array of names for displayed sequences");}	  
	  strcpy(q.sname[nn],sname[k]);
	  q.seq[nn]=new(char[strlen(seq[k])+1]);
	  if (!q.seq[nn]) {MemoryError("array of names for displayed sequences");}	  
	  strcpy(q.seq[nn],seq[k]);
	}
    }
  q.n_display=n; // how many sequences to be displayed in alignments?
  
  // Copy secondary structure information into HMM
  if (kss_dssp>=0) 
    for (i=1; i<=L; i++) q.ss_dssp[i]=X[kss_dssp][i];
  if (ksa_dssp>=0) 
    for (i=1; i<=L; i++) q.sa_dssp[i]=X[ksa_dssp][i];
  if (kss_pred>=0) 
    {
      for (i=1; i<=L; i++) q.ss_pred[i]=X[kss_pred][i];
      if (kss_conf>=0) 
	for (i=1; i<=L; i++) q.ss_conf[i]=X[kss_conf][i];
      else
	for (i=1; i<=L; i++) q.ss_conf[i]=5;
    }

  q.lamda=0.0;
  q.mu=0.0;

  // Debug: print occurence of amino acids for each position i
  if (v>=2) printf("Effective number of sequences exp(entropy) = %-4.1f\n",q.Neff_HMM); //PRINT
  if (v>=3) 
    {
      cout<<"\nMatr: ";
      for (a=0; a<20; a++) printf("%4.1f ",100*pb[a]);
      cout<<"\nAmino acid frequencies without pseudocounts:\n";
      cout<<"         A    R    N    D    C    Q    E    G    H    I    L    K    M    F    P    S    T    W    Y    V\n";  
      for (i=1; i<=L; i++)
	{
	  printf("%3i:  ",i);
	  for (a=0; a<20; a++)  printf("%4.0f ",100*q.f[i][a]);
	  cout<<endl;
	}
      cout<<"\n";
 
      printf("\nListing transition probabilities without pseudocounts:\n");
      printf("   i    M->M   M->I   M->D   I->M   I->I   D->M   D->D  Neff_M Neff_I Neff_D\n");
      for (i=0; i<=L; i++) 
	{
	  printf("%4i  %6.3f %6.3f %6.3f ",i,pow(2.0,q.tr[i][M2M]),pow(2.0,q.tr[i][M2I]),pow(2.0,q.tr[i][M2D]));
	  printf("%6.3f %6.3f ",pow(2.0,q.tr[i][I2M]),pow(2.0,q.tr[i][I2I]));
	  printf("%6.3f %6.3f  ",pow(2.0,q.tr[i][D2M]),pow(2.0,q.tr[i][D2D]));
	  printf("%6.3f %6.3f %6.3f\n",q.Neff_M[i],q.Neff_I[i],q.Neff_D[i]);
	}
    }
  
  return;
}


/////////////////////////////////////////////////////////////////////////////////////
// Calculate freqs q.f[i][a] and transitions q.tr[i][a] (a=MM,MI,MD) with pos-specific subalignments
/////////////////////////////////////////////////////////////////////////////////////
void Alignment::Amino_acid_frequencies_and_transitions_from_M_state(HMM& q, char* in)
{
  // Calculate position-dependent weights wi[k] for each i.
  // For calculation of weights in column i use sub-alignment 
  // over sequences which have a *residue* in column i (no gap, no end gap) 
  // and over columns where non of these sequences has an end gap. 
  // This is done by updating the arrays n[j][a] at each step i-1->i while letting i run from 1 to L. 
  // n[j][a] = number of occurences of amino acid a at column j of the subalignment, 
  //        => only columns with n[j][ENDGAP]=0 are contained in the subalignment!
  // If no sequences enter or leave the subalignment at the step i-1 -> i (i.e. change=0)
  // then the old values wi[k], Neff[i-1], and ncol are used for the new position i.
  // Index a can be an amino acid (0-19), ANY=20, GAP=21, or ENDGAP=22
  
  int k;                      // index of sequence
  int i,j;                    // position in alignment
  int a;                      // amino acid (0..19)
  int naa;                    // number of different amino acids
  int** n;                    // n[j][a] = number of seq's with some residue at column i AND a at position j
  float wi[MAXSEQ];           // weight of sequence k in column i, calculated from subalignment i
  float Neff[MAXRES];         // diversity of subalignment i
  float Naa;                  // NEW
  int nseqi=0;                // number of sequences in subalignment i 
  int ncol=0;                 // number of columns j that contribute to Neff[i]
  char change;                // has the set of sequences in subalignment changed? 0:no  1:yes
  float fj[NAA+3];            // to calculate entropy
  float sum;                

  // Initialization
  n = new(int*[L+2]);
  for (j=1; j<=L; j++) n[j]=new(int[NAA+3]);
  Neff[0]=0.0;   // if the first column has no residues (i.e. change==0), Neff[i]=Neff[i-1]=Neff[0]
  for (j=1; j<=L; j++) 
    for (a=0; a<NAA+3; a++) n[j][a]=0;

  //////////////////////////////////////////////////////////////////////////////////////////////
  // Main loop through alignment columns
  for (i=1; i<=L; i++) // Calculate wi[k] at position i as well as Neff[i]
    {
      change=0;
      // Check all sequences k and update n[j][a] and ri[j] if necessary
      for (k=0; k<N_in; k++) 
	{
	  if (!in[k]) continue;
	  if (X[k][i-1]>=ANY && X[k][i]<ANY)
	    { // ... if sequence k was NOT included in i-1 and has to be included for column i
	      change=1;
	      nseqi++;
	      for (int j=1; j<=L; j++) n[j][ (int)X[k][j]]++; 
	    }
	  else if (X[k][i-1]<ANY && X[k][i]>=ANY)
	    { // ... if sequence k WAS included in i-1 and has to be thrown out for column i
	      change=1;
	      nseqi--;
	      for (int j=1; j<=L; j++) n[j][ (int)X[k][j]]--; 
	    }
	} //end for (k)
      nseqs[i]=nseqi;
      
      // If subalignment changed: update weights wi[k] and Neff[i]
      if (change) 
	{
	  // Initialize weights and numbers of residues for subalignment i
	  ncol=0; 
	  for (k=0; k<N_in; k++) wi[k]=0.0;
	  Naa = 0.0; // NEW

	  // sum wi[k] over all columns j and sequences k of subalignment
	  for (j=1; j<=L; j++) 
	    {
	      //  do at least a fraction MAXENDGAPFRAC of sequences in subalignment contain an end gap in j?
 	      if (n[j][ENDGAP]>MAXENDGAPFRAC*nseqi) continue; 
	      naa=0; for (a=0; a<20; a++) if(n[j][a]) naa++;
	      if (naa==0) continue;
	      ncol++; 
	      Naa += float(naa); // NEW
	      for (k=0; k<N_in; k++) 
		{
		  if (in[k] && X[k][i]<ANY && X[k][j]<ANY) 
		    {
// 		      if (!n[j][ (int)X[k][j]]) {fprintf(stderr,"Error: Mi=%i: n[%i][X[%i]]=0! (X[%i]=%i)\n",i,j,k,k,X[k][j]);} 
		      wi[k]+=1.0/float(n[j][ (int)X[k][j] ]*naa);
		    }
		}
	    }
	  Naa/=(float)(ncol); // NEW

// NEW	  // Check whether number of columns in subalignment is sufficient
// 	  if (ncol<NCOLMIN) 
// 	    // Take global weights
// 	    for (k=0; k<N_in; k++)
// 	      if(in[k] && X[k][i]<ANY) wi[k]=wg[k]; else wi[k]=0.0;

// 	  // Calculate Neff[i]
// 	  Neff[i]=0.0;
// 	  for (j=1; j<=L; j++) 
// 	    {
// 	      //  do at least a fraction MAXENDGAPFRA of sequences in subalignment contain an end gap in j?
// 	      if (n[j][ENDGAP]>MAXENDGAPFRAC*nseqi) continue;  
// 	      for (a=0; a<20; a++) fj[a]=0;  
// 	      for (k=0; k<N_in; k++)
// 		if (in[k] && X[k][i]<ANY && X[k][j]<ANY)
// 		  fj[ (int)X[k][j] ]+=wi[k];        
// 	      NormalizeTo1(fj,NAA);
// 	      for (a=0; a<20; a++) 
// 		if (fj[a]>1E-10) Neff[i]-=fj[a]*fast_log2(fj[a]); 
// 	    }
// NEW	  if (ncol>0) Neff[i]=pow(2.0,Neff[i]/ncol); else Neff[i]=1.0;

	  Neff[i]=Naa; // NEW
	}

      else //no update was necessary; copy values for i-1
	{
	  Neff[i]=Neff[i-1];
	}
      ncols[i]=ncol;  /////DEBUG
      
      // Calculate amino acid frequencies q.f[i][a] from weights wi[k]
      for (a=0; a<20; a++) q.f[i][a]=0;
      for (k=0; k<N_in; k++) if (in[k]) q.f[i][ (int)X[k][i] ]+=wi[k];
      NormalizeTo1(q.f[i],NAA,pb);

      // Calculate transition probabilities from M state
      q.tr[i][M2M]=q.tr[i][M2D]=q.tr[i][M2I]=0.0;
      for (k=0; k<N_in; k++) //for all sequences
	{
	  if (!in[k]) continue;  
	  //if input alignment is local ignore transitions from and to end gaps
	  if (X[k][i]<ANY)            //current state is M
	    {
	      if (I[k][i])             //next state is I
		  q.tr[i][M2I]+=wi[k];
	      else if (X[k][i+1]<=ANY) //next state is M
		  q.tr[i][M2M]+=wi[k];
	      else if (X[k][i+1]==GAP) //next state is D
		  q.tr[i][M2D]+=wi[k];
	    }
	} // end for(k)
      // Normalize and take log
      sum = q.tr[i][M2M]+q.tr[i][M2I]+q.tr[i][M2D]+FLT_MIN;     
      q.tr[i][M2M]=log2(q.tr[i][M2M]/sum);
      q.tr[i][M2I]=log2(q.tr[i][M2I]/sum);
      q.tr[i][M2D]=log2(q.tr[i][M2D]/sum);
      
//       for (k=0; k<N_in; k++) if (in[k]) w[k][i]=wi[k];  
   } 
  // end loop through alignment columns i
  //////////////////////////////////////////////////////////////////////////////////////////////
  

  q.tr[0][M2M]=0;
  q.tr[0][M2I]=-100000;
  q.tr[0][M2D]=-100000;
  q.tr[L][M2M]=0;
  q.tr[L][M2I]=-100000;
  q.tr[L][M2D]=-100000;
  q.Neff_M[0]=99.999; // Neff_av[0] is used for calculation of transition pseudocounts for the start state

  // Set emission probabilities of zero'th (begin) state and L+1st (end) state to background probabilities
  for (a=0; a<20; a++) q.f[0][a]=q.f[L+1][a]=pb[a];

  // Assign Neff_M[i] and calculate average over alignment, Neff_M[0] 
  q.Neff_HMM=0.0;
  for (i=1; i<=L; i++) 
    {
      q.Neff_HMM+=Neff[i];
      q.Neff_M[i]=Neff[i];
    }
  q.Neff_HMM/=L;       

  // delete n[][]
  for (j=1; j<=L; j++) delete[](n[j]);
  delete[](n);
  return;
}


/////////////////////////////////////////////////////////////////////////////////////
// Calculate transitions q.tr[i][a] (a=DM,DD) with pos-specific subalignments
/////////////////////////////////////////////////////////////////////////////////////
void Alignment::Transitions_from_D_state(HMM& q, char* in)
{
  // Calculate position-dependent weights wi[k] for each i.
  // For calculation of weights in column i use sub-alignment 
  // over sequences which have a DELETE in column i 
  // and over columns where non of these sequences has an end gap. 
  // This is done by updating the arrays n[j][a] and rj[j] at each step i-1->i while letting i run from 1 to L. 
  // n[j][a] = number of occurences of index a at column j of the subalignment, 
  //        => only columns with n[j][ENDGAP]=0 are contained in the subalignment!
  // If no sequences enter or leave the subalignment at the step i-1 -> i (i.e. change=0)
  // then the old values wi[k], Neff[i-1], and ncol are used for the new position i.
  // Index a can be an amino acid (0-19), ANY=20, GAP=21, or ENDGAP=22
  
  int k;                      // index of sequence
  int i,j;                    // position in alignment
  int a;                      // amino acid (0..19)
  int naa;                    // number of different amino acids
  int** n;                    // n[j][a] = number of seq's with some residue at column i AND a at position j
  float wi[MAXSEQ];           // weight of sequence k in column i, calculated from subalignment i
  float Neff[MAXRES];         // diversity of subalignment i 
  int nseqi=0;                // number of sequences in subalignment i (for DEBUGGING)
  int ncol=0;                 // number of columns j that contribute to Neff[i]
  char change;                // has the set of sequences in subalignment changed? 0:no  1:yes
  float fj[NAA+3];            // to calculate entropy
  float sum;                

  // Initialization
  n = new(int*[L+2]);
  for (j=1; j<=L; j++) n[j]=new(int[NAA+3]);
  for (j=1; j<=L; j++) 
    for (a=0; a<NAA+3; a++) n[j][a]=0;
  
  //////////////////////////////////////////////////////////////////////////////////////////////
  // Main loop through alignment columns
  for (i=1; i<=L; i++) // Calculate wi[k] at position i as well as Neff[i]
    {
      change=0;
      // Check all sequences k and update n[j][a] and ri[j] if necessary
      for (k=0; k<N_in; k++) 
	{
	  if (!in[k]) continue;
	  if (X[k][i-1]!=GAP && X[k][i]==GAP)
	    { // ... if sequence k was NOT included in i-1 and has to be included for column i
	      change=1;
	      nseqi++;
	      for (int j=1; j<=L; j++) n[j][ (int)X[k][j]]++; 
	    }
	  else if (X[k][i-1]==GAP && X[k][i]!=GAP)
	    { // ... if sequence k WAS included in i-1 and has to be thrown out for column i
	      change=1;
	      nseqi--;
	      for (int j=1; j<=L; j++) n[j][ (int)X[k][j]]--; 
	    }
	} //end for (k)
      nseqs[i]=nseqi;
      
      // If there is no sequence in subalignment j ...
      if (nseqi==0) 
	{
	  ncols[i]=ncol=0;
	  Neff[i]=0.0;  // effective number of sequences = 0!
	  q.tr[i][D2M]=-100000;
	  q.tr[i][D2D]=-100000;
	  continue;
	}



      // If subalignment changed: update weights wi[k] and Neff[i]
      if (change) 
	{
	  // Initialize weights and numbers of residues for subalignment i
	  ncol=0; 
	  for (k=0; k<N_in; k++) wi[k]=0.0;

	  // sum wg[k][i] over all columns j and sequences k of subalignment
	  for (j=1; j<=L; j++) 
	    {
	      if (n[j][ENDGAP]>MAXENDGAPFRAC*nseqi) continue;  
	      naa=0; for (a=0; a<20; a++) if(n[j][a]) naa++;
	      if (naa==0) continue;
	      ncol++; 
	      for (k=0; k<N_in; k++) 
		{
		  if (in[k] && X[k][i]==GAP && X[k][j]<ANY) 
		    {
		      if (!n[j][ (int)X[k][j]]) {fprintf(stderr,"Error: Di=%i: n[%i][X[%i]]=0! (X[%i]=%i)\n",i,j,k,k,X[k][j]);} 
		      wi[k]+=1.0/float(n[j][ (int)X[k][j] ]*naa);
		    }
		}
	    }
	  
	  // Check whether number of columns in subalignment is sufficient
	  if (ncol<NCOLMIN) 
	    // Take global weights
	    for (k=0; k<N_in; k++) 
	      if(in[k] && X[k][i]==GAP) wi[k]=wg[k]; else wi[k]=0.0;

	  // Calculate Neff[i]
	  Neff[i]=0.0;
	  for (j=1; j<=L; j++) 
	    {
	      if (n[j][ENDGAP]>MAXENDGAPFRAC*nseqi) continue;  
	      for (a=0; a<20; a++) fj[a]=0;  
	      for (k=0; k<N_in; k++)
		if (in[k] && X[k][i]==GAP && X[k][j]<ANY)
		  fj[ (int)X[k][j] ]+=wi[k];        
	      NormalizeTo1(fj,NAA);
	      for (a=0; a<20; a++) 
		if (fj[a]>1E-10) Neff[i]-=fj[a]*fast_log2(fj[a]); 
	    }
	  if (ncol>0) Neff[i]=pow(2.0,Neff[i]/ncol); else Neff[i]=1.0;

	}

      else //no update was necessary; copy values for i-1
	{
	  Neff[i]=Neff[i-1];
	}

      ncols[i]=ncol;  /////DEBUG
      
      // Calculate transition probabilities from D state
      q.tr[i][D2M]=q.tr[i][D2D]=0.0;
      for (k=0; k<N_in; k++) //for all sequences
	{
	  if (!in[k]) continue;  
	  //if input alignment is local ignore transitions from and to end gaps
	  if (X[k][i]==GAP)            //current state is D
	    {
	      if (X[k][i+1]==GAP)      //next state is D
		  q.tr[i][D2D]+=wi[k];
	      else if (X[k][i+1]<=ANY) //next state is M
		  q.tr[i][D2M]+=wi[k];
	    }
	} // end for(k)
      // Normalize and take log
      sum = q.tr[i][D2M]+q.tr[i][D2D];     
      q.tr[i][D2M]=log2(q.tr[i][D2M]/sum);
      q.tr[i][D2D]=log2(q.tr[i][D2D]/sum);
      
   } 
  // end loop through alignment columns i
  //////////////////////////////////////////////////////////////////////////////////////////////
  
  q.tr[0][D2M]=0;
  q.tr[0][D2D]=-100000;
  q.Neff_D[0]=99.999;
  
  // Assign Neff_D[i]
  for (i=1; i<=L; i++)  
      q.Neff_D[i]=Neff[i];

  // delete n[][]
  for (j=1; j<=L; j++) delete[](n[j]);
  delete[](n);
  return;
}


/////////////////////////////////////////////////////////////////////////////////////
// Calculate transitions q.tr[i][a] (a=DM,DD) with pos-specific subalignments
/////////////////////////////////////////////////////////////////////////////////////
void Alignment::Transitions_from_I_state(HMM& q, char* in)
{
  // Calculate position-dependent weights wi[k] for each i.
  // For calculation of weights in column i use sub-alignment 
  // over sequences which have a INSERT in column i 
  // and over columns where non of these sequences has an end gap. 
  // This is done by calculating the arrays n[j][a] and rj[j] at each step i-1->i while letting i run from 1 to L. 
  // n[j][a] = number of occurences of amino acid a at column j of the subalignment, 
  //        => only columns with n[j][ENDGAP]=0 are contained in the subalignment!
  // If no sequences enter or leave the subalignment at the step i-1 -> i (i.e. change=0)
  // then the old values wi[k], Neff[i-1], and ncol are used for the new position i.
  // Index a can be an amino acid (0-19), ANY=20, GAP=21, or ENDGAP=22
  
  int k;                      // index of sequence
  int i,j;                    // position in alignment
  int a;                      // amino acid (0..19)
  int naa;                    // number of different amino acids
  int** n;                    // n[j][a] = number of seq's with some residue at column i AND a at position j
  float wi[MAXSEQ];           // weight of sequence k in column i, calculated from subalignment i
  float Neff[MAXRES];         // diversity of subalignment i
  int nseqi;                  // number of sequences in subalignment i 
  int ncol;                   // number of columns j that contribute to Neff[i]
  float fj[NAA+3];            // to calculate entropy
  float sum;                

  // Initialization
  n = new(int*[L+2]);
  for (j=1; j<=L; j++) n[j]=new(int[NAA+3]);
  
  //////////////////////////////////////////////////////////////////////////////////////////////
  // Main loop through alignment columns
  for (i=1; i<=L; i++) // Calculate wi[k] at position i as well as Neff[i]
    {

      // Initialization of n[j][a]
      for (j=1; j<=L; j++) 
	for (a=0; a<NAA+3; a++) n[j][a]=0;
  
      // Calculate n[j][a] and ri[j] 
      nseqi=0;
      for (k=0; k<N_in; k++) 
	{
	  if (in[k] && I[k][i]>0)
	    {
	      nseqi++;
	      for (int j=1; j<=L; j++) n[j][ (int)X[k][j]]++; 
	    }
	} //end for (k)
      nseqs[i]=nseqi;
      
      // If there is no sequence in subalignment j ...
      if (nseqi==0) 
	{
	  ncols[i]=ncol=0;
	  Neff[i]=0.0;  // effective number of sequence = 0!
	  q.tr[i][I2M]=-100000;
	  q.tr[i][I2I]=-100000;
	  continue;
	}

      // update weights wi[k] and Neff[i]
//    if (1) 
	{
	  // Initialize weights and numbers of residues for subalignment i
	  ncol=0; 
	  for (k=0; k<N_in; k++) wi[k]=0.0;

	  // sum wg[k][i] over all columns j and sequences k of subalignment
	  for (j=1; j<=L; j++) 
	    {
	      if (n[j][ENDGAP]>MAXENDGAPFRAC*nseqi) continue;  
	      naa=0; for (a=0; a<20; a++) if(n[j][a]) naa++;
	      if (naa==0) continue;
	      ncol++; 
	      for (k=0; k<N_in; k++) 
		{
		  if (in[k] && I[k][i]>0 && X[k][j]<ANY) 
		    {
		      if (!n[j][ (int)X[k][j]]) {fprintf(stderr,"Error: Ii=%i: n[%i][X[%i]]=0! (X[%i]=%i)\n",i,j,k,k,X[k][j]);} 
		      wi[k]+=1.0/float(n[j][ (int)X[k][j] ]*naa);
		    }
		}
	    }
	  
	  // Check whether number of columns in subalignment is sufficient
	  if (ncol>=NCOLMIN) 
	    // Take global weights
	    for (k=0; k<N_in; k++)
	      if(in[k] && I[k][i]>0) wi[k]=wg[k]; else wi[k]=0.0;

	  // Calculate Neff[i]
	  Neff[i]=0.0;
	  for (j=1; j<=L; j++) 
	    {
	      if (n[j][ENDGAP]>MAXENDGAPFRAC*nseqi) continue;  
	      for (a=0; a<20; a++) fj[a]=0;  
	      for (k=0; k<N_in; k++)
		if (in[k] && I[k][i]>0 && X[k][j]<ANY)
		  fj[ (int)X[k][j] ]+=wi[k];        
	      NormalizeTo1(fj,NAA);
	      for (a=0; a<20; a++) 
		if (fj[a]>1E-10) Neff[i]-=fj[a]*fast_log2(fj[a]); 
	    }
	  if (ncol>0) Neff[i]=pow(2.0,Neff[i]/ncol); else Neff[i]=1.0;

	}

      ncols[i]=ncol;  /////DEBUG
      
      // Calculate transition probabilities from I state
      q.tr[i][I2M]=q.tr[i][I2I]=0.0;
      for (k=0; k<N_in; k++) //for all sequences
	{
	  if (in[k] && I[k][i]>0)  //current state is I
	    {
	      q.tr[i][I2M]+=wi[k];
	      q.tr[i][I2I]+=wi[k]*(I[k][i]-1);
	    }
	} // end for(k)
      // Normalize and take log
      sum = q.tr[i][I2M]+q.tr[i][I2I];     
      q.tr[i][I2M]=log2(q.tr[i][I2M]/sum);
      q.tr[i][I2I]=log2(q.tr[i][I2I]/sum);
      
   } 
  // end loop through alignment columns i
  //////////////////////////////////////////////////////////////////////////////////////////////

  q.tr[0][I2M]=0;
  q.tr[0][I2I]=-100000;
  q.tr[L][I2M]=0;
  q.tr[L][I2I]=-100000;
  q.Neff_I[0]=99.999;

  // Assign Neff_M[i] and calculate average over alignment, Neff_M[0] 
  for (i=1; i<=L; i++)   // Calculate wi[k] at position i as well as Neff[i] and Neff[i]
      q.Neff_I[i]=Neff[i];

  // delete n[][]
  for (j=1; j<=L; j++) delete[](n[j]);
  delete[](n);

  return;
}


// Write alignment without insert states (lower case) to alignment file?
void Alignment::WriteWithoutInsertsToFile(char* alnfile)
{
  if (v>=2) cout<<"Writing alignment to "<<alnfile<<"\n";
  FILE* alnf;
  if (!par.append) alnf = fopen(alnfile,"w"); else alnf = fopen(alnfile,"a");
  if (!alnf) OpenFileError(alnfile);
  // If alignment name is different from that of query: write name into commentary line
  if (strncmp(longname,sname[kfirst],DESCLEN-1)) fprintf(alnf,"#%s\n",longname);
  if (v>=2) cout<<"Writing alignment to "<<alnfile<<"\n";
  for (int k=0; k<N_in; k++)
    if (keep[k] || display[k]==2) // print if either in profile (keep[k]>0) or display is obligatory (display[k]==2)
      {
	fprintf(alnf,">%s\n",sname[k]);
	for (int i=1; i<=L; i++) fprintf(alnf,"%c",i2aa(X[k][i]));
	fprintf(alnf,"\n");
      }
  fclose(alnf);  
}

// Write stored,filtered sequences WITH insert states (lower case) to alignment file?
void Alignment::WriteToFile(char* alnfile, char* format)
{
  FILE* alnf;
  if (!par.append) alnf = fopen(alnfile,"w"); else alnf = fopen(alnfile,"a");
  if (!alnf) OpenFileError(alnfile);
  // If alignment name is different from that of query: write name into commentary line
  if (strncmp(longname,sname[kfirst],DESCLEN-1)) fprintf(alnf,"#%s\n",longname);
  if (!format || !strcmp(format,"a3m"))
    {
      if (v>=2) cout<<"Writing A3M alignment to "<<alnfile<<"\n";
      for (int k=0; k<N_in; k++)
	if (keep[k] || display[k]==2) // print if either in profile (keep[k]>0) or display obligatory (display[k]==2)
	  fprintf(alnf,">%s\n%s\n",sname[k],seq[k]+1);
    }
  else // PSI-BLAST format 
    {  
      if (v>=2) cout<<"Writing PSI-BLAST-formatted alignment to "<<alnfile<<"\n";
      for (int k=0; k<N_in; k++)
	if (keep[k] || display[k]==2) // print if either in profile (keep[k]>0) or display obligatory (display[k]==2)
	  {
	    strcut(sname[k]);
	    fprintf(alnf,"%-20.20s ",sname[k]);
	    for (int i=1; i<=L; i++) fprintf(alnf,"%c",i2aa(X[k][i]));
	    fprintf(alnf,"\n");
	  }
    }
  
  fclose(alnf);  
}

// Read a3m slave alignment of hit from file and merge into (query) master alignment
void Alignment::MergeMasterSlave(Hit& hit, char ta3mfile[])
{
  Alignment Tali;
  char* cur_seq = new(char[MAXCOL]);   // Sequence currently read in
  int maxcol=MAXCOL;
  int l,ll;           // position in unaligned template (T) sequence Tali.seq[l]
  int i;              // counts match states in query (Q) HMM
  int j;              // counts match states in T sequence Tali.seq[l]
  int h;              // position in aligned T sequence cur_seq[h]
  int k;              // sequence index
  char c;             // 

  // If par.append==1 do not print query alignment
  if (par.append) for (k=0; k<N_in; k++) keep[k]=display[k]=0;

  // Read template alignment into Tali
  FILE* ta3mf=fopen(ta3mfile,"r");
  if (!ta3mf) OpenFileError(ta3mfile);
  Tali.Read(ta3mf,ta3mfile);
  fclose(ta3mf);
  if (v>=3) printf("Merging %s to query alignment\n",ta3mfile);
  
  // Record imatch[j]
  int* imatch=new(int[hit.j2+1]);
  int step = hit.nsteps;
  for (j=hit.j1; j<=hit.j2; j++) 
    {
      // Advance to position of next T match state j
      while (hit.j[step]<j) step--;
      imatch[j] = hit.i[step];
//       printf("step=%-3i  i=%-3i j=%-3i\n",step,imatch[j],j);
    }

  // Determine number of match states of Qali
  for (L=0,l=1; seq[kfirst][l]>'\0'; l++) 
    if ((seq[kfirst][l]>='A' && seq[kfirst][l]<='Z') || seq[kfirst][l]=='-') L++;

  // For each sequence in T alignment: align to Qali
  for (k=0; k<Tali.N_in; k++)
    {
      if (!Tali.keep[k]) continue;
      cur_seq[0]=' ';     // 0'th position not used
      
      // Add the hit.i1-1 left end gaps to aligned sequence
      for (h=1; h<hit.i1; h++) cur_seq[h]='-';   
      
      // Advance to match state hit.j1 of Tali.seq[k]
      for (j=0, l=1; (c=Tali.seq[k][l])>'\0'; l++)
	if ((c>='A' && c<='Z') || c=='-') // match state at position l?
	  if ((++j)==hit.j1) break;       // yes: increment j. Reached hit,j1? yes: break
      
      if (j<hit.j1) 
	{printf("Error: did not find %i match states in sequence %i of %s. Sequence:\n%s\n",hit.j1,k,Tali.name,Tali.seq[k]); exit(1);}
      
      // Write first match state to cur_seq
      int iprev=hit.i1; // index of previous query match state 
      int lprev=l;      // previous T match state in Tali.seq[k][l]
      cur_seq[h++] = Tali.seq[k][l]; // first column of alignment is Match-Match state
    
      // For each further match state j in alignment
      step = hit.nsteps;
      for (j=hit.j1+1; j<=hit.j2; j++) 
	{
	  // Advance to position of next T match state j
	  i=imatch[j];

	  // Advance to position of next T match state j
	  while ((c=Tali.seq[k][++l])>'\0' && ((c>='a' && c<='z') || c=='.')) ;

	  int di=i-iprev; // number of Match states in Q between T match state j-1 and j
	  int dl=l-lprev; // 1 + number of inserted residues in T sequence between T match state j-1 and j
	  if (di==1) 
	    {
	      // One Q match state for one T match state (treated as special case for speed reasons)
	      // i:       i-1   i         di=1
	      // Q:  XXXXXX.....XXXXXX   
	      // T:  YYYYYYyyyyyYYYYYY   
	      // j:       j-1   j     
	      // l:       lprev l         dl=6

	      // Inserts in lower case
	      for (ll=lprev+1; ll<l; ll++) 
		if (Tali.seq[k][ll]!='-' && Tali.seq[k][ll]!='.') cur_seq[h++] = lwrchr(Tali.seq[k][ll]);

	      // Template Match state -> upper case  
	      cur_seq[h++] = Tali.seq[k][ll];
	    }
	  else if (di==0) 
	    {
	      // Gap in query: no Q match state for on T match state (special case for speed reasons)
	      // i:       i-1   i-1       di=0
	      // Q:  XXXXXX.....~~~XXX   
	      // T:  YYYYYYyyyyyYYYYYY   
	      // j:       j-1   j     
	      // l:       lprev l         dl=6

	      // All T residues (including T match state) in lower case
	      for (ll=lprev+1; ll<=l; ll++) 
		if (Tali.seq[k][ll]!='-' && Tali.seq[k][ll]!='.') cur_seq[h++] = lwrchr(Tali.seq[k][ll]);
	    }
	  else if (di>=dl)
	    {
	      // More Match states in Q than Inserts in the T sequence
	      // => half T inserts y left, half right-aligned in uc, gaps to fill up
	      // Number of T insert residues to be left-aligned: (int)(dl/2)
	      // i:        iprev  i       di=7
	      // Q:  XXXXXXXXXXXXXXXXXX   
	      // T:  YYYYYYYyyy-yyYYYYY   
	      // j:        j-1    j
	      // l:        lprev  l       dl=6

	      // Add left-bounded template residues 
	      for (ll=lprev+1; ll<=lprev+(int)(dl/2); ll++) 
		cur_seq[h++]=uprchr(Tali.seq[k][ll]);

	      // Add central gaps
	      for (int gap=1; gap<=di-dl; gap++) cur_seq[h++]='-';

	      // Add right-bounded residues
	      for (; ll<=l; ll++)
		cur_seq[h++]=uprchr(Tali.seq[k][ll]);
	    }
	  else if (di<dl)
	    {
	      // Fewer Match states in Q than inserts in T sequence 
	      // => half of available space di for left- half for right-aligned T inserts, rest in lc 
	      // number of T inserts to be left-aligned in uc: (int)(di/2), 
	      // i:        iprev i       di=5 
	      // Q:  XXXXXXXXX.XXXXXXX   
	      // T:  YYYYYYYyyyyyYYYYY   
	      // j:        j-1   j
	      // l:        lprev l       dl=6

	      // Add left-bounded template residues 
	      for (ll=lprev+1; ll<=lprev+(int)(di/2); ll++) 
		cur_seq[h++]=uprchr(Tali.seq[k][ll]);

	      // Add central inserts
	      for (int ins=1; ins<=dl-di; ins++,ll++) 
		if (Tali.seq[k][ll]!='-' && Tali.seq[k][ll]!='.') cur_seq[h++] = lwrchr(Tali.seq[k][ll]);

	      // Add right-bounded residues
	      for (; ll<=l; ll++)
		cur_seq[h++]=uprchr(Tali.seq[k][ll]);
	    }
//  	  printf("i=%-3i j=%-3i l=%-3i cur_seq=%s\n",i,j,l,cur_seq);

	  iprev=i; lprev=l;
	  if (h>=maxcol-1000) // too few columns? Reserve double space
	    {
	      char* new_seq=new(char[2*maxcol]);
	      strncpy(new_seq,cur_seq,maxcol); //////// check: maxcol-1 ????
	      delete[](cur_seq);
	      cur_seq=new_seq;
	      maxcol*=2;
	    }
	}

      // Add the remaining gaps '-' to the end of the template sequence
      for (i=hit.i2+1; i<=L; i++) cur_seq[h++]='-';   
      cur_seq[h++]='\0';
     
      keep[N_in] = display[N_in] = 1;
      seq[N_in]=new(char[h]);
      if (!seq[N_in]) MemoryError("array for input sequences");
      strcpy(seq[N_in],cur_seq);
      X[N_in]=new(char[h]);
      if (!X[N_in]) MemoryError("array for input sequences");
      I[N_in]=new(short unsigned int[h]);
      if (!I[N_in]) MemoryError("array for input sequences");
      sname[N_in]=new(char[strlen(Tali.sname[k])+1]);
      if (!sname[N_in]) MemoryError("array for input sequences");
      strcpy(sname[N_in],Tali.sname[k]);
      N_in++;
      if (N_in>=MAXSEQ) 
	{
	  fprintf(stderr,"WARNING in %s: maximum number of %i sequences exceeded. Skipping following sequences\n",program_name,MAXSEQ); 
	  break;
	}

//       printf("k=%-3i %s\n",k,Tali.seq[k]);
//       printf("Query %s\n",seq[kfirst]);
//       printf("k=%-3i %s\n\n",k,cur_seq);
    
    } // end for (k)

  delete[] cur_seq;
  delete[] imatch;
}
