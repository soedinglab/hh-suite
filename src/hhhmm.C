// hhhmm.C

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
#endif

// #ifndef WNLIB
// #define WNLIB
// #include "wnconj.h"   // Will Naylor's wnlib for optimization in C
// #endif

/////////////////////////////////////////////////////////////////////////////////////
//// Class HMM
/////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////
// Object constructor
/////////////////////////////////////////////////////////////////////////////////////
HMM::HMM(int maxseqdis, int maxres)
{
  sname = new char*[maxseqdis];   // names of stored sequences
  seq = new char*[maxseqdis];     // residues of stored sequences (first at pos 1!)
  Neff_M = new float[maxres];     // Neff_M[i] = diversity of subalignment of seqs that have residue in col i
  Neff_I = new float[maxres];     // Neff_I[i] = diversity of subalignment of seqs that have insert in col i
  Neff_D = new float[maxres];     // Neff_D[i] = diversity of subalignment of seqs that have delete in col i
  longname = new char[DESCLEN];   // Full name of first sequence of original alignment (NAME field)
  ss_dssp = new char[maxres];     // secondary structure determined by dssp 0:-  1:H  2:E  3:C  4:S  5:T  6:G  7:B
  sa_dssp = new char[maxres];     // solvent accessibility state determined by dssp 0:-  1:A (absolutely buried) 2:B  3:C  4:D  5:E (exposed)
  ss_pred = new char[maxres];     // predicted secondary structure          0:-  1:H  2:E  3:C
  ss_conf = new char[maxres];     // confidence value of prediction         0:-  1:0 ... 10:9
  Xcons   = NULL;                 // create only when needed: consensus sequence in internal representation (A=0 R=1 N=2 D=3 ...)
  l = new int[maxres];            // l[i] = pos. of j'th match state in aligment
  f = new float*[maxres];         // f[i][a] = prob of finding amino acid a in column i WITHOUT pseudocounts
  g = new float*[maxres];         // f[i][a] = prob of finding amino acid a in column i WITH pseudocounts
  p = new float*[maxres];         // p[i][a] = prob of finding amino acid a in column i WITH OPTIMUM pseudocounts
  tr = new float*[maxres];        // log2 of transition probabilities M2M M2I M2D I2M I2I D2M D2D
//   tr_lin = new float*[maxres];    // linear transition probabilities M2M M2I M2D I2M I2I D2M D2D
  for (int i=0; i<maxres; i++) f[i]=new(float[NAA+3]);
  for (int i=0; i<maxres; i++) g[i]=new(float[NAA]);
  // for (int i=0; i<maxres; i++) p[i]=new(float[NAA]);
  // for (int i=0; i<maxres; i++) tr[i]=new(float[NTRANS]);
  for (int i=0; i<maxres; i++) p[i]=(float*) memalign(16,NAA*sizeof(float));  // align memory on 16b boundaries for SSE2
  for (int i=0; i<maxres; i++) tr[i]=(float*) memalign(16,NTRANS*sizeof(float));
//   for (int i=0; i<maxres; i++) tr_lin[i]=new(float[NTRANS]);
  L=0;
  Neff_HMM=0;
  n_display=N_in=N_filtered=0;
  nss_dssp=nsa_dssp=nss_pred=nss_conf=nfirst=ncons=-1;
//   lamda_hash.New(37,0.0); // Set size and NULL element for hash
//   mu_hash.New(37,0.0);    // Set size and NULL element for hash
  lamda=0.0; mu=0.0;
  name[0]=longname[0]=fam[0]='\0';
  trans_lin=0; // transition probs in log space
  dont_delete_seqs=0;
  has_pseudocounts=false;
}


/////////////////////////////////////////////////////////////////////////////////////
// Object destructor
/////////////////////////////////////////////////////////////////////////////////////
HMM::~HMM()
{
  //Delete name and seq matrices
  if (!dont_delete_seqs) // don't delete sname and seq if flat copy to hit object has been made
    {
      for (int k=0; k<n_display; k++) delete [] sname[k];
      for (int k=0; k<n_display; k++) delete [] seq[k];
    }
  delete[] sname;
  delete[] seq;
  delete[] Neff_M;
  delete[] Neff_D;
  delete[] Neff_I;
  delete[] longname;
  delete[] ss_dssp;
  delete[] sa_dssp;
  delete[] ss_pred;
  delete[] ss_conf;
  delete[] Xcons;
  delete[] l;
  for (int i=0; i<MAXRES; i++) if (f[i]) delete[] f[i]; else break;
  for (int i=0; i<MAXRES; i++) if (g[i]) delete[] g[i]; else break;
  for (int i=0; i<MAXRES; i++) if (p[i])  free(p[i]);  else break;
  for (int i=0; i<MAXRES; i++) if (tr[i]) free(tr[i]); else break;
  delete[] f;
  delete[] g;
  delete[] p;
  delete[] tr;
}

/////////////////////////////////////////////////////////////////////////////////////
// Deep-copy constructor
/////////////////////////////////////////////////////////////////////////////////////
HMM& HMM::operator=(HMM& q)
{
  L=q.L;
  for (int i=0; i<=L+1; ++i)
    {
      for (int a=0; a<NAA; ++a)
        {
          f[i][a]=q.f[i][a];
          g[i][a]=q.g[i][a];
          p[i][a]=q.p[i][a];
        }
      for (int a=0; a<NTRANS; ++a)
        tr[i][a]=q.tr[i][a];
      ss_dssp[i]=q.ss_dssp[i];
      sa_dssp[i]=q.sa_dssp[i];
      ss_pred[i]=q.ss_pred[i];
      ss_conf[i]=q.ss_conf[i];
      l[i]=q.l[i];
    }
  if (q.Xcons)
    for (int i=0; i<=L+1; ++i)
      Xcons[i]  =q.Xcons[i];

  n_display=q.n_display;
  for (int k=0; k<n_display; k++) {
    sname[k]=new(char[strlen(q.sname[k])+1]);
    if (!sname[k]) MemoryError("array of names for sequences to display");
    strcpy(sname[k],q.sname[k]);
  }
  for (int k=0; k<n_display; k++) {
    seq[k]=new(char[strlen(q.seq[k])+1]);
    if (!seq[k]) MemoryError("array of names for sequences to display");
    strcpy(seq[k],q.seq[k]);
  }
  ncons=q.ncons;
  nfirst=q.nfirst;
  nss_dssp=q.nss_dssp;
  nsa_dssp=q.nsa_dssp;
  nss_pred=q.nss_pred;
  nss_conf=q.nss_conf;

  for (int i=0; i<=L+1; ++i) Neff_M[i]=q.Neff_M[i];
  for (int i=0; i<=L+1; ++i) Neff_I[i]=q.Neff_I[i];
  for (int i=0; i<=L+1; ++i) Neff_D[i]=q.Neff_D[i];
  Neff_HMM=q.Neff_HMM;

  strcpy(longname,q.longname);
  strcpy(name,q.name);
  strcpy(fam,q.fam);
  strcpy(sfam,q.sfam);
  strcpy(fold,q.fold);
  strcpy(cl,q.cl);
  strcpy(file,q.file);

  lamda=q.lamda;
  mu=q.mu;
  has_pseudocounts=q.has_pseudocounts;

  for (int a=0; a<NAA; ++a) pav[a]=q.pav[a];
  N_in=q.N_in;
  N_filtered=q.N_filtered;
  trans_lin=q.trans_lin;
  return (HMM&) (*this);
}


/////////////////////////////////////////////////////////////////////////////////////
//// Read an HMM from an HHsearch .hhm file; return 0 at end of file
/////////////////////////////////////////////////////////////////////////////////////
int HMM::Read(FILE* dbf, char* path)
{
  char line[LINELEN]="";    // input line
  char str3[8]="",str4[8]=""; // first 3 and 4 letters of input line
  char* ptr;                // pointer for string manipulation
  int i=0;                  // index for match state (first=1)
  int a;                    // amino acid index
  static int warn=0;

  trans_lin=0;
  L=0;
  Neff_HMM=0;
  n_display=N_in=N_filtered=0;
  nss_dssp=nsa_dssp=nss_pred=nss_conf=nfirst=ncons=-1;
  lamda=mu=0.0;
  trans_lin=0; // transition probs in log space
  name[0]=longname[0]=fam[0]='\0';
  has_pseudocounts=false;
  //If at the end of while-loop L is still 0 then we have reached end of db file

  //Do not delete name and seq vectors because their adresses are transferred to hitlist as part of a hit!!

  while (fgetline(line,LINELEN-1,dbf) && !(line[0]=='/' && line[1]=='/'))
    {

      if (strscn(line)==NULL) continue;    // skip lines that contain only white space
      substr(str3,line,0,2);               // copy the first three characters into str3
      substr(str4,line,0,3);               // copy the first four characters into str4

      if (!strncmp("HH",line,2)) continue;

      if (!strcmp("NAME",str4))
        {
          ptr=strscn(line+4);              //advance to first non-white-space character
          if (ptr)
            {
              strncpy(longname,ptr,DESCLEN-1); //copy full name to longname
              longname[DESCLEN-1]='\0';
              strncpy(name,ptr,NAMELEN-1);     //copy longname to name...
              strcut(name);                    //...cut after first word...
            }
          else
            {
              strcpy(longname,"undefined");
              strcpy(name,"undefined");
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
          if (path) strncpy(file,path,NAMELEN-1); else *file='\0'; // copy path to file variable
          ptr=strscn(line+4);              //advance to first non-white-space character
          if (ptr)
            strncat(file,ptr,NAMELEN-1-strlen(file));   // append file name read from file to path
          else strcat(file,"*");
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
//        char key[IDLEN];
          sscanf(line+6,"%f %f",&lamda,&mu);
//        sscanf(line+22,"%s",key);
//        lamda_hash.Add(key,lamda);
//        mu_hash.Add(key,mu);
        }

      else if (!strcmp("PCT",str3)) { has_pseudocounts=true; }
      else if (!strcmp("DESC",str4)) continue;
      else if (!strcmp("COM",str3))  continue;
      else if (!strcmp("DATE",str4)) continue;

      /////////////////////////////////////////////////////////////////////////////////////
      // Read template sequences that should get displayed in output alignments
      else if (!strcmp("SEQ",str3))
        {
          char cur_seq[MAXCOL]=""; //Sequence currently read in
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
                  if (k>=MAXSEQDIS-1) //maximum number of allowable sequences exceeded
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
                      if (n_seq>=par.nseqdis)
                        {while (fgetline(line,LINELEN-1,dbf) && line[0]!='#'); k--; break;}
                      n_seq++;
                    }

                  //If this is not the first sequence then store residues of previous sequence
                  if (k>0) {
                    seq[k-1]=new(char[strlen(cur_seq)+1]);
                    if (!seq[k-1]) MemoryError("array of sequences to display");
                    strcpy(seq[k-1],cur_seq);
                  }

                  // store sequence name
                  strcut(line+1); //find next white-space character and overwrite it with end-of-string character
                  sname[k] = new (char[strlen(line+1)+1]); //+1 for terminating '\0'
                  if (!sname[k]) MemoryError("array of names for sequences to display");
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

                } //end else
            } //while(getline)
          //If this is not the first sequence some residues have already been read in
          if (k>=0) {
            seq[k]=new(char[strlen(cur_seq)+1]);
            if (!seq[k]) MemoryError("array of sequences to display");
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
          exit(1);
        }

      else if (!strcmp("AVER",str4)) {} // AVER line scrapped
      else if (!strcmp("NULL",str4))
        {
          ptr=line+4;
          for (a=0; a<20 && ptr; ++a)
            //s2[a]: transform amino acids Sorted by alphabet -> internal numbers for amino acids
            pb[s2a[a]] = (float) fpow2(float(-strinta(ptr))/HMMSCALE);
          if (!ptr) return Warning(dbf,line,name);
          if (v>=4)
            {
              printf("\nNULL  ");
              for (a=0; a<20; ++a) printf("%5.1f ",100.*pb[s2a[a]]);
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

          for (a=0; a<=D2D && ptr; ++a)
            tr[0][a] = float(-strinta(ptr))/HMMSCALE; //store transition probabilites as log2 values
            // strinta returns next integer in string and puts ptr to first char
            // after the integer. Returns -99999 if '*' is found.
            // ptr is set to 0 if no integer is found after ptr.
          Neff_M[0] = float(strinta(ptr))/HMMSCALE;  // Read eff. number of sequences with M->? transition
          Neff_I[0] = float(strinta(ptr))/HMMSCALE;  // Read eff. number of sequences with I->? transition
          Neff_D[0] = float(strinta(ptr))/HMMSCALE;  // Read eff. number of sequences with D->? transition
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
              next_i = strint(ptr); ++i;
              if (v && next_i!=prev_i+1)
                if (++warn<=5)
                  {
                    cerr<<endl<<"WARNING: in HMM "<<name<<" state "<<prev_i<<" is followed by state "<<next_i<<"\n";
                    if (warn==5) cerr<<endl<<"WARNING: further warnings while reading HMMs will be suppressed.\n";
                  }
              if (i>L)
                {
                  cerr<<endl<<"WARNING: in HMM "<<name<<" there are more columns than the stated length "<<L<<". Skipping HMM\n";
                  return 2;
                }
              if (i>MAXRES-2)
                {
                  fgetline(line,LINELEN-1,dbf); // Skip line
                  continue;
                }

              for (a=0; a<20 && ptr; ++a)
//              f[i][s2a[a]] = (float)pow(2.,float(-strinta(ptr))/HMMSCALE);
                f[i][s2a[a]] = fpow2(float(-strinta(ptr))/HMMSCALE);       // speed-up ~5 s for 10000 SCOP domains

              //s2a[a]: transform amino acids Sorted by alphabet -> internal numbers for amino acids
              l[i]=strint(ptr);
              if (!ptr) return Warning(dbf,line,name);
              if (v>=4)
                {
                  printf("%s",line);
                  printf("%6i ",i);
                  for (a=0; a<20; ++a) printf("%5.1f ",100*f[i][s2a[a]]);
                  printf("%5i",l[i]);
                  printf("\n");
                }

              // Read transition probabilities
              fgetline(line,LINELEN-1,dbf); // Skip line with amino acid labels
              if (line[0]!=' ' && line[0]!='\t') return Warning(dbf,line,name);
              ptr=line;
              for (a=0; a<=D2D && ptr; ++a)
                tr[i][a] = float(-strinta(ptr))/HMMSCALE; //store transition prob's as log2-values
              Neff_M[i] = float(strinta(ptr))/HMMSCALE;  // Read eff. number of sequences with M->? transition
              Neff_I[i] = float(strinta(ptr))/HMMSCALE;  // Read eff. number of sequences with I->? transition
              Neff_D[i] = float(strinta(ptr))/HMMSCALE;  // Read eff. number of sequences with D->? transition
              if (!ptr) return Warning(dbf,line,name);
              if (v>=4)
                {
                  printf("       ");
                  for (a=0; a<=D2D; ++a) printf("%5.1f ",100*fpow2(tr[i][a]));
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
  if (v && i>MAXRES-2) {i=MAXRES-2; cerr<<endl<<"WARNING: maximum number "<<MAXRES-2<<" of residues exceeded while reading HMM "<<name<<"\n";}
  if (v && !i)  cerr<<endl<<"WARNING: HMM "<<name<<" contains no match states. Check the alignment that gave rise to this HMM.\n";
  if (v>=2) cout<<"Read in HMM "<<name<<" with "<<L<<" match states and effective number of sequences = "<<Neff_HMM<<"\n";
  L = i;

  // Set emission probabilities of zero'th (begin) state and L+1st (end) state to background probabilities
  for (a=0; a<20; ++a) f[0][a]=f[L+1][a]=pb[a];
  Neff_M[L+1]=1.0f;
  Neff_I[L+1]=Neff_D[L+1]=0.0f;

  return 1; //return status: ok
}


/////////////////////////////////////////////////////////////////////////////////////
//// Read an HMM from a HMMer .hmm file; return 0 at end of file
/////////////////////////////////////////////////////////////////////////////////////
int HMM::ReadHMMer(FILE* dbf, char* filestr)
{
  char line[LINELEN]="";    // input line
  char desc[DESCLEN]="";    // description of family
  char str4[5]="";          // first 4 letters of input line
  char* ptr;                // pointer for string manipulation
  int i=0;                  // index for match state (first=1)
  int a;                    // amino acid index
  char dssp=0;              // 1 if a consensus SS has been found in the transition prob lines
  char annot=0;             // 1 if at least one annotation character in insert lines is ne '-' or ' '
  int k=0;                  // index for seq[k]
  static char ignore_hmmer_cal = 0;
  char* annotchr;           // consensus amino acids in ASCII format, or, in HMMER format, the reference annotation character in insert line
  annotchr = new char[MAXRES]; // consensus amino acids in ASCII format, or, in HMMER format, the reference annotation character in insert line
  static int warn=0;

  trans_lin=0;
  L=0;
  Neff_HMM=0;
  n_display=N_in=N_filtered=0;
  nss_dssp=nsa_dssp=nss_pred=nss_conf=nfirst=ncons=-1;
  lamda=mu=0.0;
  trans_lin=0; // transition probs in log space
  name[0]=longname[0]=desc[0]=fam[0]='\0';
  //If at the end of while-loop L is still 0 then we have reached end of db file

  // Do not delete name and seq vectors because their adresses are transferred to hitlist as part of a hit!!

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
              sname[nsa_dssp] = new(char[15]);
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
              sname[nss_pred] = new(char[15]);
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
              sname[nss_conf] = new(char[15]);
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
          for (a=0; a<20 && ptr; ++a)
            //s2a[a]: transform amino acids Sorted by alphabet -> internal numbers for amino acids
            pb[s2a[a]] = (float) 0.05 * fpow2(float(strinta(ptr,-99999))/HMMSCALE);
          if (!ptr) return Warning(dbf,line,name);
          if (v>=4)
            {
              printf("\nNULL  ");
              for (a=0; a<20; ++a) printf("%6.3g ",100.*pb[s2a[a]]);
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
              if (v>=2 && ignore_hmmer_cal==0)
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
          for (a=0; a<=M2D && ptr; ++a)
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
              for (a=0; a<=D2D && ptr; ++a) printf("%6.3g ",100*fpow2(tr[i][a]));
              printf("\n");
            }

          // Prepare to store DSSP states (if there are none, delete afterwards)
          nss_dssp=k++;
          seq[nss_dssp] = new(char[MAXRES+2]);
          sname[nss_dssp] = new(char[15]);
          strcpy(sname[nss_dssp],"ss_dssp");

          /////////////////////////////////////////////////////////////////////////////////////
          // Read columns of HMM
          int next_i=0;  // index of next column
          while (fgetline(line,LINELEN-1,dbf) &&  !(line[0]=='/' && line[1]=='/') && line[0]!='#')
            {
              if (strscn(line)==NULL) continue; // skip lines that contain only white space

              // Read in AA probabilities
              ptr=line;
              int prev_i = next_i;
              next_i = strint(ptr); ++i;
              if (v && next_i!=prev_i+1)
                if (++warn<5)
                  {
                    cerr<<endl<<"WARNING: in HMM "<<name<<" state "<<prev_i<<" is followed by state "<<next_i<<"\n";
                    if (warn==5) cerr<<endl<<"WARNING: further warnings while reading HMMs will be suppressed.\n";
                  }
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

              for (a=0; a<20 && ptr; ++a)
                f[i][s2a[a]] = (float) pb[s2a[a]]*fpow2(float(strinta(ptr,-99999))/HMMSCALE);
              //s2a[a]: transform amino acids Sorted by alphabet -> internal numbers for amino acids
              if (!ptr) return Warning(dbf,line,name);
              if (v>=4)
                {
                  printf("%6i ",i);
                  for (a=0; a<20; ++a) printf("%6.3g ",100*f[i][s2a[a]]);
                  printf("\n");
                }

              // Read insert emission line
              fgetline(line,LINELEN-1,dbf);
              ptr = strscn(line);
              if (!ptr) return Warning(dbf,line,name);
              annotchr[i]=uprchr(*ptr);
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
		case '-': // no SS available from any template
		case '.': // no clear consensus SS structure
		case 'X': // no clear consensus SS structure
		  ss_dssp[i]=0;
		  seq[nss_dssp][i]='-';
		  break;
		default:
		  ss_dssp[i]=0;
		  seq[nss_dssp][i]=*ptr;
		  break;
               }

              ptr+=2;
              for (a=0; a<=D2D && ptr; ++a)
                tr[i][a] = float(strinta(ptr,-99999))/HMMSCALE; //store transition prob's as log2-values
              if (!ptr) return Warning(dbf,line,name);
              if (v>=4)
                {
                  printf("       ");
                  for (a=0; a<=D2D; ++a) printf("%6.3g ",100*fpow2(tr[i][a]));
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
  RemoveExtension(file,filestr); // copy name of dbfile without extension into 'file'

  // Secondary structure
  if (!dssp)
    {
      // remove dssp sequence
      delete[] seq[nss_dssp];    // memory that had been allocated in case ss_dssp was given needs to be freed
      delete[] sname[nss_dssp];  // memory that had been allocated in case ss_dssp was given needs to be freed
      nss_dssp=-1;
      k--;
    }
  else { seq[nss_dssp][0]='-'; seq[nss_dssp][L+1]='\0'; }

  if (nss_pred>=0)
    {
      for (i=1; i<=L; ++i) ss_pred[i] = ss2i(seq[nss_pred][i]);
      if (nss_conf>=0)
        for (i=1; i<=L; ++i) ss_conf[i] = cf2i(seq[nss_conf][i]);
      else
        for (i=1; i<=L; ++i) ss_conf[i] = 5;
    }

  // Copy query (first sequence) and consensus  residues?
  if (par.showcons)
    {
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
      for (i=1; i<=L; ++i)
        {
          float pmax=0.0;
          int amax=0;
          for (a=0; a<NAA; ++a)
            if (f[i][a]>pmax) {amax=a; pmax=f[i][a];}
          if (pmax>0.6) seq[k][i]=i2aa(amax);
          else if (pmax>0.4) seq[k][i]=lwrchr(i2aa(amax));
          else seq[k][i]='x';
          seq[k+1][i]=i2aa(amax);
        }
      ncons=k++; // nfirst is set later!
    }
  else
    {
      sname[k]=new(char[strlen(longname)+1]);
      strcpy(sname[k],longname);
      seq[k]=new(char[L+2]);
      seq[k][0]=' ';
      seq[k][L+1]='\0';
    }

  if (annot) // read in some annotation characters?
    {
      annotchr[0]=' ';
      annotchr[L+1]='\0';
      strcpy(seq[k],annotchr); // overwrite the consensus sequence with the annotation characters
    }
  else if (!par.showcons)  // we have not yet calculated the consensus, but we need it now as query (first sequence)
    {
      for (i=1; i<=L; ++i)
        {
          float pmax=0.0;
          int amax=0;
          for (a=0; a<NAA; ++a)
            if (f[i][a]>pmax) {amax=a; pmax=f[i][a];}
          seq[k][i]=i2aa(amax);
        }
    }
//   printf("%i query name=%s  seq=%s\n",n,sname[n],seq[n]);
  nfirst=k++;

  n_display=k;

  // Calculate overall Neff_HMM
  Neff_HMM=0;
  for (i=1; i<=L; ++i)
    {
      float S=0.0;
      for (a=0; a<20; ++a)
        if (f[i][a]>1E-10) S-=f[i][a]*fast_log2(f[i][a]);
      Neff_HMM+=(float) fpow2(S);
    }
  Neff_HMM/=L;
  for (i=0; i<=L; ++i) Neff_M[i] = Neff_I[i] = Neff_D[i] = 10.0; // to add only little additional pseudocounts!
  Neff_M[L+1]=1.0f;
  Neff_I[L+1]=Neff_D[L+1]=0.0f;
  if (v>=2)
    cout<<"Read in HMM "<<name<<" with "<<L<<" match states and effective number of sequences = "<<Neff_HMM<<"\n";

  // Set emission probabilities of zero'th (begin) state and L+1st (end) state to background probabilities
  for (a=0; a<20; ++a) f[0][a]=f[L+1][a]=pb[a];
  delete[] annotchr;

  return 1; //return status: ok
}

/////////////////////////////////////////////////////////////////////////////////////
//// Read an HMM from a HMMER3 .hmm file; return 0 at end of file
/////////////////////////////////////////////////////////////////////////////////////
int HMM::ReadHMMer3(FILE* dbf, char* filestr)
{
  char line[LINELEN]="";    // input line
  char desc[DESCLEN]="";    // description of family
  char str4[5]="";          // first 4 letters of input line
  char* ptr;                // pointer for string manipulation
  int i=0;                  // index for match state (first=1)
  int a;                    // amino acid index
  char dssp=0;              // 1 if a consensus SS has been found in the transition prob lines
  char annot=0;             // 1 if at least one annotation character in insert lines is ne '-' or ' '
  int k=0;                  // index for seq[k]
  char* annotchr;           // consensus amino acids in ASCII format, or, in HMMER format, the reference annotation character in insert line
  annotchr = new char[MAXRES]; // consensus amino acids in ASCII format, or, in HMMER format, the reference annotation character in insert line
  static int warn=0;

  trans_lin=0;
  L=0;
  Neff_HMM=0;
  n_display=N_in=N_filtered=0;
  nss_dssp=nsa_dssp=nss_pred=nss_conf=nfirst=ncons=-1;
  lamda=mu=0.0;
  trans_lin=0; // transition probs in log space
  name[0]=longname[0]=desc[0]=fam[0]='\0';
  //If at the end of while-loop L is still 0 then we have reached end of db file

  // Do not delete name and seq vectors because their adresses are transferred to hitlist as part of a hit!!

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
      else if (!strcmp("EFFN",str4)) continue;
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

      //////////////////////////////////////////////////////////////////////////////////////////////////////
      // Still needed???

      else if (!strncmp("SADSS",line,5))
        {
          if (nsa_dssp<0)
            {
              nsa_dssp=k++;
              seq[nsa_dssp] = new(char[MAXRES+2]);
              sname[nsa_dssp] = new(char[15]);
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
              sname[nss_pred] = new(char[15]);
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
              sname[nss_conf] = new(char[15]);
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
      //////////////////////////////////////////////////////////////////////////////////////////////////////

      else if (!strncmp("STATS LOCAL",line,11)) continue;

      /////////////////////////////////////////////////////////////////////////////////////
      // Read transition probabilities from start state
      else if (!strncmp("HMM",line,3))
        {
          fgetline(line,LINELEN-1,dbf); // Skip line with amino acid labels
          fgetline(line,LINELEN-1,dbf); // Skip line with transition labels
          ptr=strscn(line);

	  if (!strncmp("COMPO",ptr,5))
	    {
	      ptr=ptr+5;
	      for (a=0; a<20 && ptr; ++a)
		//s2a[a]: transform amino acids Sorted by alphabet -> internal numbers for amino acids
		pb[s2a[a]] = (float) exp(-1.0*strflta(ptr,99999));
	      if (!ptr) return Warning(dbf,line,name);
	      if (v>=4)
		{
		  printf("\nNULL ");
		  for (a=0; a<20; ++a) printf("%6.3g ",100.*pb[s2a[a]]);
		  printf("\n");
		}
	      fgetline(line,LINELEN-1,dbf); // Read next line
	    }
	      
	  fgetline(line,LINELEN-1,dbf); // Skip line with 0-states insert probabilities

	  ptr = strscn(line);
	  for (a=0; a<=D2D && ptr; ++a)
	    tr[0][a] = log2((float) exp(-1.0*strflta(ptr,99999))); //store transition probabilites as log2 values
	  // strinta returns next integer in string and puts ptr to first char
	  // after the integer. Returns -99999 if '*' is found.
	  // ptr is set to 0 if no integer is found after ptr.
          if (!ptr) return Warning(dbf,line,name);
          if (v>=4)
            {
              printf("       ");
              for (a=0; a<=D2D && ptr; ++a) printf("%6.3g ",100*fpow2(tr[i][a]));
              printf("\n");
            }
	  
          // Prepare to store DSSP states (if there are none, delete afterwards)
          nss_dssp=k++;
          seq[nss_dssp] = new(char[MAXRES+2]);
          sname[nss_dssp] = new(char[15]);
          strcpy(sname[nss_dssp],"ss_dssp");

	  /////////////////////////////////////////////////////////////////////////////////////
          // Read columns of HMM
          int next_i=0;  // index of next column
          while (fgetline(line,LINELEN-1,dbf) &&  !(line[0]=='/' && line[1]=='/') && line[0]!='#')
            {
              if (strscn(line)==NULL) continue; // skip lines that contain only white space

              // Read in AA probabilities
              ptr=line;
              int prev_i = next_i;
              next_i = strint(ptr); ++i;
              if (v && next_i!=prev_i+1)
                if (++warn<5)
                  {
                    cerr<<endl<<"WARNING: in HMM "<<name<<" state "<<prev_i<<" is followed by state "<<next_i<<"\n";
                    if (warn==5) cerr<<endl<<"WARNING: further warnings while reading HMMs will be suppressed.\n";
                  }
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

              for (a=0; a<20 && ptr; ++a)
                f[i][s2a[a]] = (float) exp(-1.0*strflta(ptr,99999));
              //s2a[a]: transform amino acids Sorted by alphabet -> internal numbers for amino acids
              if (!ptr) return Warning(dbf,line,name);
              if (v>=4)
                {
                  printf("%6i ",i);
                  for (a=0; a<20; ++a) printf("%6.3g ",100*f[i][s2a[a]]);
                  printf("\n");
                }

	      // Ignore MAP annotation
	      ptr = strscn(line); //find next word
	      ptr = strscn_ws(line); // ignore word

	      // Read RF and CS annotation
	      ptr = strscn(line);
              if (!ptr) return Warning(dbf,line,name);
              annotchr[i]=uprchr(*ptr);
              if (*ptr!='-' && *ptr!=' ') annot=1;

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
		case '-': // no SS available from any template
		case '.': // no clear consensus SS structure
		case 'X': // no clear consensus SS structure
		  ss_dssp[i]=0;
		  seq[nss_dssp][i]='-';
		  break;
		default:
		  ss_dssp[i]=0;
		  seq[nss_dssp][i]=*ptr;
		  break;
               }

              // Read insert emission line
              fgetline(line,LINELEN-1,dbf);

              // Read seven transition probabilities
              fgetline(line,LINELEN-1,dbf);
             
              ptr+=2;
              for (a=0; a<=D2D && ptr; ++a)
                tr[i][a] = log2((float) exp(-1.0*strflta(ptr,99999))); //store transition prob's as log2-values
              if (!ptr) return Warning(dbf,line,name);
              if (v>=4)
                {
                  printf("       ");
                  for (a=0; a<=D2D; ++a) printf("%6.3g ",100*fpow2(tr[i][a]));
                  printf("\n");
                }
            }

          if (line[0]=='/' && line[1]=='/') break;

        }

    } //while(getline)

  if (L==0) return 0; //End of db file -> stop reading in

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
  RemoveExtension(file,filestr); // copy name of dbfile without extension into 'file'

  // Secondary structure
  if (!dssp)
    {
      // remove dssp sequence
      delete[] seq[nss_dssp];    // memory that had been allocated in case ss_dssp was given needs to be freed
      delete[] sname[nss_dssp];  // memory that had been allocated in case ss_dssp was given needs to be freed
      nss_dssp=-1;
      k--;
    }
  else { seq[nss_dssp][0]='-'; seq[nss_dssp][L+1]='\0'; }

  if (nss_pred>=0)
    {
      for (i=1; i<=L; ++i) ss_pred[i] = ss2i(seq[nss_pred][i]);
      if (nss_conf>=0)
        for (i=1; i<=L; ++i) ss_conf[i] = cf2i(seq[nss_conf][i]);
      else
        for (i=1; i<=L; ++i) ss_conf[i] = 5;
    }

  // Copy query (first sequence) and consensus  residues?
  if (par.showcons)
    {
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
      for (i=1; i<=L; ++i)
        {
          float pmax=0.0;
          int amax=0;
          for (a=0; a<NAA; ++a)
            if (f[i][a]>pmax) {amax=a; pmax=f[i][a];}
          if (pmax>0.6) seq[k][i]=i2aa(amax);
          else if (pmax>0.4) seq[k][i]=lwrchr(i2aa(amax));
          else seq[k][i]='x';
          seq[k+1][i]=i2aa(amax);
        }
      ncons=k++; // nfirst is set later!
    }
  else
    {
      sname[k]=new(char[strlen(longname)+1]);
      strcpy(sname[k],longname);
      seq[k]=new(char[L+2]);
      seq[k][0]=' ';
      seq[k][L+1]='\0';
    }

  if (annot) // read in some annotation characters?
    {
      annotchr[0]=' ';
      annotchr[L+1]='\0';
      strcpy(seq[k],annotchr); // overwrite the consensus sequence with the annotation characters
    }
  else if (!par.showcons)  // we have not yet calculated the consensus, but we need it now as query (first sequence)
    {
      for (i=1; i<=L; ++i)
        {
          float pmax=0.0;
          int amax=0;
          for (a=0; a<NAA; ++a)
            if (f[i][a]>pmax) {amax=a; pmax=f[i][a];}
          seq[k][i]=i2aa(amax);
        }
    }
//   printf("%i query name=%s  seq=%s\n",n,sname[n],seq[n]);
  nfirst=k++;

  n_display=k;

  ///////////////////////////////////////////////////////////////////
  // TODO

  // Calculate overall Neff_HMM
  Neff_HMM=0;
  for (i=1; i<=L; ++i)
    {
      float S=0.0;
      for (a=0; a<20; ++a)
        if (f[i][a]>1E-10) S-=f[i][a]*fast_log2(f[i][a]);
      Neff_HMM+=(float) fpow2(S);
    }
  Neff_HMM/=L;
  for (i=0; i<=L; ++i) Neff_M[i] = Neff_I[i] = Neff_D[i] = 10.0; // to add only little additional pseudocounts!
  Neff_M[L+1]=1.0f;
  Neff_I[L+1]=Neff_D[L+1]=0.0f;
  if (v>=2)
    cout<<"Read in HMM "<<name<<" with "<<L<<" match states and effective number of sequences = "<<Neff_HMM<<"\n";

  ///////////////////////////////////////////////////////////////////

  // Set emission probabilities of zero'th (begin) state and L+1st (end) state to background probabilities
  for (a=0; a<20; ++a) f[0][a]=f[L+1][a]=pb[a];
  delete[] annotchr;

  return 1; //return status: ok
}


/////////////////////////////////////////////////////////////////////////////////////
// Add transition pseudocounts to HMM (and calculate lin-space transition probs)
/////////////////////////////////////////////////////////////////////////////////////
void HMM::AddTransitionPseudocounts(float gapd, float gape, float gapf, float gapg, float gaph, float gapi, float gapb)
{
  int i;               //position in alignment
  float sum;
  float pM2M, pM2I, pM2D, pI2I, pI2M, pD2D, pD2M;
  float p0,p1,p2;
  if (par.gapb<=0) return;
  if (trans_lin==1) {fprintf(stderr,"Error: Adding transition pseudocounts to linear representation of %s not allowed. Please report this error to the HHsearch developers.\n",name); exit(6);}
  if (trans_lin==2) {fprintf(stderr,"Error: Adding transition pseudocounts twice is %s not allowed. Please report this error to the HHsearch developers.\n",name); exit(6);}
  trans_lin=2;

  // Calculate pseudocount transition probabilities
  pM2D=pM2I=gapd*0.0286;     //a-priori probability for inserts and deletions
  pM2M=1-pM2D-pM2I;
  // gape=0 -> pI2I=0   gape=1 -> pI2I=0.75    gape=inf -> pI2I=1.
  pI2I=1.0*gape/(gape-1+1.0/0.75);
  pI2M=1-pI2I;
  // gape=0 -> pD2D=0   gape=1 -> pD2D=0.75    gape=inf -> pD2D=1.
  pD2D=1.0*gape/(gape-1+1.0/0.75);
  pD2M=1-pD2D;

  for (i=0; i<=L; ++i) //for all columns in HMM
    {
      // Transitions from M state
      p0 = (Neff_M[i]-1)*fpow2(tr[i][M2M]) + gapb*pM2M;
      p1 = (Neff_M[i]-1)*fpow2(tr[i][M2D]) + gapb*pM2D;
      p2 = (Neff_M[i]-1)*fpow2(tr[i][M2I]) + gapb*pM2I;
      if (i==0) p1=p2=0;       //from M(0) no transition to D(1) and I(0) possible
      if (i==L) p1=p2=0;       //from M(L) no transition to D(L+1) and I(L+1) possible
      sum = p0+p1+p2+FLT_MIN;

//       p0 = p0/sum ;
//       p1 = pow(p1/sum,gapf);
//       p2 = pow(p2/sum,gapg);
//       sum = p0+p1+p2+FLT_MIN;
//       tr[i][M2M] = fast_log2(p0/sum);
//       tr[i][M2D] = fast_log2(p1/sum);
//       tr[i][M2I] = fast_log2(p2/sum);

      tr[i][M2M] = fast_log2(p0/sum);
      tr[i][M2D] = fast_log2(p1/sum)*gapf;
      tr[i][M2I] = fast_log2(p2/sum)*gapg;

      // Transitions from I state
      p0 = Neff_I[i]*fpow2(tr[i][I2M]) + gapb*pI2M;
      p1 = Neff_I[i]*fpow2(tr[i][I2I]) + gapb*pI2I;
      sum = p0+p1+FLT_MIN;

//       p0 = pow(p0/sum,gapg);
//       p1 = pow(p1/sum,gapi);
//       sum = p0+p1+FLT_MIN;
//       tr[i][I2M] = fast_log2(p0/sum);
//       tr[i][I2I] = fast_log2(p1/sum);

      tr[i][I2M] = fast_log2(p0/sum);
      tr[i][I2I] = fast_log2(p1/sum)*gapi;

      // Transitions from D state
      p0 = Neff_D[i]*fpow2(tr[i][D2M]) + gapb*pD2M;
      p1 = Neff_D[i]*fpow2(tr[i][D2D]) + gapb*pD2D;
      if (i==L) p1=0;          //from D(L) no transition to D(L+1) possible
      sum = p0+p1+FLT_MIN;

//       p0 = pow(p0/sum,gapf);
//       p1 = pow(p1/sum,gaph);
//       sum = p0+p1+FLT_MIN;
//       tr[i][D2M] = fast_log2(p0/sum);
//       tr[i][D2D] = fast_log2(p1/sum);

      tr[i][D2M] = fast_log2(p0/sum);
      tr[i][D2D] = fast_log2(p1/sum)*gaph;

    }

  if (v>=4)
    {
      printf("\nPseudocount transition probabilities:\n");
      printf("pM2M=%4.1f%%, pM2I=%4.1f%%, pM2D=%4.1f%%, ",100*pM2M,100*pM2I,100*pM2D);
      printf("pI2M=%4.1f%%, pI2I=%4.1f%%, ",100*pI2M,100*pI2I);
      printf("pD2M=%4.1f%%, pD2D=%4.1f%% ",100*pD2M,100*pD2D);
      printf("tau = %4.1f%%\n\n",100.*gapb/(Neff_HMM-1+gapb));
      printf("Listing transition probabilities WITH pseudocounts:\n");
      printf("   i dssp pred sacc     M->M   M->I   M->D   I->M   I->I   D->M   D->D\n");

      for (i=1; i<=L; ++i) //for all columns in HMM
        {
          printf("%4i  %1c    %1c    %1c    %6.3f %6.3f %6.3f ",i,i2ss(ss_dssp[i]),i2ss(ss_pred[i]),i2sa(sa_dssp[i]),fpow2(tr[i][M2M]),fpow2(tr[i][M2I]),fpow2(tr[i][M2D]));
          printf("%6.3f %6.3f ",fpow2(tr[i][I2M]),fpow2(tr[i][I2I]));
          printf("%6.3f %6.3f ",fpow2(tr[i][D2M]),fpow2(tr[i][D2D]));
          printf("%1i %2i  %1i\n",ss_pred[i],ss_conf[i],ss_dssp[i]);
        }
      printf("\n");
      printf("nss_dssp=%i  nss_pred=%i\n",nss_dssp,nss_pred);
    }
  return;
}

/////////////////////////////////////////////////////////////////////////////////////
// Generate an amino acid frequency matrix g[][] with full pseudocount admixture (tau=1)
/////////////////////////////////////////////////////////////////////////////////////
void HMM::PreparePseudocounts()
{
  for (int i=0; i<=L+1; ++i)
    for (int a=0; a<20; ++a)
      g[i][a] = // produces fast code
       R[a][0]*f[i][0]  +R[a][1]*f[i][1]  +R[a][2]*f[i][2]  +R[a][3]*f[i][3]  +R[a][4]*f[i][4]
      +R[a][5]*f[i][5]  +R[a][6]*f[i][6]  +R[a][7]*f[i][7]  +R[a][8]*f[i][8]  +R[a][9]*f[i][9]
      +R[a][10]*f[i][10]+R[a][11]*f[i][11]+R[a][12]*f[i][12]+R[a][13]*f[i][13]+R[a][14]*f[i][14]
      +R[a][15]*f[i][15]+R[a][16]*f[i][16]+R[a][17]*f[i][17]+R[a][18]*f[i][18]+R[a][19]*f[i][19];
}

/////////////////////////////////////////////////////////////////////////////////////
// Generate an AA frequency matrix g[][] with full context specific pseudocount admixture (tau=1)
/////////////////////////////////////////////////////////////////////////////////////

void HMM::AddContextSpecificPseudocounts(char pcm, float pca, float pcb, float pcc)
{
  int i;               //position in HMM
  int a;               //amino acid (0..19)
  
  cs::CountProfile<cs::AA> ali_profile(L);
  fillCountProfile(&ali_profile);

  cs::Admix *admix = NULL;

  switch (pcm)
    {
    case 1:
      admix =  new cs::ConstantAdmix(pca);
      break;
    case 2:
    case 4:
      admix =  new cs::HHsearchAdmix(pca, pcb, pcc);
      break;
    case 3:
      // TODO
      break;
    }      

  if (admix == NULL)
    {
      // write cs-pseudocounts in HMM.p
      for (i=1; i<=L; ++i)
	for (a=0; a<20; ++a)
	  p[i][a] = f[i][a];

    }
  else
    {

      cs::Profile<cs::AA> profile(lib_pc->AddTo(ali_profile, *admix));
      delete admix;
      
      // write cs-pseudocounts in HMM.p
      for (i=1; i<=L; ++i)
	for (a=0; a<20; ++a)
	  p[i][a] = profile[i-1][a];
    }

  // DEBUGGING output
  if (v>=3)
    {
      float sum;
      
      cout<<"Context specific pseudocounts added!\n\n";
      switch (pcm)
        {
        case 0:
          cout<<"No pseudocounts added (-pcm 0)\n";
          return;
        case 1:
          cout<<"Adding constant AA pseudocount admixture of "<<pca<<" to HMM "<<name<<"\n";
          break;
        case 2:
          cout<<"Adding divergence-dependent AA pseudocounts (-pcm 2) with admixture of "
              <<fmin(1.0, pca/(1. + Neff_HMM/pcb ) )<<" to HMM "<<name<<"\n";
          break;
        } //end switch (pcm)
      if (v>=4)
        {
          cout<<"\nAmino acid frequencies WITHOUT pseudocounts:\n       A    R    N    D    C    Q    E    G    H    I    L    K    M    F    P    S    T    W    Y    V\n";
          for (i=1; i<=L; ++i)
            {
              printf("%3i:  ",i);
              sum=0;
              for (a=0; a<20; ++a)
                {
                  sum+=f[i][a];
                  printf("%4.1f ",100*f[i][a]);
                }
              printf("  sum=%5.3f\n",sum);
            }
          cout<<"\nAmino acid frequencies WITH pseudocounts:\n       A    R    N    D    C    Q    E    G    H    I    L    K    M    F    P    S    T    W    Y    V\n";
          for (i=1; i<=L; ++i)
            {
              printf("%3i:  ",i);
              sum=0;
              for (a=0; a<20; ++a)
                {
                  sum+=p[i][a];
                  printf("%4.1f ",100*p[i][a]);
                }
              printf("  sum=%5.3f\n",sum);
            }
        }
   }
}

void HMM::fillCountProfile(cs::CountProfile<cs::AA> *csProfile)
{
  for (int i=0; i<L; ++i)
    {
      csProfile->neff[i] = Neff_M[i+1];
      for (int a=0; a<20; ++a)
	csProfile->counts[i][a] = f[i+1][a] * Neff_M[i+1];
    }
}

/////////////////////////////////////////////////////////////////////////////////////
// Calculate amino acid background frequencies in HMM
/////////////////////////////////////////////////////////////////////////////////////
void HMM::CalculateAminoAcidBackground()
{
    int a, i;
  // initialize vector of average aa freqs with pseudocounts
  for (a=0; a<20; ++a) pav[a]=pb[a]*100.0f/Neff_HMM;
  // calculate averages
  for (i=1; i<=L; ++i)
      for (a=0; a<20; ++a)
          pav[a] += p[i][a];
  // Normalize vector of average aa frequencies pav[a]
  NormalizeTo1(pav,NAA);
  for (a=0; a<20; ++a) p[0][a] = p[L+1][a] = pav[a];

}

/////////////////////////////////////////////////////////////////////////////////////
// Add amino acid pseudocounts to HMM and calculate average protein aa probabilities pav[a]
// Pseudocounts: t.p[i][a] = (1-tau)*f[i][a] + tau*g[i][a]
/////////////////////////////////////////////////////////////////////////////////////
void HMM::AddAminoAcidPseudocounts(char pcm, float pca, float pcb, float pcc)
{
  int i;               //position in HMM
  int a;               //amino acid (0..19)
  float sum;
  float tau;           //tau = pseudocount admixture

  // Calculate amino acid frequencies p[i][a] = (1-tau(i))*f[i][a] + tau(i)*g[i][a]
  switch (pcm)
    {
    case 0: //no pseudocounts whatsoever: tau=0
      for (i=1; i<=L; ++i)
        for (a=0; a<20; ++a)
          p[i][a]=f[i][a];
      break;
    case 1: //constant pseudocounts (for optimization): tau = pca
      tau = pca;
      for (i=1; i<=L; ++i)
        for (a=0; a<20; ++a)
          p[i][a] = (1.-tau)*f[i][a] + tau * g[i][a];
      break;
    case 2: //divergence-dependent pseudocounts
    case 4: //divergence-dependent pseudocounts and rate matrix rescaling
      if (par.pcc==1.0f)
        for (i=1; i<=L; ++i)
          {
            tau = fmin(1.0, pca/(1. + Neff_M[i]/pcb ) );
            for (a=0; a<20; ++a)
              p[i][a] = (1.-tau)*f[i][a] + tau * g[i][a];
          }
      else
        for (i=1; i<=L; ++i)
          {
            tau = fmin(1.0, pca/(1. + pow((Neff_M[i])/pcb,pcc)));
            for (a=0; a<20; ++a)
              p[i][a] = (1.-tau)*f[i][a] + tau * g[i][a];
          }
      break;
    case 3: // constant-divergence pseudocounts
      for (i=1; i<=L; ++i)
        {
          float x = Neff_M[i]/pcb;
          pca = 0.793 + 0.048*(pcb-10.0);
          tau = fmax(0.0, pca*(1-x + pcc*x*(1-x)) );
          for (a=0; a<20; ++a)
            p[i][a] = (1.-tau)*f[i][a] + tau * g[i][a];
        }
      if (v>=2) { printf("Divergence before / after addition of amino acid pseudocounts: %5.2f / %5.2f\n",Neff_HMM, CalcNeff()); }
     break;
    } //end switch (pcm)


  //turn on pseudocount switch to indicate that HMM contains pseudocounts
  if (pcm!=0) has_pseudocounts=true;

  // DEBUGGING output
  if (v>=3)
    {
      switch (pcm)
        {
        case 0:
          cout<<"No pseudocounts added (-pcm 0)\n";
          return;
        case 1:
          cout<<"Adding constant AA pseudocount admixture of "<<pca<<" to HMM "<<name<<"\n";
          break;
        case 2:
          cout<<"Adding divergence-dependent AA pseudocounts (-pcm 2) with admixture of "
              <<fmin(1.0, pca/(1. + Neff_HMM/pcb ) )<<" to HMM "<<name<<"\n";
          break;
        } //end switch (pcm)
      if (v>=4)
        {
          cout<<"\nAmino acid frequencies WITHOUT pseudocounts:\n       A    R    N    D    C    Q    E    G    H    I    L    K    M    F    P    S    T    W    Y    V\n";
          for (i=1; i<=L; ++i)
            {
              printf("%3i:  ",i);
              sum=0;
              for (a=0; a<20; ++a)
                {
                  sum+=f[i][a];
                  printf("%4.1f ",100*f[i][a]);
                }
              printf("  sum=%5.3f\n",sum);
            }
          cout<<"\nAmino acid frequencies WITH pseudocounts:\n       A    R    N    D    C    Q    E    G    H    I    L    K    M    F    P    S    T    W    Y    V\n";
          for (i=1; i<=L; ++i)
            {
              printf("%3i:  ",i);
              sum=0;
              for (a=0; a<20; ++a)
                {
                  sum+=p[i][a];
                  printf("%4.1f ",100*p[i][a]);
                }
              printf("  sum=%5.3f\n",sum);
            }
        }
   }
  return;
}


/////////////////////////////////////////////////////////////////////////////////////
// Factor Null model into HMM t
/////////////////////////////////////////////////////////////////////////////////////
void HMM::IncludeNullModelInHMM(HMM& q, HMM& t, int columnscore )
{

  int i,j;        //query and template match state indices
  int a;          //amino acid index

  switch (columnscore)
    {
    default:
    case 0: // Null model with background prob. from database
      for (a=0; a<20; ++a) pnul[a]=pb[a];
      break;

    case 1: // Null model with background prob. equal average from query and template
      for (a=0; a<20; ++a) pnul[a]=0.5*(q.pav[a]+t.pav[a]);
      break;

    case 2: // Null model with background prob. from template protein
      for (a=0; a<20; ++a) pnul[a]=t.pav[a];
      break;

    case 3: // Null model with background prob. from query protein
      for (a=0; a<20; ++a) pnul[a]=q.pav[a];
      break;

    case 4: // Null model with background prob. equal average from query and template
      for (a=0; a<20; ++a) pnul[a]=sqrt(q.pav[a]*t.pav[a]);
      break;

   case 10: // Separated column scoring for Stochastic Backtracing (STILL USED??)
      for (i=0; i<=q.L+1; ++i)
        {
          float sum = 0.0;
          for (a=0; a<20; ++a) sum += pb[a]*q.p[i][a];
          sum = 1.0/sqrt(sum);
          for (a=0; a<20; ++a) q.p[i][a]*=sum;
        }
      for (j=0; j<=t.L+1; j++)
        {
          float sum = 0.0;
          for (a=0; a<20; ++a) sum += pb[a]*t.p[j][a];
          sum = 1.0/sqrt(sum);
          for (a=0; a<20; ++a) t.p[j][a]*=sum;
        }
      break;

    case 11:  // log co-emission probability (no null model)
      for (a=0; a<20; ++a) pnul[a]=0.05;
      break;

   }

  // !!!!! ATTENTION!!!!!!!  after this t.p is not the same as after adding pseudocounts !!!
  //Introduce amino acid weights into template (for all but SOP scores)
  if (par.columnscore!=10)
    for (a=0; a<20; ++a)
      for (j=0; j<=t.L+1; j++)
        t.p[j][a]/=pnul[a];

  if (v>=4)
    {
      cout<<"\nAverage amino acid frequencies\n";
      cout<<"         A    R    N    D    C    Q    E    G    H    I    L    K    M    F    P    S    T    W    Y    V\n";
      cout<<"Q:    ";
      for (a=0; a<20; ++a) printf("%4.1f ",100*q.pav[a]);
      cout<<"\nT:    ";
      for (a=0; a<20; ++a) printf("%4.1f ",100*t.pav[a]);
      cout<<"\nNull: ";
      for (a=0; a<20; ++a) printf("%4.1f ",100*pnul[a]);
      cout<<"\npb:   ";
      for (a=0; a<20; ++a) printf("%4.1f ",100*pb[a]);
    }
 
  // cout<<"\nWeighted amino acid frequencies WITH pseudocounts:\n       A    R    N    D    C    Q    E    G    H    I    L    K    M    F    P    S    T    W    Y    V\n";
  // for (i=1; i<=L; ++i)
  //   {
  //     printf("%3i:  ",i);
  //     float sum=0;
  //     for (a=0; a<20; ++a)
  // 	{
  // 	  sum+=t.p[i][a];
  // 	  printf("%4.1f ",100*t.p[i][a]);
  // 	}
  //     printf("  sum=%5.3f\n",sum);
  //   }

  return;
}


/////////////////////////////////////////////////////////////////////////////////////
// Write HMM to output file
/////////////////////////////////////////////////////////////////////////////////////
void HMM::WriteToFile(char* outfile)
{
  const int SEQLEN=100;      // number of residues per line for sequences to be displayed
  int i,a;

  if (trans_lin==1) {fprintf(stderr,"Error: Writing transition pseudocounts in linear representation not allowed. Please report this error to the HHsearch developers.\n"); exit(6);}

  FILE *outf=NULL;
  if (strcmp(outfile,"stdout"))
    {
      if (par.append) outf=fopen(outfile,"a"); else outf=fopen(outfile,"w");
      if (!outf) OpenFileError(outfile);
    }
  else
    outf = stdout;
  if (v>=2) cout<<"Writing HMM to "<<outfile<<"\n";

//   fprintf(outf,"HHsearch HHM format 1.5\n");
  fprintf(outf,"HHsearch 1.6\n");         // format specification
  fprintf(outf,"NAME  %s\n",longname);    // name of first sequence
  fprintf(outf,"FAM   %s\n",fam);         // family name
  char file_nopath[NAMELEN];
  RemovePath(file_nopath,file);
  fprintf(outf,"FILE  %s\n",file_nopath); // base name of alignment file

  // Print command line
  fprintf(outf,"COM   ");
  for (int i=0; i<par.argc; i++)
    if (strlen(par.argv[i])<=100)
      fprintf(outf,"%s ",par.argv[i]);
    else
      fprintf(outf,"<%i characters> ",(int)strlen(par.argv[i]));
  fprintf(outf,"\n");

  // print out date stamp
  time_t* tp=new(time_t);
  *tp=time(NULL);
  fprintf(outf,"DATE  %s",ctime(tp));
  delete tp;

  // Print out some statistics of alignment
  fprintf(outf,"LENG  %i match states, %i columns in multiple alignment\n",L,l[L]);
  fprintf(outf,"FILT  %i out of %i sequences passed filter (-id %i -cov %i -qid %i -qsc %.2f -diff %i)\n",N_filtered,N_in,par.max_seqid,par.coverage,par.qid,par.qsc,par.Ndiff);
  fprintf(outf,"NEFF  %-4.1f\n",Neff_HMM);
  if (has_pseudocounts) { fprintf(outf,"PCT   true\n"); }

  // Print selected sequences from alignment (including secondary structure and confidence values, if known)
  fprintf(outf,"SEQ\n");
  for (int n=0; n<n_display; n++)
    {
      fprintf(outf,">%s\n",sname[n]);
      //first sequence character starts at 1; 0 not used.
      for(unsigned int j=0; j<strlen(seq[n]+1); j+=SEQLEN) fprintf(outf,"%-.*s\n",SEQLEN,seq[n]+1+j);
    }
  fprintf(outf,"#\n");

  // print null model background probabilities from substitution matrix
  fprintf(outf,"NULL   ");
  for (a=0; a<20; ++a) fout(outf,-iround(fast_log2(pb[s2a[a]])*HMMSCALE ));
  fprintf(outf,"\n");

  // print table header line with amino acids
  fprintf(outf,"HMM    ");
  for (a=0; a<20; ++a) fprintf(outf,"%1c\t",i2aa(s2a[a]));
  fprintf(outf,"\n");

  // print table header line with state transitions
  fprintf(outf,"       M->M\tM->I\tM->D\tI->M\tI->I\tD->M\tD->D\tNeff\tNeff_I\tNeff_D\n");

  // print out transition probabilities from begin state (virtual match state)
  fprintf(outf,"       ");
  for (a=0; a<=D2D; ++a) fout(outf,-iround(tr[0][a]*HMMSCALE));
  fout(outf,iround(Neff_M[0]*HMMSCALE));
  fout(outf,iround(Neff_I[0]*HMMSCALE));
  fout(outf,iround(Neff_D[0]*HMMSCALE));
  fprintf(outf,"\n");

  // Start loop for printing HMM columns
  int h=1;
  for (i=1; i<=L; ++i)
    {

      while(islower(seq[nfirst][h]) && seq[nfirst][h]) h++;
      fprintf(outf,"%1c %-4i ",seq[nfirst][h++],i);

      // Print emission probabilities for match state
      for (a=0; a<20; ++a) fout(outf,-iround(fast_log2(p[i][s2a[a]])*HMMSCALE ));
      fprintf(outf,"%-i",l[i]);
      fprintf(outf,"\n");

      // Print transition probabilities
      fprintf(outf,"       ");
      for (a=0; a<=D2D; ++a) fout(outf,-iround(tr[i][a]*HMMSCALE));
      fout(outf,iround(Neff_M[i]*HMMSCALE));
      fout(outf,iround(Neff_I[i]*HMMSCALE));
      fout(outf,iround(Neff_D[i]*HMMSCALE));
      fprintf(outf,"\n\n");
    } // end for(i)-loop for printing HMM columns

  fprintf(outf,"//\n");
  fclose(outf);
}

/////////////////////////////////////////////////////////////////////////////////////
// Write HMM to output file
/////////////////////////////////////////////////////////////////////////////////////
void HMM::InsertCalibration(char* infile)
{
  char* line =  new(char[LINELEN]);    // input line
  char** lines = new(char*[3*L+100000]);
  int nline=0;
  int l;
  char done=0;   // inserted new 'EVD mu sigma' line?

  // Read from infile all lines and insert the EVD line with lamda and mu coefficients
  ifstream inf;
  inf.open(infile, ios::in);
  if (!inf) OpenFileError(infile);
  if (v>=2) cout<<"Recording calibration coefficients in "<<infile<<"\n";

  while (inf.getline(line,LINELEN) && !(line[0]=='/' && line[1]=='/') && nline<2*MAXRES)
    {

      // Found an EVD lamda mu line? -> remove
      while (!done && !strncmp(line,"EVD ",3) && !(line[0]=='/' && line[1]=='/') && nline<2*MAXRES)
        inf.getline(line,LINELEN);
      if ((line[0]=='/' && line[1]=='/') || nline>=2*MAXRES)
        {fprintf(stderr,"Error: wrong format in %s. Expecting hhm format\n",infile); exit(1);}

      // Found the SEQ line? -> insert calibration before this line
      if (!done && (!strncmp("SEQ",line,3) || !strncmp("HMM",line,3)) && (isspace(line[3]) || line[3]=='\0'))
        {
          done=1;
          lines[nline]=new(char[128]);
          if (!lines[nline]) MemoryError("space to read in HHM file for calibration");
          sprintf(lines[nline],"EVD   %-7.4f %-7.4f",lamda,mu);
          nline++;
        }
      lines[nline]=new(char[strlen(line)+1]);
      if (!lines[nline]) MemoryError("space to read in HHM file for calibration");
      strcpy (lines[nline],line);
      nline++;
    }
  inf.close();

  // Write to infile all lines
  FILE* infout=fopen(infile,"w");
  if (!infout) {
    cerr<<endl<<"WARNING in "<<program_name<<": no calibration coefficients written to "<<infile<<":\n";
    cerr<<"Could not open file for writing.\n";
    return;
  }
  for (l=0; l<nline; l++) {fprintf(infout,"%s\n",lines[l]); delete[] lines[l];}
  fprintf(infout,"//\n");
  fclose(infout);
  delete[] line;
  delete[] lines;
  return;
}

/////////////////////////////////////////////////////////////////////////////////////
// Transform log to lin transition probs
/////////////////////////////////////////////////////////////////////////////////////
void HMM::Log2LinTransitionProbs(float beta)
{
  if (trans_lin==1) return;
  trans_lin=1;
  for (int i=0; i<=L; ++i)
    {
      for (int a=0; a<NTRANS; ++a)
        tr[i][a] = fpow2(beta*tr[i][a]);
    }
}


/////////////////////////////////////////////////////////////////////////////////////
// Set query columns in His-tags etc to Null model distribution
/////////////////////////////////////////////////////////////////////////////////////
void HMM::NeutralizeTags()
{
  char* qseq = seq[nfirst];
  char* pt;
  int a,i;

  // Neutralize His tag
  if ( (pt=strstr(qseq,"HHHHH")) )
    {
      int i0 = pt-qseq+1;
      if (v>=2) printf("Neutralized His-tag at position %i\n",i0);
      for (i=imax(i0-5,1); i<i0; ++i)   // neutralize leading 5 columns
        for (a=0; a<NAA; ++a) p[i][a]=pb[a];
      for (; (*pt)!='H'; ++i,++pt)      // neutralize His columns
        for (a=0; a<NAA; ++a) p[i][a]=pb[a];
      i0=i;
       for (; i<imin(i0+5,L+1); ++i)    // neutralize trailing 5 columns
        for (a=0; a<NAA; ++a) p[i][a]=pb[a];
       if (v>=3) printf("start:%i  end:%i\n",imax(i0-5,1),i-1);
    }

  // Neutralize C-myc tag
  if ( (pt=strstr(qseq,"EQKLISEEDL")) )
    {
      if (v>=2) printf("Neutralized C-myc-tag at position %i\n",int(pt-qseq)+1);
      for (i=pt-qseq+1; i<=pt-qseq+10; ++i)
        for (a=0; a<NAA; ++a) p[i][a]=pb[a];
    }
  // Neutralize FLAG tag
  if ( (pt=strstr(qseq,"DYKDDDDK")) )
    {
      if (v>=2) printf("Neutralized FLAG-tag at position %i\n",int(pt-qseq)+1);
      for (i=pt-qseq+1; i<=pt-qseq+8; ++i)
        for (a=0; a<NAA; ++a) p[i][a]=pb[a];
    }
}



/////////////////////////////////////////////////////////////////////////////////////
// Calculate effective number of sequences using profiles INCLUDING pseudocounts
/////////////////////////////////////////////////////////////////////////////////////
float HMM::CalcNeff()
{
  float Neff=0;
  for (int i=1; i<=L; ++i)
    for (int a=0; a<20; ++a)
      if (p[i][a]>1E-10) Neff-=p[i][a]*fast_log2(p[i][a]);
  return fpow2(Neff/L);
}

/////////////////////////////////////////////////////////////////////////////////////
// Add secondary structure prediction to alignment (overwrite existing prediction)
/////////////////////////////////////////////////////////////////////////////////////
void HMM::AddSSPrediction(char seq_pred[], char seq_conf[])
{
  unsigned int i;

  if ((int)strlen(seq_pred)!=L+1)
    {
      cerr<<"WARNING! Could not add secondary struture prediction - unequal length!\n";
      return;
    }

  if (nss_pred < 0)  // No ss prediction exists
    {
      nss_pred=n_display;
      seq[nss_pred]=new(char[L+2]);
      strcpy(seq[nss_pred],seq_pred);
      ss_pred=new(char[L+2]);
      for (i=0; i<strlen(seq_pred); i++) ss_pred[i]=ss2i(seq_pred[i]);
      sname[nss_pred]=new(char[50]);
      strcpy(sname[nss_pred],"ss_pred PSIPRED predicted secondary structure");
      n_display++;
    }
  else  // overwrite existing ss prediction
    {
      strcpy(seq[nss_pred],seq_pred);
      for (i=0; i<strlen(seq_pred); i++) ss_pred[i]=ss2i(seq_pred[i]);
    }

  if (nss_conf < 0)  // No ss prediction exists
    {
      nss_conf=n_display;
      seq[nss_conf]=new(char[L+2]);
      strcpy(seq[nss_conf],seq_conf);
      ss_conf=new(char[L+2]);
      for (i=0; i<strlen(seq_conf); i++) ss_conf[i]=cf2i(seq_conf[i]);
      sname[nss_conf]=new(char[50]);
      strcpy(sname[nss_conf],"ss_conf PSIPRED confidence values");
      n_display++;
    }
  else  // overwrite existing ss prediction
    {
      strcpy(seq[nss_conf],seq_conf);
      for (i=0; i<strlen(seq_conf); i++) ss_conf[i]=cf2i(seq_conf[i]);
    }

}


// #define Weff(Neff) (1.0+par.neffa*(Neff-1.0)+(par.neffb-4.0*par.neffa)/16.0*(Neff-1.0)*(Neff-1.0))


// /////////////////////////////////////////////////////////////////////////////////////
// // Normalize probabilities in total merged super-HMM
// /////////////////////////////////////////////////////////////////////////////////////
// void HMM::NormalizeHMMandTransitionsLin2Log()
// {
//   int i;      // position in query
//   int a;      // amino acid
//   for (i=0; i<=L+1; i++)
//     {
//       float sum=0.0;
//       for (a=0; a<20; a++) sum += f[i][a];
//       for (a=0; a<20; a++) f[i][a]/=sum;
//       sum = tr_lin[i][M2M] + tr_lin[i][M2I] + tr_lin[i][M2D];
//       tr_lin[i][M2M] /= sum;
//       tr_lin[i][M2I] /= sum;
//       tr_lin[i][M2D] /= sum;
//       tr[i][M2M] = fast_log2(tr_lin[i][M2M]);
//       tr[i][M2I] = fast_log2(tr_lin[i][M2I]);
//       tr[i][M2D] = fast_log2(tr_lin[i][M2D]);
//       sum = tr_lin[i][D2M] + tr_lin[i][D2D];
//       tr_lin[i][D2M] /= sum;
//       tr_lin[i][D2D] /= sum;
//       tr[i][D2M] = fast_log2(tr_lin[i][D2M]);
//       tr[i][D2D] = fast_log2(tr_lin[i][D2D]);
//       sum = tr_lin[i][I2M] + tr_lin[i][I2I];
//       tr_lin[i][I2M] /= sum;
//       tr_lin[i][I2I] /= sum;
//       tr[i][I2M] = fast_log2(tr_lin[i][I2M]);
//       tr[i][I2I] = fast_log2(tr_lin[i][I2I]);
//    }
// }


// UNCOMMENT TO ACTIVATE COMPOSITIONALLY BIASED PSEUDOCOUNTS BY RESCALING THE RATE MATRIX

// /////////////////////////////////////////////////////////////////////////////////////
// //// Function to minimize
// /////////////////////////////////////////////////////////////////////////////////////
// double RescaleMatrixFunc(double x[])
// {
//   double sum=0.0;
//   for (int a=0; a<20; ++a)
//     {
//       double za=0.0;
//       for (int b=0; b<20; ++b) za+=P[a][b]*x[b];
//       sum += (x[a]*za-qav[a])*(x[a]*za-qav[a]);
//     }
//   return sum;
// }

// /////////////////////////////////////////////////////////////////////////////////////
// //// Gradient of function to minimize
// /////////////////////////////////////////////////////////////////////////////////////
// void RescaleMatrixFuncGrad(double grad[], double x[])
// {
//   double z[20] = {0.0};
//   double w[20];
//   double tmp;
//   for (int a=0; a<20; ++a)
//     for (int b=0; b<20; ++b) z[a] += P[a][b]*x[b];

//   for (int a=0; a<20; ++a) w[a] = x[a]*z[a]-qav[a];
//   for (int a=0; a<20; ++a)
//     {
//       tmp = w[a]*z[a];
//       for (int b=0; b<20; ++b) tmp += P[a][b]*x[b]*w[b];
//       grad[a] = 2.0*tmp;
//     }
//   return;
// }


// /////////////////////////////////////////////////////////////////////////////////////
// //// Rescale a substitution matrix to biased aa frequencies in global vector qav[a]
// /////////////////////////////////////////////////////////////////////////////////////
// void HMM::RescaleMatrix()
// {
//   int a,b;
//   int code;
//   double x[21];     // scaling factor
//   double val_min;
//   const int len=20;
//   const int max_iterations=50;

//   if (v>=2) printf("Adjusting rate matrix to query amino acid composition ...\n");

//   // Put amino acid frequencies into global array (needed to call WNLIB's conjugate gradient method)
//   for (a=0; a<20; ++a) qav[a] = pav[a];

//   // Initialize scaling factors x[a]
//   for (a=0; a<20; ++a) x[a]=pow(qav[a]/pb[a],0.73);  // Initialize

//   // Call conjugate gradient minimization method from WNLIB
//   wn_conj_gradient_method(&code,&val_min,x,len,&RescaleMatrixFunc,&RescaleMatrixFuncGrad,max_iterations);


//   // Calculate z[a] = sum_b Pab*xb
//   float sum_err=0.0f;
//   float sum = 0.0f;
//   for (a=0; a<20; ++a)
//     {
//       float za=0.0f; // za = sum_b Pab*xb
//       for (b=0; b<20; ++b) za+=P[a][b]*x[b];
//       sum_err += (x[a]*za/qav[a]-1)*(x[a]*za/qav[a]-1);
//       sum += x[a]*za;
//    }
//   if (sum_err>1e-3 & v>=1) fprintf(stderr,"WARNING: adjusting rate matrix by CG resulted in residual error of %5.3f.\n",sum_err);

//   // Rescale rate matrix
//   for (a=0; a<20; ++a)
//     for (b=0; b<20; ++b)
//       {
//      P[a][b] *= x[a]*x[b]/sum;
//      R[a][b] = P[a][b]/qav[b];
//       }

//   // How well approximated?
//   if (v>=3)
//     {
//       // Calculate z[a] = sum_b Pab*xb
//       float z[21];
//       for (a=0; a<20; ++a)
//      for (z[a]=0.0, b=0; b<20; ++b) z[a]+=P[a][b];
//       printf("Adjust   A    R    N    D    C    Q    E    G    H    I    L    K    M    F    P    S    T    W    Y    V\nErr?  ");
//       for (a=0; a<20; ++a) printf("%4.0f ",1000*z[a]/qav[a]);
//       cout<<endl<<"xa    ";
//       for (a=0; a<20; ++a) fprintf(stdout,"%4.2f ",x[a]);
//       cout<<endl;
//     }

//   // Evaluate sequence identity underlying substitution matrix
//   if (v>=3)
//     {
//       float id=0.0f;
//       float entropy=0.0f;
//       float entropy_qav=0.0f;
//       float mut_info=0.0f;
//       for (a=0; a<20; ++a) id += P[a][a];
//       for (a=0; a<20; ++a)  entropy_qav-=qav[a]*fast_log2(qav[a]);
//       for (a=0; a<20; ++a)
//        for (b=0; b<20; ++b)
//          {
//            entropy-=P[a][b]*fast_log2(R[a][b]);
//            mut_info += P[a][b]*fast_log2(P[a][b]/qav[a]/qav[b]);
//          }

//       fprintf(stdout,"Rescaling rate matrix: sequence identity = %2.0f%%; entropy per column = %4.2f bits (out of %4.2f); mutual information = %4.2f bits\n",100*id,entropy,entropy_qav,mut_info);
//     }
//   return;
// }

