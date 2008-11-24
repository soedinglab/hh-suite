// hhfullalignment.C

#ifndef MAIN
#define MAIN
#include <iostream>   // cin, cout, cerr
#include <fstream>    // ofstream, ifstream 
#include <stdio.h>    // printf
#include <stdlib.h>   // exit
#include <string>     // strcmp, strstr
#include <math.h>     // sqrt, pow
#include <limits.h>   // INT_MIN
#include <float.h>    // FLT_MIN
#include <time.h>     // clock
#include <ctype.h>    // islower, isdigit etc
using std::ios;
using std::ifstream;
using std::ofstream;
using std::cout;
using std::cerr;
using std::endl;
#include "util.C"     // imax, fmax, iround, iceil, ifloor, strint, strscn, strcut, substr, uprstr, uprchr, Basename etc.
#include "list.h"     // list data structure
#include "hash.h"     // hash data structure
#include "hhdecl.C"      // constants, class 
#include "hhutil.C"      // imax, fmax, iround, iceil, ifloor, strint, strscn, strcut, substr, uprstr, uprchr, Basename etc.
#include "hhhmm.h"       // class HMM
#include "hhalignment.h" // class Alignment
#include "hhhit.h"
#include "hhhalfalignment.h"
#endif

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
// Methods of class FullAlignment
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
  

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Output Results: use classes HalfAlignment and FullAlignment
//
// Example:
// Each column list contains at least a match state 
// Insert states between the match states are omitted
//
//  step 19 18 17 16 15 14 13 12 11 10  9  8  7  6  5  4  3  2  1
//     i  0  0  1  2  3  4  5  6  7  8  9  9  9  9 10 11 12 13 14
// 
//     Q  ~  ~  X  X  X  X  X  X  X  X  X  ~  ~  ~  X  X  X  X  X
//     T  Y  Y  Y  Y  Y  Y  Y  Y  ~  Y  Y  Y  Y  Y  Y  Y  Y  ~  ~
//
//     j  7  8  9 10 11 12 13 14 14 15 16 17 18 19 20 21 22 22 22
// state IM IM MM MM MM MM MM MM DG MM MM GD GD GD MM MM MM MI MI
//
// nsteps=19
//


/////////////////////////////////////////////////////////////////////////////////////
// Constructor
FullAlignment::FullAlignment(int maxseqdis) 
{
  qa = new HalfAlignment(maxseqdis);
  ta = new HalfAlignment(maxseqdis);
  for (int h=0; h<LINELEN-1; h++) symbol[h]=' ';
}

/////////////////////////////////////////////////////////////////////////////////////
// Destructor
FullAlignment::~FullAlignment() 
{
  delete qa;
  delete ta;
}

/////////////////////////////////////////////////////////////////////////////////////
// Free memory of arrays s[][], l[][], and m[][] in HalfAlignment
void FullAlignment::FreeMemory()
{
  qa->Unset();
  ta->Unset();
}

/////////////////////////////////////////////////////////////////////////////////////
// Add columns for match (and delete) states. 
/////////////////////////////////////////////////////////////////////////////////////
void FullAlignment::AddColumns(int i, int j, char prev_state, char state, float S)
{
  switch(state)
    {
    case MM:  //MM pair state (both query and template in Match state)
      AddGaps(); //fill up gaps until query and template parts have same length
      symbol[qa->pos] =ScoreChr(S);
      qa->AddColumn(i);
      ta->AddColumn(j);
      qa->AddInsertsAndFillUpGaps(i);
      ta->AddInsertsAndFillUpGaps(j);
      break;

    case GD: //-D state
      if (prev_state==DG) AddGaps();
      symbol[ta->pos]='Q';
      ta->AddColumn(j); //query has gap -> add nothing
      ta->AddInsertsAndFillUpGaps(j);
      break;
    case IM: //IM state
      if (prev_state==MI) AddGaps();
      symbol[ta->pos]='Q';
      ta->AddColumn(j); //query has gap -> add nothing
      ta->AddInsertsAndFillUpGaps(j);
      break;

    case DG: //D- state
      if (prev_state==GD) AddGaps();
      symbol[qa->pos]='T';
      qa->AddColumn(i);//template has gap -> add nothing
      qa->AddInsertsAndFillUpGaps(i);
      break;
    case MI: //MI state
      if (prev_state==IM) AddGaps();
      symbol[qa->pos]='T';
      qa->AddColumn(i);//template has gap -> add nothing
      qa->AddInsertsAndFillUpGaps(i);
      break;
    }
}

/////////////////////////////////////////////////////////////////////////////////////
// Fill up gaps until query and template parts have same length
/////////////////////////////////////////////////////////////////////////////////////
void FullAlignment::AddGaps()
{
  while (qa->pos<ta->pos) qa->AddChar('.');
  while (ta->pos<qa->pos) ta->AddChar('.');
}


/////////////////////////////////////////////////////////////////////////////////////
// Build full alignment -> qa->s[k][h] and ta->s[k][h]
/////////////////////////////////////////////////////////////////////////////////////
void FullAlignment::Build(HMM& q, Hit& hit)
{
  int step;
  char prev_state=MM, state=MM;  
  int n;
  int hh;
  int k;
  identities=0;       // number of identical residues in query and template sequence
  score_sim=0.0f;     // substitution matrix similarity score between query and template

  ClearSymbols();
  
  // Set up half-alignments 
  // n is the sequence index up to which sequences are prepared for display
  n = imin(  q.n_display,par.nseqdis+(  q.nss_dssp>=0)+(  q.nsa_dssp>=0)+(  q.nss_pred>=0)+(  q.nss_conf>=0)+(  q.ncons>=0));
  qa->Set(  q.name,  q.seq,  q.sname, n,  q.L,  q.nss_dssp , q.nss_pred , q.nss_conf,  q.nsa_dssp, q.ncons);
  n = imin(hit.n_display,par.nseqdis+(hit.nss_dssp>=0)+(hit.nsa_dssp>=0)+(hit.nss_pred>=0)+(hit.nss_conf>=0)+(hit.ncons>=0));
  ta->Set(hit.name,hit.seq,hit.sname, n,hit.L,hit.nss_dssp,hit.nss_pred,hit.nss_conf,hit.nsa_dssp, hit.ncons);

//   printf("HMM: %s\nstep nst   i   j state hq  ht\n",hit.name);

  for (step=hit.nsteps; step>=1; step--) 
  {
    prev_state = state;
    state = hit.states[step];
    
    // Add column to alignment and compute identities and sequence-sequence similarity score
    AddColumns(hit.i[step],hit.j[step],prev_state,state,hit.S[step]);
    if (state==MM) 
      {
	char qc=qa->seq[  q.nfirst][ qa->m[  q.nfirst][hit.i[step]] ];
	char tc=ta->seq[hit.nfirst][ ta->m[hit.nfirst][hit.j[step]] ];
	if (qc==tc) identities++;  // count identical amino acids
	score_sim += S[(int)aa2i(qc)][(int)aa2i(tc)];
 	//fprintf(stderr,"%3i %3i  %3i %3i  %3i %1c %1c %6.2f %6.2f %6.2f %6.2f  \n",step,hit.nsteps,hit.i[step],hit.j[step],int(state),qc,tc,S[(int)aa2i(qc)][(int)aa2i(tc)],score_sim,hit.P_posterior[step],hit.sum_of_probs); //DEBUG (P_posterior not defined for Viterbi!)
      }    
  }
  
  AddGaps(); //fill up gaps until query and template parts have same length
  qa->AddChar('\0');
  ta->AddChar('\0');

  // Change gap symbol '.' (gap aligned to insert) to '~' if one HMM has gap with respect to other HMM
  for (hh=1; hh<qa->pos; hh++) 
    {  
      if (symbol[hh]=='Q')
	{
	  // Gap in query (IM or GD state)
	  symbol[hh]=' ';
	  for (k=0; k<qa->n; k++) if (qa->s[k][hh]=='.') qa->s[k][hh]='-';
	} 
      else if (symbol[hh]=='T') 
	{
	  // Gap in target (MI or DG state)
	  symbol[hh]=' ';
	  for (k=0; k<ta->n; k++) if (ta->s[k][hh]=='.') ta->s[k][hh]='-';
	}
    }
}


/////////////////////////////////////////////////////////////////////////////////////
// Print out header before full alignment 
/////////////////////////////////////////////////////////////////////////////////////
void FullAlignment::PrintHeader(FILE* outf, HMM& q, Hit& hit)
{
  fprintf(outf,">%s\n",hit.longname); 
  fprintf(outf,"Probab=%-.2f  E-value=%-.2g  Score=%-.2f  Aligned_cols=%i  Identities=%i%%  Similarity=%-.3f  Sum_probs=%.1f\n\n",
     hit.Probab,hit.Eval,hit.score,hit.matched_cols,iround(100.0*identities/hit.matched_cols),score_sim/hit.matched_cols,hit.sum_of_probs);
}

/////////////////////////////////////////////////////////////////////////////////////
// Print out full alignment in HHR format
/////////////////////////////////////////////////////////////////////////////////////
void FullAlignment::PrintHHR(FILE* outf, Hit& hit)
{
  const int NLEN=14;  //Length of name field in front of multiple alignment
  int h=0;    //counts position (column) in alignment
  int hh=0;   //points to column at start of present output block of alignment
  int k;      //counts sequences in query and template 
  short unsigned int lq[MAXSEQ]; // lq[k] counts index of residue from query sequence k to be printed next;
  short unsigned int lt[MAXSEQ]; // lt[k] counts index of residue from template sequence k to be printed next;
  char namestr[NAMELEN]; //name of sequence
  int iq=hit.i1;  // match state counter for query HMM (displayed in consensus line)
  int jt=hit.j1;  // match state counter for template HMM (displayed in consensus line)

  for (k=0; k<qa->n; k++) lq[k]=qa->l[k][hit.i1];
  for (k=0; k<ta->n; k++) lt[k]=ta->l[k][hit.j1];
  
  while (hh<ta->pos-1) // print alignment block 
    {
      // Print query secondary structure sequences
      for (k=0; k<qa->n; k++)
	{
	  if (k==qa->nsa_dssp) continue;
 	  if (!(k==qa->nss_dssp || k==qa->nsa_dssp || k==qa->nss_pred || k==qa->nss_conf)) continue;
	  if (k==qa->nss_dssp && !par.showdssp) continue;
	  if ((k==qa->nss_pred || k==qa->nss_conf) && !par.showpred) continue;
	  strncpy(namestr,qa->sname[k],NAMELEN-2);
	  namestr[NAMELEN-1]='\0';
	  strcut(namestr);
	  fprintf(outf,"Q %-*.*s      ",NLEN,NLEN,namestr);
	  for (h=hh; h<imin(hh+par.aliwidth,qa->pos-1); h++) fprintf(outf,"%1c",qa->s[k][h]);
	  fprintf(outf,"\n");
	}

      // Print query sequences
      for (k=0; k<qa->n; k++)
	{
 	  if (k==qa->nss_dssp || k==qa->nsa_dssp || k==qa->nss_pred || k==qa->nss_conf || k==qa->ncons) continue;
	  strncpy(namestr,qa->sname[k],NAMELEN-2);
	  namestr[NAMELEN-1]='\0';
	  strcut(namestr);
	  fprintf(outf,"Q %-*.*s %4i ",NLEN,NLEN,namestr,lq[k]);
	  for (h=hh; h<imin(hh+par.aliwidth,qa->pos-1); h++) 
	    {fprintf(outf,"%1c",qa->s[k][h]); lq[k]+=WordChr(qa->s[k][h]);}  //WordChr() returns 1 if a-z or A-Z; 0 otherwise
	  fprintf(outf," %4i (%i)\n",lq[k]-1,qa->l[k][qa->L+1]);
	}

      // Print query consensus sequence
      if (par.showcons && qa->ncons>=0)
	{
	  k=qa->ncons; 
	  strncpy(namestr,qa->sname[k],NAMELEN-2);
	  namestr[NAMELEN-1]='\0';
	  strcut(namestr);
	  fprintf(outf,"Q %-*.*s %4i ",NLEN,NLEN,namestr,iq);
	  for (h=hh; h<imin(hh+par.aliwidth,qa->pos-1); h++) 
	    {
	      if (qa->s[k][h]=='x') qa->s[k][h]='~'; 
	      if (qa->s[k][h]!='-' && qa->s[k][h]!='.') iq++;
	      fprintf(outf,"%1c",qa->s[k][h]); 
	    }
	  fprintf(outf," %4i (%i)\n",iq-1,qa->L);
	}


      // Print symbols representing the score
      fprintf(outf,"  %*.*s      ",NLEN,NLEN," ");
      for (h=hh; h<imin(hh+par.aliwidth,qa->pos-1); h++) fprintf(outf,"%1c",symbol[h]);
      fprintf(outf,"\n");

      // Print template consensus sequence
      if (par.showcons && ta->ncons>=0)
	{
	  k=ta->ncons; 
	  strncpy(namestr,ta->sname[k],NAMELEN-2);
	  namestr[NAMELEN-1]='\0';
	  strcut(namestr);
	  fprintf(outf,"T %-*.*s %4i ",NLEN,NLEN,namestr,jt);
	  for (h=hh; h<imin(hh+par.aliwidth,ta->pos-1); h++) 
	    {
	      if (ta->s[k][h]=='x') ta->s[k][h]='~'; 
	      if (ta->s[k][h]!='-' && ta->s[k][h]!='.') jt++;
	      fprintf(outf,"%1c",ta->s[k][h]); 
	    }
	  fprintf(outf," %4i (%i)\n",jt-1,ta->L);
	}
      // Print template sequences
      for (k=0; k<ta->n; k++)
	{
	  if (k==ta->nss_dssp || k==ta->nsa_dssp || k==ta->nss_pred || k==ta->nss_conf || k==ta->ncons) continue;
	  strncpy(namestr,ta->sname[k],NAMELEN-2);
	  namestr[NAMELEN-1]='\0';
	  strcut(namestr);
	  fprintf(outf,"T %-*.*s %4i ",NLEN,NLEN,namestr,lt[k]);
	  for (h=hh; h<imin(hh+par.aliwidth,ta->pos-1); h++) 
	    {fprintf(outf,"%1c",ta->s[k][h]); lt[k]+=WordChr(ta->s[k][h]);}  //WordChr() returns 1 if a-z or A-Z; 0 otherwise
	  fprintf(outf," %4i (%i)\n",lt[k]-1,ta->l[k][ta->L+1]);
	}

      // Print template secondary structure sequences
      for (k=0; k<ta->n; k++)
	{
	  if (k==ta->nsa_dssp) continue;
	  if (!(k==ta->nss_dssp || k==ta->nss_pred || k==ta->nss_conf)) continue;
	  if (k==ta->nss_dssp && !par.showdssp) continue;
	  if ((k==ta->nss_pred || k==ta->nss_conf)&& !par.showpred) continue;
	  strncpy(namestr,ta->sname[k],NAMELEN-2);
	  namestr[NAMELEN-1]='\0';
	  strcut(namestr);
	  fprintf(outf,"T %-*.*s      ",NLEN,NLEN,namestr);
	  for (h=hh; h<imin(hh+par.aliwidth,ta->pos-1); h++) fprintf(outf,"%1c",ta->s[k][h]); 
	  fprintf(outf,"\n");
	}

      hh=h;
      fprintf(outf,"\n\n");
    } 
}

/////////////////////////////////////////////////////////////////////////////////////
// Print out full alignment in A2M format
/////////////////////////////////////////////////////////////////////////////////////
void FullAlignment::PrintA2M(FILE* outf, Hit& hit)
{
  int k;      //counts sequences in query and template 
  int h,hh; 

  // Print query sequences
  for (k=0; k<qa->n; k++)
    {
      if (k==qa->nsa_dssp) continue;
      if (k==qa->nss_dssp && !par.showdssp) continue;
      if ((k==qa->nss_pred || k==qa->nss_conf) && !par.showpred) continue;
      if (k==qa->ncons && !par.showcons) continue;
      fprintf(outf,">%s\n",qa->sname[k]);
      for (h=0,hh=-par.aliwidth; qa->s[k][h]>0; h++,hh++) 
	{
	  if (!hh) {fprintf(outf,"\n"); hh-=par.aliwidth;}
	  fprintf(outf,"%1c",qa->s[k][h]); 
	}
      fprintf(outf,"\n");
    }
  
  // Print template sequences
  for (k=0; k<ta->n; k++)
    {
      if (k==ta->nsa_dssp) continue;
      if (k==ta->nss_dssp && !par.showdssp) continue;
      if ((k==ta->nss_pred || k==ta->nss_conf) && !par.showpred) continue;
      if (k==ta->ncons && !par.showcons) continue;
      fprintf(outf,">%s\n",ta->sname[k]);
      for (h=0,hh=-par.aliwidth; ta->s[k][h]>0; h++,hh++) 
	{
	  if (!hh) {fprintf(outf,"\n"); hh-=par.aliwidth;}
	  fprintf(outf,"%1c",ta->s[k][h]); 
	}
      fprintf(outf,"\n");
    }
  fprintf(outf,"\n");
}

/////////////////////////////////////////////////////////////////////////////////////
// Print out full alignment in A2M format
/////////////////////////////////////////////////////////////////////////////////////
void FullAlignment::PrintFASTA(FILE* outf, Hit& hit)
{
  // Transform sequences to uppercase and '.' to '-'
  qa->ToFASTA();
  ta->ToFASTA();
  PrintA2M(outf,hit);
}

/////////////////////////////////////////////////////////////////////////////////////
// Print out full alignment in A2M format
/////////////////////////////////////////////////////////////////////////////////////
void FullAlignment::PrintA3M(FILE* outf, Hit& hit)
{
  // Remove all '.' from sequences
  qa->RemoveChars('.');
  ta->RemoveChars('.');
  PrintA2M(outf,hit);
}






