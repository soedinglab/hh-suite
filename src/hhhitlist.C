// hhhitlist.C

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
#include "hhfullalignment.h"
#endif


/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
//// Methods of class HitList
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////////////////////////////////
// Print summary listing of hits
/////////////////////////////////////////////////////////////////////////////////////
void HitList::PrintHitList(HMM& q, char* outfile)
{
  Hit hit;
  int nhits=0;

  FILE* outf=NULL;
  if (strcmp(outfile,"stdout"))
    {
      outf=fopen(outfile,"w");
      if (!outf) OpenFileError(outfile);
    }
  else
    outf = stdout;


  fprintf(outf,"Query         %s\n",q.longname); 
//  fprintf(outf,"Family        %s\n",q.fam); 
  fprintf(outf,"Match_columns %i\n",q.L);
  fprintf(outf,"No_of_seqs    %i out of %i\n",q.N_filtered,q.N_in);
  fprintf(outf,"Neff          %-4.1f\n",q.Neff_HMM);
  fprintf(outf,"Searched_HMMs %i\n",N_searched);
  
  // Print date stamp
  time_t* tp=new(time_t);
  *tp=time(NULL);
  fprintf(outf,"Date          %s",ctime(tp));
  delete (tp);
  
  // Print command line
  fprintf(outf,"Command       "); 
  for (int i=0; i<par.argc; i++) 
    if (strlen(par.argv[i])<=par.maxdbstrlen)
      fprintf(outf,"%s ",par.argv[i]);
    else
      fprintf(outf,"<%i characters> ",(int)strlen(par.argv[i]));
  fprintf(outf,"\n\n"); 
   
#ifdef WINDOWS   
    fprintf(outf," No Hit                             Prob  E-value  P-value  Score    SS Cols Query HMM  Template HMM\n");
#else
    fprintf(outf," No Hit                             Prob E-value P-value  Score    SS Cols Query HMM  Template HMM\n");
#endif

  Reset();
  while (!End()) // print hit list
    {
      hit = ReadNext();
      if (nhits>=par.Z) break;       //max number of lines reached?
      if (nhits>=par.z && hit.Probab < par.p) break;
      if (nhits>=par.z && hit.Eval > par.E) continue;
//       if (hit.matched_cols <=1) continue; // adding this might get to intransparent... analogous statement in PrintAlignments
      nhits++;

      char Estr[10];
      char Pstr[10];
      char str[NAMELEN];
      sprintf(str,"%3i %-30.30s    ",nhits,hit.longname);
      

// #ifdef WINDOWS      
// 	 fprintf(outf,"%-34.34s %5.1f %8.2G %8.2G %6.1f %5.1f %4i ",str,hit.Probab,hit.Eval,hit.Pval,hit.score,hit.score_ss,hit.matched_cols);
// #else
// 	 fprintf(outf,"%-34.34s %5.1f %7.2G %7.2G %6.1f %5.1f %4i ",str,hit.Probab,hit.Eval,hit.Pval,hit.score,hit.score_ss,hit.matched_cols);
// #endif
#ifdef WINDOWS      
      if (hit.Eval>=1E-99) sprintf(Estr,"%8.2G",hit.Eval); else sprintf(Estr,"%8.1G",hit.Eval);
      if (hit.Pval>=1E-99) sprintf(Pstr,"%8.2G",hit.Pval); else sprintf(Pstr,"%8.1G",hit.Pval);
      fprintf(outf,"%-34.34s %5.1f %8s %8s ",str,hit.Probab,Estr,Pstr);
#else
      if (hit.Eval>=1E-99) sprintf(Estr,"%7.2G",hit.Eval); else sprintf(Estr,"%7.1G",hit.Eval);
      if (hit.Pval>=1E-99) sprintf(Pstr,"%7.2G",hit.Pval); else sprintf(Pstr,"%7.1G",hit.Pval);
      fprintf(outf,"%-34.34s %5.1f %7s %7s ",str,hit.Probab,Estr,Pstr);
#endif

      // Needed for long sequences (more than 5 digits in length)
      sprintf(str,"%6.1f",hit.score);
      fprintf(outf,"%-6.6s %5.1f %4i %4i-%-4i %4i-%-4i(%i)\n",str,hit.score_ss,hit.matched_cols,hit.i1,hit.i2,hit.j1,hit.j2,hit.L);
   } //end print hit list 

  fprintf(outf,"\n");
  if (strcmp(outfile,"stdout")) fclose(outf);
}



/////////////////////////////////////////////////////////////////////////////////////
// Print alignments of query sequences against hit sequences 
/////////////////////////////////////////////////////////////////////////////////////
void HitList::PrintAlignments(HMM& q, char* outfile, char outformat)
{
  Hit hit;
  FullAlignment qt_ali(par.nseqdis+10); // maximum 10 annotation (pseudo) sequences (ss_dssp, sa_dssp, ss_pred, ss_conf, consens,...)
  int nhits=0;

  FILE* outf=NULL;
  if (strcmp(outfile,"stdout"))
    {
      if (outformat==0)
	outf=fopen(outfile,"a"); //append to summary hitlist
      else 
	outf=fopen(outfile,"w"); //open for writing
      if (!outf) OpenFileError(outfile);
    }
  else
    outf = stdout;

  Reset();
  while (!End()) // print hit list
    {
      if (nhits>=par.B) break;       //max number of lines reached?
      hit = ReadNext();
      if (nhits>=par.b && hit.Probab < par.p) break;
      if (nhits>=par.b && hit.Eval > par.E) continue;
//       if (hit.matched_cols <=1) continue; // adding this might get to intransparent... analogous statement in PrintHitlist and hhalign.C
      nhits++;
      
      // Build double alignment of query against template sequences
      qt_ali.Build(q,hit);


      // Print out alignment 
      if (outformat==0) // HHR format
	{
	  fprintf(outf,"No %-3i\n",nhits);
	  qt_ali.PrintHeader(outf,q,hit);
	  qt_ali.PrintHHR(outf,hit);
	}
      else if (outformat==1) // FASTA format
	{
	  fprintf(outf,"# No %-3i\n",nhits);
	  qt_ali.PrintFASTA(outf,hit);
	}
      else if(outformat==2) // A2M format
	{
	  fprintf(outf,"# No %-3i\n",nhits);
	  qt_ali.PrintA2M(outf,hit);
	}
      else // A3m format
	{
	  fprintf(outf,"# No %-3i\n",nhits);
	  qt_ali.PrintA3M(outf,hit);
	}	

      qt_ali.FreeMemory();
    }
  if (strcmp(outfile,"stdout")) fclose(outf);
}





/////////////////////////////////////////////////////////////////////////////////////
// Return the ROC_5 score for optimization (changed 28.3.08 by Michael & Johannes)
/////////////////////////////////////////////////////////////////////////////////////
void HitList::Optimize(HMM& q, char* buffer)
{
  const int NFAM =5;   // calculate ROC_5 score 
  const int NSFAM=5;   // calculate ROC_5 score 
  int roc=0;           // ROC score
  int fam=0;           // number of hits from same family (at current threshold)
  int not_fam=0;       // number of hits not from same family
  int sfam=0;          // number of hits from same suporfamily (at current threshold)
  int not_sfam=0;      // number of hits not from same superfamily
  Hit hit;
  
  SortList();
  Reset();
  while (!End()) 
    {
      hit = ReadNext();
      if (!strcmp(hit.fam,q.fam)) fam++;       // query and template from same superfamily? => positive
      else if (not_fam<NFAM) // query and template from different family? => negative
	{
	  not_fam++;
	  roc += fam;
	}
      if (!strcmp(hit.sfam,q.sfam)) sfam++;       // query and template from same superfamily? => positive
      else if (not_sfam<NSFAM) // query and template from different superfamily?   => negative
	{
	  not_sfam++;
	  roc += sfam;
	}
//       printf("qfam=%s tfam=%s qsfam=%s tsfam=%s  fam=%-2i  not_fam=%3i  sfam=%-3i  not_sfam=%-5i  roc=%-3i\n",q.fam,hit.fam,q.sfam,hit.sfam,fam,not_fam,sfam,not_sfam,roc);
    }

  // Write ROC score to file or stdout
  FILE* buf=NULL;
  if (strcmp(par.buffer,"stdout"))
    {
      buf=fopen(buffer,"w");
      if (!buf) OpenFileError(par.buffer);
    }
  else
    buf = stdout;

  fprintf(buf,"%f\n",float(roc)/float(fam*NFAM+sfam*NSFAM)); // must be between 0 and 1
  if (v>=2) printf("ROC=%f\n",float(roc)/float(fam*NFAM+sfam*NSFAM)); 
  fclose(buf);
} 



/////////////////////////////////////////////////////////////////////////////////////
// Print score distribution into file score_dist
/////////////////////////////////////////////////////////////////////////////////////
void HitList::PrintScoreFile(HMM& q)
{
  int i=0, n;
  FILE* scoref=NULL;
  Hit hit;
  Hash<int> twice(10000); // make sure only one hit per HMM is listed
  twice.Null(-1);      

  if (strcmp(par.scorefile,"stdout"))
    {
      scoref=fopen(par.scorefile,"w");
      if (!scoref)
	{cerr<<endl<<"WARNING from "<<par.argv[0]<<": could not open \'"<<par.scorefile<<"\'\n"; return;}
    }
  else
    scoref = stdout;
  Reset();
  fprintf(scoref,"NAME  %s\n",q.longname);
  fprintf(scoref,"FAM   %s\n",q.fam);
  fprintf(scoref,"FILE  %s\n",q.file);
  fprintf(scoref,"LENG  %i\n",q.L);
  fprintf(scoref,"\n");
//fprintf(scoref,"TARGET      REL LEN COL LOG-PVA   S-TOT     MS NALI\n");

//For hhformat, the PROBAB field has to start at position 41 !!
//                ----+----1----+----2----+----3----+----4----+---- 
  fprintf(scoref,"TARGET                FAMILY   REL  LEN  COL  LOG-PVA  S-AASS PROBAB  SCORE  LOG-EVAL\n");
  //              d153l__               5 185 185  287.82  464.22 100.00 
  //              d1qsaa2               3 168 124  145.55  239.22  57.36
  while (!End()) 
    {
      i++;
      hit = ReadNext();
      if (twice[hit.name]==1) continue; // better hit with same HMM has been listed already
      twice.Add(hit.name,1);
     //if template and query are from the same superfamily
      if (!strcmp(hit.name,q.name)) n=5;
      else if (!strcmp(hit.fam,q.fam)) n=4;
      else if (!strcmp(hit.sfam,q.sfam)) n=3;
      else if (!strcmp(hit.fold,q.fold)) n=2;
      else if (!strcmp(hit.cl,q.cl)) n=1;
      else n=0;
      fprintf(scoref,"%-20s %-10s %1i %5i %3i %8.3f %7.2f %6.2f %7.2f %8.3f\n",hit.name,hit.fam,n,hit.L,hit.matched_cols,-1.443*hit.logPval,-hit.score_aass,hit.Probab,hit.score,-1.443*hit.logEval);
    }      
  fclose(scoref);
}

void HitList::WriteToAlifile(HMM& q, bool scop_only)
{
  Hit hit;
  int i=0, n;
  Hash<int> twice(10000); // make sure only one hit per HMM is listed
  twice.Null(-1);      
  FILE* alitabf=NULL;
  if (strcmp(par.alitabfile,"stdout")) alitabf = fopen(par.alitabfile, "w"); else alitabf = stdout;
  if (!alitabf) OpenFileError(par.alitabfile);

  fprintf(alitabf,"NAME  %s\n",q.longname);
  fprintf(alitabf,"FAM   %s\n",q.fam);
  fprintf(alitabf,"FILE  %s\n",q.file);
  fprintf(alitabf,"LENG  %i\n",q.L);
  fprintf(alitabf,"\n");

  Reset();
  while (!End()) 
    {
      i++;
      hit = ReadNext();
      if (scop_only && (!strncmp(hit.name,"cl|",3) || !strncmp(hit.name,"UP20|",5) || !strncmp(hit.name,"NR20|",5))) continue;
      if (twice[hit.name]==1) continue; // better hit with same HMM has been listed already
      twice.Add(hit.name,1);
      //if template and query are from the same superfamily
      if (!strcmp(hit.name,q.name)) n=5;
      else if (!strcmp(hit.fam,q.fam)) n=4;
      else if (!strcmp(hit.sfam,q.sfam)) n=3;
      else if (!strcmp(hit.fold,q.fold)) n=2;
      else if (!strcmp(hit.cl,q.cl)) n=1;
      else n=0;

      if (hit.P_posterior != NULL) {
	fprintf(alitabf,"\nHit %3i (%-20s %-10s Rel: %i  LOG-PVA: %6.2f  LOG-EVAL: %6.2f  Score: %6.2f  Probab: %6.2f):\n    i     j  score     SS  probab\n",i,hit.name,hit.fam,n,-1.443*hit.logPval,-1.443*hit.logEval,hit.score,hit.Probab);
	for (int step=hit.nsteps; step>=1; step--)
	  if (hit.states[step]>=MM) 
	    fprintf(alitabf,"%5i %5i %6.2f %6.2f %7.4f\n",hit.i[step],hit.j[step],hit.S[step],hit.S_ss[step],hit.P_posterior[step]);
      } else { 
	fprintf(alitabf,"\nHit %3i (%-20s %-10s Rel: %i  LOG-PVA: %6.2f  LOG-EVAL: %6.2f  Score: %6.2f  Probab: %6.2f):\n    i     j  score     SS\n",i,hit.name,hit.fam,n,-1.443*hit.logPval,-1.443*hit.logEval,hit.score,hit.Probab);
	for (int step=hit.nsteps; step>=1; step--)
	  if (hit.states[step]>=MM) 
	    fprintf(alitabf,"%5i %5i %6.2f %6.2f\n",hit.i[step],hit.j[step],hit.S[step],hit.S_ss[step]);
      }
    }
  fclose(alitabf);
}


/////////////////////////////////////////////////////////////////////////////////////
//// Evaluate the *negative* log likelihood of the data at the vertex v = (lamda,mu)
////    p(s) = lamda * exp{ -exp[-lamda*(s-mu)] - lamda*(s-mu) } = lamda * exp( -exp(-x) - x) 
/////////////////////////////////////////////////////////////////////////////////////
double HitList::LogLikelihoodEVD(double* v)
{  
  double sum=0.0, sumw=0.0;
  for (int i=0; i<Nprof; i++)
    {
      double x = v[0]*(score[i]-v[1]);
      sum += weight[i]*(exp(-x)+x);
      sumw += weight[i];
    }
  return sum - sumw*log(v[0]);
}
// Static wrapper-function for calling the nonstatic member function LogLikelihoodEVD() 
// ( see http://www.newty.de/fpt/callback.html#member )
double HitList::LogLikelihoodEVD_static(void* pt2hitlist, double* v)
{
  HitList* mySelf = (HitList*) pt2hitlist; // explicitly cast to a pointer to Hitlist
  return mySelf->LogLikelihoodEVD(v);                // call member function
}

/////////////////////////////////////////////////////////////////////////////////////
//// Subroutine to FindMin: try new point given by highest point ihigh and fac and replace ihigh if it is lower 
/////////////////////////////////////////////////////////////////////////////////////
double HitList::TryPoint(const int ndim, double* p, double* y, double* psum, int ihigh, double fac, double (*Func)(void* pt2hitlist, double* v))
{
  // New point p_try = p_c + fac*(p_high-p_c),
  // where p_c = ( sum_i (p_i) - p_high)/ndim is the center of ndim other points
  // => p_try = fac1*sum_i(p_i) + fac2*p_high
  double fac1=(1.-fac)/ndim;
  double fac2=fac-fac1;
  double ptry[ndim];   //new point to try out
  double ytry;         //function value of new point 
  int j;               //index for the ndim parameters

  for (j=0; j<ndim; j++)
    ptry[j]=psum[j]*fac1+p[ihigh*ndim+j]*fac2;
  ytry = (*Func)(this,ptry);
  if (ytry<=y[ihigh]) 
    {
//       if (v>=4) printf("Trying:                  %-7.3f %-7.3f %-7.3f -> accept\n",ptry[0],ptry[1],ytry);
      y[ihigh]=ytry;
      for (j=0; j<ndim; j++)
	{
	  psum[j] += ptry[j]-p[ihigh*ndim+j]; //update psum[j]
	  p[ihigh*ndim+j]=ptry[j];            //replace p[ihigh] with ptry
	}	                              //Note: ihigh is now not highest point anymore!
    }
//   else if (v>=4) printf("Trying:                  %-7.3f %-7.3f %-7.3f -> reject\n",ptry[0],ptry[1],ytry);

  return ytry;
}



/////////////////////////////////////////////////////////////////////////////////////
////Find minimum with simplex method of Nelder and Mead (1965)
/////////////////////////////////////////////////////////////////////////////////////
float HitList::FindMin(const int ndim, double* p, double* y, double tol, int& nfunc, double (*Func)(void* pt2hitlist, double* v))
{
  const int MAXNFUNC=99; //maximum allowed number of function evaluations
  int ihigh;    //index of highest point on simplex
  int inext;    //index of second highest point on simplex
  int ilow;     //index of lowest point on simplex
  int i;        //index for the ndim+1 points
  int j;        //index for the ndim parameters
  double rtol;  //tolerance: difference of function value between highest and lowest point of simplex
  double temp;   //dummy
  double ytry;   //function value of trial point
  double psum[ndim]; //psum[j] = j'th coordinate of sum vector (sum over all vertex vectors)

  nfunc=0;    //number of function evaluations =0
  //Calculate sum vector psum[j]
  for (j=0; j<ndim; j++)
    {
      psum[j]=p[j];
      for (i=1; i<ndim+1; i++)
	psum[j]+=p[i*ndim+j];
    }
  
  // Repeat finding better points in simplex until rtol<tol
  while(1)
    {
      // Find indices for highest, next highest and lowest point
      ilow=0;
      if (y[0]>y[1]) {inext=1; ihigh=0;} else {inext=0; ihigh=1;}
      for (i=0; i<ndim+1; i++)
	{
	  if (y[i]<=y[ilow]) ilow=i;
	  if (y[i]>y[ihigh]) {inext=ihigh; ihigh=i;}
	  else if (y[i]>y[inext] && i!= ihigh) inext=i;
	}
      
      // If tolerance in y is smaller than tol swap lowest point to index 0 and break -> return
      rtol = 2.*fabs(y[ihigh]-y[ilow]) / (fabs(y[ihigh])+fabs(y[ilow])+1E-10);
      if (rtol<tol) 
	{
	  temp=y[ilow]; y[ilow]=y[0]; y[0]=temp;
	  for (j=0; j<ndim; j++)
	    {
	      temp=p[ilow*ndim+j]; p[ilow*ndim+j]=p[j]; p[j]=temp; 
	    }
	  break;
	}

      // Max number of function evaluations exceeded?
      if (nfunc>=MAXNFUNC ) 
	{
	  if (v) fprintf(stderr,"\nWARNING: maximum likelihood fit of score distribution did not converge.\n");
	  return 1;
	}

      nfunc+=2;
      // Point-reflect highest point on the center of gravity p_c of the other ndim points of the simplex
      if (v>=3) printf("%3i  %-7.3f  %-7.3f %-12.8f %-9.3E\n",nfunc,p[ilow*ndim],p[ilow*ndim+1],y[ilow],rtol);
//       if (v>=2) printf("           %3i %-9.3E   %-7.3f %-7.3f %-7.3f   %-7.3f %-7.3f %-7.3f   %-7.3f %-7.3f %-7.3f\n",nfunc,rtol,p[ilow*ndim],p[ilow*ndim+1],y[ilow],p[inext*ndim],p[inext*ndim+1],y[inext],p[ihigh*ndim],p[ihigh*ndim+1],y[ihigh]);
      ytry = TryPoint(ndim,p,y,psum,ihigh,-1.0,Func); //reflect highest point on p_c 

      if (ytry<=y[ilow])  
	{
	  ytry = TryPoint(ndim,p,y,psum,ihigh,2.0,Func); //expand: try new point 2x further away from p_c
//  	  if (v>=2) printf("Expanded:  %3i %-9.3E   %-7.3f %-7.3f %-7.3f   %-7.3f %-7.3f %-7.3f   %-7.3f %-7.3f %-7.3f\n",nfunc,rtol,p[ilow*ndim],p[ilow*ndim+1],y[ilow],p[inext*ndim],p[inext*ndim+1],y[inext],p[ihigh*ndim],p[ihigh*ndim+1],y[ihigh]);
	}
      else if (ytry>=y[inext]) 
	{
	  // The new point is worse than the second worst point
	  temp=y[ihigh];
	  ytry=TryPoint(ndim,p,y,psum,ihigh,0.5,Func); //contract simplex by 0.5 along (p_high-p_c
//  	  if (v>=2) printf("Compressed:%3i %-9.3E   %-7.3f %-7.3f %-7.3f   %-7.3f %-7.3f %-7.3f   %-7.3f %-7.3f %-7.3f\n",nfunc,rtol,p[ilow*ndim],p[ilow*ndim+1],y[ilow],p[inext*ndim],p[inext*ndim+1],y[inext],p[ihigh*ndim],p[ihigh*ndim+1],y[ihigh]);
	  if (ytry>=temp) 
	    {
	      // Trial point is larger than worst point => contract simplex by 0.5 towards lowest point
	      for (i=0; i<ndim+1; i++)
		{
		  if (i!=ilow)
		    {
		      for (j=0; j<ndim; j++)
			p[i*ndim+j]=0.5*(p[i*ndim+j]+p[ilow+j]);
		      y[i] = (*Func)(this,p+i*ndim);
// 		      y[i] = (*Func)(p+i*ndim);
		    }
		}
	      nfunc+=ndim;
//  	      if (v>=2) printf("Contracted:%3i %-9.3E   %-7.3f %-7.3f  %-7.3f  %-7.3f %-7.3f  %-7.3f  %-7.3f %-7.3f  %-7.3f\n",nfunc,rtol,p[ilow*ndim],p[ilow*ndim+1],y[ilow],p[inext*ndim],p[inext*ndim+1],y[inext],p[ihigh*ndim],p[ihigh*ndim+1],y[ihigh]);

	      //Calculate psum[j]
	      for (j=0; j<ndim; j++)
		{
		  psum[j]=p[j];
		  for (i=1; i<ndim+1; i++)
		    psum[j]+=p[i*ndim+j];
		}
	    }
	}
      else nfunc--;
    }
  return (float)rtol;
}



/////////////////////////////////////////////////////////////////////////////////////
//// Do a maximum likelihod fit of the scores with an EV distribution with parameters lamda and mu 
/////////////////////////////////////////////////////////////////////////////////////
void HitList::MaxLikelihoodEVD(HMM& q, int nbest)
{  
  double tol=1E-6;                 // Maximum relative tolerance when minimizing -log(P)/N (~likelihood)
  static char first_call=1;
  static Hash<int> size_fam(MAXPROF/10);  // Hash counts number of HMMs in family
  static Hash<int> size_sfam(MAXPROF/10); // Hash counts number of families in superfamily
  Hash<int> excluded(50);          // Hash containing names of superfamilies to be excluded from fit
  size_fam.Null(0);                // Set int value to return when no data can be retrieved
  size_sfam.Null(0);               // Set int value to return when no data can be retrieved
  excluded.Null(0);                // Set int value to return when no data can be retrieved
  Hit hit; 

  double mu;                       // EVD[mu,lam](x) = exp(-exp(-(x-mu)/lam)) = P(score<=x)
  double vertex[2*3];              // three vertices of the simplex in lamda-mu plane
  double yvertex[3];               // log likelihood values at the three vertices of the simplex
  int nfunc=0;                     // number of function calls
  double sum_weights=0.0;
  float sum_scores=0.0;
  float rtol;

  if (first_call==1) 
    {
      first_call=0;
      // Count how many HMMs are in each family
      if (v>=4) printf("  count number of profiles in each family and families in each superfamily ...\n");
      Reset();
      while (!End()) 
	{
	  hit = ReadNext();
	  if (!size_fam.Contains(hit.fam)) (*size_sfam(hit.sfam))++; //Add one to hash element for superfamily 
	 (*size_fam(hit.fam))++;              //Add one to hash element for family 
	  //      printf("size(%s)=%i name=%s\n",hit.fam,*size_fam(hit.fam),hit.name)
	}
      fams=size_fam.Size();
      sfams=size_sfam.Size();
      if (v>=3) 
	printf("%-3i HMMs from %i families and %i superfamilies searched. Found %i hits\n",N_searched,fams,sfams,Size());
    }

  // Query has SCOP family identifier?
  if (q.fam && q.fam[0]>='a' && q.fam[0]<='k' && q.fam[1]=='.')
    { 
      char sfamid[NAMELEN];
      char* ptr_in_fam=q.fam;
      while ((ptr_in_fam=strwrd(sfamid,ptr_in_fam,'-')))
	{
	  char* ptr=strrchr(sfamid,'.');
	  if (ptr) *ptr='\0';
	  excluded.Add(sfamid); 
// 	  fprintf(stderr,"Exclude SCOP superfamily %s  ptr_in_fam='%s'\n",sfamid,ptr_in_fam);
      }
    }
  // Exclude best superfamilies from fit
  else if (nbest>0) 
    {
      if (sfams<97+nbest) return;

      // Find the nbest best-scoring superfamilies for exclusion from first ML fit
      if (v>=4) printf("  find %i best-scoring superfamilies to exclude from first fit  ...\n",nbest);
      hit = Smallest();
      excluded.Add(hit.sfam);
//       printf("Exclude in first round: %s %8.2f %s\n",hit.name,hit.score_aass,hit.sfam);
      while (excluded.Size()<nbest)
	{
	  Reset();
	  while (!End() && excluded.Contains(ReadNext().sfam)) ;
	  hit=ReadCurrent();
	  while (!End()) 
	    {
	      if (ReadNext()<hit && !excluded.Contains(ReadCurrent().sfam)) 
		hit=ReadCurrent();
	    }
	  excluded.Add(hit.sfam);
// 	  printf("Exclude in first round: %s %8.2f %s %i %i\n",hit.name,hit.score_aass,hit.sfam,excluded.Size(),excluded.Contains(hit.sfam));
	}
      tol = 0.01/size_sfam.Size(); // tol=1/N would lead to delta(log-likelihood)~1 (where N ~ number of superfamilies) since (1+1/N)^N = e
    } 
  else 
    {
      // Find the best-scoring superfamilies from first fit for exclusion from second ML fit
      if (v>=4) printf("  find best-scoring superfamilies to exclude from second fit  ...\n");
      Reset();
      while (!End()) 
	{
	  hit = ReadNext();
	  if (hit.Eval < 0.05) excluded.Add(hit.sfam); // changed from 0.5 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	}
      tol = 0.001/size_sfam.Size(); // tol=1/N would lead to delta(log-likelihood)~1 (where N ~ number of superfamilies) since (1+1/N)^N = e
    }

  // Put scores into score[] and weights into weight[]
  if (v>=3) printf("  generate scores and weights array for ML fitting ...\n");
  Nprof=0;
  Reset();
  while (!End()) 
    {
      hit = ReadNext();
      if (hit.irep > 1) continue;                   //Use only best hit per template
      if (Nprof>=MAXPROF) break;

      char sfamid[NAMELEN];
      char* ptr_in_fam=hit.fam;
      while ((ptr_in_fam=strwrd(sfamid,ptr_in_fam,'-')))
      {
	char* ptr=strrchr(sfamid,'.');
	if (ptr) *ptr='\0';
	if (excluded.Contains(sfamid)) break;    //HMM is among superfamilies to be excluded 
      }
      if (excluded.Contains(sfamid)) {
	if (v>=3) fprintf(stderr,"Exclude hit %s (family %s contains %s)\n",hit.name,hit.fam,sfamid); 
	continue;
      }
//       ScopID(hit.cl,hit.fold,hit.sfam,hit.fam);     //Get scop superfamily code for template
//       if (*hit.sfam=='\0' || excluded.Contains(hit.sfam)) continue;    // skip HMM

      score[Nprof] = hit.score;
      weight[Nprof]=1./size_fam[hit.fam]/size_sfam[hit.sfam];
      sum_scores +=hit.score*weight[Nprof];
      sum_weights+=weight[Nprof];

      //DEBUG
//       if (v>=4) printf("%-10.10s   %-12.12s %-3i   %-12.12s %-3i   %6.4f   %6.4f  %7.1f\n",hit.name,hit.fam,size_fam[hit.fam],hit.sfam,size_sfam[hit.sfam],1./size_fam[hit.fam]/size_sfam[hit.sfam],sum,hit.score);
      Nprof++;
    }
  //DEBUG
  if (v>=3) 
    printf("%i hits used for score distribution\n",Nprof);
   // for (int i=0; i<Nprof; i++) printf("%3i  score=%8.3f  weight=%7.5f\n",i,score[i],weight[i]);

  // Set simplex vertices and function values
  mu = sum_scores/sum_weights - 0.584/LAMDA;
  if (par.loc) // fit only in local mode; in global mode use fixed value LAMDA and mu mean score
    {
      double (*Func)(void*, double*); 
      Func = HitList::LogLikelihoodEVD_static;

      if (nbest>0) {vertex[0]=LAMDA;   vertex[1]=mu;}  /////////////////////////////////////////// DEBUG
      else         {vertex[0]=q.lamda; vertex[1]=mu;}
      vertex[2]=vertex[0]+0.1; vertex[3]=vertex[1]; 
      vertex[4]=vertex[0];     vertex[5]=vertex[1]+0.2; 
      yvertex[0]=Func(this,vertex  );
      yvertex[1]=Func(this,vertex+2);
      yvertex[2]=Func(this,vertex+4);

      // Find lam and mu that minimize negative log likelihood of data
      if (v>=3) printf("Fitting to EVD by maximum likelihood...\niter lamda       mu    -log(P)/N   tol\n");
      rtol = FindMin(2,vertex,yvertex,tol,nfunc,Func);
      if (v>=3) printf("%3i  %-7.3f  %-7.2f     %-7.3f %-7.1E\n\n",nfunc,vertex[0],vertex[1],yvertex[0]-(1.5772-log(vertex[0])),rtol);
//       printf("HHsearch lamda=%-6.3f   mu=%-6.3f\n",vertex[0],vertex[1]);
    }
  else 
    {
      vertex[0]=LAMDA_GLOB; vertex[1]=mu; 
    }
  
  // Set lamda and mu of profile
  q.lamda = vertex[0]; 
  q.mu = vertex[1];

  // Set P-values and E-values
  // CHECK UPDATE FROM score=-logpval to score=-logpval+SSSCORE2NATLOG*score_ss !!!!
  Reset();
  while (!End()) 
    {
      hit = ReadNext();
      
      // Calculate total score in raw score units: P-value = 1- exp(-exp(-lamda*(Saa-mu))) 
      hit.logPval = logPvalue(hit.score,vertex);
      hit.Pval=Pvalue(hit.score,vertex);
      hit.Eval=exp(hit.logPval+log(N_searched));
      hit.logEval = hit.logPval+log(N_searched);
//    hit.score_aass = hit.logPval/0.45-3.0 - hit.score_ss;  // median(lamda)~0.45, median(mu)~4.0 in EVDs for scop20.1.63 HMMs
      hit.score_aass = -q.lamda*(hit.score-q.mu)/0.45-3.0 - fmin(hit.score_ss,fmax(0.0,0.5*hit.score-5.0)); // median(lamda)~0.45, median(mu)~3.0 in EVDs for scop20.1.63 HMMs
      hit.Probab = Probab(hit);
      hit.score_sort = hit.score_aass;
      Overwrite(hit);                     // copy hit object into current position of hitlist
    }
}

/////////////////////////////////////////////////////////////////////////////////////
// Calculate output of hidden neural network units
/////////////////////////////////////////////////////////////////////////////////////
inline float calc_hidden_output(const float* weights, const float* bias, float Lqnorm, float Ltnorm, float Nqnorm, float Ntnorm)
{
  float res;
  // Calculate activation of hidden unit = sum of all inputs * weights + bias
  res = Lqnorm*weights[0] + Ltnorm*weights[1] + Nqnorm*weights[2] + Ntnorm*weights[3] + *bias;
  res = 1.0 / (1.0 + exp(-(res ))); // logistic function
  return res;
}

////////////////////////////////////////////////////////////////////////////////////
//// Neural network regressions of lamda for EVD
/////////////////////////////////////////////////////////////////////////////////////
inline float lamda_NN(float Lqnorm, float Ltnorm, float Nqnorm, float Ntnorm)
{
  const int inputs = 4;
  const int hidden = 4;
  const float biases[] = {-0.73195, -1.43792, -1.18839, -3.01141}; // bias for all hidden units
  const float weights[] = { // Weights for the neural networks (column = start unit, row = end unit)
    -0.52356, -3.37650, 1.12984, -0.46796,
    -4.71361, 0.14166, 1.66807, 0.16383,
    -0.94895, -1.24358, -1.20293, 0.95434,
    -0.00318, 0.53022, -0.04914, -0.77046,
    2.45630, 3.02905, 2.53803, 2.64379
  };
  float lamda=0.0;
  for (int h = 0; h<hidden; h++) {
    lamda += calc_hidden_output( weights+inputs*h, biases+h, Lqnorm,Ltnorm,Nqnorm,Ntnorm ) * weights[hidden*inputs+h]; 
  }
  return lamda;
}

////////////////////////////////////////////////////////////////////////////////////
//// Neural network regressions of mu for EVD
/////////////////////////////////////////////////////////////////////////////////////
inline float mu_NN(float Lqnorm, float Ltnorm, float Nqnorm, float Ntnorm)
{
  const int inputs = 4;
  const int hidden = 6;
  const float biases[] = {-4.25264, -3.63484, -5.86653, -4.78472, -2.76356, -2.21580};  // bias for all hidden units
  const float weights[] = { // Weights for the neural networks (column = start unit, row = end unit)
    1.96172, 1.07181, -7.41256, 0.26471, 
    0.84643, 1.46777, -1.04800, -0.51425,
    1.42697, 1.99927, 0.64647, 0.27834,
    1.34216, 1.64064, 0.35538, -8.08311,
    2.30046, 1.31700, -0.46435, -0.46803,
    0.90090, -3.53067, 0.59212, 1.47503,
    -1.26036, 1.52812, 1.58413, -1.90409, 0.92803, -0.66871
  };
  float mu=0.0;
  for (int h = 0; h<hidden; h++) { 
    mu += calc_hidden_output( weights+inputs*h, biases+h, Lqnorm,Ltnorm,Nqnorm,Ntnorm ) * weights[hidden*inputs+h]; 
  }
  return 20.0*mu;
}

////////////////////////////////////////////////////////////////////////////////////
//// Neural network regressions of alpha for Evalue correction factor
/////////////////////////////////////////////////////////////////////////////////////
inline float alpha_NN(float Lqnorm, float Ltnorm, float Nqnorm, float Ntnorm)
{
  const int inputs = 4;
  const int hidden = 4;
  const float biases[] = {7.89636,3.68944,2.05448,3.69149};  // bias for all hidden units
  const float alpha_bias = 1.33439;
  const float weights[] = { // Weights for the neural networks (column = start unit, row = end unit)
    -6.72336, -4.73393, -2.15446, -4.75140,
    -14.54957, 4.05462, 0.57951, 3.55780,
    2.08289, -1.81976, -1.19936, -17.35097,
    1.53268, -8.13514, -2.50677, 1.51106,
    6.37397, -0.36254, 0.16279, -1.32174
  };
  float alpha=0.0;
  for (int h = 0; h<hidden; h++) { 
    alpha += calc_hidden_output( weights+inputs*h, biases+h, Lqnorm,Ltnorm,Nqnorm,Ntnorm ) * weights[hidden*inputs+h]; 
  }
  alpha = 1.0 / (1.0 + exp(-(alpha + alpha_bias))); // logistic function
  return alpha;
}

////////////////////////////////////////////////////////////////////////////////////
//// Neural network regressions of beta for Evalue correction factor
/////////////////////////////////////////////////////////////////////////////////////
inline float beta_NN(float Lqnorm, float Ltnorm, float Nqnorm, float Ntnorm)
{
  const int inputs = 4;
  const int hidden = 4;
  const float biases[] = {7.89636,3.68944,2.05448,3.69149};  // bias for all hidden units
  const float beta_bias = 5.43347;
  const float weights[] = { // Weights for the neural networks (column = start unit, row = end unit)
    -6.72336, -4.73393, -2.15446, -4.75140,
    -14.54957, 4.05462, 0.57951, 3.55780,
    2.08289, -1.81976, -1.19936, -17.35097,
    1.53268, -8.13514, -2.50677, 1.51106,
    -2.27841, -7.79426, -9.53092, 3.65717
  };
  float beta=0.0;
  for (int h = 0; h<hidden; h++) { 
    beta += calc_hidden_output( weights+inputs*h, biases+h, Lqnorm,Ltnorm,Nqnorm,Ntnorm ) * weights[hidden*inputs+h]; 
  }
  beta = 1.0 / (1.0 + exp(-(beta + beta_bias))); // logistic function
  return beta;
}

/////////////////////////////////////////////////////////////////////////////////////
//// Calculate HHblits composite E-values 
/////////////////////////////////////////////////////////////////////////////////////
void HitList::CalculateHHblitsEvalues(HMM& q)
{
  Hit hit; 
  // OLD!!!
  //float alpha=0.75, beta=0;  // correlation factors for HHblits Evalue (correlation factor = exp(alpha * S' - beta) )
  //const float log1000=log(1000.0);
  //int nhits = 0;

  float alpha = 0;
  float log_Pcut = log(par.prefilter_evalue_thresh / par.dbsize);
  float log_dbsize = log(par.dbsize);
  //printf("log_Pcut: %7.4f  Pcut: %7.4f DBsize: %10i   a: %4.2f  b: %4.2f  c: %4.2f\n",log_Pcut, exp(log_Pcut), par.dbsize, par.alphaa, par.alphab, par.alphac);
      
  Reset();
  while (!End()) 
    {
      hit = ReadNext();

      // if (nhits++<50) 
      // 	printf("before correction  Eval: %7.4g    logEval: %7.4f\n",hit.Eval, hit.logEval);

      alpha = par.alphaa + par.alphab * (hit.Neff_HMM - 1) * (1 - par.alphac * (q.Neff_HMM - 1));
      
      hit.Eval = exp(hit.logPval + log_dbsize + (alpha * log_Pcut)); 
      hit.logEval = hit.logPval + log_dbsize + (alpha * log_Pcut); 

      // if (nhits++<50) 
      // 	printf("                   Eval: %7.4g    logEval: %7.4f   alpha: %7.4f   Neff_T: %5.2f  Neff_Q: %5.2f\n",hit.Eval, hit.logEval, alpha, hit.Neff_HMM, q.Neff_HMM);

      Overwrite(hit);   // copy hit object into current position of hitlist
    }
  ResortList(); // use InsertSort to resort list according to sum of minus-log-Pvalues
}



/////////////////////////////////////////////////////////////////////////////////////
//// Calculate Pvalues as a function of query and template lengths and diversities
/////////////////////////////////////////////////////////////////////////////////////
void HitList::CalculatePvalues(HMM& q)
{  
  Hit hit; 
  float lamda=LAMDA_GLOB, mu=3.0;   // init for global search 
  const float log1000=log(1000.0);

  if(N_searched==0) N_searched=1;
  if (v>=3) 
    printf("Calculate Pvalues as a function of query and template lengths and diversities...\n");
  Reset();
  while (!End()) 
    {
      hit = ReadNext();

      if (par.loc)
	{
	  lamda = lamda_NN( log(q.L)/log1000, log(hit.L)/log1000, q.Neff_HMM/10.0, hit.Neff_HMM/10.0 ); 
	  mu    =    mu_NN( log(q.L)/log1000, log(hit.L)/log1000, q.Neff_HMM/10.0, hit.Neff_HMM/10.0 ); 
// 	  if (v>=3 && nhits++<20) 
// 	     printf("hit=%-10.10s Lq=%-4i  Lt=%-4i  Nq=%5.2f  Nt=%5.2f  =>  lamda=%-6.3f  mu=%-6.3f\n",hit.name,q.L,hit.L,q.Neff_HMM,hit.Neff_HMM,lamda,mu);
	}
      hit.logPval = logPvalue(hit.score,lamda,mu);
      hit.Pval    = Pvalue(hit.score,lamda,mu);
      hit.Eval=exp(hit.logPval+log(N_searched));
      hit.logEval = hit.logPval+log(N_searched);
//    hit.score_aass = hit.logPval/LAMDA-3.0 - hit.score_ss;  // median(lamda)~0.45, median(mu)~3.0 in EVDs for scop20.1.63 HMMs
      // P-value = 1- exp(-exp(-lamda*(Saa-mu))) => -lamda*(Saa-mu) = log(-log(1-Pvalue))
      hit.score_aass = (hit.logPval<-10.0? hit.logPval : log(-log(1-hit.Pval)) )/0.45 - fmin(lamda*hit.score_ss,fmax(0.0,0.2*(hit.score-8.0)))/0.45 - 3.0;
      hit.score_sort = hit.score_aass;
      hit.Probab = Probab(hit);
      Overwrite(hit);
    }
  SortList();
  Reset();
  return;
}

/////////////////////////////////////////////////////////////////////////////////////
//// Calculate Pvalues from calibration of  0: query HMM, 1:template HMMs, 2: both
/////////////////////////////////////////////////////////////////////////////////////
void HitList::GetPvalsFromCalibration(HMM& q)
{  
  Hit hit; 
  char warn=0;
  if(N_searched==0) N_searched=1;
  if (v>=2) 
    {
      switch (par.calm) 
	{
	case 0: 
	  printf("Using lamda=%-5.3f and mu=%-5.2f from calibrated query HMM %s. \n",q.lamda,q.mu,q.name);
	  printf("Note that HMMs need to be recalibrated when changing HMM-HMM alignment options.\n");
	  break;
	case 1:
	  printf("Using score distribution parameters lamda and mu from database HMMs \n");
	  break;
	case 2:
	  printf("Combining score distribution parameters lamda and mu from query and database HMMs\n");
	  printf("Note that HMMs need to be recalibrated when changing HMM-HMM alignment options.\n");
	  break;
	}
    }
  Reset();
  while (!End()) 
    {
      hit = ReadNext();
      if (par.calm==0 || (hit.logPvalt==0) )
	{
	  hit.logPval = logPvalue(hit.score,q.lamda,q.mu);
	  hit.Pval    = Pvalue(hit.score,q.lamda,q.mu);
	  if (par.calm>0 && warn++<1 && v>=1) 
	    printf("WARNING: some template HMM (e.g. %s) are not calibrated. Using query calibration.\n",hit.name);
	} 
      else if (par.calm==1) 
	{
	  hit.logPval = hit.logPvalt;
	  hit.Pval    = hit.Pvalt;
	}
      else if (par.calm==2) 
	{
	  hit.logPval = 0.5*( logPvalue(hit.score,q.lamda,q.mu) + hit.logPvalt);
	  hit.Pval    = sqrt( Pvalue(hit.score,q.lamda,q.mu) * hit.Pvalt);
	  if (v>=5) printf("Score: %7.1f  lamda: %7.1f  mu: %7.1f  P-values:  query-calibrated: %8.2G   template-calibrated: %8.2G   geometric mean: %8.2G\n",hit.score,q.lamda,q.mu,Pvalue(hit.score,q.lamda,q.mu),hit.Pvalt,hit.Pval);
	}

      hit.Eval=exp(hit.logPval+log(N_searched));
      hit.logEval = hit.logPval+log(N_searched);
//    hit.score_aass = hit.logPval/LAMDA-3.0 - hit.score_ss;  // median(lamda)~0.45, median(mu)~3.0 in EVDs for scop20.1.63 HMMs
      // P-value = 1- exp(-exp(-lamda*(Saa-mu))) => -lamda*(Saa-mu) = log(-log(1-Pvalue))
      hit.score_aass = (hit.logPval<-10.0? hit.logPval : log(-log(1-hit.Pval)) ) / 0.45-3.0 - fmin(hit.score_ss,fmax(0.0,0.5*hit.score-5.0));
      hit.score_sort = hit.score_aass;
      hit.Probab = Probab(hit);
      Overwrite(hit);
    }
  SortList();
  Reset();
  return;
}


