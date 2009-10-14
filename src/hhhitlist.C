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
  char str[NAMELEN]="";

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
  if (par.trans) 
    fprintf(outf," No Hit                             Prob  E-trans  E-value  Score    SS Cols Query HMM  Template HMM\n");
  else 
    fprintf(outf," No Hit                             Prob  E-value  P-value  Score    SS Cols Query HMM  Template HMM\n");
#else
  if (par.trans) 
    fprintf(outf," No Hit                             Prob E-trans E-value  Score    SS Cols Query HMM  Template HMM\n");
  else 
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
       sprintf(str,"%3i %-30.30s    ",nhits,hit.longname);
       
	 
#ifdef WINDOWS      
       if (par.trans) // Transitive scoring 
	 fprintf(outf,"%-34.34s %5.1f %8.2G %8.2G %6.1f %5.1f %4i ",str,hit.Probab,hit.E1val,hit.Eval,hit.score,hit.score_ss,hit.matched_cols);
       else // Normal scoring
	 fprintf(outf,"%-34.34s %5.1f %8.2G %8.2G %6.1f %5.1f %4i ",str,hit.Probab,hit.Eval,hit.Pval,hit.score,hit.score_ss,hit.matched_cols);
#else
       if (par.trans) // Transitive scoring 
         fprintf(outf,"%-34.34s %5.1f %7.2G %7.2G %6.1f %5.1f %4i ",str,hit.Probab,hit.E1val,hit.Eval,hit.score,hit.score_ss,hit.matched_cols);
       else // Normal scoring
	 fprintf(outf,"%-34.34s %5.1f %7.2G %7.2G %6.1f %5.1f %4i ",str,hit.Probab,hit.Eval,hit.Pval,hit.score,hit.score_ss,hit.matched_cols);
#endif

      sprintf(str,"%4i-%-4i ",hit.i1,hit.i2);
      fprintf(outf,"%-10.10s",str);
      sprintf(str,"%4i-%-4i",hit.j1,hit.j2);
      fprintf(outf,"%-9.9s(%i)\n",str,hit.L);
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
      // Count how many HMMs are in each family; set number of multiple hits per template nrep 
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
      hit.weight=1./size_fam[hit.fam]/size_sfam[hit.sfam];   // needed for transitive scoring
      hit.logPval = logPvalue(hit.score,vertex);
      hit.Pval=Pvalue(hit.score,vertex);
      hit.Eval=exp(hit.logPval+log(N_searched));
      hit.logEval = hit.logPval+log(N_searched);
//    hit.score_aass = hit.logPval/0.45-3.0 - hit.score_ss;  // median(lamda)~0.45, median(mu)~4.0 in EVDs for scop20.1.63 HMMs
      hit.score_aass = -q.lamda*(hit.score-q.mu)/0.45-3.0 - fmin(hit.score_ss,fmax(0.0,0.5*hit.score-5.0)); // median(lamda)~0.45, median(mu)~3.0 in EVDs for scop20.1.63 HMMs
      hit.Probab = Probab(hit);
      hit.score_sort = hit.score_aass;
      Overwrite(hit);                     // copy hit object into current position of hitlist

      if (nbest==0 && par.trans==1)       // if in transitive scoring mode (weights file given)
	TransitiveScoring();
      else if (nbest==0 && par.trans==2)  // if in transitive scoring mode (weights file given)
	TransitiveScoring2();
      else if (nbest==0 && par.trans==3)  // if in transitive scoring mode (weights file given)
	TransitiveScoring3();
      else if (nbest==0 && par.trans==4)  // if in transitive scoring mode (weights file given)
	TransitiveScoring4();
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
//// Calculate HHblast composite E-values 
/////////////////////////////////////////////////////////////////////////////////////
void HitList::CalculateHHblastEvalues(HMM& q)
{
  Hit hit; 
  float alpha=0.75, beta=0;  // correlation factors for HHblast Evalue (correlation factor = exp(alpha * S' - beta) )
  const float log1000=log(1000.0);
  int nhits = 0;

  Reset();
  while (!End()) 
    {
      hit = ReadNext();
      if (par.loc)
	{
	  alpha = alpha_NN( log(q.L)/log1000, log(hit.L)/log1000, q.Neff_HMM/10.0, hit.Neff_HMM/10.0 ); 
	  beta = beta_NN( log(q.L)/log1000, log(hit.L)/log1000, q.Neff_HMM/10.0, hit.Neff_HMM/10.0 ); 
 	  if (v>=3 && nhits++<20) 
	    printf("hit=%-10.10s Lq=%-4i  Lt=%-4i  Nq=%5.2f  Nt=%5.2f  =>  alpha=%-6.3f  beta=%-6.3f  S'=%-6.3f\n",hit.name,q.L,hit.L,q.Neff_HMM,hit.Neff_HMM,alpha,beta,par.hhblast_prefilter_logpval);
	}
      else 
	{
	  //printf("WARNING: global calibration not yet implemented!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
	}

      hit.Eval = exp(hit.logPval+log(N_searched)+(alpha*par.hhblast_prefilter_logpval - beta));     // overwrite E-value from HHsearch with composite E-value from HHblast
      hit.logEval = hit.logPval+log(N_searched)+(alpha*par.hhblast_prefilter_logpval - beta);
      //printf("      Eval: %7.4f    logEval: %7.4f   logPval: %7.4f\n",hit.Eval, hit.logEval, hit.logPval);
      Overwrite(hit);   // copy hit object into current position of hitlist
    }
  ResortList(); // use InsertSort to resort list according to sum of minus-log-Pvalues
}


/////////////////////////////////////////////////////////////////////////////////////
// Read file generated by blastpgp (default output) and store P-values in hash
/////////////////////////////////////////////////////////////////////////////////////
void HitList::ReadBlastFile(HMM& q)
{
  char line[LINELEN]="";    // input line
  int Ndb;       // number of sequences in database
  int Ldb=0;     // size of database in number of amino acids
  char* templ;
  int i;
  if (!blast_logPvals) { blast_logPvals = new(Hash<float>); blast_logPvals->New(16381,0); }
  
  FILE* blaf = NULL;
  if (!strcmp(par.blafile,"stdin")) blaf=stdin;
  else
    {
      blaf = fopen(par.blafile,"rb");
      if (!blaf) OpenFileError(par.blafile);
    }
  
  // Read number of sequences and size of database
  while (fgetline(line,LINELEN-1,blaf) && !strstr(line,"sequences;"));
  if (!strstr(line,"sequences;")) FormatError(par.blafile,"No 'Database:' string found.");
  char* ptr=line;
  Ndb = strint(ptr);
  if (Ndb==INT_MIN) FormatError(par.blafile,"No integer for number of sequences in database found.");
  while ((i=strint(ptr))>INT_MIN) Ldb = 1000*Ldb + i;
  if (Ldb==0) FormatError(par.blafile,"No integer for size of database found.");
  printf("\nNumber of sequences in database = %i    Size of database = %i\n",Ndb,Ldb);

  // Read all E-values and sequence lengths
  while (fgetline(line,LINELEN-1,blaf))
    {
      if (line[0]=='>')
	{
	  // Read template name
	  templ = new(char[255]);
	  ptr = line+1;
	  strwrd(templ,ptr); 
	  if (!blast_logPvals->Contains(templ)) // store logPval only for best HSP with template 
	    {
	      // Read length
	      while (fgetline(line,LINELEN-1,blaf) && !strstr(line,"Length ="));
	      ptr = line+18;
	      int length = strint(ptr);
	      // Read E-value
	      fgetline(line,LINELEN-1,blaf);
	      fgetline(line,LINELEN-1,blaf);
	      float EvalDB; // E-value[seq-db]  = Evalue for comparison Query vs. database, from PSI-BLAST
	      float EvalQT; // E-value[seq-seq] = Evalue for comparison Query vs. template (seq-seq)
	      double logPval;
	      ptr = strstr(line+20,"Expect =");
	      if (!ptr) FormatError(par.blafile,"No 'Expect =' string found.");
	      if (sscanf(ptr+8,"%g",&EvalDB)<1) 
		{
		  ptr[7]='1';
		  if (sscanf(ptr+7,"%g",&EvalDB)<1) 
		    FormatError(par.blafile,"No Evalue found after 'Expect ='.");
		}
	      // Calculate P-value[seq-seq] = 1 -  exp(-E-value[seq-seq]) = 1 - exp(-Lt/Ldb*E-value[seq-db])
	      EvalQT = length/double(Ldb)*double(EvalDB);
	      if (EvalQT>1E-3) logPval = log(1.0-exp(-EvalQT)); else logPval=log(double(EvalQT)+1.0E-99); 
	      blast_logPvals->Add(templ,logPval);
 	      printf("template=%-10.10s  length=%-3i  EvalDB=%8.2g  EvalQT=%8.2g  P-value=%8.2g log Pval=%8.2g\n",templ,length,EvalDB,EvalQT,exp(logPval),logPval);
	    } 
	  else
	    delete[] templ;
	}
    }
  fclose(blaf);
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
	    printf("Warning: some template HMM (e.g. %s) are not calibrated. Using query calibration.\n",hit.name);
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









/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Transitive scoring
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







/////////////////////////////////////////////////////////////////////////////////////
// Calculate P-values and Probabilities from transitive scoring over whole database
/////////////////////////////////////////////////////////////////////////////////////
void HitList::TransitiveScoring()
{
  void PrintMatrix(float** V, int N);
  void PrintMatrix(double** V, int N);

  float** Z;    // matrix of intra-db Z-scores Z_kl
  float** C;    // covariance matrix for Z_k: C_kl = sum_m=1^N (Z_km * Z_lm)
  char** fold;  // fold name of HMM k
  char** fam;   // family of HMM k
  float* Prob;  // probability of HMM k
  float* Zq;    // Zq[k] = Z-score between query and database HMM k
  float* Ztq;   // Ztq[k] = transitive Z-score from query to database HMM k: Ztq[k] = sum_l[ w_ql * Z_lk] / normalization_q
  float* Zrq;   // Zrq[k] = transitive Z-score from database HMM k to query: Zrq[k] = sum_l[ w_kl * Z_lq] / normalization_k
  float* w;     // unnormalized weight matrix; w[l] is w_ql or w_kl, respectively
  int* ll;      // ll[m] is the m'th index l for which Z_lq, Z_lk > Zmin_trans
  int N;        // dimension of weight matrix is NxN
  int M;        // number of HMMs l with Z_ql>Ztrans_min (or Z_lk>Ztrans_min, respectively)
  int k,l,m,n;  // indices for database HMMs 
  char name[NAMELEN];
  Hash<int> index(MAXPROF+7);  // index{name} = index of HMM name in {1,...,N} 
  index.Null(-1);              // Set int value to return when no data can be retrieved
  Hash<int> excluded(13);      // Hash containing names of superfamilies to be excluded from fit
  excluded.Null(0);            // Set int value to return when no data can be retrieved
  Hit hit; 
  size_t dummy;

  // Read weights matrix W with index hash and names array
  fprintf(stderr,"Reading in weights file\n");
  FILE* wfile = fopen(par.wfile,"rb");
  if (v>=1 && wfile==NULL) 
    {
      fprintf(stderr,"Error: %s could not be opened: (N_searched=%i) ",par.wfile,N_searched);
      perror("fopen");
      fprintf(stderr,"Skipping caclulation of transitive P-values\n"); 
      par.trans=0;
      return;
    }
  dummy=fread(&N,sizeof(int),1,wfile);  // read matrix dimension (i.e. number of HMMs in database)
  if (v>=1 && N!=N_searched) 
    {
      fprintf(stderr,"Error: Number %i of HMMs in weight file is different from number %i of HMMs in searched databases. \n",N,N_searched);
      fprintf(stderr,"Skipping caclulation of transitive P-values\n"); 
      par.trans=0;
      return;
    }
  if (v>=2) fprintf(stderr,"Calculating transitive P-values for %i HMMs\n",N);
  // Read names of HMMs (to specify mapping of HMM to matrix indices)
  for (k=0; k<N; k++) 
    {
      dummy=fread(name,sizeof(char),IDLEN,wfile);
      index.Add(name,k);
    }
  // Read symmetric Z-scores matrix
  Z = new(float*[N]);
  for (k=0; k<N; k++) 
    {
      Z[k] = new(float[N]);
      for (l=0; l<k; l++) Z[k][l] = Z[l][k];
      dummy=fread(Z[k]+k,sizeof(float),N-k,wfile);   
    }
  // Read symmetric covariance matrix
  C = new(float*[N]);
  for (k=0; k<N; k++) 
    {
      C[k] = new(float[N]);
      for (l=0; l<k; l++) C[k][l] = C[l][k];
      dummy=fread(C[k]+k,sizeof(float),N-k,wfile);
    }
  fclose(wfile);

  // Allocate memory
  Zq = new(float[N]);
  Ztq = new(float[N]);
  Zrq = new(float[N]);
  fold = new(char*[N]);
  fam = new(char*[N]);
  Prob = new(float[N]);
  ll = new(int[N]);
  w = new(float[N]);

  // Transform P-values to normally distributed Z-scores and store in Zq vector
  fprintf(stderr,"Transform P-values to Z-scores\n");
  float Zmax_neg   = Score2Z( -log(MINEVALEXCL) + log(N_searched) ); // calculate Z-score corresponding to E-value MINEVALEXCL
  float Zmin_trans = Score2Z( -log(par.Emax_trans) + log(N_searched) ); // calculate Z-score corresponding to E-value par.Emax_trans
  printf("Zmax = %6.2f   Zmin = %6.2f \n",Zmax_neg,Zmin_trans);

  Reset();
  while (!End()) 
    {
      hit = ReadNext();
      if (hit.irep>1) continue;
      k = index.Show(hit.name);
      if (k<0) {fprintf(stderr,"Error: no index found in weights file for domain %s\n",hit.name); exit(1);}      
      if (hit.logPvalt<0)
	Zq[k] = 0.5*Score2Z(fabs(hit.logPval)) + 0.5*Score2Z(fabs(hit.logPvalt));  // Zq[k] = 0.5*(Zkq + Zqk)
      else 
	Zq[k] = Score2Z(fabs(hit.logPval));                           // Zq[k] = Zqk 
//      printf("%4i  %-10.10s logPvalt=%9g  Zq=%9f\n",k,hit.name,hit.logPvalt,Zq[k]);
//       if (isnan(Zq[k])) {
// 	fprintf(stderr,"Error: a floating point exception occurred. Skipping transitive scoring\n"); 
// 	printf("%4i  %-10.10s logPval=%9g  logPvalt=%9g  Zq=%9f\n",k,hit.name,hit.logPval,hit.logPvalt,Zq[k]);
//   	par.trans=0;
// 	return;
//       }
      if (Zq[k]>Zmax_neg) excluded.Add(hit.fold);
      fold[k] = new(char[IDLEN]);
      fam[k] = new(char[IDLEN]);
      strcpy(fold[k],hit.fold);
      strcpy(fam[k],hit.fam);
      weight[k] = hit.weight;
      Prob[k] = hit.Probab;
   }
  
  if (v>=3) 
    {
      excluded.Reset();
      while (!excluded.End())
	{
	  excluded.ReadNext(name);
	  printf("Excluded fold %s from fitting to Ztq\n",name);
	}
    }


  ////////////////////////////////////////////////////////////////
  // Calculate transitive score (query->l) Zt[l]
  
  // Construct vector ll of indices l for which Z_lq > Zmin_trans
  m = 0;
  for (l=0; l<N; l++)
    if (Zq[l]>=Zmin_trans) ll[m++]=l;
  M = m;  // number of indices l for which Z_lq,Z_lk > Zmin_trans
  
//   for (m=0; m<M; m++)
//     fprintf(stderr,"m=%-4i l=%-4i  %-10.10s  Zq[l]=%7f\n",m,ll[m],fam[ll[m]],Zq[ll[m]]);

  if (M<=1) 
    for (k=0; k<N; k++) Ztq[k]=0.0;
  else
    {
      // Generate submatrix of C for indices l for which Z_lq,Z_lk > Zmin_trans 
      double** Csub = new(double*[M]);
      double** Cinv = new(double*[M]);
      for (m=0; m<M; m++) 
	{
	  Csub[m] = new(double[M]);
	  Cinv[m] = new(double[M]);
	  for (n=0; n<M; n++)
	    Csub[m][n] = double(C[ll[m]][ll[n]]);
	}
      
      if (v>=3) 
	{
	  fprintf(stderr,"Covariance matrix\n");
	  PrintMatrix(Csub,M);
	}
      
      // Invert Csub
      fprintf(stderr,"Calculate inverse of covariance submatrix\n");
      InvertMatrix(Cinv,Csub,M);
      
      if (v>=3) 
	{
	  fprintf(stderr,"Inverse covariance matrix\n");
	  PrintMatrix(Cinv,M);
	}
      
      // Calculate weights w[l] 
      for (m=0; m<M; m++) 
	{
	  double sum = 0.0;
	  for (n=0; n<M; n++)
	    sum += 1.0 * Cinv[m][n]; 
	  w[m] = fmax(sum,0.0);
	}
      for (l=0; l<M; l++) delete[](Cinv[l]);
      delete[](Cinv);
      
      // Calculate Ztq[k] for all HMMs k
      fprintf(stderr,"Calculate Ztq vector of transitive Z-scores\n");
      float norm = NormalizationFactor(Csub,w,M);
      for (k=0; k<N; k++) 
	{
	  double sumZ = 0.0;
	  for (m=0; m<M; m++) 
	    sumZ += w[m] * Z[ll[m]][k];
	  Ztq[k] = sumZ/norm;   
	}
      
      for (l=0; l<M; l++) delete[](Csub[l]);
      delete[](Csub);
    }

  ////////////////////////////////////////////////////////////////
  // Calculate reverse transitive score (l->query-) Zrq[l]

  fprintf(stderr,"Calculate Zrq vector of transitive Z-scores\n");
  for (k=0; k<N; k++) 
    {
      // Construct vector ll of indices l for which Z_lk > Zmin_tran
      m = 0;
      for (l=0; l<N; l++)
	  if (Z[l][k]+Z[k][l]>=2*Zmin_trans) ll[m++]=l;
      int M = m;  // number of indices l for which Z_lq,Z_lk > Zmin_tran


//    fprintf(stderr,"\nfam[k]: %s\n",fam[k]);
//    for (m=0; m<M; m++)
//    printf(stderr,"m=%-4i k=%-4i  l=%-4i  %-10.10s  Zq[l]=%7f  Z_lk=%7f  \n",m,k,ll[m],fold[ll[m]],Zq[ll[m]],Z[k][ll[m]]);
           
      if (M<=1) 
	{
	  Zrq[k] = Zq[k];
	}
      else 
	{
	  // Generate submatrix of C for indices l for which Z_lq,Z_lk > Zmin_trans 
	  double** Csub = new(double*[M]);
	  for (m=0; m<M; m++) 
	    {
	      Csub[m] = new(double[M]);
	      for (n=0; n<M; n++)
		Csub[m][n] = double(C[ll[m]][ll[n]]);
	    }
//        fprintf(stderr,"Covariance matrix\n");
//        PrintMatrix(Csub,M);
	      
	  if (M==2) 
	    {
	      for (m=0; m<M; m++) w[m] = 1.0/M;
	    }
	  else 
	    {
	      
	      double** Cinv = new(double*[M]);
	      for (m=0; m<M; m++) Cinv[m] = new(double[M]);

	      // Invert Csub
	      InvertMatrix(Cinv,Csub,M);
	      
	      //         fprintf(stderr,"Inverse covariance matrix\n");
	      //         PrintMatrix(Cinv,M);
	      
	      // Calculate weights w[l] 
	      for (m=0; m<M; m++) 
		{
		  double sum = 0.0;
		  for (n=0; n<M; n++)
		    sum += 1.0 * Cinv[m][n]; 
		  w[m] = fmax(sum,0.0);
		}

//            for (m=0; m<M; m++) fprintf(stderr,"w[%i]=%8.2g\n",m,w[m]);

	      for (l=0; l<M; l++) delete[](Cinv[l]);
	      delete[](Cinv);
	    }

	  // Calculate Zrq[k] and normalize
	  float norm = NormalizationFactor(Csub,w,M);
	  double sumZ = 0.0;
	  for (m=0; m<M; m++) 
	    sumZ += w[m] * Zq[ll[m]];
	  Zrq[k] = sumZ/norm;   
	  
	  for (l=0; l<M; l++) delete[](Csub[l]);
	  delete[](Csub);
	} 

//    fprintf(stderr,"\nZq[k]=%8.2g  Zq1[k]=%8.2g\n",Zq[k],Zrq[k]);
    }

  // Total Z-score = weighted sum over original Z-score, forward transitive and reverse transitive Z-score
  for (k=0; k<N; k++) 
    {
      float Zqtot =  Zq[k] + par.wtrans*(Ztq[k]+Zrq[k]);
//       if (isnan(Zqtot))
// 	{
// 	  fprintf(stderr,"Error: a floating point exception occurred. Skipping transitive scoring\n"); 
// 	  printf("%4i  %-10.10s Zq=%6.2f  Ztq=%6.2f  Zrq=%6.2f  Zqtot=%6.2f\n",k,fam[k],Zq[k],Ztq[k],Zrq[k],Zqtot);
// 	  par.trans=0;
// 	  return;
// 	}
      if (v>=2 &&  Zq[k] + Zqtot > 2*Zmin_trans) {
	printf("%4i  %-10.10s Zq=%6.2f  Ztq=%6.2f  Zrq=%6.2f  -> Zqtot=%6.2f\n",k,fam[k],Zq[k],Ztq[k],Zrq[k],Zqtot);
      }
      Ztq[k] = Zqtot;
    }

  // Calculate mean and standard deviation of Z1q
  fprintf(stderr,"Calculate mean and standard deviation of Ztq\n");
  double sumw=0.0;
  double sumZ=0.0;
  double sumZ2=0.0;
  for (k=0; k<N; k++) 
    {  
      if (excluded.Contains(fold[k])) continue;
      sumw  += weight[k];
      sumZ  += weight[k]*Ztq[k];
      sumZ2 += weight[k]*Ztq[k]*Ztq[k];
//       if (isnan(sumZ)) 
// 	{
// 	  fprintf(stderr,"Error: a floating point exception occurred. Skipping transitive scoring\n"); 
// 	  printf("%4i  %-10.10s Zq=%9f  Zrq=%9f  Ztq=%9f\n",k,fam[k],Zq[k],Zrq[k],Ztq[k]);
// 	  par.trans=0;
// 	  return;
// 	}
    }
  float mu = sumZ/sumw;  
  float sigma = sqrt(sumZ2/sumw-mu*mu);
  if (v>=2) printf("mu(Ztq)=%6.3f  sigma(Ztq)=%6.2f\n",mu,sigma);
  sigma *= 1.01;// correct different fitting of EVD and normal variables

  // Normalize Ztq and calculate P1-values
  fprintf(stderr,"Normalize Ztq and calculate P1-values\n");
  Reset();
  while (!End()) 
    {
      hit = ReadNext();
      hit.logPval = -Z2Score((Ztq[index.Show(hit.name)]-mu)/sigma);
      hit.E1val = N_searched*(hit.logPval<-100.0? 0.0 : exp(hit.logPval));
      // P-value = 1- exp(-exp(-lamda*(Saa-mu))) => -lamda*(Saa-mu) = log(-log(1-Pvalue))
      hit.score_aass = (hit.logPval<-10.0? hit.logPval : log(-log(1-exp(hit.logPval))) ) / 0.45-3.0 - hit.score_ss;
      hit.Probab = Probab(hit);
      hit.score_sort = hit.logPval;
      Overwrite(hit);                                        // copy hit object into current position of hitlist
   }

  for (k=0; k<N; k++) delete[](Z[k]);
  for (k=0; k<N; k++) delete[](C[k]);
  for (k=0; k<N; k++) delete[](fold[k]);
  for (k=0; k<N; k++) delete[](fam[k]);
  delete[](C);
  delete[](Z);
  delete[](fold);
  delete[](fam);
  delete[](Prob);
  delete[](ll);
  delete[](Zq);
  delete[](Ztq);
}



/////////////////////////////////////////////////////////////////////////////////////
// Calculate P-values and Probabilities from transitive scoring over whole database
/////////////////////////////////////////////////////////////////////////////////////
void HitList::TransitiveScoring2()
{
  void PrintMatrix(float** V, int N);
  void PrintMatrix(double** V, int N);

  float** Z;    // matrix of intra-db Z-scores Z_kl
  float** C;    // covariance matrix for Z_k: C_kl = sum_m=1^N (Z_km * Z_lm)
  char** fold;  // fold name of HMM k
  char** fam;   // family of HMM k
  float* Prob;  // probability of HMM k
  float* Zq;    // Zq[k] = Z-score between query and database HMM k
  float* Ztq;   // Ztq[k] = transitive Z-score from query to database HMM k: Ztq[k] = sum_l[ w_ql * Z_lk] / normalization_q
  float* Zrq;   // Zrq[k] = transitive Z-score from database HMM k to query: Zrq[k] = sum_l[ w_kl * Z_lq] / normalization_k
  float* w;     // unnormalized weight matrix; w[l] is w_ql or w_kl, respectively
  int* ll;      // ll[m] is the m'th index l for which Z_lq, Z_lk > Zmin_trans
  int N;        // dimension of weight matrix is NxN
  int M;        // number of HMMs l with Z_ql>Ztrans_min (or Z_lk>Ztrans_min, respectively)
  int k,l,m,n;  // indices for database HMMs 
  char name[NAMELEN];
  Hash<int> index(MAXPROF+7);  // index{name} = index of HMM name in {1,...,N} 
  index.Null(-1);              // Set int value to return when no data can be retrieved
  Hash<int> excluded(13);      // Hash containing names of superfamilies to be excluded from fit
  excluded.Null(0);            // Set int value to return when no data can be retrieved
  Hit hit; 
  size_t dummy;

  // Read weights matrix W with index hash and names array
  fprintf(stderr,"Reading in weights file\n");
  FILE* wfile = fopen(par.wfile,"rb");
  if (v>=1 && wfile==NULL) 
    {
      fprintf(stderr,"Error: %s could not be opened: (N_searched=%i) ",par.wfile,N_searched);
      perror("fopen");
      fprintf(stderr,"Skipping caclulation of transitive P-values\n"); 
      par.trans=0;
      return;
    }
  dummy=fread(&N,sizeof(int),1,wfile);  // read matrix dimension (i.e. number of HMMs in database)
  if (v>=1 && N!=N_searched) 
    {
      fprintf(stderr,"Error: Number %i of HMMs in weight file is different from number %i of HMMs in searched databases. \n",N,N_searched);
      fprintf(stderr,"Skipping caclulation of transitive P-values\n"); 
      par.trans=0;
      return;
    }
  if (v>=2) fprintf(stderr,"Calculating transitive P-values for %i HMMs\n",N);
  // Read names of HMMs (to specify mapping of HMM to matrix indices)
  for (k=0; k<N; k++) 
    {
      dummy=fread(name,sizeof(char),IDLEN,wfile);
      index.Add(name,k);
    }
  // Read symmetric Z-scores matrix
  Z = new(float*[N]);
  for (k=0; k<N; k++) 
    {
      Z[k] = new(float[N]);
      for (l=0; l<k; l++) Z[k][l] = Z[l][k];
      dummy=fread(Z[k]+k,sizeof(float),N-k,wfile);   
    }
  // Read symmetric covariance matrix
  C = new(float*[N]);
  for (k=0; k<N; k++) 
    {
      C[k] = new(float[N]);
      for (l=0; l<k; l++) C[k][l] = C[l][k];
      dummy=fread(C[k]+k,sizeof(float),N-k,wfile);
    }
  fclose(wfile);

  // Allocate memory
  Zq = new(float[N]);
  Ztq = new(float[N]);
  Zrq = new(float[N]);
  fold = new(char*[N]);
  fam = new(char*[N]);
  Prob = new(float[N]);
  ll = new(int[N]);
  w = new(float[N]);

  // Transform P-values to normally distributed Z-scores and store in Zq vector
  fprintf(stderr,"Transform P-values to Z-scores\n");
  float Zmax_neg   = Score2Z( -log(MINEVALEXCL) + log(N_searched) ); // calculate Z-score corresponding to E-value MINEVALEXCL
  float Zmin_trans = Score2Z( -log(par.Emax_trans) + log(N_searched) ); // calculate Z-score corresponding to E-value par.Emax_trans
  printf("Zmax = %6.2f   Zmin = %6.2f \n",Zmax_neg,Zmin_trans);

  Reset();
  while (!End()) 
    {
      hit = ReadNext();
      if (hit.irep>1) continue;
      k = index.Show(hit.name);
      if (k<0) {fprintf(stderr,"Error: no index found in weights file for domain %s\n",hit.name); exit(1);}      
      if (hit.logPvalt<0)
	Zq[k] = 0.5*Score2Z(fabs(hit.logPval)) + 0.5*Score2Z(fabs(hit.logPvalt));  // Zq[k] = 0.5*(Zkq + Zqk)
      else 
	Zq[k] = Score2Z(fabs(hit.logPval));                           // Zq[k] = Zqk 
//      printf("%4i  %-10.10s logPvalt=%9g  Zq=%9f\n",k,hit.name,hit.logPvalt,Zq[k]);
//      if (isnan(Zq[k])) 
//       {
// 	fprintf(stderr,"Error: a floating point exception occurred. Skipping transitive scoring\n"); 
// 	printf("%4i  %-10.10s logPval=%9g  logPvalt=%9g  Zq=%9f\n",k,hit.name,hit.logPval,hit.logPvalt,Zq[k]);
//   	par.trans=0;
// 	return;
//       }
      if (Zq[k]>Zmax_neg) excluded.Add(hit.fold);
      fold[k] = new(char[IDLEN]);
      fam[k] = new(char[IDLEN]);
      strcpy(fold[k],hit.fold);
      strcpy(fam[k],hit.fam);
      weight[k] = hit.weight;
      Prob[k] = hit.Probab;
   }
  
  if (v>=3) 
    {
      excluded.Reset();
      while (!excluded.End())
	{
	  excluded.ReadNext(name);
	  printf("Excluded fold %s from fitting to Ztq\n",name);
	}
    }


  ////////////////////////////////////////////////////////////////
  // Calculate transitive score (query->l) Zt[l]
  
  // Construct vector ll of indices l for which Z_lq > Zmin_trans
  m = 0;
  for (l=0; l<N; l++)
    if (Zq[l]>=Zmin_trans) ll[m++]=l;
  M = m;  // number of indices l for which Z_lq,Z_lk > Zmin_trans
  
//   for (m=0; m<M; m++)
//     fprintf(stderr,"m=%-4i l=%-4i  %-10.10s  Zq[l]=%7f\n",m,ll[m],fam[ll[m]],Zq[ll[m]]);

  if (M<=1) 
    for (k=0; k<N; k++) Ztq[k]=0.0;
  else
    {
      // Generate submatrix of C for indices l for which Z_lq,Z_lk > Zmin_trans 
      double** Csub = new(double*[M]);
      double** Cinv = new(double*[M]);
      for (m=0; m<M; m++) 
	{
	  Csub[m] = new(double[M]);
	  Cinv[m] = new(double[M]);
	  for (n=0; n<M; n++)
	    Csub[m][n] = double(C[ll[m]][ll[n]]);
	}
      
      if (v>=3) 
	{
	  fprintf(stderr,"Covariance matrix\n");
	  PrintMatrix(Csub,M);
	}
      
//       // Invert Csub
//       fprintf(stderr,"Calculate inverse of covariance submatrix\n");
//       InvertMatrix(Cinv,Csub,M);
      
//       if (v>=3) 
// 	{
// 	  fprintf(stderr,"Inverse covariance matrix\n");
// 	  PrintMatrix(Cinv,M);
// 	}
      

      // Calculate weights w[l] 
      for (m=0; m<M; m++) 
	{
	  double sum = 0.0;
	  for (n=0; n<M; n++)
	    sum += 1.0 * Csub[m][n]; 
	  printf("w[%4i] = %-8.5f\n",ll[m],1.0/sum);
	  w[m] = (sum>0? Zq[ll[m]] / sum : 0.0);
	}
      for (l=0; l<M; l++) delete[](Cinv[l]);
      delete[](Cinv);
      
      // Calculate Ztq[k] for all HMMs k
      fprintf(stderr,"Calculate Ztq vector of transitive Z-scores\n");
      float norm = NormalizationFactor(Csub,w,M);
      for (k=0; k<N; k++) 
	{
	  double sumZ = 0.0;
	  for (m=0; m<M; m++) 
	    sumZ += w[m] * Z[ll[m]][k];
	  Ztq[k] = sumZ/norm;   
	}
      
      for (l=0; l<M; l++) delete[](Csub[l]);
      delete[](Csub);
    }

  ////////////////////////////////////////////////////////////////
  // Calculate reverse transitive score (l->query-) Zrq[l]

  fprintf(stderr,"Calculate Zrq vector of transitive Z-scores\n");
  for (k=0; k<N; k++) 
    {
      // Construct vector ll of indices l for which Z_lk > Zmin_tran
      m = 0;
      for (l=0; l<N; l++)
	  if (Z[l][k]+Z[k][l]>=2*Zmin_trans) ll[m++]=l;  
      int M = m;  // number of indices l for which Z_lq,Z_lk > Zmin_tran


//    fprintf(stderr,"\nfam[k]: %s\n",fam[k]);
//    for (m=0; m<M; m++)
//    printf(stderr,"m=%-4i k=%-4i  l=%-4i  %-10.10s  Zq[l]=%7f  Z_lk=%7f  \n",m,k,ll[m],fold[ll[m]],Zq[ll[m]],Z[k][ll[m]]);
           
      if (M<=1) 
	{
	  Zrq[k] = Zq[k];
	}
      else 
	{
	  // Generate submatrix of C for indices l for which Z_lq,Z_lk > Zmin_trans 
	  double** Csub = new(double*[M]);
	  for (m=0; m<M; m++) 
	    {
	      Csub[m] = new(double[M]);
	      for (n=0; n<M; n++)
		Csub[m][n] = double(C[ll[m]][ll[n]]);
	    }
//        fprintf(stderr,"Covariance matrix\n");
//        PrintMatrix(Csub,M);
	      
	  if (M<=2) 
	    {
	      for (m=0; m<M; m++) w[m] = 1.0/M;
	    }
	  else 
	    {
	      
	      double** Cinv = new(double*[M]);
	      for (m=0; m<M; m++) Cinv[m] = new(double[M]);

// 	      // Invert Csub
// 	      InvertMatrix(Cinv,Csub,M);
	      
// //	      fprintf(stderr,"Inverse covariance matrix\n");
// //	      PrintMatrix(Cinv,M);
	      
	      // Calculate weights w[l] 
	      for (m=0; m<M; m++) 
		{
		  double sum = 0.0;
		  for (n=0; n<M; n++)
		    sum += 1.0 * Csub[m][n]; 
		  w[m] = (sum>0? Z[ll[m]][k] / sum : 0.0);
		}

//            for (m=0; m<M; m++) fprintf(stderr,"w[%i]=%8.2g\n",m,w[m]);

	      for (l=0; l<M; l++) delete[](Cinv[l]);
	      delete[](Cinv);
	    }

	  // Calculate Zrq[k] and normalize
	  float norm = NormalizationFactor(Csub,w,M);
	  double sumZ = 0.0;
	  for (m=0; m<M; m++) 
	    sumZ += w[m] * Zq[ll[m]];
	  Zrq[k] = sumZ/norm;   
	  
	  for (l=0; l<M; l++) delete[](Csub[l]);
	  delete[](Csub);
	} 

//    fprintf(stderr,"\nZq[k]=%8.2g  Zq1[k]=%8.2g\n",Zq[k],Zrq[k]);
    }

  // Total Z-score = weighted sum over original Z-score, forward transitive and reverse transitive Z-score
  for (k=0; k<N; k++) 
    {
      float Zqtot =  Zq[k] + par.wtrans*(Ztq[k]+Zrq[k]);
//        if (isnan(Zqtot))
// 	{
// 	  fprintf(stderr,"Error: a floating point exception occurred. Skipping transitive scoring\n"); 
// 	  printf("%4i  %-10.10s Zq=%6.2f  Ztq=%6.2f  Zrq=%6.2f  Zqtot=%6.2f\n",k,fam[k],Zq[k],Ztq[k],Zrq[k],Zqtot);
// 	  par.trans=0;
// 	  return;
// 	}
      if (v>=2 &&  Zq[k] + Zqtot > 2*Zmin_trans) {
	printf("%4i  %-10.10s Zq=%6.2f  Ztq=%6.2f  Zrq=%6.2f  -> Zqtot=%6.2f\n",k,fam[k],Zq[k],Ztq[k],Zrq[k],Zqtot);
      }
      Ztq[k] = Zqtot;
    }

  // Calculate mean and standard deviation of Z1q
  fprintf(stderr,"Calculate mean and standard deviation of Ztq\n");
  double sumw=0.0;
  double sumZ=0.0;
  double sumZ2=0.0;
  for (k=0; k<N; k++) 
    {  
      if (excluded.Contains(fold[k])) continue;
      sumw  += weight[k];
      sumZ  += weight[k]*Ztq[k];
      sumZ2 += weight[k]*Ztq[k]*Ztq[k];
//       if (isnan(sumZ)) 
// 	{
// 	  fprintf(stderr,"Error: a floating point exception occurred. Skipping transitive scoring\n"); 
// 	  printf("%4i  %-10.10s Zq=%9f  Zrq=%9f  Ztq=%9f\n",k,fam[k],Zq[k],Zrq[k],Ztq[k]);
// 	  par.trans=0;
// 	  return;
// 	}
    }
  float mu = sumZ/sumw;  
  float sigma = sqrt(sumZ2/sumw-mu*mu);
  if (v>=2) printf("mu(Ztq)=%6.3f  sigma(Ztq)=%6.2f\n",mu,sigma);
  sigma *= 1.01;// correct different fitting of EVD and normal variables

  // Normalize Ztq and calculate P1-values
  fprintf(stderr,"Normalize Ztq and calculate P1-values\n");
  Reset();
  while (!End()) 
    {
      hit = ReadNext();
      hit.logPval = -Z2Score((Ztq[index.Show(hit.name)]-mu)/sigma);
      hit.E1val = N_searched*(hit.logPval<-100? 0.0 : exp(hit.logPval));
      // P-value = 1- exp(-exp(-lamda*(Saa-mu))) => -lamda*(Saa-mu) = log(-log(1-Pvalue))
      hit.score_aass = (hit.logPval<-10.0? hit.logPval : log(-log(1-exp(hit.logPval))) ) / 0.45-3.0 - hit.score_ss;
      hit.Probab = Probab(hit);
      hit.score_sort = hit.logPval;
      Overwrite(hit);                                        // copy hit object into current position of hitlist
   }

  for (k=0; k<N; k++) delete[](Z[k]);
  for (k=0; k<N; k++) delete[](C[k]);
  for (k=0; k<N; k++) delete[](fold[k]);
  for (k=0; k<N; k++) delete[](fam[k]);
  delete[](C);
  delete[](Z);
  delete[](fold);
  delete[](fam);
  delete[](Prob);
  delete[](ll);
  delete[](Zq);
  delete[](Ztq);
}


/////////////////////////////////////////////////////////////////////////////////////
// Calculate P-values and Probabilities from transitive scoring over whole database
// Like TransitiveScoring(), 
// but in transitive scoring, Z1_qk = sum_l w_l*Z_lk,  use all  l:E_ql<=E_qk
// and in reverse    scoring, Z1_kr = sum_l w_l*Z_lq,  use all  l:E_kl<=E_kq
/////////////////////////////////////////////////////////////////////////////////////
void HitList::TransitiveScoring3()
{
  void PrintMatrix(float** V, int N);
  void PrintMatrix(double** V, int N);

  float** Z;    // matrix of intra-db Z-scores Z_kl
  float** C;    // covariance matrix for Z_k: C_kl = sum_m=1^N (Z_km * Z_lm)
  char** fold;  // fold name of HMM k
  char** fam;   // family of HMM k
  float* Prob;  // probability of HMM k
  float* Zq;    // Zq[k] = Z-score between query and database HMM k
  float* Ztq;   // Ztq[k] = transitive Z-score from query to database HMM k: Ztq[k] = sum_l[ w_ql * Z_lk] / normalization_q
  float* Zrq;   // Zrq[k] = transitive Z-score from database HMM k to query: Zrq[k] = sum_l[ w_kl * Z_lq] / normalization_k
  float* w;     // unnormalized weight matrix; w[l] is w_ql or w_kl, respectively
  int* ll;      // ll[m] is the m'th index l for which Z_lq, Z_lk > Zmin_trans
  int N;        // dimension of weight matrix is NxN
  int M;        // number of HMMs l with Z_ql>Ztrans_min (or Z_lk>Ztrans_min, respectively)
  int k,l,m,n;  // indices for database HMMs 
  char name[NAMELEN];
  Hash<int> index(MAXPROF+7);  // index{name} = index of HMM name in {1,...,N} 
  index.Null(-1);              // Set int value to return when no data can be retrieved
  Hash<int> excluded(13);      // Hash containing names of superfamilies to be excluded from fit
  excluded.Null(0);            // Set int value to return when no data can be retrieved
  Hit hit; 
  size_t dummy;

  // Read weights matrix W with index hash and names array
  fprintf(stderr,"Reading in weights file\n");
  FILE* wfile = fopen(par.wfile,"rb");
  if (v>=1 && wfile==NULL) 
    {
      fprintf(stderr,"Error: %s could not be opened: (N_searched=%i) ",par.wfile,N_searched);
      perror("fopen");
      fprintf(stderr,"Skipping caclulation of transitive P-values\n"); 
      par.trans=0;
      return;
    }
  dummy=fread(&N,sizeof(int),1,wfile);  // read matrix dimension (i.e. number of HMMs in database)
  if (v>=1 && N!=N_searched) 
    {
      fprintf(stderr,"Error: Number %i of HMMs in weight file is different from number %i of HMMs in searched databases. \n",N,N_searched);
      fprintf(stderr,"Skipping caclulation of transitive P-values\n"); 
      par.trans=0;
      return;
    }
  if (v>=2) fprintf(stderr,"Calculating transitive P-values for %i HMMs\n",N);
  // Read names of HMMs (to specify mapping of HMM to matrix indices)
  for (k=0; k<N; k++) 
    {
      dummy=fread(name,sizeof(char),IDLEN,wfile);
      index.Add(name,k);
    }
  // Read symmetric Z-scores matrix
  Z = new(float*[N]);
  for (k=0; k<N; k++) 
    {
      Z[k] = new(float[N]);
      for (l=0; l<k; l++) Z[k][l] = Z[l][k];
      dummy=fread(Z[k]+k,sizeof(float),N-k,wfile);   
    }
  // Read symmetric covariance matrix
  C = new(float*[N]);
  for (k=0; k<N; k++) 
    {
      C[k] = new(float[N]);
      for (l=0; l<k; l++) C[k][l] = C[l][k];
      dummy=fread(C[k]+k,sizeof(float),N-k,wfile);
    }
  fclose(wfile);

  // Allocate memory
  Zq = new(float[N]);
  Ztq = new(float[N]);
  Zrq = new(float[N]);
  fold = new(char*[N]);
  fam = new(char*[N]);
  Prob = new(float[N]);
  ll = new(int[N]);
  w = new(float[N]);

  // Transform P-values to normally distributed Z-scores and store in Zq vector
  fprintf(stderr,"Transform P-values to Z-scores\n");
  float Zmax_neg   = Score2Z( -log(MINEVALEXCL) + log(N_searched) ); // calculate Z-score corresponding to E-value MINEVALEXCL
  float Zmin_trans = Score2Z( -log(par.Emax_trans) + log(N_searched) ); // calculate Z-score corresponding to E-value par.Emax_trans
  printf("Zmax = %6.2f   Zmin = %6.2f \n",Zmax_neg,Zmin_trans);

  Reset();
  while (!End()) 
    {
      hit = ReadNext();
      if (hit.irep>1) continue;
      k = index.Show(hit.name);
      if (k<0) {fprintf(stderr,"Error: no index found in weights file for domain %s\n",hit.name); exit(1);}      
      if (hit.logPvalt<0)
	Zq[k] = 0.5*Score2Z(fabs(hit.logPval)) + 0.5*Score2Z(fabs(hit.logPvalt));  // Zq[k] = 0.5*(Zkq + Zqk)
      else 
	Zq[k] = Score2Z(fabs(hit.logPval));                           // Zq[k] = Zqk 
//      printf("%4i  %-10.10s logPvalt=%9g  Zq=%9f\n",k,hit.name,hit.logPvalt,Zq[k]);
//       if (isnan(Zq[k])) 
// 	{
// 	  fprintf(stderr,"Error: a floating point exception occurred. Skipping transitive scoring\n"); 
// 	  printf("%4i  %-10.10s logPval=%9g  logPvalt=%9g  Zq=%9f\n",k,hit.name,hit.logPval,hit.logPvalt,Zq[k]);
// 	  par.trans=0;
// 	  return;
// 	}
      if (Zq[k]>Zmax_neg) excluded.Add(hit.fold);
      fold[k] = new(char[IDLEN]);
      fam[k] = new(char[IDLEN]);
      strcpy(fold[k],hit.fold);
      strcpy(fam[k],hit.fam);
      weight[k] = hit.weight;
      Prob[k] = hit.Probab;
   }
  
  if (v>=3) 
    {
      excluded.Reset();
      while (!excluded.End())
	{
	  excluded.ReadNext(name);
	  printf("Excluded fold %s from fitting to Ztq\n",name);
	}
    }


  ////////////////////////////////////////////////////////////////
  // Calculate transitive score (query->l) Ztq[l]

  fprintf(stderr,"Calculate Ztq vector of transitive Z-scores\n");
  for (k=0; k<N; k++) 
    {
      // Construct vector ll of indices l for which Z_lq OR Z_lk >= max(Z_kq,Zmin_trans)
      float Zmink = fmax(Zq[k],Zmin_trans);
      for (m=l=0; l<N; l++)
	if (Zq[l]>=Zmink) ll[m++]=l;
      M = m;  // number of indices l for which  Z_lq OR Z_lk >= max(Z_kq,Zmin_trans)
  
//    for (m=0; m<M; m++)
//    fprintf(stderr,"m=%-4i l=%-4i  %-10.10s  Zq[l]=%7f\n",m,ll[m],fam[ll[m]],Zq[ll[m]]);

      if (M<=1) 
	{
	  Ztq[k]=Zq[k];
	}
      else
	{
	  // Generate submatrix of C for indices l for which Z_lq,Z_lk > Zmin_trans 
	  double** Csub = new(double*[M]);
	  double** Cinv = new(double*[M]);
	  for (m=0; m<M; m++) 
	    {
	      Csub[m] = new(double[M]);
	      Cinv[m] = new(double[M]);
	      for (n=0; n<M; n++)
		Csub[m][n] = double(C[ll[m]][ll[n]]);
	    }
	  
// 	  fprintf(stderr,"Covariance matrix\n");
// 	  PrintMatrix(Csub,M);
	  
	  // Invert Csub
// 	  fprintf(stderr,"Calculate inverse of covariance submatrix\n");
	  InvertMatrix(Cinv,Csub,M);  
	  
// 	  fprintf(stderr,"Inverse covariance matrix\n");
// 	  PrintMatrix(Cinv,M);
	  
	  // Calculate weights w[l] 
	  for (m=0; m<M; m++) 
	    {
	      double sum = 0.0;
	      for (n=0; n<M; n++)
		sum += 1.0 * Cinv[m][n]; // signal ~ sum_l w_l*Z_lq !
	      w[m] = fmax(sum,0.0);
	    }
	  for (l=0; l<M; l++) delete[](Cinv[l]);
	  delete[](Cinv);
	  
	  // Calculate Ztq[k]
	  float norm = NormalizationFactor(Csub,w,M);
	  double sumZ = 0.0;
	  for (m=0; m<M; m++) 
	    sumZ += w[m] * fmin(Zq[ll[m]],Z[ll[m]][k]);
// 	    sumZ += w[m] * Z[ll[m]][k];
	  Ztq[k] = sumZ/norm;   
      
	  for (l=0; l<M; l++) delete[](Csub[l]);
	  delete[](Csub);
	}
    }

  ////////////////////////////////////////////////////////////////
  // Calculate reverse transitive score (l->query-) Zrq[l]

  fprintf(stderr,"Calculate Zrq vector of transitive Z-scores\n");
  for (k=0; k<N; k++) 
    {
      // Construct vector ll of indices l for which Z_lk > Zmin_tran
      float Zmink = fmax(Zq[k],Zmin_trans);
      for (m=l=0; l<N; l++)
	if (Z[l][k]>=Zmink) ll[m++]=l;
      int M = m;  // number of indices l for which Z_lq,Z_lk > Zmin_tran


//    fprintf(stderr,"\nfam[k]: %s\n",fam[k]);
//    for (m=0; m<M; m++)
//    printf(stderr,"m=%-4i k=%-4i  l=%-4i  %-10.10s  Zq[l]=%7f  Z_lk=%7f  \n",m,k,ll[m],fold[ll[m]],Zq[ll[m]],Z[k][ll[m]]);
           
     if (M<=1) 
	{
	  Zrq[k] = Zq[k];
	}
      else 
	{
	  // Generate submatrix of C for indices l for which Z_lq,Z_lk > Zmin_trans 
	  double** Csub = new(double*[M]);
	  for (m=0; m<M; m++) 
	    {
	      Csub[m] = new(double[M]);
	      for (n=0; n<M; n++)
		Csub[m][n] = double(C[ll[m]][ll[n]]);
	    }
//        fprintf(stderr,"Covariance matrix\n");
//        PrintMatrix(Csub,M);
	      
	  if (M==2) 
	    {
	      for (m=0; m<M; m++) w[m] = 1.0/M;
	    }
	  else 
	    {
	      
	      double** Cinv = new(double*[M]);
	      for (m=0; m<M; m++) Cinv[m] = new(double[M]);

	      // Invert Csub
 	      InvertMatrix(Cinv,Csub,M); 
	      
	      //         fprintf(stderr,"Inverse covariance matrix\n");
	      //         PrintMatrix(Cinv,M);
	      
	      // Calculate weights w[l] 
	      for (m=0; m<M; m++) 
		{
		  double sum = 0.0;
		  for (n=0; n<M; n++)
		    sum += 1.0 * Cinv[m][n]; // signal ~ sum_l w_l*Z_lq !
		  w[m] = fmax(sum,0.0);
		}
//            for (m=0; m<M; m++) fprintf(stderr,"w[%i]=%8.2g\n",m,w[m]);
	      for (l=0; l<M; l++) delete[](Cinv[l]);
	      delete[](Cinv);
	    }

	  // Calculate Zrq[k] and normalize
	  float norm = NormalizationFactor(Csub,w,M);
	  double sumZ = 0.0;
	  for (m=0; m<M; m++) 
	    sumZ += w[m] * fmin(Zq[ll[m]],Z[ll[m]][k]);
// 	    sumZ += w[m] * Zq[ll[m]];
	  Zrq[k] = sumZ/norm;   
	  
	  for (l=0; l<M; l++) delete[](Csub[l]);
	  delete[](Csub);
	} 

//    fprintf(stderr,"\nZq[k]=%8.2g  Zq1[k]=%8.2g\n",Zq[k],Zrq[k]);
    }

  // Total Z-score = weighted sum over original Z-score, forward transitive and reverse transitive Z-score
  for (k=0; k<N; k++) 
    { 

      float Zqtot = Zq[k] + par.wtrans*(Ztq[k]+Zrq[k]);
//       if (isnan(Zqtot))
// 	{
// 	  fprintf(stderr,"Error: a floating point exception occurred. Skipping transitive scoring\n"); 
// 	  printf("%4i  %-10.10s Zq=%6.2f  Ztq=%6.2f  Zrq=%6.2f  ->  Zqtot=%6.2f\n",k,fam[k],Zq[k],Ztq[k],Zrq[k],Zqtot);
// 	  par.trans=0;
// 	  return;
// 	}
      if (v>=3 && Zqtot > 2*Zmin_trans) {
	printf("%4i  %-10.10s Zq=%6.2f  Ztq=%6.2f  Zrq=%6.2f  ->  Zqtot=%6.2f\n",k,fam[k],Zq[k],Ztq[k],Zrq[k],Zqtot);
      }
      Ztq[k] = Zqtot;
    }

  // Calculate mean and standard deviation of Z1q
  fprintf(stderr,"Calculate mean and standard deviation of Ztq\n");
  double sumw=0.0;
  double sumZ=0.0;
  double sumZ2=0.0;
  for (k=0; k<N; k++) 
    {  
      if (excluded.Contains(fold[k])) continue;
      sumw  += weight[k];
      sumZ  += weight[k]*Ztq[k];
      sumZ2 += weight[k]*Ztq[k]*Ztq[k];
//       if (isnan(sumZ)) 
// 	{
// 	  fprintf(stderr,"Error: a floating point exception occurred. Skipping transitive scoring\n"); 
// 	  printf("%4i  %-10.10s Zq=%9f  Zrq=%9f  Ztq=%9f\n",k,fam[k],Zq[k],Zrq[k],Ztq[k]);
// 	  par.trans=0;
// 	  return;
// 	}
    }
  float mu = sumZ/sumw;  
  float sigma = sqrt(sumZ2/sumw-mu*mu);
  if (v>=2) printf("mu(Ztq)=%6.3f  sigma(Ztq)=%6.2f\n",mu,sigma);
  sigma *= 1.01;// correct different fitting of EVD and normal variables

  // Normalize Ztq and calculate P1-values
  fprintf(stderr,"Normalize Ztq and calculate P1-values\n");
  Reset();
  while (!End()) 
    {
      hit = ReadNext();
      hit.logPval = -Z2Score((Ztq[index.Show(hit.name)]-mu)/sigma);
      hit.E1val = N_searched*(hit.logPval<-100? 0.0 : exp(hit.logPval));
      // P-value = 1- exp(-exp(-lamda*(Saa-mu))) => -lamda*(Saa-mu) = log(-log(1-Pvalue))
      hit.score_aass = (hit.logPval<-10.0? hit.logPval : log(-log(1-exp(hit.logPval))) ) / 0.45-3.0 - hit.score_ss;
      hit.Probab = Probab(hit);
      hit.score_sort = hit.logPval;
      Overwrite(hit);                                        // copy hit object into current position of hitlist
   }

  for (k=0; k<N; k++) delete[](Z[k]);
  for (k=0; k<N; k++) delete[](C[k]);
  for (k=0; k<N; k++) delete[](fold[k]);
  for (k=0; k<N; k++) delete[](fam[k]);
  delete[](C);
  delete[](Z);
  delete[](fold);
  delete[](fam);
  delete[](Prob);
  delete[](ll);
  delete[](Zq);
  delete[](Ztq);

}


/////////////////////////////////////////////////////////////////////////////////////
// Calculate P-values and Probabilities from transitive scoring over whole database
// Best tested scheme. Use fmin(Zq[ll[m]],Z[ll[m]][k]) 
// and fast approximation for weights (not inverse covariance matrix)
/////////////////////////////////////////////////////////////////////////////////////
void HitList::TransitiveScoring4()
{
  void PrintMatrix(float** V, int N);
  void PrintMatrix(double** V, int N);

  float** Z;    // matrix of intra-db Z-scores Z_kl
  float** C;    // covariance matrix for Z_k: C_kl = sum_m=1^N (Z_km * Z_lm)
  char** fold;  // fold name of HMM k
  char** fam;   // family of HMM k
  float* Prob;  // probability of HMM k
  float* Zq;    // Zq[k] = Z-score between query and database HMM k
  float* Ztq;   // Ztq[k] = transitive Z-score from query to database HMM k: Ztq[k] = sum_l[ w_ql * Z_lk] / normalization_q
  float* Zrq;   // Zrq[k] = transitive Z-score from database HMM k to query: Zrq[k] = sum_l[ w_kl * Z_lq] / normalization_k
  float* w;     // unnormalized weight matrix; w[l] is w_ql or w_kl, respectively
  int* ll;      // ll[m] is the m'th index l for which Z_lq, Z_lk > Zmin_trans
  int N;        // dimension of weight matrix is NxN
  int M;        // number of HMMs l with Z_ql>Ztrans_min (or Z_lk>Ztrans_min, respectively)
  int k,l,m,n;  // indices for database HMMs 
  char name[NAMELEN];
  Hash<int> index(MAXPROF+7);  // index{name} = index of HMM name in {1,...,N} 
  index.Null(-1);              // Set int value to return when no data can be retrieved
  Hash<int> excluded(13);      // Hash containing names of superfamilies to be excluded from fit
  excluded.Null(0);            // Set int value to return when no data can be retrieved
  Hit hit; 
  size_t dummy;

  // Read weights matrix W with index hash and names array
  fprintf(stderr,"Reading in weights file\n");
  FILE* wfile = fopen(par.wfile,"rb");
  if (v>=1 && wfile==NULL) 
    {
      fprintf(stderr,"Error: %s could not be opened: (N_searched=%i) ",par.wfile,N_searched);
      perror("fopen");
      fprintf(stderr,"Skipping caclulation of transitive P-values\n"); 
      par.trans=0;
      return;
    }
  dummy=fread(&N,sizeof(int),1,wfile);  // read matrix dimension (i.e. number of HMMs in database)
  if (v>=1 && N!=N_searched) 
    {
      fprintf(stderr,"Error: Number %i of HMMs in weight file is different from number %i of HMMs in searched databases. \n",N,N_searched);
      fprintf(stderr,"Skipping caclulation of transitive P-values\n"); 
      par.trans=0;
      return;
    }
  if (v>=2) fprintf(stderr,"Calculating transitive P-values for %i HMMs\n",N);
  // Read names of HMMs (to specify mapping of HMM to matrix indices)
  for (k=0; k<N; k++) 
    {
      dummy=fread(name,sizeof(char),IDLEN,wfile);
      index.Add(name,k);
    }
  // Read symmetric Z-scores matrix
  Z = new(float*[N]);
  for (k=0; k<N; k++) 
    {
      Z[k] = new(float[N]);
      for (l=0; l<k; l++) Z[k][l] = Z[l][k];
      dummy=fread(Z[k]+k,sizeof(float),N-k,wfile);   
    }
  // Read symmetric covariance matrix
  C = new(float*[N]);
  for (k=0; k<N; k++) 
    {
      C[k] = new(float[N]);
      for (l=0; l<k; l++) C[k][l] = C[l][k];
      dummy=fread(C[k]+k,sizeof(float),N-k,wfile);
    }
  fclose(wfile);

  // Allocate memory
  Zq = new(float[N]);
  Ztq = new(float[N]);
  Zrq = new(float[N]);
  fold = new(char*[N]);
  fam = new(char*[N]);
  Prob = new(float[N]);
  ll = new(int[N]);
  w = new(float[N]);

  // Transform P-values to normally distributed Z-scores and store in Zq vector
  fprintf(stderr,"Transform P-values to Z-scores\n");
  float Zmax_neg   = Score2Z( -log(MINEVALEXCL) + log(N_searched) ); // calculate Z-score corresponding to E-value MINEVALEXCL
  float Zmin_trans = Score2Z( -log(par.Emax_trans) + log(N_searched) ); // calculate Z-score corresponding to E-value par.Emax_trans
  printf("Zmax = %6.2f   Zmin = %6.2f \n",Zmax_neg,Zmin_trans);

  Reset();
  while (!End()) 
    {
      hit = ReadNext();
      if (hit.irep>1) continue;
      k = index.Show(hit.name);
      if (k<0) {fprintf(stderr,"Error: no index found in weights file for domain %s\n",hit.name); exit(1);}      
      if (hit.logPvalt<0)
	Zq[k] = 0.5*Score2Z(fabs(hit.logPval)) + 0.5*Score2Z(fabs(hit.logPvalt));  // Zq[k] = 0.5*(Zkq + Zqk)
      else 
	Zq[k] = Score2Z(fabs(hit.logPval));                           // Zq[k] = Zqk 
//      printf("%4i  %-10.10s logPvalt=%9g  Zq=%9f\n",k,hit.name,hit.logPvalt,Zq[k]);
//      if (isnan(Zq[k])) {
// 	fprintf(stderr,"Error: a floating point exception occurred. Skipping transitive scoring\n"); 
// 	printf("%4i  %-10.10s logPval=%9g  logPvalt=%9g  Zq=%9f\n",k,hit.name,hit.logPval,hit.logPvalt,Zq[k]);
//   	par.trans=0;
// 	return;
//       }
      if (Zq[k]>Zmax_neg) excluded.Add(hit.fold);
      fold[k] = new(char[IDLEN]);
      fam[k] = new(char[IDLEN]);
      strcpy(fold[k],hit.fold);
      strcpy(fam[k],hit.fam);
      weight[k] = hit.weight;
      Prob[k] = hit.Probab;
   }
  
  if (v>=3) 
    {
      excluded.Reset();
      while (!excluded.End())
	{
	  excluded.ReadNext(name);
	  printf("Excluded fold %s from fitting to Ztq\n",name);
	}
    }

  ////////////////////////////////////////////////////////////////
  // Calculate transitive score (query->l) Zt[l]
  
  // Construct vector ll of indices l for which Z_lq > Zmin_trans
  m = 0;
  for (l=0; l<N; l++)
    if (Zq[l]>=Zmin_trans) ll[m++]=l;
  M = m;  // number of indices l for which Z_lq,Z_lk > Zmin_trans
  
//   for (m=0; m<M; m++)
//     fprintf(stderr,"m=%-4i l=%-4i  %-10.10s  Zq[l]=%7f\n",m,ll[m],fam[ll[m]],Zq[ll[m]]);

  if (M<=1) 
    for (k=0; k<N; k++) Ztq[k]=0.0;
  else
    {
      // Generate submatrix of C for indices l for which Z_lq,Z_lk > Zmin_trans 
      double** Csub = new(double*[M]);
      for (m=0; m<M; m++) 
	{
	  Csub[m] = new(double[M]);
	  for (n=0; n<M; n++)
	    Csub[m][n] = double(C[ll[m]][ll[n]]);
	}
      
      if (v>=3) 
	{
	  fprintf(stderr,"Covariance matrix\n");
	  PrintMatrix(Csub,M);
	}
      

      // Calculate weights w[l] 
      for (m=0; m<M; m++) 
	{
	  double sum = 0.0;
	  for (n=0; n<M; n++)
	    sum += fmax(0.0,Csub[m][n]);
	  printf("w[%4i] = %-8.5f\n",ll[m],1.0/sum);
 	  w[m] = 1.0/sum;
	}
      
      // Calculate Ztq[k] for all HMMs k
      fprintf(stderr,"Calculate Ztq vector of transitive Z-scores\n");
      float norm = NormalizationFactor(Csub,w,M);
      for (k=0; k<N; k++) 
	{
	  double sumZ = 0.0;
	  for (m=0; m<M; m++) 
	    sumZ += w[m] * fmin(Zq[ll[m]],Z[ll[m]][k]);
	  Ztq[k] = sumZ/norm;   
	}
      
      for (l=0; l<M; l++) delete[](Csub[l]);
      delete[](Csub);
    }

  ////////////////////////////////////////////////////////////////
  // Calculate reverse transitive score (l->query-) Zrq[l]

  fprintf(stderr,"Calculate Zrq vector of transitive Z-scores\n");
  for (k=0; k<N; k++) 
    {
      // Construct vector ll of indices l for which Z_lk > Zmin_tran
      m = 0;
      for (l=0; l<N; l++)
	  if (Z[k][l]>=Zmin_trans) ll[m++]=l;  
      int M = m;  // number of indices l for which Z_lq,Z_lk > Zmin_tran


//    fprintf(stderr,"\nfam[k]: %s\n",fam[k]);
//    for (m=0; m<M; m++)
//    printf(stderr,"m=%-4i k=%-4i  l=%-4i  %-10.10s  Zq[l]=%7f  Z_lk=%7f  \n",m,k,ll[m],fold[ll[m]],Zq[ll[m]],Z[k][ll[m]]);
           
      if (M<=1) 
	{
	  Zrq[k] = Zq[k];
	}
      else 
	{
	  // Generate submatrix of C for indices l for which Z_lq,Z_lk > Zmin_trans 
	  double** Csub = new(double*[M]);
	  for (m=0; m<M; m++) 
	    {
	      Csub[m] = new(double[M]);
	      for (n=0; n<M; n++)
		Csub[m][n] = double(C[ll[m]][ll[n]]);
	    }
//        fprintf(stderr,"Covariance matrix\n");
//        PrintMatrix(Csub,M);
	      
	  // Calculate weights w[l] 
	  for (m=0; m<M; m++) 
	    {
	      double sum = 0.0;
	      for (n=0; n<M; n++)
		sum += fmax(0.0,Csub[m][n]);
	      w[m] = 1.0/sum; 
	    }
	  
//        for (m=0; m<M; m++) fprintf(stderr,"w[%i]=%8.2g\n",m,w[m]);


	  // Calculate Zrq[k] and normalize
	  float norm = NormalizationFactor(Csub,w,M);
	  double sumZ = 0.0;
	  for (m=0; m<M; m++) 
	    sumZ += w[m] * fmin(Zq[ll[m]],Z[ll[m]][k]);
	  Zrq[k] = sumZ/norm;   
	  
	  for (l=0; l<M; l++) delete[](Csub[l]);
	  delete[](Csub);
	} 

//    fprintf(stderr,"\nZq[k]=%8.2g  Zq1[k]=%8.2g\n",Zq[k],Zrq[k]);
    }

  // Total Z-score = weighted sum over original Z-score, forward transitive and reverse transitive Z-score
  for (k=0; k<N; k++) 
    {
      float Zqtot =  Zq[k] + par.wtrans*(Ztq[k]+Zrq[k]);
//        if (isnan(Zqtot))
// 	{
// 	  fprintf(stderr,"Error: a floating point exception occurred. Skipping transitive scoring\n"); 
// 	  printf("%4i  %-10.10s Zq=%6.2f  Ztq=%6.2f  Zrq=%6.2f  Zqtot=%6.2f\n",k,fam[k],Zq[k],Ztq[k],Zrq[k],Zqtot);
// 	  par.trans=0;
// 	  return;
// 	}
      if (v>=3 &&  Zq[k] + Zqtot > 2*Zmin_trans) {
	printf("%4i  %-10.10s Zq=%6.2f  Ztq=%6.2f  Zrq=%6.2f  -> Zqtot=%6.2f\n",k,fam[k],Zq[k],Ztq[k],Zrq[k],Zqtot);
      }
      Ztq[k] = Zqtot;
    }

  // Calculate mean and standard deviation of Z1q
  fprintf(stderr,"Calculate mean and standard deviation of Ztq\n");
  double sumw=0.0;
  double sumZ=0.0;
  double sumZ2=0.0;
  for (k=0; k<N; k++) 
    {  
      if (excluded.Contains(fold[k])) continue;
      sumw  += weight[k];
      sumZ  += weight[k]*Ztq[k];
      sumZ2 += weight[k]*Ztq[k]*Ztq[k];
//       if (isnan(sumZ)) 
// 	{
// 	  fprintf(stderr,"Error: a floating point exception occurred. Skipping transitive scoring\n"); 
// 	  printf("%4i  %-10.10s Zq=%9f  Zrq=%9f  Ztq=%9f\n",k,fam[k],Zq[k],Zrq[k],Ztq[k]);
// 	  par.trans=0;
// 	  return;
// 	}
    }
  float mu = sumZ/sumw;  
  float sigma = sqrt(sumZ2/sumw-mu*mu);
  if (v>=2) printf("mu(Ztq)=%6.3f  sigma(Ztq)=%6.2f\n",mu,sigma);
  sigma *= 1.01;// correct different fitting of EVD and normal variables

  // Normalize Ztq and calculate P1-values
  fprintf(stderr,"Normalize Ztq and calculate P1-values\n");
  Reset();
  while (!End()) 
    {
      hit = ReadNext();
      hit.logPval = -Z2Score((Ztq[index.Show(hit.name)]-mu)/sigma);
      hit.E1val = N_searched*(hit.logPval<-100? 0.0 : exp(hit.logPval));
      // P-value = 1- exp(-exp(-lamda*(Saa-mu))) => -lamda*(Saa-mu) = log(-log(1-Pvalue))
      hit.score_aass = (hit.logPval<-10.0? hit.logPval : log(-log(1-exp(hit.logPval))) ) / 0.45-3.0 - hit.score_ss;
      hit.Probab = Probab(hit);
      hit.score_sort = hit.logPval;
      Overwrite(hit);                                        // copy hit object into current position of hitlist
   }

  for (k=0; k<N; k++) delete[](Z[k]);
  for (k=0; k<N; k++) delete[](C[k]);
  for (k=0; k<N; k++) delete[](fold[k]);
  for (k=0; k<N; k++) delete[](fam[k]);
  delete[](C);
  delete[](Z);
  delete[](fold);
  delete[](fam);
  delete[](Prob);
  delete[](ll);
  delete[](Zq);
  delete[](Ztq);
}


/////////////////////////////////////////////////////////////////////////////////////
// Score2Z transforms the -log(P-value) score into a Z-score for 0 < S
// Score2Z(S) = sqrt(2)*dierfc(2*e^(-S)), where dierfc is the inverse of the complementary error function
/////////////////////////////////////////////////////////////////////////////////////
double HitList::Score2Z(double S)
{
  double s, t, u, w, x, y, z;
  if (S<=0) return double(-100000);
  y = ( S>200 ? 0.0 : 2.0*exp(-S) );
  if (y > 1) 
    {
      z =  (S<1e-6? 2*S : 2-y);
      w = 0.916461398268964 - log(z);
    }
  else 
    {
      z = y; 
      w = 0.916461398268964 - (0.69314718056-S);
    }

  u = sqrt(w);
  s = (log(u) + 0.488826640273108) / w;
  t = 1 / (u + 0.231729200323405);

  x = u * (1 - s * (s * 0.124610454613712 + 0.5)) - 
    ((((-0.0728846765585675 * t + 0.269999308670029) * t + 
       0.150689047360223) * t + 0.116065025341614) * t + 
     0.499999303439796) * t;
  t = 3.97886080735226 / (x + 3.97886080735226);
  u = t - 0.5;
  s = (((((((((0.00112648096188977922 * u + 
	1.05739299623423047e-4) * u - 0.00351287146129100025) * u - 
	7.71708358954120939e-4) * u + 0.00685649426074558612) * u + 
	0.00339721910367775861) * u - 0.011274916933250487) * u - 
	0.0118598117047771104)  * u + 0.0142961988697898018) * u + 
        0.0346494207789099922)  * u + 0.00220995927012179067;
  s = ((((((((((((s * u - 0.0743424357241784861) * u - 
	0.105872177941595488) * u + 0.0147297938331485121) * u + 
	0.316847638520135944) * u + 0.713657635868730364) * u + 
	1.05375024970847138)  * u + 1.21448730779995237) * u + 
	1.16374581931560831)  * u + 0.956464974744799006) * u + 
	0.686265948274097816) * u + 0.434397492331430115) * u + 
        0.244044510593190935) * t - 
    (z==0? 0: z * exp(x * x - 0.120782237635245222));
  x += s * (x * s + 1);
  if (y > 1) {
    x = -x;
  }
  return double (1.41421356237*x);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// Z2Score transforms the Z-score into a -log(P-value) value
// Z2Score(Z) = log(2) - log( erfc(Z/sqrt(2)) ) , where derfc is the complementary error function
/////////////////////////////////////////////////////////////////////////////////////////////////////////
double HitList::Z2Score(double Z)
{
    double t, u, x, y;
    x = 0.707106781188*Z;
    if (x>10) return 0.69314718056 - (-x*x - log( (1-0.5/x/x)/x/1.772453851) );
    t = 3.97886080735226 / (fabs(x) + 3.97886080735226);
    u = t - 0.5;
    y = (((((((((0.00127109764952614092 * u + 1.19314022838340944e-4) * u - 
        0.003963850973605135)   * u - 8.70779635317295828e-4) * u + 
        0.00773672528313526668) * u + 0.00383335126264887303) * u - 
        0.0127223813782122755)  * u - 0.0133823644533460069) * u + 
        0.0161315329733252248)  * u + 0.0390976845588484035) * u + 
        0.00249367200053503304;
    y = ((((((((((((y * u - 0.0838864557023001992) * u - 
        0.119463959964325415) * u + 0.0166207924969367356) * u + 
        0.357524274449531043) * u + 0.805276408752910567) * u + 
        1.18902982909273333)  * u + 1.37040217682338167) * u + 
        1.31314653831023098)  * u + 1.07925515155856677) * u + 
        0.774368199119538609) * u + 0.490165080585318424) * u + 
        0.275374741597376782) * t * (x>10? 0.0 : exp(-x * x));
    return 0.69314718056 - log( x < 0 ? 2 - y : y );
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
void PrintMatrix(float** V, int N)
{
  int k,l;
  for (k=0; k<N; k++)
    {
      fprintf(stderr,"k=%4i \n",k);
      for (l=0; l<N; l++)
	{
	  fprintf(stderr,"%4i:%6.3f ",l,V[k][l]);
	  if ((l+1)%10==0) fprintf(stderr,"\n");
	}
      fprintf(stderr,"\n");
    }
  fprintf(stderr,"\n");
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
void PrintMatrix(double** V, int N)
{
  int k,l;
  for (k=0; k<N; k++)
    {
      fprintf(stderr,"k=%4i \n",k);
      for (l=0; l<N; l++)
	{
	  fprintf(stderr,"%4i:%6.3f ",l,V[k][l]);
	  if ((l+1)%10==0) fprintf(stderr,"\n");
	}
      fprintf(stderr,"\n");
    }
  fprintf(stderr,"\n");
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
void HitList::Normalize(float* Ztq, char** fold, Hash<int>& excluded) 
{ 
  double sumw=0.0;
  double sumZ=0.0;
  double sumZ2=0.0;
  for (int k=0; k<N_searched; k++) 
    {  
      if (excluded.Contains(fold[k])) continue;
      sumw  += weight[k];
      sumZ  += weight[k]*Ztq[k];
      sumZ2 += weight[k]*Ztq[k]*Ztq[k];
    }
  float mu = sumZ/sumw;  
  float sigma = sqrt(sumZ2/sumw-mu*mu);
  printf("Transitive score Ztq: mu=%8.3g  sigma=%8.3g\n",mu,sigma);
  for (int k=0; k<N_searched; k++) Ztq[k] = (Ztq[k]-mu)/sigma;
  return;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// Calculate standard deviation of Z1 = sum_m [ w_m * Z_m ], where Csub_mn = cov(Z_m,Z_n)  
float HitList::NormalizationFactor(double** Csub, float* w, int M)
  {
    double sum=0.0;
    for (int m=0; m<M; m++) 
      {
	double summ=0.0;
	for (int n=0; n<M; n++) summ += Csub[m][n]*w[n];
	sum += w[m]*summ;
      }
    return sqrt(sum);  
  }

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// Calculate inverse of matrix A and store result in B
void HitList::InvertMatrix(double** B, double** A, int N)
{
  if (N==0) 
    {
      printf("Error: InvertMatrix called with matrix of dimension 0\n");
      exit(6);
    }
  if (N==1) 
    {
      B[0][0] = (A[0][0]==0.0? 0 :1.0/A[0][0]);
      return;
    }

  int k,l,m;
  double** V = new(double*[N]);
  double* s  = new(double[N]);
  for (k=0; k<N; k++) V[k] = new(double[N]);

  // Copy original matrix A into B since B will be overwritten by SVD()
  for (k=0; k<N; k++) 
    for (l=0; l<N; l++) 
      B[k][l] = A[k][l];

  SVD(B, N, s, V);  // U replaces B on output; s[] contains singluar values
  
  // Calculate inverse of A: A^-1 = V * diag(1/s) * U^t
  double** U = B;
  // Calculate V[k][m] -> V[k][m] *diag(1/s)
  for (k=0; k<N; k++) 
    for (m=0; m<N; m++) 
      if (s[m]!=0.0) V[k][m] /= s[m]; else V[k][m] = 0.0;
  // Calculate V[k][l] -> (V * U^t)_kl
  for (k=0; k<N; k++) 
    {
      if (v>=4 && k%100==0) printf("%i\n",k); 
      for (l=0; l<N; l++) 
	{
	  s[l] = 0.0; // use s[] as temporary memory to avoid overwriting B[k][] as long as it is needed
	  for (m=0; m<N; m++) 
	    s[l] += V[k][m]*U[l][m];
	}
      for (l=0; l<N; l++) V[k][l]=s[l];
    }  
  for (k=0; k<N; k++) 
    for (l=0; l<N; l++) 
      B[k][l] = V[k][l];

  for (k=0; k<N; k++) delete[](V[k]);
  delete[](V);
  return;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
void HitList::TransposeMatrix(double** V, int N)
{
  int k,l;
  for (k=0; k<N; k++) // transpose Z for efficiency of ensuing matrix multiplication
    for (l=0; l<k; l++) 
      {
	double buf = V[k][l];
	V[k][l] = V[l][k];
	V[l][k] = buf;
      }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
static double sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg) 
static double maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ? (maxarg1) : (maxarg2)) 
static int iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ? (iminarg1) : (iminarg2)) 
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a)) 

// This is a version of the Golub and Reinsch algorithm for singular value decomposition for a quadratic 
// (n x n) matrix A. It is sped up by transposing A amd V matrices at various places in the algorithm.
// On a 400x400 matrix it runs in 1.6 s or 2.3 times faster than the original (n x m) version.
// On a 4993x4993 matrix it runs in 2h03 or 4.5 times faster than the original (n x m) version.

// Given a matrix a[0..n-1][0..n-1], this routine computes its singular value decomposition, A = U · W · V^t . 
// The matrix U replaces a on output. The diagonal matrix of singular values W is out-put as a vector w[0..n-1]. 
// The matrix V (not the transpose V^t) is output as V[0..n-1][0..n-1] ./
void HitList::SVD(double **A, int n, double w[], double **V)
{
  int m=n; // in general algorithm A is an (m x n) matrix instead of (n x n)
 
  double pythag(double a, double b);
  int flag,i,its,j,jj,k,l=1,nm=1;
  double anorm,c,f,g,h,s,scale,x,y,z,*rv1;
  rv1=new(double[n]);
  g=scale=anorm=0.0;    
  
  // Householder reduction to bidiagonal form.
  if (v>=5) printf("\nHouseholder reduction to bidiagonal form\n");
  for (i=0;i<n;i++) {
    if (v>=4 && i%100==0) printf("i=%i\n",i);
    if (v>=4) fprintf(stderr,".");
    l=i+1;
    rv1[i]=scale*g;
    g=s=scale=0.0;
    if (i < m) {
      for (k=i;k<m;k++) scale += fabs(A[k][i]);
      if (scale) {
	for (k=i;k<m;k++) {
	  A[k][i] /= scale;
	  s += A[k][i]*A[k][i];
	}
	f=A[i][i];
	g = -SIGN(sqrt(s),f);
	h=f*g-s;
	A[i][i]=f-g;
	for (j=l;j<n;j++) {
	  for (s=0.0,k=i;k<m;k++) s += A[k][i]*A[k][j];
	  f=s/h;
	  for (k=i;k<m;k++) A[k][j] += f*A[k][i];
	}
	for (k=i;k<m;k++) A[k][i] *= scale;
      }
    }
    w[i]=scale *g;
    g=s=scale=0.0;
    if (i < m && i != n-1) {
      for (k=l;k<n;k++) scale += fabs(A[i][k]);
      if (scale) {
	for (k=l;k<n;k++) {
	  A[i][k] /= scale;
	  s += A[i][k]*A[i][k];
	}
	f=A[i][l];
	g = -SIGN(sqrt(s),f);
	h=f*g-s;
	A[i][l]=f-g;
	for (k=l;k<n;k++) rv1[k]=A[i][k]/h;
	for (j=l;j<m;j++) {
	  for (s=0.0,k=l;k<n;k++) s += A[j][k]*A[i][k];
	  for (k=l;k<n;k++) A[j][k] += s*rv1[k];
	}
	for (k=l;k<n;k++) A[i][k] *= scale;
      }
    }
    anorm=FMAX(anorm,(fabs(w[i])+fabs(rv1[i])));
  }
  // Accumulation of right-hand transformations.
  if (v>=5) printf("\nAccumulation of right-hand transformations\n");
  TransposeMatrix(V,n);
  for (i=n-1;i>=0;i--) {
    if (v>=4 && i%100==0) printf("i=%i\n",i);
    if (v>=4) fprintf(stderr,".");
    if (i < n-1) {
      if (g) {
	// Double division to avoid possible underflow.
	for (j=l;j<n;j++)
	  V[i][j]=(A[i][j]/A[i][l])/g;
	for (j=l;j<n;j++) {
	  for (s=0.0,k=l;k<n;k++) s += A[i][k]*V[j][k];
	  for (k=l;k<n;k++) V[j][k] += s*V[i][k];
	}
      }
      for (j=l;j<n;j++) V[j][i]=V[i][j]=0.0;
    }
    V[i][i]=1.0;
    g=rv1[i];
    l=i;
  }
  // Accumulation of left-hand transformations.
  if (v>=5) printf("\nAccumulation of left-hand transformations\n");
  TransposeMatrix(A,n);
  for (i=IMIN(m,n)-1;i>=0;i--) {
    if (v>=4 && i%100==0) printf("i=%i\n",i);
    if (v>=4) fprintf(stderr,".");
    l=i+1;
    g=w[i];
    for (j=l;j<n;j++) A[j][i]=0.0;
    if (g) {
      g=1.0/g;
      for (j=l;j<n;j++) {
	for (s=0.0,k=l;k<m;k++) s += A[i][k]*A[j][k];
	f=(s/A[i][i])*g;
	for (k=i;k<m;k++) A[j][k] += f*A[i][k];
      }
      for (j=i;j<m;j++) A[i][j] *= g;
    } else for (j=i;j<m;j++) A[i][j]=0.0;
    ++A[i][i];
  }

  // Diagonalization of the bidiagonal form: Loop over singular values, and over allowed iterations.
  if (v>=5) printf("\nDiagonalization of the bidiagonal form\n");
  for (k=n-1;k>=0;k--) {
    if (v>=4 && k%100==0) printf("k=%i\n",k);
    if (v>=4) fprintf(stderr,".");
    for (its=1;its<=30;its++) {
      flag=1;
      // Test for splitting. Note that rv1[1] is always zero.
      for (l=k;l>=0;l--) {
	nm=l-1;
	if ((double)(fabs(rv1[l])+anorm) == anorm) {
	  flag=0;
	  break;
	}
	if ((double)(fabs(w[nm])+anorm) == anorm) break;
      }
      if (flag) {
	// Cancellation of rv1[l], if l > 1.
	c=0.0;
	s=1.0;
	for (i=l;i<=k;i++) {
	  f=s*rv1[i];
	  rv1[i]=c*rv1[i];
	  if ((double)(fabs(f)+anorm) == anorm) break;
	  g=w[i];
	  h=pythag(f,g);
	  w[i]=h;
	  h=1.0/h;
	  c=g*h;
	  s = -f*h;
	  for (j=0;j<m;j++) {
	    y=A[nm][j];
	    z=A[i][j];
	    A[nm][j]=y*c+z*s;
	    A[i][j]=z*c-y*s;
	  }
	}
      }
      z=w[k];
      // Convergence.
      if (l == k) {
	// Singular value is made nonnegative.
	if (z < 0.0) {
	  w[k] = -z;
	  for (j=0;j<n;j++) V[k][j] = -V[k][j];
	}
	break;
      }
      if (its == 30) {printf("Error in SVD: no convergence in 30 iterations\n"); exit(7);}
      // Shift from bottom 2-by-2 minor.
      x=w[l];
      nm=k-1;
      y=w[nm];
      g=rv1[nm];
      h=rv1[k];
      f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
      g=pythag(f,1.0);
      f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
      // Next QR transformation:
      c=s=1.0;
      for (j=l;j<=nm;j++) {
	i=j+1;
	g=rv1[i];
	y=w[i];
	h=s*g;
	g=c*g;
	z=pythag(f,h);
	rv1[j]=z;
	c=f/z;
	s=h/z;
	f=x*c+g*s;
	g = g*c-x*s;
	h=y*s;
	y *= c;
	for (jj=0;jj<n;jj++) {
	  x=V[j][jj];
	  z=V[i][jj];
	  V[j][jj]=x*c+z*s;
	  V[i][jj]=z*c-x*s;
	}
	z=pythag(f,h);
	// Rotation can be arbitrary if z = 0.
	w[j]=z;
	if (z) {
	  z=1.0/z;
	  c=f*z;
	  s=h*z;
	}
	f=c*g+s*y;
	x=c*y-s*g;
	
	for (jj=0;jj<m;jj++) {
	  y=A[j][jj];
	  z=A[i][jj];
	  A[j][jj]=y*c+z*s;
	  A[i][jj]=z*c-y*s;
	}
      }
      rv1[l]=0.0;
      rv1[k]=f;
      w[k]=x;
    }
  }
  TransposeMatrix(V,n);
  TransposeMatrix(A,n);
  delete[](rv1);
}

// Computes (a2 + b2 )^1/2 without destructive underflow or overflow.
double pythag(double a, double b)
{
  double absa,absb;
  absa=fabs(a);
  absb=fabs(b);
  if (absa > absb) 
    return absa*sqrt(1.0+SQR(absb/absa));
  else 
    return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
}



