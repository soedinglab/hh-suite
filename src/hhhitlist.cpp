// hhhitlist.C

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
//// Methods of class HitList
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

#include "hhhitlist.h"

using std::ios;
using std::ifstream;
using std::ofstream;
using std::cout;
using std::cerr;
using std::endl;

/////////////////////////////////////////////////////////////////////////////////////
// Print summary listing of hits
/////////////////////////////////////////////////////////////////////////////////////
void HitList::PrintHitList(HMM* q, std::stringstream& out) {
  int nhits = 0;

  out << "Query         " << q->longname << std::endl;
  out << "Match_columns " << q->L << std::endl;
  out << "No_of_seqs    " << q->N_filtered << " out of " << q->N_in << std::endl;
  out << "Neff          " << q->Neff_HMM << std::endl;
  out << "Searched_HMMs " << N_searched << std::endl;

  // Print date stamp
  time_t* tp = new (time_t);
  *tp = time(NULL);

  out << "Date          " << ctime(tp);
  delete (tp);

  // Print command line
  out << "Command       ";
  for (int i = 0; i < par.argc; i++)
    if (strlen(par.argv[i]) <= par.maxdbstrlen)
      out << par.argv[i] << " ";
    else
      out << "<" << strlen(par.argv[i]) << "characters> ";
  out << std::endl << std::endl;

#ifdef WINDOWS
  out << " No Hit                             Prob  E-value  P-value  Score    SS Cols Query HMM  Template HMM" << std::endl;
#else
  out << " No Hit                             Prob E-value P-value  Score    SS Cols Query HMM  Template HMM" << std::endl;
#endif

  char line[LINELEN];
  Reset();
  while (!End()) {
    Hit hit = ReadNext();
    if (nhits >= par.Z)
      break;       //max number of lines reached?
    if (nhits >= par.z && hit.Probab < par.p)
      break;
    if (nhits >= par.z && hit.Eval > par.E)
      continue;
    nhits++;

    char Estr[10];
    char Pstr[10];
    char str[NAMELEN];
    sprintf(str, "%3i %-30.30s    ", nhits, hit.longname);

#ifdef WINDOWS
    if (hit.Eval>=1E-99) sprintf(Estr,"%8.2G",hit.Eval);
    else sprintf(Estr,"%8.1E",hit.Eval);
    if (hit.Pval>=1E-99) sprintf(Pstr,"%8.2G",hit.Pval);
    else sprintf(Pstr,"%8.1E",hit.Pval);
    sprintf(line,"%-34.34s %5.1f %8s %8s ",str,hit.Probab,Estr,Pstr);
    out << line;
#else
    if (hit.Eval >= 1E-99)
      sprintf(Estr, "%7.2G", hit.Eval);
    else
      sprintf(Estr, "%7.0E", hit.Eval);
    if (hit.Pval >= 1E-99)
      sprintf(Pstr, "%7.2G", hit.Pval);
    else
      sprintf(Pstr, "%7.0E", hit.Pval);
    sprintf(line, "%-34.34s %5.1f %7s %7s ", str, hit.Probab, Estr, Pstr);
    out << line;
#endif

    // Needed for long sequences (more than 5 digits in length)
    sprintf(str, "%6.1f", hit.score);
    sprintf(line, "%-6.6s %5.1f %4i %4i-%-4i %4i-%-4i(%i)\n", str, hit.score_ss,
        hit.matched_cols, hit.i1, hit.i2, hit.j1, hit.j2, hit.L);
    out << line;
  }

  out << std::endl;
}

void HitList::PrintHitList(HMM* q, char* outfile) {
  std::stringstream out;
  PrintHitList(q, out);

  if (strcmp(outfile, "stdout") == 0) {
    std::cout << out.str();
  }
  else {
    std::ofstream outf(outfile);
    if (!outf)
      OpenFileError(outfile);

    outf << out.str();

    outf.close();
  }
}

void HitList::PrintHHR(HMM* q, char* outfile) {
  std::stringstream out;
  PrintHHR(q, out);

  if (strcmp(outfile, "stdout") == 0) {
    std::cout << out.str();
  }
  else {
    std::ofstream outf(outfile);

    if (!outf)
      OpenFileError(outfile);

    outf << out.str();
    outf.close();
  }
}

void HitList::PrintHHR(HMM* q, std::stringstream& out) {
  PrintHitList(q, out);
  PrintAlignments(q, out, 0);
}

/////////////////////////////////////////////////////////////////////////////////////
// Print alignments of query sequences against hit sequences 
/////////////////////////////////////////////////////////////////////////////////////
void HitList::PrintAlignments(HMM* q, char* outfile, char outformat) {
  std::stringstream out;
  PrintAlignments(q, out, outformat);

  if (strcmp(outfile, "stdout") == 0) {
    std::cout << out.str();
  }
  else {
    std::fstream outf;
    if (outformat == 0) {
      outf.open(outfile, std::ios::out | std::ios::app);
    }
    else
      outf.open(outfile, std::ios::out);

    if (!outf.good())
      OpenFileError(outfile);

    outf << out.str();
    outf.close();
  }
}

void HitList::PrintAlignments(HMM* q, std::stringstream& out, char outformat) {
  FullAlignment qt_ali(par.nseqdis + 10); // maximum 10 annotation (pseudo) sequences (ss_dssp, sa_dssp, ss_pred, ss_conf, consens,...)
  int nhits = 0;

  Reset();
  while (!End()) {
    if (nhits >= par.B)
      break;
    Hit hit = ReadNext();
    if (nhits >= par.b && hit.Probab < par.p)
      break;
    if (nhits >= par.b && hit.Eval > par.E)
      continue;
    nhits++;

    // Build double alignment of query against template sequences
    qt_ali.Build(q, hit);

    // Print out alignment
    // HHR format
    out << "No " << nhits << std::endl;
    if (outformat == 0)
    {
      qt_ali.PrintHeader(out, q, hit);
      qt_ali.PrintHHR(out, hit);
    }
    // FASTA format
    else if (outformat == 1)
    {
      qt_ali.PrintFASTA(out, hit);
    }
    // A2M format
    else if (outformat == 2)
    {
      qt_ali.PrintA2M(out, hit);
    }
    // A3m format
    else
    {
      qt_ali.PrintA3M(out, hit);
    }

    qt_ali.FreeMemory();
  }
}


/////////////////////////////////////////////////////////////////////////////////////
// Return the ROC_5 score for optimization (changed 28.3.08 by Michael & Johannes)
/////////////////////////////////////////////////////////////////////////////////////
void HitList::Optimize(HMM* q) {
  const int NFAM = 5;   // calculate ROC_5 score
  const int NSFAM = 5;   // calculate ROC_5 score
  int roc = 0;           // ROC score
  int fam = 0;         // number of hits from same family (at current threshold)
  int not_fam = 0;       // number of hits not from same family
  int sfam = 0;   // number of hits from same suporfamily (at current threshold)
  int not_sfam = 0;      // number of hits not from same superfamily
  Hit hit;
  
  SortList();
  Reset();
  while (!End()) {
    hit = ReadNext();
    if (!strcmp(hit.fam, q->fam))
      fam++;       // query and template from same superfamily? => positive
    else if (not_fam < NFAM) // query and template from different family? => negative
        {
      not_fam++;
      roc += fam;
    }
    if (!strcmp(hit.sfam, q->sfam))
      sfam++;       // query and template from same superfamily? => positive
    else if (not_sfam < NSFAM) // query and template from different superfamily?   => negative
        {
      not_sfam++;
      roc += sfam;
    }
//       printf("qfam=%s tfam=%s qsfam=%s tsfam=%s  fam=%-2i  not_fam=%3i  sfam=%-3i  not_sfam=%-5i  roc=%-3i\n",q->fam,hit.fam,q->sfam,hit.sfam,fam,not_fam,sfam,not_sfam,roc);
  }

  // Write ROC score to file or stdout
  printf("%f\n", float(roc) / float(fam * NFAM + sfam * NSFAM)); // must be between 0 and 1
  if (v >= 2)
    printf("ROC=%f\n", float(roc) / float(fam * NFAM + sfam * NSFAM));
}

/////////////////////////////////////////////////////////////////////////////////////
// Print score distribution into file score_dist
/////////////////////////////////////////////////////////////////////////////////////
void HitList::PrintScoreFile(HMM* q, char* outputfile) {
  std::stringstream outbuffer;
  PrintScoreFile(q, outbuffer);

  if (strcmp(outputfile, "stdout") == 0) {
    std::cout << outbuffer.str();
  }
  else {
    std::ofstream out(outputfile);
    if (!out.good()) {
      std::cerr << "WARNING from " << par.argv[0] << ": could not open \'" << outputfile << std::endl;
      return;
    }

    out << outbuffer.str();

    out.close();
  }
}

void HitList::PrintScoreFile(HMM* q, std::stringstream& outbuffer) {
  Hash<int> twice(10000); // make sure only one hit per HMM is listed
  twice.Null(-1);

  Reset();
  outbuffer << "NAME  " << q->longname << std::endl;
  outbuffer << "FAM   " << q->fam << std::endl;
  outbuffer << "FILE  " << q->file << std::endl;
  outbuffer << "LENG  " << q->L << std::endl;

//For hhformat, the PROBAB field has to start at position 41 !!
//                ----+----1----+----2----+----3----+----4----+----
  outbuffer << "TARGET                FAMILY   REL  LEN  COL  LOG-PVA  S-AASS PROBAB  SCORE  LOG-EVAL\n";
  //              d153l__               5 185 185  287.82  464.22 100.00 
  //              d1qsaa2               3 168 124  145.55  239.22  57.36

  int i = 0;
  while (!End()) {
    i++;
    Hit hit = ReadNext();
    if (twice[hit.name] == 1)
      continue; // better hit with same HMM has been listed already
    twice.Add(hit.name, 1);

    int n;
    //if template and query are from the same superfamily
    if (!strcmp(hit.name, q->name))
      n = 5;
    else if (!strcmp(hit.fam, q->fam))
      n = 4;
    else if (!strcmp(hit.sfam, q->sfam))
      n = 3;
    else if (!strcmp(hit.fold, q->fold))
      n = 2;
    else if (!strcmp(hit.cl, q->cl))
      n = 1;
    else
      n = 0;

    char line[LINELEN];

    sprintf(line, "%-20s %-10s %1i %5i %3i %8.3f %7.2f %6.2f %7.2f %8.3f\n",
            hit.name, hit.fam, n, hit.L, hit.matched_cols, -1.443 * hit.logPval,
            -hit.score_aass, hit.Probab, hit.score, -1.443 * hit.logEval);

    outbuffer << line;
  }
}

void HitList::WriteToAlifile(HMM* q, bool scop_only) {
  Hit hit;
  int i = 0, n;
  Hash<int> twice(10000); // make sure only one hit per HMM is listed
  twice.Null(-1);
  FILE* alitabf = NULL;
  if (strcmp(par.alitabfile, "stdout"))
    alitabf = fopen(par.alitabfile, "w");
  else
    alitabf = stdout;
  if (!alitabf)
    OpenFileError(par.alitabfile);

  fprintf(alitabf, "NAME  %s\n", q->longname);
  fprintf(alitabf, "FAM   %s\n", q->fam);
  fprintf(alitabf, "FILE  %s\n", q->file);
  fprintf(alitabf, "LENG  %i\n", q->L);
  fprintf(alitabf, "\n");

  Reset();
  while (!End()) {
    i++;
    hit = ReadNext();
    if (scop_only
        && (!strncmp(hit.name, "cl|", 3) || !strncmp(hit.name, "UP20|", 5)
            || !strncmp(hit.name, "NR20|", 5)))
      continue;
    if (twice[hit.name] == 1)
      continue; // better hit with same HMM has been listed already
    twice.Add(hit.name, 1);
    //if template and query are from the same superfamily
    if (!strcmp(hit.name, q->name))
      n = 5;
    else if (!strcmp(hit.fam, q->fam))
      n = 4;
    else if (!strcmp(hit.sfam, q->sfam))
      n = 3;
    else if (!strcmp(hit.fold, q->fold))
      n = 2;
    else if (!strcmp(hit.cl, q->cl))
      n = 1;
    else
      n = 0;

    if (hit.P_posterior != NULL) {
      fprintf(alitabf,
          "\nHit %3i (%-20s %-10s Rel: %i  LOG-PVA: %6.2f  LOG-EVAL: %6.2f  Score: %6.2f  Probab: %6.2f):\n    i     j  score     SS  probab\n",
          i, hit.name, hit.fam, n, -1.443 * hit.logPval, -1.443 * hit.logEval,
          hit.score, hit.Probab);
      for (int step = hit.nsteps; step >= 1; step--)
        if (hit.states[step] >= MM)
          fprintf(alitabf, "%5i %5i %6.2f %6.2f %7.4f\n", hit.i[step],
              hit.j[step], hit.S[step], hit.S_ss[step], hit.P_posterior[step]);
    }
    else {
      fprintf(alitabf,
          "\nHit %3i (%-20s %-10s Rel: %i  LOG-PVA: %6.2f  LOG-EVAL: %6.2f  Score: %6.2f  Probab: %6.2f):\n    i     j  score     SS\n",
          i, hit.name, hit.fam, n, -1.443 * hit.logPval, -1.443 * hit.logEval,
          hit.score, hit.Probab);
      for (int step = hit.nsteps; step >= 1; step--)
        if (hit.states[step] >= MM)
          fprintf(alitabf, "%5i %5i %6.2f %6.2f\n", hit.i[step], hit.j[step],
              hit.S[step], hit.S_ss[step]);
    }
  }
  fclose(alitabf);
}

/////////////////////////////////////////////////////////////////////////////////////
//// Evaluate the *negative* log likelihood of the data at the vertex v = (lamda,mu)
////    p(s) = lamda * exp{ -exp[-lamda*(s-mu)] - lamda*(s-mu) } = lamda * exp( -exp(-x) - x) 
/////////////////////////////////////////////////////////////////////////////////////
double HitList::LogLikelihoodEVD(double* v) {
  double sum = 0.0, sumw = 0.0;
  for (int i = 0; i < Nprof; i++) {
    double x = v[0] * (score[i] - v[1]);
    sum += weight[i] * (exp(-x) + x);
    sumw += weight[i];
  }
  return sum - sumw * log(v[0]);
}
// Static wrapper-function for calling the nonstatic member function LogLikelihoodEVD() 
// ( see http://www.newty.de/fpt/callback.html#member )
double HitList::LogLikelihoodEVD_static(void* pt2hitlist, double* v) {
  HitList* mySelf = (HitList*) pt2hitlist; // explicitly cast to a pointer to Hitlist
  return mySelf->LogLikelihoodEVD(v);                // call member function
}

/////////////////////////////////////////////////////////////////////////////////////
//// Subroutine to FindMin: try new point given by highest point ihigh and fac and replace ihigh if it is lower 
/////////////////////////////////////////////////////////////////////////////////////
double HitList::TryPoint(const int ndim, double* p, double* y, double* psum,
    int ihigh, double fac, double (*Func)(void* pt2hitlist, double* v)) {
  // New point p_try = p_c + fac*(p_high-p_c),
  // where p_c = ( sum_i (p_i) - p_high)/ndim is the center of ndim other points
  // => p_try = fac1*sum_i(p_i) + fac2*p_high
  double fac1 = (1. - fac) / ndim;
  double fac2 = fac - fac1;
  double ptry[ndim];   //new point to try out
  double ytry;         //function value of new point 
  int j;               //index for the ndim parameters

  for (j = 0; j < ndim; j++)
    ptry[j] = psum[j] * fac1 + p[ihigh * ndim + j] * fac2;
  ytry = (*Func)(this, ptry);
  if (ytry <= y[ihigh]) {
//       if (v>=4) printf("Trying:                  %-7.3f %-7.3f %-7.3f -> accept\n",ptry[0],ptry[1],ytry);
    y[ihigh] = ytry;
    for (j = 0; j < ndim; j++) {
      psum[j] += ptry[j] - p[ihigh * ndim + j]; //update psum[j]
      p[ihigh * ndim + j] = ptry[j];            //replace p[ihigh] with ptry
    }	                           //Note: ihigh is now not highest point anymore!
  }
//   else if (v>=4) printf("Trying:                  %-7.3f %-7.3f %-7.3f -> reject\n",ptry[0],ptry[1],ytry);

  return ytry;
}

/////////////////////////////////////////////////////////////////////////////////////
////Find minimum with simplex method of Nelder and Mead (1965)
/////////////////////////////////////////////////////////////////////////////////////
float HitList::FindMin(const int ndim, double* p, double* y, double tol,
    int& nfunc, double (*Func)(void* pt2hitlist, double* v)) {
  const int MAXNFUNC = 99; //maximum allowed number of function evaluations
  int ihigh;    //index of highest point on simplex
  int inext;    //index of second highest point on simplex
  int ilow;     //index of lowest point on simplex
  int i;        //index for the ndim+1 points
  int j;        //index for the ndim parameters
  double rtol; //tolerance: difference of function value between highest and lowest point of simplex
  double temp;   //dummy
  double ytry;   //function value of trial point
  double psum[ndim]; //psum[j] = j'th coordinate of sum vector (sum over all vertex vectors)

  nfunc = 0;    //number of function evaluations =0
  //Calculate sum vector psum[j]
  for (j = 0; j < ndim; j++) {
    psum[j] = p[j];
    for (i = 1; i < ndim + 1; i++)
      psum[j] += p[i * ndim + j];
  }
  
  // Repeat finding better points in simplex until rtol<tol
  while (1) {
    // Find indices for highest, next highest and lowest point
    ilow = 0;
    if (y[0] > y[1]) {
      inext = 1;
      ihigh = 0;
    }
    else {
      inext = 0;
      ihigh = 1;
    }
    for (i = 0; i < ndim + 1; i++) {
      if (y[i] <= y[ilow])
        ilow = i;
      if (y[i] > y[ihigh]) {
        inext = ihigh;
        ihigh = i;
      }
      else if (y[i] > y[inext] && i != ihigh)
        inext = i;
    }

    // If tolerance in y is smaller than tol swap lowest point to index 0 and break -> return
    rtol = 2. * fabs(y[ihigh] - y[ilow])
        / (fabs(y[ihigh]) + fabs(y[ilow]) + 1E-10);
    if (rtol < tol) {
      temp = y[ilow];
      y[ilow] = y[0];
      y[0] = temp;
      for (j = 0; j < ndim; j++) {
        temp = p[ilow * ndim + j];
        p[ilow * ndim + j] = p[j];
        p[j] = temp;
      }
      break;
    }

    // Max number of function evaluations exceeded?
    if (nfunc >= MAXNFUNC) {
      if (v)
        fprintf(stderr,
            "\nWARNING: maximum likelihood fit of score distribution did not converge.\n");
      return 1;
    }

    nfunc += 2;
    // Point-reflect highest point on the center of gravity p_c of the other ndim points of the simplex
    if (v >= 3)
      printf("%3i  %-7.3f  %-7.3f %-12.8f %-9.3E\n", nfunc, p[ilow * ndim],
          p[ilow * ndim + 1], y[ilow], rtol);
//       if (v>=2) printf("           %3i %-9.3E   %-7.3f %-7.3f %-7.3f   %-7.3f %-7.3f %-7.3f   %-7.3f %-7.3f %-7.3f\n",nfunc,rtol,p[ilow*ndim],p[ilow*ndim+1],y[ilow],p[inext*ndim],p[inext*ndim+1],y[inext],p[ihigh*ndim],p[ihigh*ndim+1],y[ihigh]);
    ytry = TryPoint(ndim, p, y, psum, ihigh, -1.0, Func); //reflect highest point on p_c

    if (ytry <= y[ilow]) {
      ytry = TryPoint(ndim, p, y, psum, ihigh, 2.0, Func); //expand: try new point 2x further away from p_c
//  	  if (v>=2) printf("Expanded:  %3i %-9.3E   %-7.3f %-7.3f %-7.3f   %-7.3f %-7.3f %-7.3f   %-7.3f %-7.3f %-7.3f\n",nfunc,rtol,p[ilow*ndim],p[ilow*ndim+1],y[ilow],p[inext*ndim],p[inext*ndim+1],y[inext],p[ihigh*ndim],p[ihigh*ndim+1],y[ihigh]);
    }
    else if (ytry >= y[inext]) {
      // The new point is worse than the second worst point
      temp = y[ihigh];
      ytry = TryPoint(ndim, p, y, psum, ihigh, 0.5, Func); //contract simplex by 0.5 along (p_high-p_c
//  	  if (v>=2) printf("Compressed:%3i %-9.3E   %-7.3f %-7.3f %-7.3f   %-7.3f %-7.3f %-7.3f   %-7.3f %-7.3f %-7.3f\n",nfunc,rtol,p[ilow*ndim],p[ilow*ndim+1],y[ilow],p[inext*ndim],p[inext*ndim+1],y[inext],p[ihigh*ndim],p[ihigh*ndim+1],y[ihigh]);
      if (ytry >= temp) {
        // Trial point is larger than worst point => contract simplex by 0.5 towards lowest point
        for (i = 0; i < ndim + 1; i++) {
          if (i != ilow) {
            for (j = 0; j < ndim; j++)
              p[i * ndim + j] = 0.5 * (p[i * ndim + j] + p[ilow + j]);
            y[i] = (*Func)(this, p + i * ndim);
// 		      y[i] = (*Func)(p+i*ndim);
          }
        }
        nfunc += ndim;
//  	      if (v>=2) printf("Contracted:%3i %-9.3E   %-7.3f %-7.3f  %-7.3f  %-7.3f %-7.3f  %-7.3f  %-7.3f %-7.3f  %-7.3f\n",nfunc,rtol,p[ilow*ndim],p[ilow*ndim+1],y[ilow],p[inext*ndim],p[inext*ndim+1],y[inext],p[ihigh*ndim],p[ihigh*ndim+1],y[ihigh]);

        //Calculate psum[j]
        for (j = 0; j < ndim; j++) {
          psum[j] = p[j];
          for (i = 1; i < ndim + 1; i++)
            psum[j] += p[i * ndim + j];
        }
      }
    }
    else
      nfunc--;
  }
  return (float) rtol;
}

/////////////////////////////////////////////////////////////////////////////////////
//// Do a maximum likelihod fit of the scores with an EV distribution with parameters lamda and mu 
/////////////////////////////////////////////////////////////////////////////////////
void HitList::MaxLikelihoodEVD(HMM* q, int nbest) {
  double tol = 1E-6; // Maximum relative tolerance when minimizing -log(P)/N (~likelihood)
  static char first_call = 1;
  static Hash<int> size_fam(MAXPROF / 10); // Hash counts number of HMMs in family
  static Hash<int> size_sfam(MAXPROF / 10); // Hash counts number of families in superfamily
  Hash<int> excluded(50); // Hash containing names of superfamilies to be excluded from fit
  size_fam.Null(0);     // Set int value to return when no data can be retrieved
  size_sfam.Null(0);    // Set int value to return when no data can be retrieved
  excluded.Null(0);     // Set int value to return when no data can be retrieved
  Hit hit;

  double mu;            // EVD[mu,lam](x) = exp(-exp(-(x-mu)/lam)) = P(score<=x)
  double vertex[2 * 3];       // three vertices of the simplex in lamda-mu plane
  double yvertex[3]; // log likelihood values at the three vertices of the simplex
  int nfunc = 0;                     // number of function calls
  double sum_weights = 0.0;
  float sum_scores = 0.0;
  float rtol;

  if (first_call == 1) {
    first_call = 0;
    // Count how many HMMs are in each family
    if (v >= 4)
      printf(
          "  count number of profiles in each family and families in each superfamily ...\n");
    Reset();
    while (!End()) {
      hit = ReadNext();
      if (!size_fam.Contains(hit.fam))
        (*size_sfam(hit.sfam))++; //Add one to hash element for superfamily
      (*size_fam(hit.fam))++;              //Add one to hash element for family
      //      printf("size(%s)=%i name=%s\n",hit.fam,*size_fam(hit.fam),hit.name)
    }
    fams = size_fam.Size();
    sfams = size_sfam.Size();
    if (v >= 3)
      printf(
          "%-3i HMMs from %i families and %i superfamilies searched. Found %i hits\n",
          N_searched, fams, sfams, Size());
  }

  // Query has SCOP family identifier?
  if (q->fam && q->fam[0] >= 'a' && q->fam[0] <= 'k' && q->fam[1] == '.') {
    char sfamid[NAMELEN];
    char* ptr_in_fam = q->fam;
    while ((ptr_in_fam = strwrd(sfamid, ptr_in_fam, '-'))) {
      char* ptr = strrchr(sfamid, '.');
      if (ptr)
        *ptr = '\0';
      excluded.Add(sfamid);
// 	  fprintf(stderr,"Exclude SCOP superfamily %s  ptr_in_fam='%s'\n",sfamid,ptr_in_fam);
    }
  }
  // Exclude best superfamilies from fit
  else if (nbest > 0) {
    if (sfams < 97 + nbest)
      return;

    // Find the nbest best-scoring superfamilies for exclusion from first ML fit
    if (v >= 4)
      printf(
          "  find %i best-scoring superfamilies to exclude from first fit  ...\n",
          nbest);
    hit = Smallest();
    excluded.Add(hit.sfam);
//       printf("Exclude in first round: %s %8.2f %s\n",hit.name,hit.score_aass,hit.sfam);
    while (excluded.Size() < nbest) {
      Reset();
      while (!End() && excluded.Contains(ReadNext().sfam))
        ;
      hit = ReadCurrent();
      while (!End()) {
        if (ReadNext() < hit && !excluded.Contains(ReadCurrent().sfam))
          hit = ReadCurrent();
      }
      excluded.Add(hit.sfam);
// 	  printf("Exclude in first round: %s %8.2f %s %i %i\n",hit.name,hit.score_aass,hit.sfam,excluded.Size(),excluded.Contains(hit.sfam));
    }
    tol = 0.01 / size_sfam.Size(); // tol=1/N would lead to delta(log-likelihood)~1 (where N ~ number of superfamilies) since (1+1/N)^N = e
  }
  else {
    // Find the best-scoring superfamilies from first fit for exclusion from second ML fit
    if (v >= 4)
      printf(
          "  find best-scoring superfamilies to exclude from second fit  ...\n");
    Reset();
    while (!End()) {
      hit = ReadNext();
      if (hit.Eval < 0.05)
        excluded.Add(hit.sfam); // changed from 0.5 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    }
    tol = 0.001 / size_sfam.Size(); // tol=1/N would lead to delta(log-likelihood)~1 (where N ~ number of superfamilies) since (1+1/N)^N = e
  }

  // Put scores into score[] and weights into weight[]
  if (v >= 3)
    printf("  generate scores and weights array for ML fitting ...\n");
  Nprof = 0;
  Reset();
  while (!End()) {
    hit = ReadNext();
    if (hit.irep > 1)
      continue;                   //Use only best hit per template
    if (Nprof >= MAXPROF)
      break;

    char sfamid[NAMELEN];
    char* ptr_in_fam = hit.fam;
    while ((ptr_in_fam = strwrd(sfamid, ptr_in_fam, '-'))) {
      char* ptr = strrchr(sfamid, '.');
      if (ptr)
        *ptr = '\0';
      if (excluded.Contains(sfamid))
        break;    //HMM is among superfamilies to be excluded
    }
    if (excluded.Contains(sfamid)) {
      if (v >= 3)
        fprintf(stderr, "Exclude hit %s (family %s contains %s)\n", hit.name,
            hit.fam, sfamid);
      continue;
    }
//       ScopID(hit.cl,hit.fold,hit.sfam,hit.fam);     //Get scop superfamily code for template
//       if (*hit.sfam=='\0' || excluded.Contains(hit.sfam)) continue;    // skip HMM

    score[Nprof] = hit.score;
    weight[Nprof] = 1. / size_fam[hit.fam] / size_sfam[hit.sfam];
    sum_scores += hit.score * weight[Nprof];
    sum_weights += weight[Nprof];

    //DEBUG
//       if (v>=4) printf("%-10.10s   %-12.12s %-3i   %-12.12s %-3i   %6.4f   %6.4f  %7.1f\n",hit.name,hit.fam,size_fam[hit.fam],hit.sfam,size_sfam[hit.sfam],1./size_fam[hit.fam]/size_sfam[hit.sfam],sum,hit.score);
    Nprof++;
  }
  //DEBUG
  if (v >= 3)
    printf("%i hits used for score distribution\n", Nprof);
  // for (int i=0; i<Nprof; i++) printf("%3i  score=%8.3f  weight=%7.5f\n",i,score[i],weight[i]);

  // Set simplex vertices and function values
  mu = sum_scores / sum_weights - 0.584 / LAMDA;
  if (par.loc) // fit only in local mode; in global mode use fixed value LAMDA and mu mean score
  {
    double (*Func)(void*, double*);
    Func = HitList::LogLikelihoodEVD_static;

    if (nbest > 0) {
      vertex[0] = LAMDA;
      vertex[1] = mu;
    }  /////////////////////////////////////////// DEBUG
    else {
      vertex[0] = q->lamda;
      vertex[1] = mu;
    }
    vertex[2] = vertex[0] + 0.1;
    vertex[3] = vertex[1];
    vertex[4] = vertex[0];
    vertex[5] = vertex[1] + 0.2;
    yvertex[0] = Func(this, vertex);
    yvertex[1] = Func(this, vertex + 2);
    yvertex[2] = Func(this, vertex + 4);

    // Find lam and mu that minimize negative log likelihood of data
    if (v >= 3)
      printf(
          "Fitting to EVD by maximum likelihood...\niter lamda       mu    -log(P)/N   tol\n");
    rtol = FindMin(2, vertex, yvertex, tol, nfunc, Func);
    if (v >= 3)
      printf("%3i  %-7.3f  %-7.2f     %-7.3f %-7.1E\n\n", nfunc, vertex[0],
          vertex[1], yvertex[0] - (1.5772 - log(vertex[0])), rtol);
//       printf("HHsearch lamda=%-6.3f   mu=%-6.3f\n",vertex[0],vertex[1]);
  }
  else {
    vertex[0] = LAMDA_GLOB;
    vertex[1] = mu;
  }
  
  // Set lamda and mu of profile
  q->lamda = vertex[0];
  q->mu = vertex[1];

  // Set P-values and E-values
  // CHECK UPDATE FROM score=-logpval to score=-logpval+SSSCORE2NATLOG*score_ss !!!!
  Reset();
  while (!End()) {
    hit = ReadNext();
    hit.logPval = logPvalue(hit.score, q->lamda, q->mu);
    hit.Pval = Pvalue(hit.score, q->lamda, q->mu);
    hit.CalcEvalScoreProbab(N_searched, q->lamda); // calculate Evalue, score_aass, Proba from logPval and score_ss

    Overwrite(hit);
  }
}


/////////////////////////////////////////////////////////////////////////////////////
//// Calculate HHblits composite E-values 
/////////////////////////////////////////////////////////////////////////////////////
void HitList::CalculateHHblitsEvalues(HMM* q) {
  Hit hit;
  double alpha = 0;
  double log_Pcut = log(par.prefilter_evalue_thresh / par.dbsize);
  double log_dbsize = log((double) par.dbsize);
  //printf("log_Pcut: %7.4f  Pcut: %7.4f DBsize: %10i   a: %4.2f  b: %4.2f  c: %4.2f\n",log_Pcut, exp(log_Pcut), par.dbsize, par.alphaa, par.alphab, par.alphac);
  // int nhits=0;?

  Reset();
  while (!End()) {
    hit = ReadNext();

    // if (nhits++<50)
    // 	printf("before correction  Eval: %7.4g    logEval: %7.4f\n",hit.Eval, hit.logEval);

    alpha = par.alphaa
        + par.alphab * (hit.Neff_HMM - 1)
            * (1 - par.alphac * (q->Neff_HMM - 1));

    hit.Eval = exp(hit.logPval + log_dbsize + (alpha * log_Pcut));
    hit.logEval = hit.logPval + log_dbsize + (alpha * log_Pcut);

    // if (nhits++<50) //DEBUG?????
    //  	printf("                   Eval: %11.3E    logEval: %7.4f   logPval: %7.4f   alpha: %7.4f   Neff_T: %5.2f  Neff_Q: %5.2f\n",hit.Eval, hit.logEval, hit.logPval, alpha, hit.Neff_HMM, q->Neff_HMM);//DEBUG?????

    Overwrite(hit);   // copy hit object into current position of hitlist
  }
  ResortList(); // use InsertSort to resort list according to sum of minus-log-Pvalues
}

/////////////////////////////////////////////////////////////////////////////////////
//// Calculate Pvalues as a function of query and template lengths and diversities
/////////////////////////////////////////////////////////////////////////////////////
void HitList::CalculatePvalues(HMM* q) {
  Hit hit;
  float lamda = LAMDA_GLOB, mu = 3.0;   // init for global search
  const float log1000 = log(1000.0);

  if (N_searched == 0)
    N_searched = 1;
  if (v >= 3)
    printf(
        "Calculate Pvalues as a function of query and template lengths and diversities...\n");
  Reset();
  while (!End()) {
    hit = ReadNext();

    if (par.loc) {
      lamda = lamda_NN(log(q->L) / log1000, log(hit.L) / log1000,
          q->Neff_HMM / 10.0, hit.Neff_HMM / 10.0);
      mu = mu_NN(log(q->L) / log1000, log(hit.L) / log1000, q->Neff_HMM / 10.0,
          hit.Neff_HMM / 10.0);
// 	  if (v>=3 && nhits++<20) 
// 	     printf("hit=%-10.10s Lq=%-4i  Lt=%-4i  Nq=%5.2f  Nt=%5.2f  =>  lamda=%-6.3f  mu=%-6.3f\n",hit.name,q->L,hit.L,q->Neff_HMM,hit.Neff_HMM,lamda,mu);
    }
    hit.logPval = logPvalue(hit.score, lamda, mu);
    hit.Pval = Pvalue(hit.score, lamda, mu);
    hit.CalcEvalScoreProbab(N_searched, lamda); // calculate Evalue, score_aass, Proba from logPval and score_ss

    Overwrite(hit);
  }

  SortList();
  Reset();
  return;
}

/////////////////////////////////////////////////////////////////////////////////////
//// Calculate Pvalues from calibration of  0: query HMM, 1:template HMMs, 2: both
/////////////////////////////////////////////////////////////////////////////////////
void HitList::GetPvalsFromCalibration(HMM* q) {
  Hit hit;
  char warn = 0;
  if (N_searched == 0)
    N_searched = 1;
  if (v >= 2) {
    switch (par.calm) {
      case 0:
        printf(
            "Using lamda=%-5.3f and mu=%-5.2f from calibrated query HMM %s. \n",
            q->lamda, q->mu, q->name);
        printf(
            "Note that HMMs need to be recalibrated when changing HMM-HMM alignment options.\n");
        break;
      case 1:
        printf(
            "Using score distribution parameters lamda and mu from database HMMs \n");
        break;
      case 2:
        printf(
            "Combining score distribution parameters lamda and mu from query and database HMMs\n");
        printf(
            "Note that HMMs need to be recalibrated when changing HMM-HMM alignment options.\n");
        break;
    }
  }
  Reset();
  while (!End()) {
    hit = ReadNext();
    if (par.calm == 0 || (hit.logPvalt == 0)) {
      hit.logPval = logPvalue(hit.score, q->lamda, q->mu);
      hit.Pval = Pvalue(hit.score, q->lamda, q->mu);
      if (par.calm > 0 && warn++ < 1 && v >= 1)
        printf(
            "WARNING: some template HMM (e.g. %s) are not calibrated. Using query calibration.\n",
            hit.name);
    }
    else if (par.calm == 1) {
      hit.logPval = hit.logPvalt;
      hit.Pval = hit.Pvalt;
    }
    else if (par.calm == 2) {
      hit.logPval = 0.5
          * (logPvalue(hit.score, q->lamda, q->mu) + hit.logPvalt);
      hit.Pval = sqrt(Pvalue(hit.score, q->lamda, q->mu) * hit.Pvalt);
      if (v >= 5)
        printf(
            "Score: %7.1f  lamda: %7.1f  mu: %7.1f  P-values:  query-calibrated: %8.2G   template-calibrated: %8.2G   geometric mean: %8.2G\n",
            hit.score, q->lamda, q->mu, Pvalue(hit.score, q->lamda, q->mu),
            hit.Pvalt, hit.Pval);
    }

    hit.CalcEvalScoreProbab(N_searched, q->lamda); // calculate Evalue, score_aass, Proba from logPval and score_ss

    Overwrite(hit);
  }
  SortList();
  Reset();
  return;
}

/////////////////////////////////////////////////////////////////////////////////////
//// Generate new, very short hit list with only non-overlapping match alignments
/////////////////////////////////////////////////////////////////////////////////////
// void HitList::GetPvalsFromCalibration(Hitlist &novlap_hitlist)
// {

// }
