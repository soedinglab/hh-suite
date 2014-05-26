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
void HitList::PrintHitList(HMM* q, std::stringstream& out, const unsigned int maxdbstrlen, const int z, const int Z, const float p, const double E, const int argc, char** argv) {
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
  for (int i = 0; i < argc; i++)
    if (strlen(argv[i]) <= maxdbstrlen)
      out << argv[i] << " ";
    else
      out << "<" << strlen(argv[i]) << "characters> ";
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
    if (nhits >= Z)
      break;       //max number of lines reached?
    if (nhits >= z && hit.Probab < p)
      break;
    if (nhits >= z && hit.Eval > E)
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

void HitList::PrintHitList(HMM* q, char* outfile, const unsigned int maxdbstrlen, const int z, const int Z, const float p, const double E, const int argc, char** argv) {
  std::stringstream out;
  PrintHitList(q, out, maxdbstrlen, z, Z, p, E, argc, argv);

  if (strcmp(outfile, "stdout") == 0) {
    std::cout << out.str();
  }
  else {
    std::ofstream outf(outfile);
    if (!outf)
      OpenFileError(outfile, __FILE__, __LINE__, __func__);

    outf << out.str();

    outf.close();
  }
}

void HitList::PrintHHR(HMM* q, char* outfile, const unsigned int maxdbstrlen,
		const char showconf, const char showcons, const char showdssp, const char showpred,
		const int b, const int B, const int z, const int Z, const int aliwidth, const int nseqdis,
		const float p, const double E, const int argc, char** argv, const float S[20][20]) {
  std::stringstream out;
  PrintHHR(q, out, maxdbstrlen, showconf, showcons, showdssp, showpred, b, B, z, Z, aliwidth, nseqdis, p, E, argc, argv, S);

  if (strcmp(outfile, "stdout") == 0) {
    std::cout << out.str();
  }
  else {
    std::ofstream outf(outfile);

    if (!outf)
      OpenFileError(outfile, __FILE__, __LINE__, __func__);

    outf << out.str();
    outf.close();
  }
}

void HitList::PrintHHR(HMM* q, std::stringstream& out, const unsigned int maxdbstrlen,
		const char showconf, const char showcons, const char showdssp, const char showpred,
		const int b, const int B, const int z, const int Z, const int aliwidth, const int nseqdis,
		const float p, const double E, const int argc, char** argv, const float S[20][20]) {
  PrintHitList(q, out, maxdbstrlen, z, Z, p, E, argc, argv);
  PrintAlignments(q, out, showconf, showcons, showdssp, showpred, p, aliwidth, nseqdis, b, B, E, S, 0);
}

/////////////////////////////////////////////////////////////////////////////////////
// Print alignments of query sequences against hit sequences
/////////////////////////////////////////////////////////////////////////////////////
void HitList::PrintAlignments(HMM* q, char* outfile, const char showconf, const char showcons,
		const char showdssp, const char showpred, const float p, const int aliwidth, const int nseqdis,
		const int b, const int B, const double E, const float S[20][20], char outformat) {
  std::stringstream out;
  PrintAlignments(q, out, showconf, showcons, showdssp, showpred, p, aliwidth, nseqdis, b, B, E, S, outformat);

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
      OpenFileError(outfile, __FILE__, __LINE__, __func__);

    outf << out.str();
    outf.close();
  }
}

void HitList::PrintAlignments(HMM* q, std::stringstream& out, const char showconf, const char showcons,
		const char showdssp, const char showpred, const float p, const int aliwidth, const int nseqdis,
		const int b, const int B, const double E, const float S[20][20], char outformat) {
  FullAlignment qt_ali(nseqdis + 10); // maximum 10 annotation (pseudo) sequences (ss_dssp, sa_dssp, ss_pred, ss_conf, consens,...)
  int nhits = 0;

  Reset();
  while (!End()) {
    if (nhits >= B)
      break;
    Hit hit = ReadNext();
    if (nhits >= b && hit.Probab < p)
      break;
    if (nhits >= b && hit.Eval > E)
      continue;
    nhits++;

    // Build double alignment of query against template sequences
    qt_ali.Build(q, hit, nseqdis, S);

    // Print out alignment
    // HHR format
    out << "No " << nhits << std::endl;
    if (outformat == 0)
    {
      qt_ali.PrintHeader(out, q, hit);
      qt_ali.PrintHHR(out, hit, showconf, showcons, showdssp, showpred, aliwidth);
    }
    // FASTA format
    else if (outformat == 1)
    {
      qt_ali.PrintFASTA(out, hit, showcons, showdssp, showpred, aliwidth);
    }
    // A2M format
    else if (outformat == 2)
    {
      qt_ali.PrintA2M(out, hit, showcons, showdssp, showpred, aliwidth);
    }
    // A3m format
    else
    {
      qt_ali.PrintA3M(out, hit, showcons, showdssp, showpred, aliwidth);
    }

    qt_ali.FreeMemory();
  }
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
    	std::cerr << "Warning in " << __FILE__ << ":" << __LINE__ << ": " << __func__ << ":" << std::endl;
      std::cerr << "\tcould not open \'" << outputfile << std::endl;
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

void HitList::WriteToAlifile(HMM* q, char* alitabfile, bool scop_only) {
	std::stringstream out;
	WriteToAlifile(q, out, scop_only);

  if (strcmp(alitabfile, "stdout") == 0) {
	  std::cout << out.str();
  }
  else {
    std::ofstream alitabf(alitabfile);

    if (!alitabf.good())
    	OpenFileError(alitabfile, __FILE__, __LINE__, __func__);

    alitabf << out.str();

    alitabf.close();
  }
}

void HitList::WriteToAlifile(HMM* q, std::stringstream& out, bool scop_only) {
  Hash<int> twice(10000); // make sure only one hit per HMM is listed
  twice.Null(-1);

  char line[LINELEN];

  out << "NAME  " <<  q->longname << std::endl;
  out << "FAM   " << q->fam << std::endl;
  out << "FILE  " << q->file << std::endl;
  out << "LENG  " << q->L << std::endl;
  out << std::endl;

  int i = 0;
  Reset();
  while (!End()) {
    i++;
    Hit hit = ReadNext();
    if (scop_only
        && (!strncmp(hit.name, "cl|", 3) || !strncmp(hit.name, "UP20|", 5)
            || !strncmp(hit.name, "NR20|", 5)))
      continue;

    if (twice[hit.name] == 1)
      continue; // better hit with same HMM has been listed already
    twice.Add(hit.name, 1);

    //if template and query are from the same superfamily
    int n;
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
      sprintf(line,
          "\nHit %3i (%-20s %-10s Rel: %i  LOG-PVA: %6.2f  LOG-EVAL: %6.2f  Score: %6.2f  Probab: %6.2f):\n    i     j  score     SS  probab\n",
          i, hit.name, hit.fam, n, -1.443 * hit.logPval, -1.443 * hit.logEval,
          hit.score, hit.Probab);
      out << line;

      for (int step = hit.nsteps; step >= 1; step--) {
        if (hit.states[step] >= MM) {
          sprintf(line, "%5i %5i %6.2f %6.2f %7.4f\n", hit.i[step],
              hit.j[step], hit.S[step], hit.S_ss[step], hit.P_posterior[step]);
          out << line;
        }
      }
    }
    else {
      sprintf(line,
          "\nHit %3i (%-20s %-10s Rel: %i  LOG-PVA: %6.2f  LOG-EVAL: %6.2f  Score: %6.2f  Probab: %6.2f):\n    i     j  score     SS\n",
          i, hit.name, hit.fam, n, -1.443 * hit.logPval, -1.443 * hit.logEval,
          hit.score, hit.Probab);
      out << line;

      for (int step = hit.nsteps; step >= 1; step--) {
        if (hit.states[step] >= MM) {
          sprintf(line, "%5i %5i %6.2f %6.2f\n", hit.i[step], hit.j[step],
              hit.S[step], hit.S_ss[step]);
          out << line;
        }
      }
    }
  }
}


/////////////////////////////////////////////////////////////////////////////////////
//// Calculate HHblits composite E-values
/////////////////////////////////////////////////////////////////////////////////////
void HitList::CalculateHHblitsEvalues(HMM* q, const int dbsize,
		const float alphaa, const float alphab, const float alphac, const double prefilter_evalue_thresh) {
  Hit hit;
  double alpha = 0;
  double log_Pcut = log(prefilter_evalue_thresh / dbsize);
  double log_dbsize = log((double) dbsize);
  //printf("log_Pcut: %7.4f  Pcut: %7.4f DBsize: %10i   a: %4.2f  b: %4.2f  c: %4.2f\n",log_Pcut, exp(log_Pcut), par.dbsize, par.alphaa, par.alphab, par.alphac);
  // int nhits=0;?

  Reset();
  while (!End()) {
    hit = ReadNext();

    // if (nhits++<50)
    // 	printf("before correction  Eval: %7.4g    logEval: %7.4f\n",hit.Eval, hit.logEval);

    alpha = alphaa
        + alphab * (hit.Neff_HMM - 1)
            * (1 - alphac * (q->Neff_HMM - 1));

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
void HitList::CalculatePvalues(HMM* q, const char loc, const char ssm, const float ssw) {
  Hit hit;
  float lamda = LAMDA_GLOB, mu = 3.0;   // init for global search
  const float log1000 = log(1000.0);

  if (N_searched == 0)
    N_searched = 1;

  HH_LOG(LogLevel::DEBUG) << "Calculate Pvalues as a function of query and template lengths and diversities..." << std::endl;

  Reset();
  while (!End()) {
    hit = ReadNext();

    if (loc) {
      lamda = lamda_NN(log(q->L) / log1000, log(hit.L) / log1000,
          q->Neff_HMM / 10.0, hit.Neff_HMM / 10.0);
      mu = mu_NN(log(q->L) / log1000, log(hit.L) / log1000, q->Neff_HMM / 10.0,
          hit.Neff_HMM / 10.0);
    }
    hit.logPval = logPvalue(hit.score, lamda, mu);
    hit.Pval = Pvalue(hit.score, lamda, mu);
    hit.CalcEvalScoreProbab(N_searched, lamda, loc, ssm, ssw); // calculate Evalue, score_aass, Proba from logPval and score_ss

    Overwrite(hit);
  }

  SortList();
  Reset();
}
