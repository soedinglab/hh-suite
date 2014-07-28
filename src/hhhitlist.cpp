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

SearchCounter::SearchCounter() {
}

SearchCounter::~SearchCounter() {
}

int SearchCounter::getCounter() {
  return already_seen.size();
}

void SearchCounter::append(std::string id) {
  already_seen.insert(id);
}

/////////////////////////////////////////////////////////////////////////////////////
// Print summary listing of hits
/////////////////////////////////////////////////////////////////////////////////////
void HitList::PrintHitList(HMM* q, std::stringstream& out,
    const unsigned int maxdbstrlen, const int z, const int Z, const float p,
    const double E, const int argc, char** argv) {
  int nhits = 0;

  out << "Query         " << q->longname << std::endl;
  out << "Match_columns " << q->L << std::endl;
  out << "No_of_seqs    " << q->N_filtered << " out of " << q->N_in
      << std::endl;
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
  out
      << " No Hit                             Prob E-value P-value  Score    SS Cols Query HMM  Template HMM"
      << std::endl;
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

void HitList::PrintHitList(HMM* q, char* outfile,
    const unsigned int maxdbstrlen, const int z, const int Z, const float p,
    const double E, const int argc, char** argv) {
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
    const char showconf, const char showcons, const char showdssp,
    const char showpred, const int b, const int B, const int z, const int Z,
    const int aliwidth, const int nseqdis, const float p, const double E,
    const int argc, char** argv, const float S[20][20]) {
  std::stringstream out;
  PrintHHR(q, out, maxdbstrlen, showconf, showcons, showdssp, showpred, b, B, z,
      Z, aliwidth, nseqdis, p, E, argc, argv, S);

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

void HitList::PrintHHR(HMM* q, std::stringstream& out,
    const unsigned int maxdbstrlen, const char showconf, const char showcons,
    const char showdssp, const char showpred, const int b, const int B,
    const int z, const int Z, const int aliwidth, const int nseqdis,
    const float p, const double E, const int argc, char** argv,
    const float S[20][20]) {
  PrintHitList(q, out, maxdbstrlen, z, Z, p, E, argc, argv);
  PrintAlignments(q, out, showconf, showcons, showdssp, showpred, p, aliwidth,
      nseqdis, b, B, E, S, 0);
}

/////////////////////////////////////////////////////////////////////////////////////
// Print alignments of query sequences against hit sequences
/////////////////////////////////////////////////////////////////////////////////////
void HitList::PrintAlignments(HMM* q, char* outfile, const char showconf,
    const char showcons, const char showdssp, const char showpred,
    const float p, const int aliwidth, const int nseqdis, const int b,
    const int B, const double E, const float S[20][20], char outformat) {
  std::stringstream out;
  PrintAlignments(q, out, showconf, showcons, showdssp, showpred, p, aliwidth,
      nseqdis, b, B, E, S, outformat);

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

void HitList::PrintAlignments(HMM* q, std::stringstream& out,
    const char showconf, const char showcons, const char showdssp,
    const char showpred, const float p, const int aliwidth, const int nseqdis,
    const int b, const int B, const double E, const float S[20][20],
    char outformat) {
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
    if (outformat == 0) {
      qt_ali.PrintHeader(out, q, hit);
      qt_ali.PrintHHR(out, hit, showconf, showcons, showdssp, showpred,
          aliwidth);
    }
    // FASTA format
    else if (outformat == 1) {
      qt_ali.PrintFASTA(out, hit, showcons, showdssp, showpred, aliwidth);
    }
    // A2M format
    else if (outformat == 2) {
      qt_ali.PrintA2M(out, hit, showcons, showdssp, showpred, aliwidth);
    }
    // A3m format
    else {
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
      std::cerr << "Warning in " << __FILE__ << ":" << __LINE__ << ": "
          << __func__ << ":" << std::endl;
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
  outbuffer
      << "TARGET                FAMILY   REL  LEN  COL  LOG-PVA  S-AASS PROBAB  SCORE  LOG-EVAL\n";
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

void HitList::WriteToAlifile(HMM* q, char* alitabfile,
    const int b, const int B, const int z, const int Z,
    const float p, const double E) {

  std::stringstream out;
  WriteToAlifile(q, out, b, B, z, Z, p, E);

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

void HitList::WriteToAlifile(HMM* q, std::stringstream& out,
    const int b, const int B, const int z, const int Z,
    const float p, const double E) {
  char line[LINELEN];

  // Store all dbfiles and ftell positions of templates to be displayed and realigned
  int nhits = 0;
  Reset();
  while (!End()) {
    Hit hit = ReadNext();
    if (nhits >= imax(B, Z))
      break;
    if (nhits >= imax(b, z) && hit.Probab < p)
      break;
    if (nhits >= imax(b, z) && hit.Eval > E)
      continue;

    sprintf(line, ">%s\n", hit.longname);
    out << line;

    if (hit.P_posterior != NULL) {
      if (hit.nss_dssp >= 0) {
        // secondary structure determined by dssp 0:-  1:H  2:E  3:C  4:S  5:T  6:G  7:B
        sprintf(line, "    i     j  score     SS  probab  dssp\n");
        out << line;
        for (int step = hit.nsteps; step >= 1; step--)
          if (hit.states[step] == MM) {
            sprintf(line, "%5i %5i %6.2f %6.2f %7.4f %5c\n", hit.i[step],
                hit.j[step], hit.S[step], hit.S_ss[step], hit.P_posterior[step],
                hit.seq[hit.nss_dssp][hit.j[step]]);
            out << line;
          }
      }
      else {
        sprintf(line, "missing dssp\n");
        out << line;
        sprintf(line, "    i     j  score     SS  probab\n");
        out << line;
        for (int step = hit.nsteps; step >= 1; step--)
          if (hit.states[step] == MM) {
            sprintf(line, "%5i %5i %6.2f %6.2f %7.4f\n", hit.i[step],
                hit.j[step], hit.S[step], hit.S_ss[step],
                hit.P_posterior[step]);
            out << line;
          }
      }
    }
    else {
      sprintf(line, "    i     j  score     SS\n");
      out << line;
      for (int step = hit.nsteps; step >= 1; step--)
        if (hit.states[step] == MM) {
          sprintf(line, "%5i %5i %6.2f %6.2f\n", hit.i[step], hit.j[step],
              hit.S[step], hit.S_ss[step]);
          out << line;
        }
    }

    nhits++;
  }
}

/////////////////////////////////////////////////////////////////////////////////////
//// Calculate HHblits composite E-values
/////////////////////////////////////////////////////////////////////////////////////
void HitList::CalculateHHblitsEvalues(HMM* q, const int dbsize,
    const float alphaa, const float alphab, const float alphac,
    const double prefilter_evalue_thresh) {
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
        + alphab * (hit.Neff_HMM - 1) * (1 - alphac * (q->Neff_HMM - 1));

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
void HitList::CalculatePvalues(HMM* q, const char loc, const char ssm,
    const float ssw) {
  Hit hit;
  float lamda = LAMDA_GLOB, mu = 3.0;   // init for global search
  const float log1000 = log(1000.0);

  if (N_searched == 0)
    N_searched = 1;

  HH_LOG(LogLevel::DEBUG)
      << "Calculate Pvalues as a function of query and template lengths and diversities..."
      << std::endl;

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

void HitList::PrintMatrices(HMM* q, const char* matricesOutputFileName,
    const size_t max_number_matrices, const float S[20][20]) {
  std::stringstream out_ss(
      std::stringstream::in | std::stringstream::out
          | std::stringstream::binary);
  PrintMatrices(q, out_ss, max_number_matrices, S);

  if (strcmp(matricesOutputFileName, "stdout") == 0) {
    std::cout << out_ss.str();
  }
  else {
    std::ofstream out(matricesOutputFileName, std::ios::out | std::ios::binary);
    if (!out.good()) {
      std::cerr << "Warning in " << __FILE__ << ":" << __LINE__ << ": "
          << __func__ << ":" << std::endl;
      std::cerr << "\tcould not open \'" << matricesOutputFileName << std::endl;
      return;
    }

    out << out_ss.str();

    out.close();
  }
}

void HitList::PrintMatrices(HMM* q, std::stringstream& out,
    const size_t max_number_matrices, const float S[20][20]) {
  //limit matrices to par.max_number_matrices
  std::vector<Hit> hits;

  const float tolerance = 0.01;

  Reset();
  while (!End()) {
    Hit hit_cur = ReadNext();

    if (!hit_cur.forward_profile || !hit_cur.backward_profile) {
      continue;
    }

    float forward_profile_sum = 0.0;
    float backward_profile_sum = 0.0;

    for (int i = 1; i <= q->L; i++) {
      forward_profile_sum += hit_cur.forward_profile[i];
      backward_profile_sum += hit_cur.backward_profile[i];
    }

    if (forward_profile_sum < 1.0 + tolerance
        && forward_profile_sum > 1.0 - tolerance
        && backward_profile_sum < 1.0 + tolerance
        && backward_profile_sum > 1.0 - tolerance) {
      hits.push_back(hit_cur);
    }
  }

  std::vector<bool> picked_alignments(hits.size(), true);
  unsigned int chosen = hits.size();
  float matix_probability_threshold = 20;

  if (hits.size() == 0) {
    HH_LOG(LogLevel::WARNING) << "There are no alignment matrices to print!"
        << std::endl;
  }

  //remove duplicate alignments (mostly alignments which were realigned afterwards)
  for (int index1 = hits.size() - 1; index1 >= 0; index1--) {
    Hit it = hits[index1];
    if (it.Probab < matix_probability_threshold) {
      picked_alignments[index1] = false;
      chosen--;
    }
    else if (picked_alignments[index1]) {
      for (int index2 = index1 - 1; index2 >= 0; index2--) {
        Hit it_comp = hits[index2];

        if ((picked_alignments[index2] && strcmp(it.name, it_comp.name) == 0
            && it.irep == it_comp.irep)
            || it.Probab < matix_probability_threshold) {
          picked_alignments[index2] = false;
          chosen--;
        }
      }
    }
  }

  float** similarity_scores = new float*[hits.size()];
  for (unsigned int i = 0; i < hits.size(); i++) {
    similarity_scores[i] = new float[hits.size()];
  }

  for (unsigned int k = 0; k < hits.size(); k++) {
    similarity_scores[k][k] = 1.0;
    for (unsigned int k_comp = k + 1; k_comp < hits.size(); k_comp++) {
      similarity_scores[k][k_comp] = 0;

      Hit it = hits[k];
      Hit it_comp = hits[k_comp];

      for (int i = 1; i <= q->L; i++) {
        similarity_scores[k][k_comp] += pow(
            it.forward_profile[i] * it_comp.forward_profile[i], 0.5)
            + pow(it.backward_profile[i] * it_comp.backward_profile[i], 0.5);
      }

      similarity_scores[k][k_comp] /= 2.0;
      similarity_scores[k_comp][k] = similarity_scores[k][k_comp];
    }
  }

  //remove too similar alignments
  while (chosen > max_number_matrices) {
    float max_value = 0.0;
    int max_index = 0;

    for (unsigned int k = 0; k < hits.size(); k++) {
      float summed_similarity = 0;
      for (unsigned int k_prim = 0; k_prim < hits.size(); k_prim++) {
        if (picked_alignments[k_prim] && picked_alignments[k]) {
          summed_similarity += similarity_scores[k][k_prim];
        }
      }

      if (summed_similarity > max_value) {
        max_value = summed_similarity;
        max_index = k;
      }
    }

    picked_alignments[max_index] = false;
    chosen--;
  }

  for (unsigned int k = 0; k < hits.size(); k++) {
    delete[] similarity_scores[k];
  }
  delete[] similarity_scores;

  HH_LOG(LogLevel::INFO) << "Printing alignment matrices..." << std::endl;
  HH_LOG(LogLevel::INFO) << "Total number of alignments    : " << hits.size()
      << std::endl;
  HH_LOG(LogLevel::INFO) << "Number of accepted alignments : " << chosen
      << std::endl;

  if (chosen == 0) {
    HH_LOG(LogLevel::WARNING)
        << "Warning: No homologs found for printing matrix!" << std::endl;
    return;
  }

  //set outputstream
  const unsigned short int delimiter_16_bit = 0;
  const unsigned char delimiter_8_bit = 0;

  out.write(q->name, strlen(q->name));
  out.write(reinterpret_cast<const char*>(&delimiter_8_bit),
      sizeof(delimiter_8_bit));

  unsigned short int query_length = q->L;
  writeU16(out, query_length);

  for (size_t index = 0; index < hits.size(); index++) {
    if (!picked_alignments[index]) {
      continue;
    }

    Hit it = hits[index];

    const char* name = it.name;

    out.write(name, strlen(it.name));
    out.write(reinterpret_cast<const char*>(&delimiter_8_bit),
        sizeof(delimiter_8_bit));
    unsigned short int template_length = it.L;
    writeU16(out, template_length);

    unsigned char ali_probability = it.Probab;
    out.write(reinterpret_cast<const char*>(&ali_probability),
        sizeof(unsigned char));

    unsigned short int alignment_similarity = it.calculateSimilarity(q, S) * 10;
    writeS16(out, alignment_similarity);

    unsigned short int forwardProbability;
    float printForwardThreshold;
    float_to_16_bit(0.0, forwardProbability);
    bit_16_to_float(forwardProbability, printForwardThreshold);

    unsigned short int startPrintForward = q->L;
    unsigned short int endPrintForward = 0;

    for (unsigned short int i = 1; i <= q->L; i++) {
      float_to_16_bit(it.forward_profile[i], forwardProbability);
      if (it.forward_profile[i] > printForwardThreshold) {
        startPrintForward = std::min(i, startPrintForward);
        endPrintForward = std::max(i, endPrintForward);
      }
    }

    writeU16(out, startPrintForward);

    for (unsigned short int i = startPrintForward; i <= endPrintForward; i++) {
      if (it.forward_profile[i] > 2 * printForwardThreshold) {
        float_to_16_bit(it.forward_profile[i], forwardProbability);
        writeU16(out, forwardProbability);
      }
      else {
        float_to_16_bit(2 * printForwardThreshold, forwardProbability);

        writeU16(out, forwardProbability);
      }
    }

    writeU16(out, delimiter_16_bit);

    unsigned short int backwardProbability;
    float printBackwardThreshold;
    float_to_16_bit(0.0, backwardProbability);
    bit_16_to_float(backwardProbability, printBackwardThreshold);

    unsigned short int startPrintBackward = q->L;
    unsigned short int endPrintBackward = 0;

    for (unsigned short int i = 1; i <= q->L; i++) {
      float_to_16_bit(it.backward_profile[i], backwardProbability);
      if (it.backward_profile[i] > printBackwardThreshold) {
        startPrintBackward = std::min(i, startPrintBackward);
        endPrintBackward = std::max(i, endPrintBackward);
      }
    }

    writeU16(out, startPrintBackward);

    for (unsigned short int i = startPrintBackward; i <= endPrintBackward;
        i++) {
      if (it.backward_profile[i] > 2 * printBackwardThreshold) {
        float_to_16_bit(it.backward_profile[i], backwardProbability);
        writeU16(out, backwardProbability);
      }
      else {
        float_to_16_bit(2 * printBackwardThreshold, backwardProbability);
        writeU16(out, backwardProbability);
      }
    }

    writeU16(out, delimiter_16_bit);

    unsigned char posteriorProbability;

    int last_i = -1;
    int last_j = -1;
    for (size_t posterior_index = 0;
        posterior_index < it.posterior_probabilities.size();
        posterior_index++) {
      Posterior_Triple* triple = it.posterior_probabilities[posterior_index];

      float_to_8_bit(triple->posterior_probability, posteriorProbability);

      if (last_i != triple->query_pos || last_j + 1 != triple->template_pos) {
        if (last_i != -1 && last_j != -1) {
          out.write(reinterpret_cast<const char*>(&delimiter_8_bit),
              sizeof(delimiter_8_bit));
        }
        writeU16(out, triple->query_pos);
        writeU16(out, triple->template_pos);
      }

      out.write(reinterpret_cast<const char*>(&posteriorProbability),
          sizeof(unsigned char));

      last_i = triple->query_pos;
      last_j = triple->template_pos;
    }
    out.write(reinterpret_cast<const char*>(&delimiter_8_bit),
        sizeof(delimiter_8_bit));
    writeU16(out, delimiter_16_bit);
  }
  out.flush();
}
