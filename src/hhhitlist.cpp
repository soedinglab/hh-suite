// hhhitlist.C

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
//// Methods of class HitList
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

#include "hhhitlist.h"


/////////////////////////////////////////////////////////////////////////////////////
// Print summary listing of hits
/////////////////////////////////////////////////////////////////////////////////////
void HitList::PrintHitList(HMM* q, std::stringstream& out,
    const unsigned int maxdbstrlen, const int z, const int Z, const float p,
    const double E, const int argc, const char** argv) {
  int nhits = 0;

  out << "Query         " << q->longname << std::endl;
  out << "Match_columns " << q->L << std::endl;
  out << "No_of_seqs    " << q->N_filtered << " out of " << q->N_in
      << std::endl;
  out << "Neff          " << q->Neff_HMM << std::endl;
  out << "Searched_HMMs " << N_searched << std::endl;

  // Print date stamp
  time_t* tp = new time_t;
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
    const double E, const int argc, const char** argv) {
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
    const int argc, const char** argv, const float S[20][20], const int maxseq) {
  std::stringstream out;
  PrintHHR(q, out, maxdbstrlen, showconf, showcons, showdssp, showpred, b, B, z, Z, aliwidth, nseqdis, p, E, argc, argv, S, maxseq);

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
    const float p, const double E, const int argc, const char** argv,
    const float S[20][20], const int maxseq) {
  PrintHitList(q, out, maxdbstrlen, z, Z, p, E, argc, argv);
  PrintAlignments(q, out, showconf, showcons, showdssp, showpred, p, aliwidth, nseqdis, b, B, E, S, maxseq, 0);
}

/////////////////////////////////////////////////////////////////////////////////////
// Print alignments of query sequences against hit sequences
/////////////////////////////////////////////////////////////////////////////////////
void HitList::PrintAlignments(HMM* q, char* outfile, const char showconf,
    const char showcons, const char showdssp, const char showpred,
    const float p, const int aliwidth, const int nseqdis, const int b,
    const int B, const double E, const float S[20][20], const int maxseq, char outformat) {
  std::stringstream out;
  PrintAlignments(q, out, showconf, showcons, showdssp, showpred, p, aliwidth, nseqdis, b, B, E, S, maxseq, outformat);

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
    const int b, const int B, const double E, const float S[20][20], const int maxseq,
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
      qt_ali.PrintHHR(out, hit, showconf, showcons, showdssp, showpred, aliwidth, maxseq);
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
      HH_LOG(WARNING) << "In " << __FILE__ << ":" << __LINE__ << ": " << __func__ << ":" << std::endl;
      HH_LOG(WARNING) << "\tCould not open \'" << outputfile << std::endl;
      return;
    }

    out << outbuffer.str();

    out.close();
  }
}


/////////////////////////////////////////////////////////////////////////////////////
// Print score distribution into a blast tab file
/////////////////////////////////////////////////////////////////////////////////////
void HitList::PrintM8File(HMM* q, char* outputfile) {
    std::stringstream outbuffer;
    PrintM8File(q, outbuffer);
    
    if (strcmp(outputfile, "stdout") == 0) {
        std::cout << outbuffer.str();
    }
    else {
        std::ofstream out(outputfile);
        if (!out.good()) {
            HH_LOG(WARNING) << "In " << __FILE__ << ":" << __LINE__ << ": " << __func__ << ":" << std::endl;
            HH_LOG(WARNING) << "\tCould not open \'" << outputfile << std::endl;
            return;
        }
        
        out << outbuffer.str();
        
        out.close();
    }
}


void HitList::PrintM8File(HMM* q, std::stringstream& outbuffer) {
    Hash<int> twice(10000); // make sure only one hit per HMM is listed
    twice.Null(-1);
    
    Reset();

    
    //Blast tab format
    // query target evalue score
    // d1c8da_	Q32Z53	0.618	547	334	0	1	548	1	542	3.11E-185	646
    
    int i = 0;
    while (!End()) {
        i++;
        Hit hit = ReadNext();
        char line[LINELEN];
        int gapOpenCount = 0;
        int missMatchCount = 0;
        int matchCount  = 0;
        bool isGapOpen = false;
        for (int step = hit.nsteps; step >= 1; step--){
            if (hit.states[step] == GD || hit.states[step] == DG) {
                if(isGapOpen == false){
                gapOpenCount++;
                }
                isGapOpen = true;
            }else if (hit.states[step] == MM){
                if(hit.seq[hit.nfirst][hit.j[step]] == q->seq[q->nfirst][hit.i[step]]){
                    matchCount++;
                }else{
                    missMatchCount++;
                }
                isGapOpen = false;
            }else{
                isGapOpen = false;
            }
        }
        sprintf(line, "%s\t%s\t%1.3f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.2E\t%.1f\n",
                q->name, hit.file, static_cast<float>(matchCount)/static_cast<float>(hit.L), hit.L, missMatchCount, gapOpenCount,
                hit.i1, hit.i2, hit.j1, hit.j2, hit.Eval, -hit.score_aass);
        outbuffer << line;
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

  HH_LOG(DEBUG)
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
    const bool filter_matrices, const size_t max_number_matrices, const float S[20][20]) {
  std::stringstream out_ss(
      std::stringstream::in | std::stringstream::out
          | std::stringstream::binary);
  PrintMatrices(q, out_ss, filter_matrices, max_number_matrices, S);

  if (strcmp(matricesOutputFileName, "stdout") == 0) {
    std::cout << out_ss.str();
  }
  else {
    std::ofstream out(matricesOutputFileName, std::ios::out | std::ios::binary);
    if (!out.good()) {
      HH_LOG(WARNING) << "Iin " << __FILE__ << ":" << __LINE__ << ": "
          << __func__ << ":" << std::endl;
      HH_LOG(WARNING) << "\tcould not open \'" << matricesOutputFileName << std::endl;
      return;
    }

    out << out_ss.str();

    out.close();
  }
}

void HitList::PrintMatrices(HMM* q, std::stringstream& out,
     const bool filter_matrices, const size_t max_number_matrices, const float S[20][20]) {
  //limit matrices to par.max_number_matrices
  std::vector<Hit> hits;
  int protein_max_length = 4000;

  if(q->L >= protein_max_length) {
    return;
  }

  //remove invalid alignments
  const float tolerance = 0.10;

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
        && backward_profile_sum > 1.0 - tolerance
        && hit_cur.forward_entries > 0 
        && hit_cur.backward_entries > 0 
        && hit_cur.posterior_entries > 0) {
      hits.push_back(hit_cur);
    }
  }

  std::vector<bool> picked_alignments(hits.size(), true);
  unsigned int chosen = hits.size();
  float matix_probability_threshold = 20;

  if (hits.size() == 0) {
    HH_LOG(WARNING) << "There are no alignment matrices to print!"
        << std::endl;
  }

  //remove duplicate alignments (mostly alignments which were realigned afterwards)
  for (int index1 = hits.size() - 1; index1 >= 0; index1--) {
    Hit it = hits[index1];
    
    if (it.Probab < matix_probability_threshold || it.L >= protein_max_length) {
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

  if(filter_matrices) {
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
  }

  HH_LOG(INFO) << "Printing alignment matrices..." << std::endl;
  HH_LOG(INFO) << "Total number of alignments    : " << hits.size()
      << std::endl;
  HH_LOG(INFO) << "Number of accepted alignments : " << chosen
      << std::endl;

  if (chosen == 0) {
    HH_LOG(WARNING)
        << "No homologs found for printing matrix!" << std::endl;
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

    int last_i;
    int last_j;
    unsigned char probability;

    //write backward matrix
    last_i = -1;
    last_j = -1;
    for (size_t backward_index = 0;
        backward_index < it.backward_entries;
        backward_index++) {

      float_to_8_bit(it.backward_matrix[backward_index][2], probability);

      if (last_i != it.backward_matrix[backward_index][0] || last_j + 1 != it.backward_matrix[backward_index][1]) {
        if (last_i != -1 && last_j != -1) {
          out.write(reinterpret_cast<const char*>(&delimiter_8_bit),
              sizeof(delimiter_8_bit));
        }
        writeU16(out, it.backward_matrix[backward_index][0]);
        writeU16(out, it.backward_matrix[backward_index][1]);
      }

      out.write(reinterpret_cast<const char*>(&probability),
          sizeof(unsigned char));

      last_i = it.backward_matrix[backward_index][0];
      last_j = it.backward_matrix[backward_index][1];
    }
    out.write(reinterpret_cast<const char*>(&delimiter_8_bit),
        sizeof(delimiter_8_bit));
    writeU16(out, delimiter_16_bit);

    //write forward matrix
    last_i = -1;
    last_j = -1;
    for (size_t forward_index = 0;
        forward_index < it.forward_entries;
        forward_index++) {

      float_to_8_bit(it.forward_matrix[forward_index][2], probability);

      if (last_i != it.forward_matrix[forward_index][0] || last_j + 1 != it.forward_matrix[forward_index][1]) {
        if (last_i != -1 && last_j != -1) {
          out.write(reinterpret_cast<const char*>(&delimiter_8_bit),
              sizeof(delimiter_8_bit));
        }
        writeU16(out, it.forward_matrix[forward_index][0]);
        writeU16(out, it.forward_matrix[forward_index][1]);
      }

      out.write(reinterpret_cast<const char*>(&probability),
          sizeof(unsigned char));

      last_i = it.forward_matrix[forward_index][0];
      last_j = it.forward_matrix[forward_index][1];
    }
    out.write(reinterpret_cast<const char*>(&delimiter_8_bit),
        sizeof(delimiter_8_bit));
    writeU16(out, delimiter_16_bit);



    //write posterior matrix
    last_i = -1;
    last_j = -1;
    for (size_t posterior_index = 0;
        posterior_index < it.posterior_entries;
        posterior_index++) {

      float_to_8_bit(it.posterior_matrix[posterior_index][2], probability);

      if (last_i != it.posterior_matrix[posterior_index][0] || last_j + 1 != it.posterior_matrix[posterior_index][1]) {
        if (last_i != -1 && last_j != -1) {
          out.write(reinterpret_cast<const char*>(&delimiter_8_bit),
              sizeof(delimiter_8_bit));
        }
        writeU16(out, it.posterior_matrix[posterior_index][0]);
        writeU16(out, it.posterior_matrix[posterior_index][1]);
      }

      out.write(reinterpret_cast<const char*>(&probability),
          sizeof(unsigned char));

      last_i = it.posterior_matrix[posterior_index][0];
      last_j = it.posterior_matrix[posterior_index][1];
    }
    out.write(reinterpret_cast<const char*>(&delimiter_8_bit),
        sizeof(delimiter_8_bit));
    writeU16(out, delimiter_16_bit);
  }
  out.flush();
}
