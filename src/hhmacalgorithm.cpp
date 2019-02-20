/*
 * hhmacalgorithmscalar.C
 *
 *  Created on: Apr 16, 2014
 *      Author: stefan
 */

#include "hhposteriordecoder.h"

#define CALCULATE_MAX4(max, var1, var2, var3, var4, varb) \
if (var1>var2) { max=var1; varb=ViterbiMatrix::STOP;}     \
else           { max=var2; varb=ViterbiMatrix::MM;};      \
if (var3>max)  { max=var3; varb=ViterbiMatrix::MI;};      \
if (var4>max)  { max=var4; varb=ViterbiMatrix::IM;};


void PosteriorDecoder::macAlgorithm(HMM & q, HMM & t, Hit & hit,
        PosteriorMatrix & p_mm, ViterbiMatrix & viterbi_matrix, float par_mact,  const int elem) {

    // Use Forward and Backward matrices to find that alignment which
    // maximizes the expected number of correctly aligned pairs of residues (mact=0)
    // or, more generally, which maximizes the expectation value of the number of
    // correctly aligned pairs minus (mact x number of aligned pairs)
    // "Correctly aligned" can be based on posterior probabilities calculated with
    // a local or a global version of the Forward-Backward algorithm.

    int i,j;           // query and template match state indices
    int jmin,jmax;     // range of dynamic programming for j
    float S_prev[t.L + 1];    // scores
    float S_curr[t.L + 1];    // scores
    float score_MAC;   // score of the best MAC alignment

    //	fprintf(stdout, "%s\n", t.name);
    //	char * nam = "PF07714";
    //	bool eq = false;
    //	if (!strcmp(t.name, nam)) {
    //		eq = true;
    //	}

    float term1, term2, term3, term4 = 0.0f;

    // Initialization of top row, i.e. cells (0,j)
    for (j=0; j <= t.L; ++j)
        S_prev[j] = 0.0;

    score_MAC = -FLT_MAX;
    hit.i2 = hit.j2 = 0;
    //	hit.bMM[0][0] = STOP;	//
    viterbi_matrix.setMatMat(0, 0, elem, ViterbiMatrix::STOP);
    //	char * c_ptr = new char;
    char val;
    hit.min_overlap = 0;
    // Dynamic programming
    for (i = 1; i <= q.L; ++i) { // Loop through query positions i
        // If q is compared to t, exclude regions where overlap of q with t < min_overlap residues
        jmin = imax( 1, i + hit.min_overlap - q.L);  // Lq-i+j>=Ovlap => j>=i+Ovlap-Lq => jmin=max{1, i+Ovlap-Lq}
        jmax = imin(t.L, i - hit.min_overlap + t.L);  // Lt-j+i>=Ovlap => j<=i-Ovlap+Lt => jmax=min{Lt,i-Ovlap+Lt}

        // Initialize cells
        S_curr[jmin-1] = 0.0;
        if (jmax < t.L)
            S_prev[jmax] = 0.0; // initialize at (i-1,jmax) if upper right triagonal is excluded due to min overlap

        //		for (j = jmin; j <= jmax; ++j) { // Loop through template positions j
        for (j = 1; j <= t.L; ++j) { // Loop through template positions j
            //			hit.bMM[i][j] = 0x00;	//
            viterbi_matrix.setMatMat(i, j, elem, 0x00);
            //			printf("hit.bMM[%i][%i]: %i; m: %i\n", i, j, hit.bMM[i][j], viterbi_matrix.getMatMat(i, j, elem));
            //			if (hit.bMM[i][j] != viterbi_matrix.getMatMat(i, j, elem)) {printf("VALUE WROOOOOONG!!!"); return; }
            if (viterbi_matrix.getCellOff(i, j, elem)) {
                S_curr[j] = -FLT_MIN;
                viterbi_matrix.setMatMat(i, j, elem, ViterbiMatrix::STOP);	// hit.bMM[i][j] = STOP;
                //	      if (i>135 && i<140)
                // 		printf("Cell off   i=%i  j=%i b:%i\n",i,j,bMM[i][j]);
            } else {
                // Recursion
                //				unsigned char c = viterbi_matrix.getMatMat(i, j, elem);
                // NOT the state before the first MM state)
                //				term1 = fpow2(p_mm.getPosteriorValue(i, j, elem)) - par.mact;
                term1 = p_mm.getPosteriorValue(i, j) - par_mact;
                //				term2 = S_prev[j-1] + fpow2(p_mm.getPosteriorValue(i, j, elem)) - par.mact;
                term2 = S_prev[j-1] + p_mm.getPosteriorValue(i, j) - par_mact;
                term3 = S_prev[j] - 0.5 * par_mact;
                term4 = S_curr[j-1] - 0.5 * par_mact;

                // check p_mm values -> OK
                //				if (eq) {
                //					fprintf(stdout, "co: %i,%2.20f\n", viterbi_matrix.getCellOff(i, j, elem), p_mm.getPosteriorValue(i, j, elem));
                //				}

                CALCULATE_MAX4(
                        S_curr[j],
                        term1,
                        term2,
                        term3,
                        term4,
                        val	//  hit.bMM[i][j]   backtracing matrix
                );
                //				CALCULATE_MAX4(
                //								S_curr[j],
                ////											fpow2(hit.P_MM[i][j]) - par.mact,  // STOP signifies the first MM state, NOT the state before the first MM state (as in Viterbi)
                ////								powf(2, p_mm.getPosteriorValue(i, j, elem)) - par.mact,  // STOP signifies the first MM state, NOT the state before the first MM state (as in Viterbi)
                //								fpow2(p_mm.getPosteriorValue(i, j, elem)) - par.mact,  // STOP signifies the first MM state, NOT the state before the first MM state (as in Viterbi)
                ////											S_prev[j-1] + fpow2(hit.P_MM[i][j]) - par.mact, // hit.P_MM[i][j] contains log2-posterior probability
                ////								S_prev[j-1] + powf(2, p_mm.getPosteriorValue(i, j, elem)) - par.mact, // hit.P_MM[i][j] contains log2-posterior probability
                //								S_prev[j-1] + fpow2(p_mm.getPosteriorValue(i, j, elem)) - par.mact, // hit.P_MM[i][j] contains log2-posterior probability
                //								S_prev[j] - 0.5 * par.mact,  // gap penalty prevents alignments such as this: XX--xxXX
                //								S_curr[j-1] - 0.5 * par.mact,  //                                               YYyy--YY
                //								val	//  hit.bMM[i][j]   backtracing matrix
                //				);

                //				if (eq) {
                // check terms ->
                //					fprintf(stdout, "%2.20f,%2.20f,%2.20f,%2.20f\n", term1, term2, term3, term4);
                //					fprintf(stdout, "%2.20f\n", S_curr[j]);
                //				}
                viterbi_matrix.setMatMat(i, j, elem, (unsigned char) val);
                //				viterbi_matrix.setMatMat(i, j, elem, ViterbiMatrix::MM);
                //				printf("hit.bMM[%i][%i]: %i; m: %i\n", i, j, hit.bMM[i][j], viterbi_matrix.getMatMat(i, j, elem));
                //				printf("value: %i\n", viterbi_matrix.getMatMat(i, j, elem));


                //if (i>36 && i<40 && j>2200 && j<2230)
                //printf("i=%i  j=%i  S[i][j]=%8.3f  MM:%7.3f  MI:%7.3f  IM:%7.3f  b:%i\n",i,j,S[i][j],S[i-1][j-1]+B_MM[i][j]-par.mact,S[i-1][j],S[i][j-1],bMM[i][j]);

                // Find maximum score; global alignment: maximize only over last row and last column
                if(S_curr[j] > score_MAC && (m_local || i == q.L)) {
                    hit.i2 = i; hit.j2 = j;
                    score_MAC = S_curr[j];
                }

            } // end if

            //			if (eq) {
            //				// check score_mac ->
            ////				fprintf(stdout, "%2.20f\n", score_MAC);
            //				fprintf(stdout, "state:%i\n", viterbi_matrix.getMatMat(i, j, elem));
            //			}

        } //end for j

        // if global alignment: look for best cell in last column
        if (!m_local && S_curr[jmax] > score_MAC) {
            hit.i2 = i;
            hit.j2 = jmax;
            score_MAC = S_curr[jmax];
        }

        for (j = 0; j <= t.L; ++j)
            S_prev[j] = S_curr[j];
    } // end for i

    /*
     // DEBUG
     if (v>=5)
     {
     printf("\nScore  ");
     for (j=0; j<=t.L; ++j) printf("%3i   ",j);
     printf("\n");
     for (i=0; i<=q.L; ++i)
     {
     printf("%2i:    ",i);
     for (j=0; j<=t.L; ++j)
     printf("%5.2f ",S[i][j]);
     printf("\n");
     }
     printf("\n");
     printf("Template=%-12.12s  i=%-4i j=%-4i score=%6.3f\n",t.name,hit.i2,hit.j2,score_MAC);
     }
     */

    //	if (eq) {
    //		// check i2, j2 ->
    //		fprintf(stdout, "[i,j] values: [%i,%i]\n", hit.i2, hit.j2);
    //	}

    return;

}

#undef CALCULATE_MAX4
