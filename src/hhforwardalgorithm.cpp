/*
 * hhforwardalgorithmscalar.C
 *
 *  Created on: Apr 16, 2014
 *      Author: stefan
 */

#include "hhposteriordecoder.h"

void PosteriorDecoder::forwardAlgorithm(HMM & q, HMM & t, Hit & hit,
		PosteriorMatrix & p_mm, ViterbiMatrix & celloff_matrix,
		float shift, const int elem) {

	int i, j, j_min;

	float Pmax_i;                        // maximum of F_MM in row i

	float F_MM_prev[t.L + 1];
	float F_GD_prev[t.L + 1];
	float F_DG_prev[t.L + 1];
	float F_IM_prev[t.L + 1];
	float F_MI_prev[t.L + 1];

	float F_MM_curr[t.L + 1];
	float F_GD_curr[t.L + 1];
	float F_DG_curr[t.L + 1];
	float F_IM_curr[t.L + 1];
	float F_MI_curr[t.L + 1];

	//	fprintf(stdout, "%s\n", t.name);
	//	char * nam = "PF07714";
	//	bool eq = false;
	//	if (!strcmp(t.name, nam)) {
	//  	eq = true;
	//	}

	// Initialize F_XX_prev (representing i=1) andhit.P_MM[1][j]
	F_MM_prev[0] = F_IM_prev[0] = F_GD_prev[0] = F_MI_prev[0] = F_DG_prev[0] = -FLT_MAX;


	for (j=1; j <= t.L; ++j) {
		if (celloff_matrix.getCellOff(1, j, elem)) {
			F_MM_prev[j] = F_MI_prev[j] = F_DG_prev[j] = F_IM_prev[j] = F_GD_prev[j] = -FLT_MAX;
			//			if (eq) {
			//				fprintf(stdout, "1-C,%2.20f,%2.20f,%2.20f,%2.20f,%2.20f\n", F_MM_prev[j], F_GD_prev[j], F_IM_prev[j], F_DG_prev[j], F_MI_prev[j]);
			//			}
		} else {
			//			F_MM_prev[j] = flog2(hit.ProbFwd(q.p[1],t.p[j])) + hit.ScoreSS(&q,&t,1,j) + par.shift;
			//			F_MM_prev[j] = flog2(hit.ProbFwd(q.p[1],t.p[j])) + par.shift;
			F_MM_prev[j] = log2f(ScalarProd20(q.p[1],t.p[j])) + shift;
			F_MI_prev[j] = F_DG_prev[j] = -FLT_MAX;
			F_IM_prev[j] = flog2_sum_fpow2(
					F_MM_prev[j-1] + q.tr[1][M2I] + t.tr[j-1][M2M],	// MM -> IM
					F_IM_prev[j-1] + q.tr[1][I2I] + t.tr[j-1][M2M]	// IM -> IM
			);
			F_GD_prev[j] = flog2_sum_fpow2(
					F_MM_prev[j-1] + t.tr[j-1][M2D],	// DG -> MM
					F_GD_prev[j-1] + t.tr[j-1][D2D]	// GD -> GD
			);
			//			if (eq) {
			//				fprintf(stdout, "1-NC,%2.20f,%2.20f,%2.20f,%2.20f,%2.20f\n", F_MM_prev[j], F_GD_prev[j], F_IM_prev[j], F_DG_prev[j], F_MI_prev[j]);
			//			}
		}
		//		if (eq && celloff_matrix.getCellOff(1, j, 2) == 0) {
		//		if (eq) {
		//			fprintf(stdout, "1,co:%i,%2.20f,%2.20f,%2.20f,%2.20f,%2.20f\n", celloff_matrix.getCellOff(1, j, elem), F_MM_prev[j], F_GD_prev[j], F_IM_prev[j], F_DG_prev[j], F_MI_prev[j]);
		//		}
		p_mm.setSingleValue(1, j, elem, F_MM_prev[j]);
		//		if (eq) {
		//			fprintf(stdout, "1,%20.20f\n", p_mm.getSingleValue(1, j, elem));
		//		}
	}

	// Forward algorithm

	// Loop through query positions i
	for (i = 2; i <= q.L; ++i) {
		if (hit.self)
			j_min = imin(i + SELFEXCL + 1, t.L);
		else
			j_min = 1;

		// Initialize cells at (i,0)
		if (celloff_matrix.getCellOff(i, j_min, elem)) {
			F_MM_curr[j_min] = F_MI_curr[j_min] = F_DG_curr[j_min] = F_IM_curr[j_min] = F_GD_curr[j_min] = -FLT_MAX;
			//			if (eq) {
			//				fprintf(stdout, "2-C,%2.20f,%2.20f,%2.20f,%2.20f,%2.20f\n", F_MM_curr[j_min], F_GD_curr[j_min], F_IM_curr[j_min], F_DG_curr[j_min], F_MI_curr[j_min]);
			//			}
		} else {
			//			F_MM_curr[j_min] = flog2(hit.ProbFwd(q.p[i],t.p[j_min])) + hit.ScoreSS(&q,&t,i,j_min) + par.shift;
			//			F_MM_curr[j_min] = flog2(hit.ProbFwd(q.p[i],t.p[j_min])) + par.shift;
			F_MM_curr[j_min] = log2f(ScalarProd20(q.p[i],t.p[j_min])) + shift;
			F_IM_curr[j_min] = F_GD_curr[j_min] = -FLT_MAX;
			F_MI_curr[j_min] = flog2_sum_fpow2(
					F_MM_prev[j_min] + q.tr[i-1][M2M] + t.tr[j_min][M2I], // MI -> MM
					F_MI_prev[j_min] + q.tr[i-1][M2M] + t.tr[j_min][I2I]  // MI -> MI
			);
			F_DG_curr[j_min] = flog2_sum_fpow2(
					F_MM_prev[j_min] + q.tr[i-1][M2D], // DG -> MM
					F_DG_prev[j_min] + q.tr[i-1][D2D]	// DG -> DG
			);
		}
		//		if (eq && celloff_matrix.getCellOff(i, j_min, 2) == 0) {
		//		if (eq) {
		//			fprintf(stdout, "2,co: %i,%2.20f,%2.20f,%2.20f,%2.20f,%2.20f\n", celloff_matrix.getCellOff(i, j_min, elem), F_MM_curr[j_min], F_GD_curr[j_min], F_IM_curr[j_min], F_DG_curr[j_min], F_MI_curr[j_min]);
		//		}

		/* copy back */
		p_mm.setSingleValue(i, j_min, elem, F_MM_curr[j_min]);

		//		if (eq) {
		//			fprintf(stdout, "2,%20.20f\n", p_mm.getSingleValue(i, j_min, elem));
		//		}

		Pmax_i = -FLT_MAX;

		// Loop through template positions j
		for (j = j_min + 1; j <= t.L; ++j) {

			// Recursion relations
			if (celloff_matrix.getCellOff(i, j, elem)) {
				F_MM_curr[j] = F_MI_curr[j] = F_DG_curr[j] = F_IM_curr[j] = F_GD_curr[j] = -FLT_MAX;
				//				if (eq) {
				//					fprintf(stdout, "3-C,%2.20f,%2.20f,%2.20f,%2.20f,%2.20f\n", F_MM_curr[j], F_GD_curr[j], F_IM_curr[j], F_DG_curr[j], F_MI_curr[j]);
				//				}
			} else {
				//				F_MM_curr[j] = flog2(hit.ProbFwd(q.p[i],t.p[j])) + hit.ScoreSS(&q,&t,i,j) + par.shift +
				//				F_MM_curr[j] = flog2(hit.ProbFwd(q.p[i],t.p[j])) + par.shift +
				F_MM_curr[j] = log2f(ScalarProd20(q.p[i],t.p[j])) + shift +
						flog2_sum_fpow2(
								m_p_min_scalar,
								F_MM_prev[j-1] + q.tr[i-1][M2M] + t.tr[j-1][M2M],	// BB -> MM (BB = Begin/Begin, for local alignment)
								F_GD_prev[j-1] + q.tr[i-1][M2M] + t.tr[j-1][D2M],	// GD -> MM
								F_IM_prev[j-1] + q.tr[i-1][I2M] + t.tr[j-1][M2M],	// IM -> MM
								F_DG_prev[j-1] + q.tr[i-1][D2M] + t.tr[j-1][M2M],	// DG -> MM
								F_MI_prev[j-1] + q.tr[i-1][M2M] + t.tr[j-1][I2M]	// MI -> MM
						);
				F_GD_curr[j] = flog2_sum_fpow2(
						F_MM_curr[j-1] + t.tr[j-1][M2D],										// GD -> MM
						F_GD_curr[j-1] + t.tr[j-1][D2D]										// GD -> GD
						//												F_DG_curr[j-1] + t.tr[j-1][M2D] + q.tr[i][D2M] 	  // DG -> GD (only when structure scores given)
				);
				F_IM_curr[j] = flog2_sum_fpow2(
						F_MM_curr[j-1] + q.tr[i][M2I] + t.tr[j-1][M2M],	// MM -> IM
						F_IM_curr[j-1] + q.tr[i][I2I] + t.tr[j-1][M2M]	// IM -> IM
						//                        F_MI_curr[j-1] + q.tr[i][M2I] + t.tr[j-1][I2M]  // MI -> IM (only when structure scores given)
				);
				F_DG_curr[j] = flog2_sum_fpow2(
						F_MM_prev[j] + q.tr[i-1][M2D],										// DG -> MM
						F_DG_prev[j] + q.tr[i-1][D2D]										// DG -> DG
				);
				F_MI_curr[j] = flog2_sum_fpow2(
						F_MM_prev[j] + q.tr[i-1][M2M] + t.tr[j][M2I],		// MI -> MM
						F_MI_prev[j] + q.tr[i-1][M2M] + t.tr[j][I2I]		// MI -> MI
				);
				//				fprintf(stdout, "%i,%2.20f,%2.20f,%2.20f,%2.20f,%2.20,%2.20f,%2.20f,%2.20f,%2.20f,%2.20f\n", counter, q.tr[i-1][M2M], q.tr[i-1][I2M], q.tr[i-1][D2M], q.tr[i-1][M2D], q.tr[i-1][D2D], q.tr[i-1][M2M], t.tr[i-1][I2M], t.tr[i-1][D2M], t.tr[i-1][M2D], t.tr[i-1][D2D]);
				//				if (eq) {
				//					fprintf(stdout, "3-NC,%2.20f,%2.20f,%2.20f,%2.20f,%2.20f\n", F_MM_curr[j], F_GD_curr[j], F_IM_curr[j], F_DG_curr[j], F_MI_curr[j]);
				//				}
				if(F_MM_curr[j] > Pmax_i)
					Pmax_i = F_MM_curr[j];

			} // end else
			//			if (eq && celloff_matrix.getCellOff(i, j, 2) == 0) {
			//			if (eq) {
			//				fprintf(stdout, "3,co: %i,%2.20f,%2.20f,%2.20f,%2.20f,%2.20f\n", celloff_matrix.getCellOff(i, j, elem), F_MM_curr[j], F_GD_curr[j], F_IM_curr[j], F_DG_curr[j], F_MI_curr[j]);
			//			}

		} //end for j

		/* F_MM_prev = F_MM_curr */
		for (int jj = 0; jj <= t.L; jj++) {
			F_MM_prev[jj] = F_MM_curr[jj];
			F_MI_prev[jj] = F_MI_curr[jj];
			F_IM_prev[jj] = F_IM_curr[jj];
			F_DG_prev[jj] = F_DG_curr[jj];
			F_GD_prev[jj] = F_GD_curr[jj];

			/* and fill matrix because it is reused elsewhere */
			p_mm.setSingleValue(i, jj, elem, F_MM_curr[jj]);
			//			if (eq) {
			//				fprintf(stdout, "3,%20.20f\n", p_mm.getSingleValue(i, jj, elem));
			//			}
		}

	} // end for i

	// Calculate the sum of the log-probabilities
	if (m_local) { // local alignment
		hit.Pforward = 0.0; // alignment contains no residues (see Mueckstein, Stadler et al.)
		// Loop through query positions i
		for (i = 1; i <= q.L; ++i) {
			if (hit.self)
				j_min = imin(i + SELFEXCL + 1, t.L);
			else
				j_min = 1;

			// Loop through template positions j
			for (j = j_min; j <= t.L; ++j) {
				hit.Pforward = flog2_sum_fpow2(hit.Pforward, p_mm.getSingleValue(i, j, elem)); //hit.P_MM already contains log-values
				//				fprintf(stdout, "%2.20f,",hit.P_MM[i][j]);
			}
			//			fprintf(stdout, "\n");
		}
	} else { // global alignment
		hit.Pforward = -FLT_MAX;
		for (i = 1; i < q.L; ++i)
			hit.Pforward = flog2_sum_fpow2(hit.Pforward, p_mm.getSingleValue(i, t.L, elem));

		for (j = 1; j <= t.L; ++j)
			hit.Pforward = flog2_sum_fpow2(hit.Pforward, p_mm.getSingleValue(q.L, j, elem));
	}

	// Calculate log2(P_forward)
	//  score = log2(hit.Pforward)-10.0f;
	//  score = Pforward - 10.0f;
	hit.score = hit.Pforward;

	if (m_local) {
		if (hit.self)
			hit.score -= log(0.5 * t.L * q.L) / LAMDA + 14.; // +14.0 to get approx same mean as for -global
		else
			hit.score -= log(t.L * q.L) / LAMDA + 14.; // +14.0 to get approx same mean as for -global
	}

	// Debugging output
	/*
     if (v>=4)
     {
     const int i0=0, i1=q.L;
     const int j0=0, j1=t.L;
     scale_prod=1;
     fprintf(stderr,"\nFwd      scale     ");
     for (j=j0; j<=j1; ++j) fprintf(stderr,"%3i     ",j);
     fprintf(stderr,"\n");
     for (i=i0; i<=i1; ++i)
     {
     scale_prod *= scale[i];
     fprintf(stderr,"%3i: %9.3G ",i,1/scale_prod);
     for (j=j0; j<=j1; ++j)
     fprintf(stderr,"%7.4f ",(F_MM[i][j]+F_MI[i][j]+F_IM[i][j]+F_DG[i][j]+F_GD[i][j]));
     fprintf(stderr,"\n");
     // 	  printf(" MM  %9.5f ",1/scale[i]);
     // 	  for (j=j0; j<=j1; ++j)
     // 	    printf("%7.4f ",F_MM[i][j]);
     // 	  printf("\n");
     }
     fprintf(stderr,"Template=%-12.12s  score=%6.3f i2=%i  j2=%i \n",t.name,score,i2,j2);
     fprintf(stderr,"\nForward total probability ratio: %8.3G\n",Pforward);
     }
     */
	return;

}


