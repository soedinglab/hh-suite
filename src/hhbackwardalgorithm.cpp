/*
 * hhbackwardalgorithmscalar.C
 *
 *  Created on: Apr 10, 2014
 *      Author: stefan
 */

#include "hhposteriordecoder.h"

void PosteriorDecoder::backwardAlgorithm(HMM & q, HMM & t, Hit & hit,
		PosteriorMatrix & p_mm, ViterbiMatrix & celloff_matrix, float shift, const int elem) {

	// Variable declarations
	int i, j;      // query and template match state indices
	// this is the scaled 1 in the SW algorithm that represents a starting alignment
	int jmin;

	//	char * nam = "PF07714";
	//	char * nam = "PF01275";
	//	bool eq = false;
	//	if (!strcmp(t.name, nam)) {
	//		eq = true;
	//	}

	float B_MM_prev[t.L + 1];
	float B_DG_prev[t.L + 1];
	float B_MI_prev[t.L + 1];

	float B_MM_curr[t.L + 1];
	float B_GD_curr[t.L + 1];
	float B_DG_curr[t.L + 1];
	float B_IM_curr[t.L + 1];
	float B_MI_curr[t.L + 1];

	// Initialization of top row, i.e. cells (0,j)
	for (int j = t.L; j >= 1; j--) {
		if (celloff_matrix.getCellOff(q.L, j, elem)) {
			p_mm.setSingleValue(q.L, j, elem, -FLT_MAX);
			B_MM_prev[j] = -FLT_MAX;
		} else {
			B_MM_prev[j] = 0.0;
			p_mm.setSingleValue(q.L, j, elem,
					p_mm.getSingleValue(q.L, j, elem) - hit.Pforward);
		}
		B_MI_prev[j] = B_DG_prev[j] = -FLT_MAX;
		//		if (eq) {
		//			fprintf(stdout, "1,co:%i,%2.20f,%2.20f,%2.20f\n", celloff_matrix.getCellOff(q.L, j, elem), B_MM_prev[j], B_DG_prev[j], B_MI_prev[j]);
		//		}
	}

	//backward probability calculation
	// Backward algorithm
	// Loop through query positions i
	for (i = q.L - 1; i >= 1; i--) {
		//       if (v>=5) printf("\n");
		if (hit.self)
			jmin = imin(i + SELFEXCL, t.L);
		else
			jmin = 1; // jmin = i+SELFEXCL and not (i+SELFEXCL+1) to set matrix element at boundary to zero

		// Initialize cells at (i,t.L+1)
		//		scale_prod *= scale[i + 1];
		//		if (scale_prod < DBL_MIN * 100)
		//			scale_prod = 0.0;

		if (celloff_matrix.getCellOff(i, t.L, elem)) {
			p_mm.setSingleValue(i, t.L, elem, -FLT_MAX);
			B_MM_curr[t.L] = -FLT_MAX;
		} else {
			B_MM_curr[t.L] = 0.0;
			p_mm.setSingleValue(i, t.L, elem,
					p_mm.getSingleValue(i, t.L, elem) - hit.Pforward);
		}



		B_IM_curr[t.L] = B_MI_curr[t.L] = B_DG_curr[t.L] = B_GD_curr[t.L] =
				-FLT_MAX;

		//		if (eq) {
		//			fprintf(stdout, "p_mm: %2.20f\n", p_mm.getSingleValue(i, t.L, elem));
		////				fprintf(stdout, "2,co: %i,%2.20f,%2.20f,%2.20f,%2.20f,%2.20f\n", celloff_matrix.getCellOff(i, j+1, elem), B_MM_curr[j+1], B_GD_curr[j+1], B_IM_curr[j+1], B_DG_curr[j+1], B_MI_curr[j+1]);
		//		}

		// Loop through template positions j
		for (j = t.L - 1; j >= jmin; j--) {
			// Recursion relations
			//	      printf("S[%i][%i]=%4.1f  ",i,j,Score(q.p[i],t.p[j]));

			//			if (eq) {
			//				fprintf(stdout, "p_mm: %2.20f\n", p_mm.getSingleValue(i, j+1, elem));
			////				fprintf(stdout, "2,co: %i,%2.20f,%2.20f,%2.20f,%2.20f,%2.20f\n", celloff_matrix.getCellOff(i, j+1, elem), B_MM_curr[j+1], B_GD_curr[j+1], B_IM_curr[j+1], B_DG_curr[j+1], B_MI_curr[j+1]);
			//			}

			if (celloff_matrix.getCellOff(i, j, elem)) {
				B_MM_curr[j] = B_GD_curr[j] = B_IM_curr[j] = B_DG_curr[j] =
						B_MI_curr[j] = -FLT_MAX;
			} else {
				//                float pmatch = B_MM_prev[j+1] + flog2(ProbFwd(q.p[i+1],t.p[j+1])) + ScoreSS(q,t,i+1,j+1);
				//				float pmatch = B_MM_prev[j + 1]
				//						+ log2f(hit.ProbFwd(q.p[i + 1], t.p[j + 1]))
				////						+ hit.ScoreSS(&q, &t, i + 1, j + 1)
				//						+ par.shift;
				float pmatch = B_MM_prev[j + 1]
						//						+ flog2(hit.ProbFwd(q.p[i + 1], t.p[j + 1]))
						+ log2f(ScalarProd20(q.p[i + 1], t.p[j + 1]))
						+ shift;
				B_MM_curr[j] = flog2_sum_fpow2(
						m_p_min_scalar, // MM -> EE (End/End, for local alignment)
						pmatch           + q.tr[i][M2M] + t.tr[j][M2M], // MM -> MM
						B_GD_curr[j + 1]                + t.tr[j][M2D], // MM -> GD (q.tr[i][M2M] is already contained in GD->MM)
						B_IM_curr[j + 1] + q.tr[i][M2I] + t.tr[j][M2M], // MM -> IM
						B_DG_prev[j]     + q.tr[i][M2D], // MM -> DG (t.tr[j][M2M] is already contained in DG->MM)
						B_MI_prev[j]     + q.tr[i][M2M] + t.tr[j][M2I]  // MM -> MI
				);
				B_GD_curr[j] = flog2_sum_fpow2(
						pmatch   				+ q.tr[i][M2M] + t.tr[j][D2M], // GD -> MM
						B_GD_curr[j + 1] 							 + t.tr[j][D2D]  // DG -> DG
				);
				B_IM_curr[j] = flog2_sum_fpow2(
						pmatch 					 + q.tr[i][I2M] + t.tr[j][M2M], // IM -> MM
						B_IM_curr[j + 1] + q.tr[i][I2I] + t.tr[j][M2M]	 // IM -> IM
				);
				B_DG_curr[j] = flog2_sum_fpow2(
						pmatch 			 + q.tr[i][D2M] + t.tr[j][M2M], // DG -> MM
						B_DG_prev[j] + q.tr[i][D2D]					        // DG -> DG
						//					B_GD[i][j+1] + q.tr[i][D2M] + t.tr[j][M2D])   // DG -> GD
				);
				B_MI_curr[j] = flog2_sum_fpow2(
						pmatch 			 + q.tr[i][M2M] + t.tr[j][I2M],	// MI -> MM
						B_MI_prev[j] + q.tr[i][M2M] + t.tr[j][I2I] // MI -> MI
						//			    B_IM[i][j+1] + q.tr[i][M2I] + t.tr[j][I2M])   // MI -> IM
				);

			} // end else

			/* Copy back to matrix */
			p_mm.setSingleValue(i, j, elem,
					p_mm.getSingleValue(i, j, elem) + B_MM_curr[j] - hit.Pforward);

			//			if (eq) {
			//				fprintf(stdout, "p_mm [i,j]: %i,%i,%2.20f\n", i, j, p_mm.getSingleValue(i, j, elem));
			//			}
		} //end for j

		for (int jj = 0; jj <= t.L; jj++) {
			B_MM_prev[jj] = B_MM_curr[jj];
			B_DG_prev[jj] = B_DG_curr[jj];
			B_MI_prev[jj] = B_MI_curr[jj];
		}
	} // end for i

	return;

}

