/*
 * hhforwardalgorithmscalar.C
 *
 *  Created on: Apr 16, 2014
 *      Author: stefan
 */

#include <Python/Python.h>
#include "hhposteriordecoder.h"

void PosteriorDecoder::forwardAlgorithm(HMM & q, HMM & t, Hit & hit,
		PosteriorMatrix & p_mm, ViterbiMatrix & celloff_matrix,
		float shift, const int elem) {
	int i, j;      // query and template match state indices
	double pmin = (m_local ? 1.0 : 0.0); // used to distinguish between SW and NW algorithms in maximization
	double Cshift = pow(2.0, shift); // score offset transformed into factor in lin-space
	double Pmax_i;                        // maximum of F_MM in row i
	double scale_prod = 1.0;                // Prod_i=1^i (scale[i])
	int jmin;

	double F_MM_prev[t.L + 1];
	double F_GD_prev[t.L + 1];
	double F_DG_prev[t.L + 1];
	double F_IM_prev[t.L + 1];
	double F_MI_prev[t.L + 1];

	double F_MM_curr[t.L + 1];
	double F_GD_curr[t.L + 1];
	double F_DG_curr[t.L + 1];
	double F_IM_curr[t.L + 1];
	double F_MI_curr[t.L + 1];

	//	fprintf(stdout, "%s\n", t.name);
	//	char * nam = "PF07714";
	//	bool eq = false;
	//	if (!strcmp(t.name, nam)) {
	//  	eq = true;
	//	}

	// Initialize F_XX_prev (representing i=1) andhit.P_MM[1][j]
	// Initialize F_XX_prev (representing i=1) and P_MM[1][j]
	for (j=1; j<=t.L; ++j)
		F_MM_curr[j] = 0.0;

	F_MM_curr[0] = 0.0;
	F_IM_curr[0] = 0.0;
	F_GD_curr[0] = 0.0;
	for (j=1; j<=t.L; ++j)
	{
		if (celloff_matrix.getCellOff(1,j,elem))
			F_MM_curr[j] = F_MI_curr[j] = F_DG_curr[j] = F_IM_curr[j] = F_GD_curr[j] = 0.0;
		else
		{
			F_MM_curr[j] = ProbFwd(q.p[1], t.p[j]) * Cshift ;
			F_MI_curr[j] = F_DG_curr[j] = 0.0;
			F_IM_curr[j] = F_MM_curr[j-1] * q.tr[1][M2I] * t.tr[j-1][M2M] + F_IM_curr[j-1] * q.tr[1][I2I] * t.tr[j-1][M2M];
			F_GD_curr[j] = F_MM_curr[j-1] * t.tr[j-1][M2D]                + F_GD_curr[j-1] * t.tr[j-1][D2D];
			//printf("%f %f %f %f %f\n",F_MM_curr[j], F_GD_curr[j], F_IM_curr[j], F_DG_curr[j], F_MI_curr[j]);
			//printf("%f %f %f %f %f\n",q.tr[1][M2I], t.tr[j-1][M2M], q.tr[1][I2I], t.tr[j-1][M2D] , t.tr[j-1][D2D]);
			//printf("%f %f %f\n",ProbFwd(q.p[1], t.p[j]), F_IM_curr[j-1], F_GD_curr[j-1]);
		}
	}

	for (int j = 0; j <= t.L; j++)
	{
		p_mm.setSingleValue(0,j,F_MM_prev[j]);
		p_mm.setSingleValue(1,j,F_MM_curr[j]);
		F_MM_prev[j] = F_MM_curr[j];
		F_MI_prev[j] = F_MI_curr[j];
		F_IM_prev[j] = F_IM_curr[j];
		F_DG_prev[j] = F_DG_curr[j];
		F_GD_prev[j] = F_GD_curr[j];
	}



	scale[0] = scale[1] = scale[2] = 1.0;

	// Forward algorithm

	// Loop through query positions i
	for (i = 2; i <= q.L; ++i) {
		if (hit.self)
			jmin = imin(i + SELFEXCL + 1, t.L);
		else
			jmin = 1;

		if (scale_prod < DBL_MIN * 100)
			scale_prod = 0.0;
		else
			scale_prod *= scale[i];

		// Initialize cells at (i,0)
		if (celloff_matrix.getCellOff(i, jmin, elem))
			F_MM_curr[jmin] = F_MI_curr[jmin] = F_DG_curr[jmin] = F_IM_curr[jmin] =
					F_GD_curr[jmin] = 0.0;
		else {
			F_MM_curr[jmin] = scale_prod * ProbFwd(q.p[i], t.p[jmin]) * Cshift;
			F_IM_curr[jmin] = F_GD_curr[jmin] = 0.0;
			F_MI_curr[jmin] = scale[i]
					* (F_MM_prev[jmin] * q.tr[i - 1][M2M] * t.tr[jmin][M2I]
					+ F_MI_prev[jmin] * q.tr[i - 1][M2M] * t.tr[jmin][I2I]);
			F_DG_curr[jmin] = scale[i]
					* (F_MM_prev[jmin] * q.tr[i - 1][M2D]
					+ F_DG_prev[jmin] * q.tr[i - 1][D2D]);
		}

		/* copy back */
		p_mm.setSingleValue(i, jmin, F_MM_curr[jmin]);

		Pmax_i = 0;

		// Loop through template positions j
		for (j = jmin + 1; j <= t.L; ++j) {

			// Recursion relations
			if (celloff_matrix.getCellOff(i, j, elem))
				F_MM_curr[j] = F_MI_curr[j] = F_DG_curr[j] = F_IM_curr[j] =
						F_GD_curr[j] = 0.0;
			else {
				F_MM_curr[j] = ProbFwd(q.p[i], t.p[j]) * Cshift
						* scale[i]
						* (pmin
						+ F_MM_prev[j - 1] * q.tr[i - 1][M2M] * t.tr[j - 1][M2M] // BB -> MM (BB = Begin/Begin, for local alignment)
						+ F_GD_prev[j - 1] * q.tr[i - 1][M2M] * t.tr[j - 1][D2M] // GD -> MM
						+ F_IM_prev[j - 1] * q.tr[i - 1][I2M] * t.tr[j - 1][M2M] // IM -> MM
						+ F_DG_prev[j - 1] * q.tr[i - 1][D2M] * t.tr[j - 1][M2M] // DG -> MM
						+ F_MI_prev[j - 1] * q.tr[i - 1][M2M] * t.tr[j - 1][I2M] // MI -> MM
				);
				F_GD_curr[j] = (
						F_MM_curr[j - 1] * t.tr[j - 1][M2D]         // GD -> MM
						+ F_GD_curr[j - 1] * t.tr[j - 1][D2D]                    // GD -> GD
				);
				F_IM_curr[j] = (F_MM_curr[j - 1] * q.tr[i][M2I] * t.tr[j - 1][M2M] // MM -> IM
						+ F_IM_curr[j - 1] * q.tr[i][I2I] * t.tr[j - 1][M2M]     // IM -> IM
				);
				F_DG_curr[j] = scale[i] * (F_MM_prev[j] * q.tr[i - 1][M2D]  // DG -> MM
						+ F_DG_prev[j] * q.tr[i - 1][D2D]                    // DG -> DG
				);
				F_MI_curr[j] = scale[i]
						* (F_MM_prev[j] * q.tr[i - 1][M2M] * t.tr[j][M2I]     // MI -> MM
						+ F_MI_prev[j] * q.tr[i - 1][M2M] * t.tr[j][I2I]     // MI -> MI
				);

				Pmax_i = fmax(Pmax_i, F_MM_curr[j]);
				//printf("%f %f %f %f %f\n",F_MM_prev[j], F_GD_prev[j], F_IM_prev[j], F_DG_prev[j], F_MI_prev[j]);
				//printf("%f %f %f %f %f\n",F_MM_curr[j], F_GD_curr[j], F_IM_curr[j], F_DG_curr[j], F_MI_curr[j]);
				//printf("%f %f %f %f %f\n",q.tr[i - 1][M2M],  q.tr[i - 1][I2M], q.tr[i - 1][D2M], t.tr[j - 1][M2M], t.tr[j - 1][D2M]);

			} // end else

		} //end for j

		/* F_MM_prev = F_MM_curr */
		for (int jj = 0; jj <= t.L; jj++) {
			F_MM_prev[jj] = F_MM_curr[jj];
			F_MI_prev[jj] = F_MI_curr[jj];
			F_IM_prev[jj] = F_IM_curr[jj];
			F_DG_prev[jj] = F_DG_curr[jj];
			F_GD_prev[jj] = F_GD_curr[jj];

			// Fill posterior probability matrix with forward score
			p_mm.setSingleValue(i, jj, F_MM_curr[jj]);
		}
		pmin *= scale[i];
		if (pmin < DBL_MIN * 100)
			pmin = 0.0;

	  scale[i + 1] = 1.0 / (Pmax_i + 1.0);
//    scale[i+1] = 1.0;   // to debug scaling

	} // end for i

// Calculate P_forward * Product_{i=1}^{Lq+1}(scale[i])
	if (m_local) {
		hit.Pforward  = 1.0; // alignment contains no residues (see Mueckstein, Stadler et al.)
		// Loop through query positions i
		for (i = 1; i <= q.L; ++i) {
			if (hit.self)
				jmin = imin(i + SELFEXCL + 1, t.L);
			else
				jmin = 1;

			for (j = jmin; j <= t.L; ++j) // Loop through template positions j
				hit.Pforward  += p_mm.getSingleValue(i, j);

			hit.Pforward *= scale[i + 1];
		}
	} else { // global alignment
		hit.Pforward  = 0.0;
		for (i = 1; i < q.L; ++i)
			hit.Pforward  = (hit.Pforward  + p_mm.getSingleValue(i, t.L) * scale[i + 1]);
		for (j = 1; j <= t.L; ++j)
			hit.Pforward  += p_mm.getSingleValue(q.L, j);
		hit.Pforward  *= scale[q.L + 1];
	}

	// Calculate log2(P_forward)
	hit.score = log2(hit.Pforward) - 10.0f;
	for (i = 1; i <= q.L + 1; ++i)
		hit.score -= log2(scale[i]);

	if (m_local) {
		if (hit.self)
			hit.score -= log(0.5 * t.L * q.L) / LAMDA + 14.; // +14.0 to get approx same mean as for -global
		else
			hit.score -= log(t.L * q.L) / LAMDA + 14.; // +14.0 to get approx same mean as for -global
	}
	return;

}


