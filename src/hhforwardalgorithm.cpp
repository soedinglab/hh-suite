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
	int i, j;      // query and template match state indices
	double pmin = (m_local ? 1.0 : 0.0); // used to distinguish between SW and NW algorithms in maximization
	double Cshift = pow(2.0, shift); // score offset transformed into factor in lin-space
	double Pmax_i;                        // maximum of F_MM in row i
	double scale_prod = 1.0;                // Prod_i=1^i (scale[i])
	int jmin;

	// Initialize F_XX_prev (representing i=1) andhit.P_MM[1]
	// Initialize F_XX_prev (representing i=1) and P_MM[1][j]
	for (j=1; j<=t.L; ++j)
		m_fwd[j].mm_curr = 0.0;

	m_fwd[0].mm_curr = 0.0;
	m_fwd[0].im_curr = 0.0;
	m_fwd[0].gd_curr = 0.0;
	for (j=1; j<=t.L; ++j)
	{
		if (celloff_matrix.getCellOff(1,j,elem))
			m_fwd[j].mm_curr = m_fwd[j].mi_curr = m_fwd[j].dg_curr = m_fwd[j].im_curr = m_fwd[j].gd_curr = 0.0;
		else
		{
			m_fwd[j].mm_curr = ProbFwd(q.p[1], t.p[j]) * Cshift ;
			m_fwd[j].mi_curr = m_fwd[j].dg_curr = 0.0;
			m_fwd[j].im_curr = m_fwd[j-1].mm_curr * q.tr[1][M2I] * t.tr[j-1][M2M] + m_fwd[j-1].im_curr * q.tr[1][I2I] * t.tr[j-1][M2M];
			m_fwd[j].gd_curr = m_fwd[j-1].mm_curr * t.tr[j-1][M2D]                + m_fwd[j-1].gd_curr * t.tr[j-1][D2D];
			//printf("%f %f %f %f %f\n",m_mm_curr[j], m_gd_curr[j], m_im_curr[j], m_dg_curr[j], m_mi_curr[j]);
			//printf("%f %f %f %f %f\n",q.tr[1][M2I], t.tr[j-1][M2M], q.tr[1][I2I], t.tr[j-1][M2D] , t.tr[j-1][D2D]);
			//printf("%f %f %f\n",ProbFwd(q.p[1], t.p[j]), m_im_curr[j-1], m_gd_curr[j-1]);
		}
	}

	for (int j = 0; j <= t.L; j++)
	{
		p_mm.setPosteriorValue(0,j, m_fwd[j].mm_prev);
		p_mm.setPosteriorValue(1,j, m_fwd[j].mm_curr);
		m_fwd[j].mm_prev = m_fwd[j].mm_curr;
		m_fwd[j].mi_prev = m_fwd[j].mi_curr;
		m_fwd[j].im_prev = m_fwd[j].im_curr;
		m_fwd[j].dg_prev = m_fwd[j].dg_curr;
		m_fwd[j].gd_prev = m_fwd[j].gd_curr;
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
			m_fwd[jmin].mm_curr = m_fwd[jmin].mi_curr = m_fwd[jmin].dg_curr = m_fwd[jmin].im_curr =
					m_fwd[jmin].gd_curr = 0.0;
		else {
			m_fwd[jmin].mm_curr = scale_prod * ProbFwd(q.p[i], t.p[jmin]) * Cshift;
			m_fwd[jmin].im_curr = m_fwd[jmin].gd_curr = 0.0;
			m_fwd[jmin].mi_curr = scale[i]
					* (m_fwd[jmin].mm_prev * q.tr[i - 1][M2M] * t.tr[jmin][M2I]
					+ m_fwd[jmin].mi_prev * q.tr[i - 1][M2M] * t.tr[jmin][I2I]);
			m_fwd[jmin].dg_curr = scale[i]
					* (m_fwd[jmin].mm_prev * q.tr[i - 1][M2D]
					+ m_fwd[jmin].dg_prev * q.tr[i - 1][D2D]);
		}

		/* copy back */
		p_mm.setPosteriorValue(i, jmin, m_fwd[jmin].mm_curr);

		Pmax_i = 0;

		// Loop through template positions j
		for (j = jmin + 1; j <= t.L; ++j) {

			// Recursion relations
			if (celloff_matrix.getCellOff(i, j, elem))
				m_fwd[j].mm_curr = m_fwd[j].mi_curr = m_fwd[j].dg_curr = m_fwd[j].im_curr =
						m_fwd[j].gd_curr = 0.0;
			else {
				m_fwd[j].mm_curr = ProbFwd(q.p[i], t.p[j]) * Cshift
						* scale[i]
						* (pmin
						+ m_fwd[j - 1].mm_prev * q.tr[i - 1][M2M] * t.tr[j - 1][M2M] // BB -> MM (BB = Begin/Begin, for local alignment)
						+ m_fwd[j - 1].gd_prev * q.tr[i - 1][M2M] * t.tr[j - 1][D2M] // GD -> MM
						+ m_fwd[j - 1].im_prev * q.tr[i - 1][I2M] * t.tr[j - 1][M2M] // IM -> MM
						+ m_fwd[j - 1].dg_prev * q.tr[i - 1][D2M] * t.tr[j - 1][M2M] // DG -> MM
						+ m_fwd[j - 1].mi_prev * q.tr[i - 1][M2M] * t.tr[j - 1][I2M] // MI -> MM
				);
				m_fwd[j].gd_curr = (
						m_fwd[j - 1].mm_curr * t.tr[j - 1][M2D] +         // GD -> MM
						m_fwd[j - 1].gd_curr * t.tr[j - 1][D2D]                    // GD -> GD
				);
				m_fwd[j].im_curr = (m_fwd[j - 1].mm_curr * q.tr[i][M2I] * t.tr[j - 1][M2M] // MM -> IM
						+ m_fwd[j - 1].im_curr * q.tr[i][I2I] * t.tr[j - 1][M2M]     // IM -> IM
				);
				m_fwd[j].dg_curr = scale[i] * (m_fwd[j].mm_prev * q.tr[i - 1][M2D]  // DG -> MM
						+ m_fwd[j].dg_prev * q.tr[i - 1][D2D]                    // DG -> DG
				);
				m_fwd[j].mi_curr = scale[i]
						* (m_fwd[j].mm_prev * q.tr[i - 1][M2M] * t.tr[j][M2I]   // MI -> MM
						+ m_fwd[j].mi_prev * q.tr[i - 1][M2M] * t.tr[j][I2I]     // MI -> MI
				);

				Pmax_i = fmax(Pmax_i, m_fwd[j].mm_curr);
				//printf("%f %f %f %f %f\n",F_MM_prev[j], F_GD_prev[j], m_im_prev[j], m_dg_prev[j], m_mi_prev[j]);
				//printf("%f %f %f %f %f\n",m_mm_curr[j], m_gd_curr[j], m_im_curr[j], m_dg_curr[j], m_mi_curr[j]);
				//printf("%f %f %f %f %f\n",q.tr[i - 1][M2M],  q.tr[i - 1][I2M], q.tr[i - 1][D2M], t.tr[j - 1][M2M], t.tr[j - 1][D2M]);

			} // end else

		} //end for j
		for (int jj = 0; jj <= t.L; jj++) {
			// Fill posterior probability matrix with forward score
			p_mm.setPosteriorValue(i, jj, m_fwd[jj].mm_curr);
			m_fwd[jj].mm_prev = m_fwd[jj].mm_curr;
			m_fwd[jj].mi_prev = m_fwd[jj].mi_curr;
			m_fwd[jj].im_prev = m_fwd[jj].im_curr;
			m_fwd[jj].dg_prev = m_fwd[jj].dg_curr;
			m_fwd[jj].gd_prev = m_fwd[jj].gd_curr;

		}
		/* F_MM_prev = m_mm_curr */


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
				hit.Pforward  += p_mm.getPosteriorValue(i, j);

			hit.Pforward *= scale[i + 1];
		}
	} else { // global alignment
		hit.Pforward  = 0.0;
		for (i = 1; i < q.L; ++i)
			hit.Pforward  = (hit.Pforward  + p_mm.getPosteriorValue(i, t.L) * scale[i + 1]);
		for (j = 1; j <= t.L; ++j)
			hit.Pforward  += p_mm.getPosteriorValue(q.L, j);
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


