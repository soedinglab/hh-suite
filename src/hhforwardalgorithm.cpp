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
		m_curr[j].mm = 0.0;

	m_curr[0].mm = 0.0;
	m_curr[0].im = 0.0;
	m_curr[0].gd = 0.0;
	for (j=1; j<=t.L; ++j)
	{
		if (celloff_matrix.getCellOff(1,j,elem))
			m_curr[j].mm = m_curr[j].mi = m_curr[j].dg = m_curr[j].im = m_curr[j].gd = 0.0;
		else
		{
			m_curr[j].mm = ProbFwd(q.p[1], t.p[j]) * Cshift ;
			m_curr[j].mi = m_curr[j].dg = 0.0;
			m_curr[j].im = m_curr[j-1].mm * q.tr[1][M2I] * t.tr[j-1][M2M] + m_curr[j-1].im * q.tr[1][I2I] * t.tr[j-1][M2M];
			m_curr[j].gd = m_curr[j-1].mm * t.tr[j-1][M2D]                + m_curr[j-1].gd * t.tr[j-1][D2D];
			//printf("%f %f %f %f %f\n",m_mm_curr[j], m_gd_curr[j], m_im_curr[j], m_dg_curr[j], m_mi_curr[j]);
			//printf("%f %f %f %f %f\n",q.tr[1][M2I], t.tr[j-1][M2M], q.tr[1][I2I], t.tr[j-1][M2D] , t.tr[j-1][D2D]);
			//printf("%f %f %f\n",ProbFwd(q.p[1], t.p[j]), m_im_curr[j-1], m_gd_curr[j-1]);
		}
	}

	for (int j = 0; j <= t.L; j++)
	{
		p_mm.setPosteriorValue(0,j, m_prev[j].mm);
		p_mm.setPosteriorValue(1,j, m_curr[j].mm);
		m_prev[j].mm = m_curr[j].mm;
		m_prev[j].mi = m_curr[j].mi;
		m_prev[j].im = m_curr[j].im;
		m_prev[j].dg = m_curr[j].dg;
		m_prev[j].gd = m_curr[j].gd;
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
			m_curr[jmin].mm = m_curr[jmin].mi = m_curr[jmin].dg = m_curr[jmin].im =
					m_curr[jmin].gd = 0.0;
		else {
			m_curr[jmin].mm = scale_prod * ProbFwd(q.p[i], t.p[jmin]) * Cshift;
			m_curr[jmin].im = m_curr[jmin].gd = 0.0;
			m_curr[jmin].mi = scale[i]
					* (m_prev[jmin].mm * q.tr[i - 1][M2M] * t.tr[jmin][M2I]
					+ m_prev[jmin].mi * q.tr[i - 1][M2M] * t.tr[jmin][I2I]);
			m_curr[jmin].dg = scale[i]
					* (m_prev[jmin].mm * q.tr[i - 1][M2D]
					+ m_prev[jmin].dg * q.tr[i - 1][D2D]);
		}

		/* copy back */
		p_mm.setPosteriorValue(i, jmin, m_curr[jmin].mm);

		Pmax_i = 0;
		memset(m_curr+(jmin + 1), 0, (t.L ) * sizeof(PosteriorMatrixCol));
		// Loop through template positions j
		for (j = jmin + 1; j <= t.L; ++j) {

			// Recursion relations
			if (!(celloff_matrix.getCellOff(i, j, elem)))
			{
				m_curr[j].mm = ProbFwd(q.p[i], t.p[j]) * Cshift
						* scale[i]
						* (pmin
						+ m_prev[j - 1].mm * q.tr[i - 1][M2M] * t.tr[j - 1][M2M] // BB -> MM (BB = Begin/Begin, for local alignment)
						+ m_prev[j - 1].gd * q.tr[i - 1][M2M] * t.tr[j - 1][D2M] // GD -> MM
						+ m_prev[j - 1].im * q.tr[i - 1][I2M] * t.tr[j - 1][M2M] // IM -> MM
						+ m_prev[j - 1].dg * q.tr[i - 1][D2M] * t.tr[j - 1][M2M] // DG -> MM
						+ m_prev[j - 1].mi * q.tr[i - 1][M2M] * t.tr[j - 1][I2M] // MI -> MM
				);
				m_curr[j].gd = (
						m_curr[j - 1].mm * t.tr[j - 1][M2D] +         // GD -> MM
								m_curr[j - 1].gd * t.tr[j - 1][D2D]                    // GD -> GD
				);
				m_curr[j].im = (m_curr[j - 1].mm * q.tr[i][M2I] * t.tr[j - 1][M2M] // MM -> IM
						+ m_curr[j - 1].im * q.tr[i][I2I] * t.tr[j - 1][M2M]     // IM -> IM
				);
				m_curr[j].dg = scale[i] * (m_prev[j].mm * q.tr[i - 1][M2D]  // DG -> MM
						+ m_prev[j].dg * q.tr[i - 1][D2D]                    // DG -> DG
				);
				m_curr[j].mi = scale[i]
						* (m_prev[j].mm * q.tr[i - 1][M2M] * t.tr[j][M2I]   // MI -> MM
						+ m_prev[j].mi * q.tr[i - 1][M2M] * t.tr[j][I2I]     // MI -> MI
				);

				Pmax_i = fmax(Pmax_i, m_curr[j].mm);
				//printf("%f %f %f %f %f\n",F_MM_prev[j], F_GD_prev[j], m_im_prev[j], m_dg_prev[j], m_mi_prev[j]);
				//printf("%f %f %f %f %f\n",m_mm_curr[j], m_gd_curr[j], m_im_curr[j], m_dg_curr[j], m_mi_curr[j]);
				//printf("%f %f %f %f %f\n",q.tr[i - 1][M2M],  q.tr[i - 1][I2M], q.tr[i - 1][D2M], t.tr[j - 1][M2M], t.tr[j - 1][D2M]);

			} // end else

		} //end for j
		for (int jj = 0; jj <= t.L; jj++) {
			// Fill posterior probability matrix with forward score
			p_mm.setPosteriorValue(i, jj, m_curr[jj].mm);
		}
		std::swap(m_prev, m_curr);
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


