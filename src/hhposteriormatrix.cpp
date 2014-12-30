/*
 * hhposteriormatrix.cpp
 *
 *  Created on: Mar 2, 2014
 *      Author: stefan
 */

#include "hhposteriormatrix.h"

PosteriorMatrix::PosteriorMatrix() {
		m_probabilities = NULL;
		m_q_max_length = 0;
		m_t_max_length = 0;
}

PosteriorMatrix::~PosteriorMatrix() {
  DeleteProbabilityMatrix();
}


///////////////////////////////////////////////////////////////////////////////////////////////
// Allocate memory in 2nd dimension
///////////////////////////////////////////////////////////////////////////////////////////////
void PosteriorMatrix::allocateMatrix(int q_length_max, int t_length_max) {

    // Allocate posterior probability matrix
    q_length_max += 1;
    t_length_max += 1;
    
    // If the already allocated matrix is sufficiently large, we are done
    if(q_length_max < m_q_max_length && t_length_max < m_t_max_length)
        return;
    if (m_q_max_length>0)
        DeleteProbabilityMatrix();

    m_q_max_length = ICEIL(q_length_max,VECSIZE_FLOAT);
    m_t_max_length = ICEIL(t_length_max,VECSIZE_FLOAT);
   
    
////TODO make just one allocation
//    m_probabilities = new float * [q_length_max];
////	printf("p_mm: Allocate 2nd dimension of m_p_mm\n");
//
//	for (int i = 0; i < q_length_max; i++) {
//		m_probabilities[i] = new float[t_length_max];
//		if (!m_probabilities[i]) {
//			fprintf(stderr,"Error: out of memory while allocating row %i (out of %i) for dynamic programming matrices \n", i+1, q_length_max);
//			fprintf(stderr,"Please decrease your memory requirements to the available memory using option -maxmem <GBs>\n");
//			fprintf(stderr,"You may want to check and increase your stack size limit (Linux: ulimit -a)\n");
//			exit(3);
//		}
//	}

    // Allocate posterior prob matrix (matrix rows are padded to make them aligned to multiples of ALIGN_FLOAT)
    m_probabilities = malloc_matrix<float>(m_q_max_length, m_t_max_length);
    if (!m_probabilities)
        MemoryError("m_probabilities", "hhposteriormatrix.cpp", 55, "PosteriorMatrix::allocateMatrix");

};


void PosteriorMatrix::DeleteProbabilityMatrix() {
    if(m_q_max_length == 0)
        return;

//  for (int i = 0; i < m_q_max_length; i++) {
//          delete [] m_probabilities[i];
//  }
//  delete[] m_probabilities;
 
    free(m_probabilities);
    m_probabilities = NULL;
    m_q_max_length = 0;
    m_t_max_length = 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////
// Return a simd_float pointer to a single row
///////////////////////////////////////////////////////////////////////////////////////////////
float * PosteriorMatrix::getRow(const int row) const {
    return m_probabilities[row];
}



