/*
 * hhposteriormatrix.cpp
 *
 *  Created on: Mar 2, 2014
 *      Author: stefan
 */

#include "hhposteriormatrix.h"

PosteriorMatrix::PosteriorMatrix() {
		m_allocated = false;
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
  q_length_max += 1;
  t_length_max += 1;

  if(q_length_max > m_q_max_length || t_length_max > m_t_max_length) {
    DeleteProbabilityMatrix();
  }
  else {
    return;
  }
//TODO make just one allocation
    m_probabilities = new float * [q_length_max];
//	printf("p_mm: Allocate 2nd dimension of m_p_mm\n");

	for (int i = 0; i < q_length_max; i++) {
		m_probabilities[i] = new float[t_length_max];
		if (!m_probabilities[i]) {
			fprintf(stderr,"Error: out of memory while allocating row %i (out of %i) for dynamic programming matrices \n", i+1, q_length_max);
			fprintf(stderr,"Please decrease your memory requirements to the available memory using option -maxmem <GBs>\n");
			fprintf(stderr,"You may want to check and increase your stack size limit (Linux: ulimit -a)\n");
			exit(3);
		}
	}
	m_allocated = true;

	m_q_max_length = q_length_max;
	m_t_max_length = t_length_max;
};

void PosteriorMatrix::DeleteProbabilityMatrix() {
  if(m_q_max_length == 0) {
    return;
  }

  for (int i = 0; i < m_q_max_length; i++) {
          delete [] m_probabilities[i];
  }
  delete[] m_probabilities;

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
///////////////////////////////////////////////////////////////////////////////////////////////
// Return a float value of a selected element of a vector
///////////////////////////////////////////////////////////////////////////////////////////////
float PosteriorMatrix::getSingleValue(const int row, const int col) const {
	return m_probabilities[row][col];
//	return ptr[IDX_CORR - elem];
}
///////////////////////////////////////////////////////////////////////////////////////////////
// Set a single float value to an element of a vector
///////////////////////////////////////////////////////////////////////////////////////////////
void PosteriorMatrix::setSingleValue(const int row, const int col, const float value) {
	m_probabilities[row][col] = value;
}

///////////////////////////////////////////////////////////////////////////////////////////////
// Return a simd_float pointer to a vector
///////////////////////////////////////////////////////////////////////////////////////////////
float PosteriorMatrix::getValue(const int row, const int col) const {
	return m_probabilities[row][col];
}


///////////////////////////////////////////////////////////////////////////////////////////////
// Returns true when the posterior matrix has been allocated
///////////////////////////////////////////////////////////////////////////////////////////////
bool PosteriorMatrix::isAllocated() const {
	return m_allocated;
}
///////////////////////////////////////////////////////////////////////////////////////////////
// Set to true when the posterior matrix has been allocated
///////////////////////////////////////////////////////////////////////////////////////////////
void PosteriorMatrix::setAllocated(bool allocated) {
	m_allocated = allocated;
}
