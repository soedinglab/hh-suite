/*
 * hhposteriormatrix.h
 *
 *  Created on: Mar 2, 2014
 *      Author: stefan
 */

#ifndef HHPOSTERIORMATRIX_H_
#define HHPOSTERIORMATRIX_H_

#include "simd.h"
#include "hhhmmsimd.h"

static const int VEC_SIZE = HMMSimd::VEC_SIZE;
static const int IDX_CORR = HMMSimd::VEC_SIZE - 1;

class PosteriorMatrix {
public:
	PosteriorMatrix();
	virtual ~PosteriorMatrix();

	void allocateMatrix(const int q_length_max, const int t_length_max);

	float * getRow(const int row) const;
    float * getColScoreRow(const int row) const;
	///////////////////////////////////////////////////////////////////////////////////////////////
	// Return a float value of a selected element of matrix
	///////////////////////////////////////////////////////////////////////////////////////////////
	inline float getPosteriorValue(const int row, const int col) const {
		return m_probabilities[row][col];
	}

	///////////////////////////////////////////////////////////////////////////////////////////////
	// Set a single float value to an element of matrix
	///////////////////////////////////////////////////////////////////////////////////////////////
	inline void setPosteriorValue(const int row, const int col, const float value) {
		m_probabilities[row][col] = value;
	}

	///////////////////////////////////////////////////////////////////////////////////////////////
	// Multiply a single float value to an element of matrix
	///////////////////////////////////////////////////////////////////////////////////////////////
	inline void multiplyPosteriorValue(const int row, const int col, const float value) {
		m_probabilities[row][col] *= value;
	}

	void DeleteProbabilityMatrix();

private:
	int m_q_max_length;
	int m_t_max_length;
	float ** m_probabilities;

};




#endif /* HHPOSTERIORMATRIX_H_ */
