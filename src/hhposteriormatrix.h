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

#include "hhposteriormatrix-inl.h"

class PosteriorMatrix {
public:
	PosteriorMatrix(const int q_length);
	virtual ~PosteriorMatrix();

	void allocateMatrix(const int t_length_max);

	simd_float * getRow(const int row) const;

	void setSingleValue(const int row, const int col, const int elem, const float value);
	float getSingleValue(const int row, const int col, const int elem) const;

	simd_float getValue(const int row, const int col) const;

	bool isAllocated() const;

	void setAllocated(bool allocated);

private:

	const int m_q_length;
	simd_float ** m_probabilities;
	bool m_allocated;
//	float ** m_p_mm;

};




#endif /* HHPOSTERIORMATRIX_H_ */
