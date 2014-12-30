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

	void setSingleValue(const int row, const int col, const float value);
	float getSingleValue(const int row, const int col) const;
    void setColScoreValue(const int row, const int col, const float value);
	float getColScoreValue(const int row, const int col) const;
    

	bool isAllocated() const;

	void setAllocated(bool allocated);
	void DeleteProbabilityMatrix();

private:
	int m_q_max_length;
	int m_t_max_length;
    
	float ** m_probabilities;
    float ** m_column_scores;
//	bool m_allocated;
//	float ** m_p_mm;

};




#endif /* HHPOSTERIORMATRIX_H_ */
