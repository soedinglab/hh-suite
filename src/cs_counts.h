#ifndef AB_CS_COUNTS_H
#define AB_CS_COUNTS_H
/***************************************************************************
 *   Copyright (C) 2006 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <utility> // pair

#include "amino_acid.h"
#include "sequence.h"
#include "profile.h"
#include "simple_cluster.h"
#include "hhhmm.h"

//#define __AB_DEBUG_PRINT_CLUSTERS 1

class CSCounts {
public:

    CSCounts(const char* fn, AminoAcid *a, Matrix *m) throw (std::exception);
    ~CSCounts();

    inline const std::vector<SimpleCluster*>& get_clusters() const { return clusters; }
    inline const float get_number_of_clusters() const { return K; }
    inline const float get_window_length() const { return W; }
    inline const float get_weight_center() const { return w_center; }
    inline const float get_beta() const { return beta; }
    void set_beta(const float);
    void set_weight_center(const float);

    void prepare_pseudocounts(HMM& q) throw (std::exception);

private:
    AminoAcid *aa;
    Matrix *mat;

    int NAA;
    const char* i2aa;

    std::vector<SimpleCluster*> clusters;
    size_t K;
    size_t W;
    size_t center;
    float w_center;
    float beta;

    float *w;

    void init(const char* fn) throw (std::exception);
};

#endif
