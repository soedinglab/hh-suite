/***************************************************************************
 *   Copyright (C) 2006 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/
#include "cs_counts.h"

CSCounts::CSCounts( const char* fn, AminoAcid* a, Matrix* m ) throw (std::exception):
    aa(a),
    mat(m),
    NAA(a->get_naa()),
    i2aa(a->get_i2aa()),
    K(0),
    W(0),
    center(0),
    w_center(1.5),
    beta(0.85),
    w(NULL)
{
    init(fn);
}

CSCounts::~CSCounts()
{
    for(std::vector<SimpleCluster*>::iterator it=clusters.begin(); it!=clusters.end(); ++it) delete *it;
    if (w!=NULL) delete [] w;
}

void CSCounts::init(const char* fn) throw (std::exception)
{
    std::ifstream fin(fn);
    size_t i=0;
    while (fin.good() && !fin.eof() && fin.peek()!=-1) {
        SimpleCluster *c = new SimpleCluster(fin, aa, mat);
        clusters.push_back(c);
        ++i;
    }
    fin.close();
    K = i;
    if (K==0) throw MyException("No clusters provided in '%s'!", fn);
    W = clusters[0]->get_length();
    center = clusters[0]->get_center();

    //init weights
    w = new float[W];
    w[center] = w_center;
    for(size_t i=1; i<=center; ++i) {
        float weight = w_center * pow(beta, i);
        w[center-i] = weight;
        w[center+i] = weight;
    }
}

void CSCounts::set_weight_center(const float wcenter)
{
    w_center = wcenter;
    w[center] = w_center;
    for(size_t i=1; i<=center; ++i) {
        float weight = w_center * pow(beta, i);
        w[center-i] = weight;
        w[center+i] = weight;
    }
}

void CSCounts::set_beta(const float b) { beta = b; }

void CSCounts::prepare_pseudocounts(HMM& q) throw (std::exception)
{
    float** pq   = q.f;
    float* neffi = q.Neff_M;
    double** pc  = new double*[q.L+1];

    //init arrays
    for(size_t i=1; i<=static_cast<size_t>(q.L); ++i) {
        pc[i] = new double[NAA];
        for(int a=0; a<NAA; ++a) {
            pc[i][a] =0.0;
        }
    }

    //compute pseudocount vectors
    for(size_t k=0; k<K; ++k) {
        float** log_pk = clusters[k]->get_log_profile();

        for(size_t i=1; i<=static_cast<size_t>(q.L); ++i) {
            double log_pki = clusters[k]->get_log_alpha();
            size_t beg = std::max(1, static_cast<int>(i)-static_cast<int>(center));
            size_t end = std::min(q.L,static_cast<int>(i+center));

            for(size_t l=beg; l<=end; ++l) {
                size_t j = l-i+center;

                double sum = 3.0*neffi[l]; // 3.0*neffi[l] to keep log_pki within range of -127 to +127
                sum += neffi[l]*pq[l][0]*log_pk[j][0];
                sum += neffi[l]*pq[l][1]*log_pk[j][1];
                sum += neffi[l]*pq[l][2]*log_pk[j][2];
                sum += neffi[l]*pq[l][3]*log_pk[j][3];
                sum += neffi[l]*pq[l][4]*log_pk[j][4];
                sum += neffi[l]*pq[l][5]*log_pk[j][5];
                sum += neffi[l]*pq[l][6]*log_pk[j][6];
                sum += neffi[l]*pq[l][7]*log_pk[j][7];
                sum += neffi[l]*pq[l][8]*log_pk[j][8];
                sum += neffi[l]*pq[l][9]*log_pk[j][9];
                sum += neffi[l]*pq[l][10]*log_pk[j][10];
                sum += neffi[l]*pq[l][11]*log_pk[j][11];
                sum += neffi[l]*pq[l][12]*log_pk[j][12];
                sum += neffi[l]*pq[l][13]*log_pk[j][13];
                sum += neffi[l]*pq[l][14]*log_pk[j][14];
                sum += neffi[l]*pq[l][15]*log_pk[j][15];
                sum += neffi[l]*pq[l][16]*log_pk[j][16];
                sum += neffi[l]*pq[l][17]*log_pk[j][17];
                sum += neffi[l]*pq[l][18]*log_pk[j][18];
                sum += neffi[l]*pq[l][19]*log_pk[j][19];
                log_pki += w[j] * sum;
            }

//            if (log_pki>127) std::cerr << "log_pki=" << log_pki << std::endl;

            for(int a=0; a<NAA; ++a) {
                float tmp = 1.442695041 * (log_pki + log_pk[center][a]);
                if (tmp<=127)
                    pc[i][a] += fast_pow2(tmp);
                else
                    pc[i][a] += exp(log_pki + log_pk[center][a]);
            }

//             pc[i][0]  += fast_pow2(1.442695041 * (log_pki + log_pk[center][0]));
//             pc[i][1]  += fast_pow2(1.442695041 * (log_pki + log_pk[center][1]));
//             pc[i][2]  += fast_pow2(1.442695041 * (log_pki + log_pk[center][2]));
//             pc[i][3]  += fast_pow2(1.442695041 * (log_pki + log_pk[center][3]));
//             pc[i][4]  += fast_pow2(1.442695041 * (log_pki + log_pk[center][4]));
//             pc[i][5]  += fast_pow2(1.442695041 * (log_pki + log_pk[center][5]));
//             pc[i][6]  += fast_pow2(1.442695041 * (log_pki + log_pk[center][6]));
//             pc[i][7]  += fast_pow2(1.442695041 * (log_pki + log_pk[center][7]));
//             pc[i][8]  += fast_pow2(1.442695041 * (log_pki + log_pk[center][8]));
//             pc[i][9]  += fast_pow2(1.442695041 * (log_pki + log_pk[center][9]));
//             pc[i][10] += fast_pow2(1.442695041 * (log_pki + log_pk[center][10]));
//             pc[i][11] += fast_pow2(1.442695041 * (log_pki + log_pk[center][11]));
//             pc[i][12] += fast_pow2(1.442695041 * (log_pki + log_pk[center][12]));
//             pc[i][13] += fast_pow2(1.442695041 * (log_pki + log_pk[center][13]));
//             pc[i][14] += fast_pow2(1.442695041 * (log_pki + log_pk[center][14]));
//             pc[i][15] += fast_pow2(1.442695041 * (log_pki + log_pk[center][15]));
//             pc[i][16] += fast_pow2(1.442695041 * (log_pki + log_pk[center][16]));
//             pc[i][17] += fast_pow2(1.442695041 * (log_pki + log_pk[center][17]));
//             pc[i][18] += fast_pow2(1.442695041 * (log_pki + log_pk[center][18]));
//             pc[i][19] += fast_pow2(1.442695041 * (log_pki + log_pk[center][19]));
        }
    }

    //add pseudocounts to profile
    for(size_t i=1; i<=static_cast<size_t>(q.L); ++i) {
        normalize_to_one(pc[i], NAA);
        for(int a=0; a<NAA; ++a) {
            assert(pc[i][a]>0);
            q.g[i][a]=pc[i][a];
        }

    }

    //free memory
    for(size_t i=1; i<=static_cast<size_t>(q.L); ++i) delete [] pc[i];
    delete [] pc;
}


