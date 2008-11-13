#ifndef AB_CLUSTER_H
#define AB_CLUSTER_H
/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

//DESCRIPTION:
//A class for clusters in EM-clustering

#include <string>
#include <fstream>
#include <iostream>
#include <limits>   // infinity
#include <iomanip>
#include <vector>
#include <utility>  //pair
#include <cstdlib>
#include <cstdio>

#include "profile.h"
#include "training_profile.h"
#include "my_exception.h"
#include "amino_acid.h"

class Cluster : public Profile {
public:

    //constructor for dummy cluster object
    Cluster() throw (std::exception);
    //constructor for empty profile object
    Cluster( AminoAcid* aa, Matrix* m ) throw (std::exception);
    //constructor for cluster object from serialized profile object
    Cluster( std::istream& in,
             AminoAcid* aa,
             Matrix* m ) throw (std::exception);
    // copy constructor
    Cluster(const Cluster& rhs);
    // copy constructor from profile object
    Cluster(const Profile& rhs);
    // assignment operator
    virtual Cluster& operator=(const Cluster& rhs);

    inline const size_t get_iteration()  const { return iter; }
    inline const float get_log_alpha()  const { return log_alpha; }
    inline const float get_alpha()  const { return exp(log_alpha); }
    inline const float get_cost()  const { return cost; }
    inline const std::vector<TrainingProfile*>& get_members() { return members; }
    inline const size_t size()  const { return members.size(); }

    inline void set_log_alpha(const float a) { log_alpha = a; }
    inline void set_alpha(const float a) { log_alpha = log(a); }
    inline void set_iteration(const size_t i) { iter = i; }
    inline void set_cost(const float c) { cost = c; }
    inline size_t increment_iteration() { return ++iter; }
    inline void clear_members() { members.clear(); }
    void add_member(TrainingProfile *tp);

protected:

    size_t iter;                                              //EM iterations performed on this profile
    float  log_alpha;                                         //log value of cluster size alpha
    float  cost;                                              //contribution to overall cost function
    std::vector<TrainingProfile*> members;                    //members of this cluster (alll n's with pc[n][k] maximal for this cluster k)

    virtual void copy_internals(const Cluster& rhs) throw (std::exception);
    virtual void read_header(std::istream& in) throw (std::exception);
    virtual std::ostream& write_header(std::ostream& out);
};

#endif

