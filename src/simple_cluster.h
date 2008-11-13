#ifndef AB_SIMPLE_CLUSTER_H
#define AB_SIMPLE_CLUSTER_H
/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

//DESCRIPTION:
//A class for simple clusters specifically tailored to memoriy efficiency

#include <string>
#include <fstream>
#include <iostream>
#include <limits>   // infinity
#include <iomanip>
#include <vector>
#include <cstdlib>
#include <cstdio>

#include "cluster.h"
#include "my_exception.h"
#include "amino_acid.h"

class SimpleCluster : public Cluster {
public:
    //constructor for cluster object from serialized profile object
    SimpleCluster( std::istream& in,
                   AminoAcid* aa,
                   Matrix* m ) throw (std::exception);
    // copy constructor
    SimpleCluster(const SimpleCluster& rhs);
    // assignment operator
    virtual SimpleCluster& operator=(const SimpleCluster& rhs);

protected:

    virtual void read_columns(std::istream& in) throw (std::exception);
};

#endif

