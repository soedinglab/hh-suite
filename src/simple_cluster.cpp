/***************************************************************************
 *   Copyright (C) 2006 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/
#include "simple_cluster.h"

// constructor for cluster object from input stream
SimpleCluster::SimpleCluster( std::istream& in,
                              AminoAcid* aa,
                              Matrix* m ) throw (std::exception) :
    Cluster::Cluster(aa, m) {
    read(in);
}

SimpleCluster::SimpleCluster(const SimpleCluster& rhs) { copy_internals(rhs); }

SimpleCluster& SimpleCluster::operator=(const SimpleCluster& rhs) {
    // identity test, if self-assignment, do nothing
    if (this == &rhs) return *this;
    copy_internals(rhs);
    return *this;
}

void SimpleCluster::read_columns(std::istream &in) throw (std::exception) {
    std::string line;            // line read from input stream
    const char* ptr;             // pointer for string manipulation
    size_t i_curr = 0;           // current index for profile column
    size_t i_prev = 0;           // previous index for profile column
    size_t line_mode = 0;        // mode of next line; 0: line with probabilities, 1: line with additional info
    std::vector<float*> p_tmp;   // temporary float matrix
    std::vector<float> neff_tmp; // temporary vector with number of effective sequences
    std::string seq_tmp;         // temporary sequence string

    getline(in, line); //skip first two description lines
    while( getline(in, line) ) {
        if (i_curr==0 && i_prev==0 && isspace(line[0])) getline(in, line); // read next line if Neff label is present
        if (strscn_c(line.c_str())==NULL) continue;                          // skip lines that contain only white space
        if (line[0]=='#') continue;                                        // skip comment lines
        if (line[0]=='/' && line[1]=='/') break;                           // read separator line

        if (line_mode==1 && isalpha(line[0])) {
            line_mode=0;
            neff_tmp.push_back(0);
        }
        if (line_mode==0) { // read probability line
            // read sequence character
            ptr = line.c_str();
            seq_tmp.push_back(aa2i[ptr[0]]);
            // read column index
            i_prev = i_curr;
            i_curr = strtoi_(++ptr);
            if( i_curr!=i_prev+1 )
                throw MyException("In profile %s column %i is followed by column %i!", name.c_str(), i_prev, i_curr);
            // read amino acid probabilities of current column and store them in col
            float *col = new float[NAA_ANY];
            for(int i=0; i<NAA; ++i)
                col[al2i[i]] = pow(2.0, static_cast<float>(-strtoi_(ptr)) / SCALE_FAC);
            col[ANY]=1.0;
            p_tmp.push_back(col);
            line_mode = 1;
        } else { // read line with additional column information
            // read number of effective sequences in column i_curr
            ptr = line.c_str();
            neff_tmp.push_back( static_cast<float>(strtoi(++ptr))/SCALE_FAC );
            line_mode = 0;
        }
    }
    if (neff_tmp.size()==p_tmp.size()-1) neff_tmp.push_back(0.0);
    if( p_tmp.size()!=len )
        throw MyException("Profile %s has length %i but has %i column records!", name.c_str(), len, p_tmp.size());
    if( !neff_tmp.empty() && neff_tmp.size()!=len )
        throw MyException("Profile %s has length %i but has %i Neff records!", name.c_str(), len, neff_tmp.size());

    //put columns into p
    logp = new float*[len];
    for(size_t i=0; i<len; ++i) {
        logp[i] = new float[NAA_ANY];
        for(int a=0; a<NAA_ANY; ++a) logp[i][a] = log( p_tmp[i][a] );
    }
    if (!neff_tmp.empty()) {
        //put neff_tmp values into neff
        neffi = new float[len];
        for(size_t i=0; i<len; ++i) neffi[i] = neff_tmp[i];
    }

    for(std::vector<float*>::iterator it=p_tmp.begin(); it!=p_tmp.end(); ++it) delete [] *it;
    p_tmp.clear();

    //put sequence into seq
    seq = new char[len];
    for(size_t i=0; i<len; ++i) seq[i]=seq_tmp[i];
    seq_tmp.clear();
}

void SimpleCluster::free_memory()
{
    // no need to free matrix p since it's never allocated
    if (logp) {
        for(size_t i=0; i<len; ++i) delete [] logp[i];
        delete [] logp;
    }
    if (neffi)
        delete [] neffi;
    if (seq)
        delete [] seq;
    logp  = 0;
    neffi = 0;
    seq   = 0;
}
