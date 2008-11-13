/***************************************************************************
 *   Copyright (C) 2006 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/
#include "cluster.h"

// constructor for empty cluster object
Cluster::Cluster() throw (std::exception) :
    Profile::Profile(),
    iter(0),
    log_alpha(0.0),
    cost(0.0) {
}

// constructor for empty cluster object
Cluster::Cluster( AminoAcid* aa,
                  Matrix* m ) throw (std::exception) :
    Profile::Profile(aa, m),
    iter(0),
    log_alpha(0.0),
    cost(0.0) {
}

// constructor for cluster object from input stream
Cluster::Cluster( std::istream& in,
                  AminoAcid* aa,
                  Matrix* m ) throw (std::exception) :
    Profile::Profile(aa, m),
    iter(0),
    log_alpha(0.0),
    cost(0.0) {

    read(in);
}

Cluster::Cluster(const Cluster& rhs) { copy_internals(rhs); }

Cluster::Cluster(const Profile& prof) : iter(0), log_alpha(0.0) {
    Profile::copy_internals(prof);
}

Cluster& Cluster::operator=(const Cluster& rhs) {
    // identity test, if self-assignment, do nothing
    if (this == &rhs) return *this;
    copy_internals(rhs);
    return *this;
}

void Cluster::copy_internals(const Cluster& rhs) throw (std::exception)  {
    Profile::copy_internals(rhs);
    iter     = rhs.iter;
    log_alpha = rhs.log_alpha;
    for(std::vector<TrainingProfile*>::const_iterator it=rhs.members.begin(); it!=rhs.members.end(); ++it)
        members.push_back(*it);
}

void Cluster::add_member(TrainingProfile *tp) {
    members.push_back(tp);
}

void Cluster::read_header(std::istream &in) throw (std::exception) {
    std::string line;
    while( getline(in, line) ) {
        if (strscn_c(line.c_str())==NULL) continue; //skip lines that contain only white space
        if (line[0]=='#') break;                  //end of header section
        std::string str3( line.substr(0, 3) );
        std::string str4( line.substr(0, 4) );

        if (str4 == "NAME")
            name = line.substr( line.find_first_not_of(" \t\n", 4) );
        else if (str4 == "FILE")
            file = line.substr( line.find_first_not_of(" \t\n", 4) );
        else if (str4 == "DATE")
            date = line.substr( line.find_first_not_of(" \t\n", 4) );
        else if (str4 == "LENG") {
            len = atoi( line.substr(line.find_first_not_of(" \t\n", 4)).c_str() );
            center = (len-1)/2;
        } else if (str3 == "IDX")
            index = atoi( line.substr(line.find_first_not_of(" \t\n", 3)).c_str() );
        else if (str4 == "ITER")
            iter = atoi( line.substr(line.find_first_not_of(" \t\n", 4)).c_str() );
        else if (str4 == "ALPH")
            log_alpha = log( atof( line.substr(line.find_first_not_of(" \t\n", 4)).c_str() ) );
    }
}

std::ostream& Cluster::write_header(std::ostream &out) {
    // print profile header
    out << "NAME  " << name << std::endl;
    out << "FILE  " << file << std::endl;
    out << "DATE  " << date << std::endl;
    out << "LENG  " << len << std::endl;
    out << "IDX   " << index << std::endl;
    out << "ITER  " << iter << std::endl;
    out << "ALPH  " << std::setprecision(10) << exp(log_alpha) << std::endl;
    out << "#" << std::endl;

    return out;
}
