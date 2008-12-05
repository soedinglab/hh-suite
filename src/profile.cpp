/***************************************************************************
 *   Copyright (C) 2006 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/
#include "profile.h"

// constructor for dummy profile object
Profile::Profile() throw (std::exception):
    NAA(20),
    NAA_ANY(21),
    ANY(20),
    GAP(21),
    aa2i(0),
    i2aa(0),
    al2i(0),
    al2aa(0),
    p_cond(0),
    p_background(0),
    p(0),
    len(0),
    index(0),
    center(0),
    name(),
    file(),
    date(),
    seq(0),
    neffi(0),
    neff(0.0)
{
}

// constructor for empty profile object
Profile::Profile( AminoAcid* aa, Matrix* m ) throw (std::exception):
    NAA(aa->get_naa()),
    NAA_ANY(aa->get_naa_any()),
    ANY(aa->get_any()),
    GAP(aa->get_gap()),
    aa2i(aa->get_aa2i()),
    i2aa(aa->get_i2aa()),
    al2i(aa->get_al2i()),
    al2aa(aa->get_al2aa()),
    p_cond(m->get_p_cond()),
    p_background(m->get_p_back()),
    p(0),
    logp(0),
    len(0),
    index(0),
    center(0),
    name(),
    file(),
    date(),
    seq(0),
    neffi(0),
    neff(0.0)
{
}

// constructor for profile of given length with probabilities initialized to zero
Profile::Profile( size_t l, AminoAcid* aa, Matrix* m ) throw (std::exception):
    NAA(aa->get_naa()),
    NAA_ANY(aa->get_naa_any()),
    ANY(aa->get_any()),
    GAP(aa->get_gap()),
    aa2i(aa->get_aa2i()),
    i2aa(aa->get_i2aa()),
    al2i(aa->get_al2i()),
    al2aa(aa->get_al2aa()),
    p_cond(m->get_p_cond()),
    p_background(m->get_p_back()),
    p(0),
    logp(0),
    len(l),
    index(0),
    center(len%2==1 ? (len-1)/2 : 0),
    name(),
    file(),
    date(),
    seq(0),
    neffi(0),
    neff(0.0)
{
    init();
}

// constructor for profile object from input stream
Profile::Profile( std::istream& in,
                  AminoAcid* aa, Matrix* m ) throw (std::exception):
    NAA(aa->get_naa()),
    NAA_ANY(aa->get_naa_any()),
    ANY(aa->get_any()),
    GAP(aa->get_gap()),
    aa2i(aa->get_aa2i()),
    i2aa(aa->get_i2aa()),
    al2i(aa->get_al2i()),
    al2aa(aa->get_al2aa()),
    p_cond(m->get_p_cond()),
    p_background(m->get_p_back()),
    p(0),
    logp(0),
    len(0),
    index(0),
    center(0),
    name(),
    file(),
    date(),
    seq(0),
    neffi(0),
    neff(0.0)
{
    read(in);
}

//constructor for profile object from sequence object
Profile::Profile( const Sequence& sequence,
         AminoAcid* aa,
         Matrix* m ) throw (std::exception):
    NAA(aa->get_naa()),
    NAA_ANY(aa->get_naa_any()),
    ANY(aa->get_any()),
    GAP(aa->get_gap()),
    aa2i(aa->get_aa2i()),
    i2aa(aa->get_i2aa()),
    al2i(aa->get_al2i()),
    al2aa(aa->get_al2aa()),
    p_cond(m->get_p_cond()),
    p_background(m->get_p_back()),
    p(0),
    logp(0),
    len(sequence.length()),
    index(0),
    center(sequence.length()%2==1 ? (sequence.length()-1)/2 : 0),
    name(sequence.get_header()),
    file(),
    date(),
    seq(0),
    neffi(0),
    neff(0.0)
{
    init();

    const char *s = sequence.get_sequence();
    for(size_t i=0; i<len; ++i) {
        seq[i]   = s[i];
        neffi[i] = 1.0;

        if (s[i]<ANY) {
            p[i][s[i]]    = 1.0;
            logp[i][s[i]] = 0.0;
        } else if (s[i]==ANY) {
            for(int a=0; a<NAA; ++a) {
                p[i][a]    = m->get_p_back()[a];
                logp[i][a] = log(p[i][a]);
            }
        } else
            throw MyException("Found invalid character at position %i of sequence %s!", i+1, name.c_str());
    }
    neff=1.0;
}

Profile::Profile(const Profile& rhs):
    NAA(rhs.NAA),
    NAA_ANY(rhs.NAA_ANY),
    ANY(rhs.ANY),
    GAP(rhs.GAP)
{
    copy_internals(rhs);
}

Profile::~Profile() { free_memory(); }

void Profile::free_memory()
{
    if (p) {
        for(size_t i=0; i<len; ++i) delete [] p[i];
        delete [] p;
    }
    if (logp) {
        for(size_t i=0; i<len; ++i) delete [] logp[i];
        delete [] logp;
    }
    if (neffi)
        delete [] neffi;
    if (seq)
        delete [] seq;

    p     = 0;
    logp  = 0;
    neffi = 0;
    seq   = 0;
}

Profile& Profile::operator=(const Profile& rhs)
{
    // identity test, if self-assignment, do nothing
    if (this == &rhs) return *this;
    free_memory();
    copy_internals(rhs);
    return *this;
}

void Profile::init() throw (std::exception)
{
    seq = new char[len];
    for(size_t i=0; i<len; ++i) seq[i]=' ';

    p = new float*[len];
    for(size_t i=0; i<len; ++i) {
        p[i] = new float[NAA_ANY];
        for(int a=0; a<NAA; ++a)
            p[i][a] = 0.0;
        p[i][ANY]=1.0;
    }

    logp = new float*[len];
    for(size_t i=0; i<len; ++i) {
        logp[i] = new float[NAA_ANY];
        for(int a=0; a<NAA; ++a)
            logp[i][a] = -std::numeric_limits<float>::infinity();
        logp[i][ANY]=0.0;
    }

    neff=0.0;
    neffi = new float[len];
    for (size_t i=0; i<len; ++i) neffi[i] = 0.0;

    time_t now = time(0);
    struct tm* timeinfo  = localtime(&now);
    date = asctime(timeinfo);
    date = date.substr(0,date.size()-1);
}

void Profile::copy_internals(const Profile& rhs) throw (std::exception)
{
    //flat copy simple data members
    aa2i = rhs.aa2i;
    i2aa = rhs.i2aa;
    al2i = rhs.al2i;
    al2aa = rhs.al2aa;
    p_cond = rhs.p_cond;
    p_background = rhs.p_background;
    len = rhs.len;
    index = rhs.index;
    center = rhs.center;
    name = rhs.name;
    file = rhs.file;
    date = rhs.date;
    neff = rhs.neff;

    seq = new char[rhs.len];
    for(size_t i=0; i<len; ++i) seq[i]=rhs.seq[i];

    p = new float*[rhs.len];
    for(size_t i=0; i<rhs.len; ++i) {
        p[i] = new float[rhs.NAA_ANY];
        for(int a=0; a<rhs.NAA_ANY; ++a)
            p[i][a] = rhs.p[i][a];
    }

    logp = new float*[rhs.len];
    for(size_t i=0; i<rhs.len; ++i) {
        logp[i] = new float[rhs.NAA_ANY];
        for(int a=0; a<rhs.NAA_ANY; ++a)
            logp[i][a] = rhs.logp[i][a];
    }

    neffi = new float[rhs.len];
    for(size_t i=0; i<rhs.len; ++i)
        neffi[i] = rhs.neffi[i];
}

void Profile::set_index(const size_t idx) { index = idx; }

void Profile::set_name(const std::string newname) { name = newname; }

void Profile::set_file(const std::string filename) { file = filename; }

void Profile::read(std::istream &in) throw (std::exception) {
    read_header(in);
    read_columns(in);
}

void Profile::calculate_number_of_effective_sequences() throw (std::exception)
{
    for (size_t i=0; i<len; ++i) {
        neffi[i]=0;
        for (int a=0; a<NAA; ++a)
            if (p[i][a]>1E-10) neffi[i]-=p[i][a]*log2(p[i][a]);
        neffi[i] = pow(2.0, neffi[i]);
    }
    neff=0.0;
    for (size_t i=0; i<len; ++i) neff+=neffi[i];
    neff/=len;
}

void Profile::reset() throw (std::exception)
{
    neff=0.0;
    for(size_t i=0; i<len; ++i) {
        for(int a=0; a<NAA_ANY; ++a) {
            p[i][a] = 0.0;
            logp[i][a] = -std::numeric_limits<float>::infinity();
        }
        if (neffi!=0) neffi[i]=0.0;
    }
}

void Profile::read_header(std::istream &in) throw (std::exception) {
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
        else if (str4 == "NEFF")
            neff = atof( line.substr(line.find_first_not_of(" \t\n", 4)).c_str() );
    }
}

void Profile::read_columns(std::istream &in) throw (std::exception) {
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
    p = new float*[len];
    logp = new float*[len];
    for(size_t i=0; i<len; ++i) {
        p[i] = p_tmp[i];
        logp[i] = new float[NAA_ANY];
        for(int a=0; a<NAA_ANY; ++a) logp[i][a] = log( p[i][a] );
    }
    if (!neff_tmp.empty()) {
        //put neff_tmp values into neff
        neffi = new float[len];
        for(size_t i=0; i<len; ++i) neffi[i] = neff_tmp[i];
    }

    p_tmp.clear();

    //put sequence into seq
    seq = new char[len];
    for(size_t i=0; i<len; ++i) seq[i]=seq_tmp[i];
    seq_tmp.clear();
}

std::ostream& Profile::write(std::ostream &out) {
    write_header(out);
    write_columns(out);
    return out;
}

std::ostream& Profile::write_hmm(std::ostream &out) {
    write_hmm_header(out);
    write_hmm_columns(out);
    return out;
}

std::ostream& Profile::write_header(std::ostream &out) {
    // print profile header
    out << "NAME  " << name << std::endl;
    out << "FILE  " << file << std::endl;
    out << "DATE  " << date << std::endl;
    out << "LENG  " << len << std::endl;
    out << "IDX   " << index << std::endl;
    out << "NEFF  " << neff << std::endl;
    out << "#" << std::endl;

    return out;
}

std::ostream& Profile::write_hmm_header(std::ostream &out)
{
    // print profile header
    out << "HHsearch 1.5" << std::endl;
    out << "NAME  " << name << std::endl;
    out << "FILE  " << file << std::endl;
    out << "DATE  " << date << std::endl;
    out << "LENG  " << len << std::endl;
    out << "NEFF  " << 1.0 << std::endl;
    out << "SEQ   " << std::endl;
    out << ">" << name << std::endl;
    for(size_t i=0; i<len; ++i) {
        if (i%100==0 &&i>0) out << std::endl;
        out << i2aa[seq[i]];
    }
    out << std::endl << "#" << std::endl;
    return out;
}

std::ostream& Profile::write_columns(std::ostream &out) {
    // print amino acid letters
    out << std::left << std::setw(7) << "PROF";
    for(int k=0; k<NAA; ++k) out << al2aa[k] << "\t";
    out << std::endl;
    out << std::left << std::setw(7) << "" << "Neff" << std::endl;
    // print aminio acid probabilities
    for(size_t i=0; i<len; ++i) {
        out << std::left << std::setw(2) << i2aa[seq[i]] << std::left << std::setw(5) << i+1;
        for(int j=0; j<NAA; ++j) {
            double logval = log2(p[i][al2i[j]]);
            if (-logval != std::numeric_limits<double>::infinity())
                out << -iround(logval*SCALE_FAC) << "\t";
            else
                out << "*\t";
        }
        out << i+1 << std::endl;
        out << std::left << std::setw(7) << "" << iround(neffi[i]*SCALE_FAC) << std::endl;
    }
    out << "//" << std::endl;

    return out;
}

std::ostream& Profile::write_hmm_columns(std::ostream &out)
{
    // print amino acid background frequencies
    out << std::left << std::setw(7) << "NULL";
    for(int k=0; k<NAA; ++k) {
        double logval = log2(p_background[al2i[k]]);
        if (-logval != std::numeric_limits<double>::infinity())
            out << -iround(logval*SCALE_FAC) << "\t";
        else
            out << "*\t";
    }
    out << std::endl;
    // print amino acid letters
    out << std::left << std::setw(7) << "HMM";
    for(int k=0; k<NAA; ++k) out << al2aa[k] << "\t";
    out << std::endl;
    // print transition labels
    out << std::left << std::setw(7) << "" << "M->M\tM->I\tM->D\tI->M\tI->I\tD->M\tD->D\tNeff\tNeff_I\tNeff_D\n";
    // print transition probabilities from begin state
    out << std::left << std::setw(7) << "" << "0\t*\t*\t0\t*\t0\t*\t*\t*\t*\n";
    // print aminio acid probabilities
    for(size_t i=0; i<len; ++i) {
        out << std::left << std::setw(2) << i2aa[seq[i]] << std::left << std::setw(5) << i+1;
        for(int j=0; j<NAA; ++j) {
            double logval = log2(p[i][al2i[j]]);
            if (-logval != std::numeric_limits<double>::infinity())
                out << -iround(logval*SCALE_FAC) << "\t";
            else
                out << "*\t";
        }
        out << i+1 << std::endl;

        // print transition probabilities
        out << std::left << std::setw(7) << "" << "0\t*\t*\t*\t*\t*\t*\t1000\t0\t0\n";
    }
    out << "//" << std::endl;

    return out;
}

std::ostream& Profile::write_binary(std::ostream &out)
{
    int len_tmp = len;
    out.write(reinterpret_cast<char*>(&len_tmp),sizeof(int));
    for(size_t i=0; i<len; ++i) out.write(&i2aa[seq[i]],sizeof(char));
    for(size_t i=0; i<len; ++i)
        for(int a=0; a<NAA; ++a) {
            double p_tmp = p[i][a];
            out.write(reinterpret_cast<char*>(&p_tmp),sizeof(double));
        }
    return out;
}

void Profile::add_pseudocounts(float pca, float pcb)
{
    float tau=0.0;
    float **f = new float*[len];

    //copy p to f
    for(size_t i=0; i<len; ++i) {
        f[i] = new float[NAA];
        for(int a=0; a<NAA; ++a)
            f[i][a] = p[i][a];
    }

    //add substitution matrix pseudocounts
    for(size_t i=0; i<len; ++i) {
        tau = std::min(1.0, pca * (1.0 + 1.0/pcb) / (1.0 + neffi[i]/pcb) );
        for(int a=0; a<NAA; ++a) {
            float sum = 0.0;
            for(int b=0; b<NAA; ++b)
                sum += p_cond[a][b]*f[i][b];
            p[i][a] = (1.0-tau)*f[i][a] + tau*sum;
        }
    }
    delete [] f;

    normalize();
}

void Profile::log2lin() throw (std::exception)
{
    for(size_t i=0; i<len; ++i)
        for(int a=0; a<NAA_ANY; ++a)
            p[i][a] = exp( logp[i][a] );
}

void Profile::lin2log() throw (std::exception)
{
    for(size_t i=0; i<len; ++i)
        for(int a=0; a<NAA_ANY; ++a) {
            logp[i][a] = p[i][a]==0 ? -FLT_MAX : log(p[i][a]);
        }
}

void Profile::normalize() throw (std::exception)
{
    for(size_t i=0; i<len; ++i) normalize_to_one(p[i], NAA, p_background);
    lin2log();
}


