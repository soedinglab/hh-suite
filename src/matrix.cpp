/***************************************************************************
 *   Copyright (C) 2006 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/
#include "matrix.h"

Matrix::Matrix(const Matrix::mtype t)throw (std::exception):AMINOACID_DIM(21) {
    matrix = new float*[AMINOACID_DIM];
    for( size_t i=0; i<AMINOACID_DIM; ++i ) matrix[i] = new float[AMINOACID_DIM];
    
    original = new float*[AMINOACID_DIM];
    for( size_t i=0; i<AMINOACID_DIM; ++i ) original[i] = new float[AMINOACID_DIM];
    
    p_background = new float[AMINOACID_DIM];
    for( size_t i=0; i<AMINOACID_DIM; ++i ) p_background[i] = 0.0f;

    p_cond = new float*[AMINOACID_DIM];
    for( size_t i=0; i<AMINOACID_DIM; ++i ) p_cond[i] = new float[AMINOACID_DIM];
    
    i2aa  = new char[AMINOACID_DIM];
    i2aa[AMINOACID_DIM-1] = 'X';
    //set x-row and -colum 
    for( size_t i=0; i<AMINOACID_DIM; ++i){
        matrix[i][AMINOACID_DIM-1] = -1.0f;
        matrix[AMINOACID_DIM-1][i] = -1.0f;
        original[i][AMINOACID_DIM-1] = 0.0f;
        original[AMINOACID_DIM-1][i] = 0.0f;
    }
    
    const size_t char_size = (int)pow(2, 8*sizeof(char));
    aa2i = new int[char_size];
    for(size_t i=0; i<char_size; ++i) aa2i[i]=-1;
    
    if( t!=static_blosum62 ) _read_blosum_matrix(_get_fn(t).c_str());
    else             _copy_dummy();
    
    for(size_t i=0; i<AMINOACID_DIM; ++i) aa2i[i2aa[i]] = i;	
    
    reset();
}

Matrix::~Matrix(){
    for( size_t i=0; i<AMINOACID_DIM; ++i ) delete [] matrix[i];
    delete [] matrix;
    for( size_t i=0; i<AMINOACID_DIM; ++i ) delete [] original[i];
    delete [] original;
    for( size_t i=0; i<AMINOACID_DIM; ++i ) delete [] p_cond[i];
    delete [] p_cond;
    delete [] p_background;
    delete [] i2aa;
    delete [] aa2i;
}

void Matrix::_copy_dummy(){
    for(size_t i=0; i<21; ++i)
        for(size_t j=0; j<21; ++j)
            original[i][j] = __DUMMY_MATRIX[i][j];
    for(size_t i=0; i<21; ++i) i2aa[i]       = __DUMMY_I2AA[i];
    for(size_t i=0; i<21; ++i) p_background[i] = __DUMMY_P_BACKGROUND[i];
    scale = __DUMMY_SCALE;
}

const float** const Matrix::get_matrix()      const { return (const float**) matrix;   }
const float** const Matrix::get_prob_matrix() const { return (const float**) original; }
const float* const Matrix::get_p_back()       const { return (const float*) p_background; }
const float** const Matrix::get_p_cond()       const { return (const float**) p_cond;   }
const size_t Matrix::get_aa_dim()             const { return AMINOACID_DIM;           }
const char* const Matrix::get_i2aa()        const { return i2aa;                  }
const int* const Matrix::get_aa2i()         const { return aa2i;                  }

void Matrix::reset() {
    for(size_t i=0; i<20; ++i)
        for(size_t j=0; j<20; ++j)
            matrix[i][j] = log2( original[i][j]/(p_background[i]*p_background[j]) ) / scale;
    for(size_t i=0; i<20; ++i)
        for(size_t j=0; j<20; ++j)
            p_cond[i][j] = original[i][j]/p_background[j];    
}

void Matrix::round_bit_scores(){
    for(size_t i=0; i<20; ++i)
        for(size_t j=0; j<20; ++j)
            matrix[i][j] = round( matrix[i][j]);
    
}

void Matrix::_copy(){
    for(size_t i=0; i<21; ++i)
        for(size_t j=0; j<21; ++j)
            matrix[i][j] = original[i][j];
}

std::ostream& Matrix::print(std::ostream &out){
    out << "Matrix in bits (float precision), scaling:" << scale << std::endl;
    char buf[21];
    out << " ";
    for(size_t i=0; i<21; ++i){
        sprintf(buf, "   %c  ", i2aa[i]);
        out << buf;
    }
    out << std::endl;
    for(size_t i=0; i<21; ++i){
        for(size_t j=0; j<21; ++j){
            sprintf(buf, "%+2.2f ", matrix[i][j]);
            out << buf;
        }
        out << std::endl;
	}
    return out;
}

std::ostream& Matrix::print_frequencies(std::ostream &out){
    out << "{" << std::endl;
    for(size_t i=0; i<21; ++i){
        out << "{ ";
        for(size_t j=0; j<21; ++j){
            out << original[i][j];
            if( j!=20 ) out << ", ";
            else out << "}";
        }
        if(i!=20) out << "," << std::endl;
        else      out << "}" << std::endl;
    }
    return out;
}

std::ostream& Matrix::print_p_back(std::ostream &out){
    for(size_t j=0; j<21; ++j)
		out << j << " " << i2aa[j] << " " << p_background[j] << std::endl;
    return out;
}

std::ostream& Matrix::print_int(std::ostream &out){
    out << "Matrix in bits (float precision), scaling:" << scale << std::endl;
    char buf[21];
    out << " ";
    for(size_t i=0; i<21; ++i){
        sprintf(buf, "   %c", i2aa[i]);
        out << buf;
    }
    out << std::endl;
    for(size_t i=0; i<21; ++i){
        sprintf(buf, "%c ", i2aa[i]);
        out << buf;
        for(size_t j=0; j<21; ++j){
            sprintf(buf, "%+3i ", (int)round(matrix[i][j]));
            out << buf;
        }
        out << std::endl;
    }
    return out;
}

const std::string Matrix::_get_fn(const Matrix::mtype t){
    std::string p="./";
    if(getenv("BLOSUM_MATRICES")) p = std::string(getenv("BLOSUM_MATRICES"));
    switch(t){
    case Matrix::blosum30: return (p+"/blosum30.out");
    case Matrix::blosum35: return (p+"/blosum35.out");
    case Matrix::blosum40: return (p+"/blosum40.out");
    case Matrix::blosum45: return (p+"/blosum45.out");
    case Matrix::blosum50: return (p+"/blosum50.out");
    case Matrix::blosum55: return (p+"/blosum55.out");
    case Matrix::blosum60: return (p+"/blosum60.out");
    case Matrix::blosum62: return (p+"/blosum62.out");
    case Matrix::blosum65: return (p+"/blosum65.out");
    case Matrix::blosum70: return (p+"/blosum70.out");
    case Matrix::blosum75: return (p+"/blosum75.out");
    case Matrix::blosum80: return (p+"/blosum80.out");
    case Matrix::blosum85: return (p+"/blosum85.out");
    case Matrix::blosum90: return (p+"/blosum90.out");
    case Matrix::blosum95: return (p+"/blosum95.out");
    case Matrix::blosum100: return (p+"/blosum100.out");
    case Matrix::static_blosum62 : return ("no-file");
    }
    return "";
}

void Matrix::_read_blosum_matrix(const char *fn) throw (std::exception){
    std::ifstream in(fn);
    if( in.fail() ) 
        throw MyException("Cannot read '%s'!\nPlease export $BLOSUM_MATRICES=/dir-where-the-blosumXX.out-files-reside/", fn);
    int c      = 0;
    int row    = 0;
    int column = 0;
    std::string line;
    bool capture = false;
    while( in.good() ){
        getline( in, line );
        if( line.length()>11 && line.substr(0, 11)!="Frequencies" && !capture ) continue;
        if( line.length()>11 && line.substr(0, 11)=="Frequencies"){
            capture=true;
            continue;
        }
        if( row==20 ) break;
        std::stringstream stream(line); std::string h; stream >> h;
        if( h=="" ) continue;
        if( isalpha(h.at(0)) ){
            i2aa[c++] = toupper( h.at(0) );
            while(	stream >> h ){
                i2aa[c++] = toupper( h.at(0) );
                if( c>20 ) throw MyException("Blosum matrix file '%s' has wrong format!\n", fn);
            }
        }else{
            column = 0;
            stream.clear();
            stream.str(line);
            float f;
            while( stream >> f ){
                original[row][column] = f;
                original[column][row] = f;
                ++column;
            }
            ++row;
        }
    }
    if( c!=20 ) throw MyException("Blosum matrix file '%s' has wrong format!\n", fn);
    in.close();
    
    float sum=0.0f;
    for(size_t i=0; i<20; ++i)
        for(size_t j=0; j<20; ++j){
            if( i==j ) p_background[i] += original[i][j];
            else       p_background[i] += (original[i][j]/2.0f);
            if( j<=i ) sum += original[i][j]; 
        }
    
    const float _2sum = 2.0*sum;	
    float pbsum = 0.0f;
    for(size_t i=0; i<20; ++i){
        pbsum += p_background[i];
        for(size_t j=0; j<20; ++j)
            if( i==j ) original[i][j] = original[i][j] / sum;
            else       original[i][j] = original[i][j] / _2sum;
    }
    
    for(size_t i=0; i<20; ++i)p_background[i] /= sum;
    
    float entropy=0.0f;
    //compute entropy	
    for(size_t i=0; i<20; ++i)
        for(size_t j=0; j<=i; ++j)
            entropy += original[i][j] * log2( original[i][j] / (p_background[i]*p_background[j]) ) ;
    //set scaling factor for blosum matrices half-bits, third-bits,...
    //std::cerr << "Entropy:" << entropy << std::endl;
    scale = std::min( 0.5, 1.0/(round(2.0/sqrt(entropy))) );
}

std::ostream& Matrix::print_debug(std::ostream &out){
    char buf[21];
    out << "Frequencies of amino acid substitutions (unscaled, float precision):" << std::endl;
    for(size_t i=0; i<21; ++i){
        sprintf(buf, "%6c ", i2aa[i]);
        out << buf;
    }
    out << std::endl;
    for(size_t i=0; i<21; ++i){
        for(size_t j=0; j<21; ++j){
            sprintf(buf, "%+2.3f ", original[i][j]);
            out << buf;
        }
        out << std::endl;
    }
    out << std::endl;
    out << "Matrix in bits (float precision), scaling:" << scale << std::endl;
    for(size_t i=0; i<21; ++i){
        sprintf(buf, "%6c ", i2aa[i]);
        out << buf;
    }
    out << std::endl;
    for(size_t i=0; i<21; ++i){
        for(size_t j=0; j<21; ++j){
            sprintf(buf, "%+2.3f ", matrix[i][j]);
            out << buf;
        }
        out << std::endl;
    }
    out << std::endl;
    out << "Matrix of conditional probabilities P(a|b)" << std::endl;
    for(size_t i=0; i<21; ++i){
        sprintf(buf, "%6c ", i2aa[i]);
        out << buf;
    }
    out << std::endl;
    for(size_t i=0; i<21; ++i){
        for(size_t j=0; j<21; ++j){
            sprintf(buf, "%+2.3f ", p_cond[i][j]);
            out << buf;
        }
        out << std::endl;
    }
    out << std::endl;
    return out;
}

void Matrix::read_matrix(const float src[][21]){
    for(size_t i=0; i<21; ++i)
        for(size_t j=0; j<21; ++j)
            matrix[i][j] = src[i][j];
}
