/***************************************************************************
 *   Copyright (C) 2006 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/
#include "sequence.h"

Sequence::Sequence(	char* header, 
			char* sequence, 
			const int*  const aa2int, 
			const char* const int2aa, 
			size_t ADIM, size_t index)throw (std::exception):
			NAA(ADIM),
			aa2i(aa2int), 
			i2aa(int2aa),
			len(0),
			idx(index)
{
    const size_t hlen = strlen(header);
    const size_t slen = strlen(sequence);
    init(header, hlen, sequence, slen);
}

Sequence::Sequence( const char* const header, 
                    const char* const sequence, 
                    AminoAcid *aa, 
                    size_t index) throw (std::exception):
    NAA(aa->get_naa()),
    aa2i(aa->get_aa2i()), 
    i2aa(aa->get_i2aa()),
    len(0),
    idx(index) 
{
    const size_t hlen = strlen(header);
    const size_t slen = strlen(sequence);
    init(header, hlen, sequence, slen);
}

Sequence::Sequence( const char* fn, 
                    const Sequence::format &f, 
                    AminoAcid *aa, 
                    size_t index) throw (std::exception):
    NAA(aa->get_naa()),
    aa2i(aa->get_aa2i()), 
    i2aa(aa->get_i2aa()),
    len(0),
    idx(index) 
{
    switch(f){
    case Sequence::fasta : read_fasta(fn); break;
    default              : throw MyException("Unknown sequence format!\n");
    }    
}

Sequence::Sequence( const size_t length,
                    AminoAcid *aa,
                    Matrix *m,  
                    const size_t index ) throw (std::exception):
    NAA(aa->get_naa()),
    aa2i(aa->get_aa2i()), 
    i2aa(aa->get_i2aa()),
    len(length),
    idx(index)
{
    header = new char[60];
    seq    = new char[len];
    sprintf(header, ">random, length: %zu", len);
    
    float aa_freq[NAA];
    float sum=0.0;
    for (size_t i=0; i<NAA; ++i ){
        sum += m->get_p_back()[i];
        aa_freq[i] = sum;
    }
    
    for(size_t i=0; i<len; ++i ){
        int aa=0;
        float r = frand();            
        for (size_t a=0; a<NAA; ++a ){
            aa=a;
            if ( r<aa_freq[a] ) break;
        }
        seq[i] = aa;
    }
}

void Sequence::set_index(const size_t index) { idx = index; }

void Sequence::read_fasta(const char *fn) throw (std::exception) 
{
    char *hbuffer = new char[READ_HBUFFER_SIZE];
    char *sbuffer = new char[READ_SBUFFER_SIZE];
    std::ifstream in(fn);
    if( in.fail() ) throw MyException("Cannot read '%s'", fn);
    // move to first header line 
    while( !in.eof() ){
        char c = in.get();
        if( c=='>' ){ in.unget(); break; }	
        if( isspace(c) ) continue;
        throw MyException("File '%s' is not in correct FastA format!", fn);
    }
    //read header and sequence
    if(!in.eof()){	
        //extracts with appended '\0'
        in.getline(hbuffer, READ_HBUFFER_SIZE, '\n');
        size_t hb = in.gcount();
        if( hb==(READ_HBUFFER_SIZE-1) ) 
            throw MyException("Error reading '%s'\nLong header line found, increase READ_HBUFFER_SIZE!", fn);	
        if( hbuffer[0] == '>' ){
            in.get(sbuffer, READ_SBUFFER_SIZE, '>');
            size_t sb = in.gcount();
            if( sb==(READ_SBUFFER_SIZE-1) ) 
                throw MyException("Error reading '%s'\nLong sequence found, increase READ_SBUFFER_SIZE!", fn);
            if( sb==0 ) 
                throw MyException("File '%s' is not in correct FastA format!", fn);
            init(hbuffer+1, hb-1, sbuffer, sb);
        }else throw MyException("File '%s' is not in correct FastA format!", fn);
    }else         throw MyException("File '%s' is not in correct FastA format!", fn);
    
    delete [] hbuffer;
    delete [] sbuffer;
}


Sequence::~Sequence()
{
    delete [] header;
    delete [] seq;
}

void Sequence::init(const char* const head, const size_t hl, const char* const sequence, const size_t sl) throw (std::exception)
{    
    header = new char[hl+1];
    memcpy(header, head, hl+1);
    hlen = hl;
    
    len = Sequence::get_length(sequence, sl);
    if (len==0)throw MyException("No sequence data found for '%s'!\n", header);
    seq = new char[len];
    
    char *ptr = seq;
    for(size_t i=0; i<sl; ++i){
        if( isalpha(sequence[i]) ){ 
            const char cur = toupper(sequence[i]);
            switch(cur){
            case 'A': *ptr = aa2i['A']; break;
            case 'B': *ptr = aa2i['N']; break;
            case 'C': *ptr = aa2i['C']; break;
            case 'D': *ptr = aa2i['D']; break;
            case 'E': *ptr = aa2i['E']; break;
            case 'F': *ptr = aa2i['F']; break;
            case 'G': *ptr = aa2i['G']; break;
            case 'H': *ptr = aa2i['H']; break;
            case 'I': *ptr = aa2i['I']; break;
            case 'J': *ptr = aa2i['L']; break;
                
            case 'K': *ptr = aa2i['K']; break;
            case 'L': *ptr = aa2i['L']; break;
            case 'M': *ptr = aa2i['M']; break;
            case 'N': *ptr = aa2i['N']; break;
            case 'P': *ptr = aa2i['P']; break;
            case 'Q': *ptr = aa2i['Q']; break;
            case 'R': *ptr = aa2i['R']; break;
            case 'S': *ptr = aa2i['S']; break;
            case 'T': *ptr = aa2i['T']; break;
            case 'U': *ptr = aa2i['C']; break;
            case 'V': *ptr = aa2i['V']; break;
            case 'W': *ptr = aa2i['W']; break;
            case 'X': *ptr = aa2i['X']; break;
            case 'Y': *ptr = aa2i['Y']; break;
            case 'Z': *ptr = aa2i['E']; break;
                                    
            default: 
                std::cerr << "Warning: Ignoring invalid character '";
                std::cerr << sequence[i] << "' at position " << i;
                std::cerr << " in sequence " << header << std::endl;
                continue;
            }
            ++ptr;
        }else if( !isspace(sequence[i]) ){
            std::cerr << "Warning: Ignoring invalid character '";
            std::cerr << sequence[i] << "' at position " << i;
            std::cerr << " in sequence " << header << std::endl;
        }
    }
}

size_t Sequence::get_length( const char* const ptr, const size_t l)
{
    size_t ret = 0;
    for( size_t i=0; i<l; ++i )if( isalpha(ptr[i]) ) ++ret;
    return ret;
}

std::ostream& Sequence::print_debug(std::ostream &out)
{
    out << "INDEX     : " << idx       << std::endl;
    out << "LENGTH    : " << len       << std::endl;
    out << "HEADER    : " << header   << "$" << std::endl;
    out << "NUMERIC   : ";
    for(size_t i=0; i<len; ++i) out << std::setw(2) << (int)seq[i] << " ";
    out << "$" << std::endl;
    out << "SEQUENCE  : ";
    for(size_t i=0; i<len; ++i) out << std::setw(2) << i2aa[(int)seq[i]] << " ";
    out << "$" << std::endl;
    return out;
}

std::ostream& Sequence::write(std::ostream &out, const size_t width) 
{
    out << '>' << header << std::endl;
    for(size_t i=0; i<len; i+=width) {
        for(size_t j=i; j<std::min(i+width, len); ++j)out << i2aa[(int)seq[j]];	
        out << std::endl;
    }    
    return out;
}

std::ostream& Sequence::write_sequence(std::ostream &out, const size_t width)
{
    for(size_t i=0; i<len; i+=width) {
        for(size_t j=i; j<std::min(i+width, len); ++j) out << i2aa[(int)seq[j]];	
        out << std::endl;
    }    
    return out;
}

void Sequence::replace_J()
{
    const size_t slen = strlen(seq);
    for( size_t i=0; i<slen; ++i){ 
        if( seq[i]=='J' )      seq[i] = 'L';
        else if( seq[i]=='j' ) seq[i] = 'l';
    }
}
