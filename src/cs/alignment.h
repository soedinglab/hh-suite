/*
  Copyright 2009 Andreas Biegert

  This file is part of the CS-BLAST package.

  The CS-BLAST package is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  The CS-BLAST package is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef CS_ALIGNMENT_H_
#define CS_ALIGNMENT_H_

#include <valarray>
#include <vector>

#include "sequence.h"
#include "globals.h"
#include "matrix.h"
#include "vector.h"

namespace cs {

// Forward declarations
template<class Abc>
class Alignment;

// Convince the compiler that operator<< is a template friend.
template<class Abc>
std::ostream& operator<< (std::ostream& out, const Alignment<Abc>& ali);

// Supported alignment formats for in- and output.
enum AlignmentFormat {
    FASTA_ALIGNMENT    = 0,
    A2M_ALIGNMENT      = 1,
    A3M_ALIGNMENT      = 2,
    CLUSTAL_ALIGNMENT  = 3,
    PSI_ALIGNMENT      = 4
};

// A container class for multiple sequence alignments.
template<class Abc>
class Alignment {
	private:
	  // Row major matrix with sequences in integer representation.
	  Matrix<uint8_t> seqs_;
	  // Array with indices of all columns [0,1,2,...,num_cols-1].
	  std::valarray<size_t> col_idx_;
	  // Array with indices of match columns.
	  std::valarray<size_t> match_idx_;
	  // Array mask indicating match and insert columns.
	  std::valarray<bool> is_match_;
	  // Headers of sequences in the alignment.
	  std::vector<std::string> headers_;
	  // Name of the alignment as given by comment line in FASTA file
	  std::string name_;

	  // Buffer size for reading
	  static const size_t kBufferSize = MB;

	  // Initializes alignment with given headers and sequences.
	  void Init(const std::vector<std::string>& headers, const std::vector<std::string>& seqs);

	  // Resize the sequence matrix and header vector to given dimensions.
	  void Resize(size_t num_seqs, size_t num_cols);

	  // Fills match_idx__ with the indices of all match columns.
	  void SetMatchIndices();

    // Removes sequences with headers indicating non-protein sequences (secondary structure predictions)
    void FilterSequencesByHeaders(std::vector<std::string>& headers, std::vector<std::string>& seqs);

	  // Reads an alignment in FASTA format.
	  void ReadFasta(FILE* fin, std::vector<std::string>& headers, std::vector<std::string>& seqs);

	  // Reads an alignment in A2M format from given stream.
	  void ReadA2M(FILE* fin, std::vector<std::string>& headers, std::vector<std::string>& seqs);

	  // Reads an alignment in A3M format from given stream.
	  void ReadA3M(FILE* fin, std::vector<std::string>& headers, std::vector<std::string>& seqs);

	  // Helper method that reads a FASTA, A2M, or A3M formatted alignment.
	  void ReadFastaFlavors(FILE* fin, std::vector<std::string>& headers, std::vector<std::string>& seqs);

	  // Reads an alignment in PSI format.
	  void ReadPsi(FILE* fin, std::vector<std::string>& headers, std::vector<std::string>& seqs);

	  // Writes the alignment in FASTA, A2M, or A3M format to output stream.
	  void WriteFastaFlavors(FILE* fout, AlignmentFormat format, size_t width = 100) const;

	  // Writes the alignment in CLUSTAL or PSI format to output stream.
	  void WriteClustalFlavors(FILE* fout, AlignmentFormat format, size_t width = 100) const;

  public:
    // Constructs alignment from multi FASTA formatted alignment read from input
    // stream.
    Alignment(FILE* fin, AlignmentFormat format);

    // Constructs an all-gaps alignment with 'ncols' columns and 'nseqs' sequences
    Alignment(size_t ncols, size_t nseqs);

    // Constructs an alignment from a single sequence.
    Alignment(const Sequence<Abc>& seq);

    // All memeber are automatically destructed
    ~Alignment() {}

    // Accessors for integer representation of character in MATCH column i of
    // sequence k.
    uint8_t* operator[](size_t i) { return seqs_[match_idx_[i]]; }
    const uint8_t* operator[](size_t i) const { return seqs_[match_idx_[i]]; }
    uint8_t& match(size_t i, size_t k) { return seqs_[match_idx_[i]][k]; }
    const uint8_t& match(size_t i, size_t k) const { return seqs_[match_idx_[i]][k]; }

    // Accessors for integer representation of the character at position i
    // (NOT match column i) of sequence k.
    uint8_t& operator() (size_t k, size_t i) { return seqs_[i][k]; }
    const uint8_t& operator() (size_t k, size_t i) const { return seqs_[i][k]; }
    uint8_t& seq(size_t k, size_t i) { return seqs_[i][k]; }
    const uint8_t& seq(size_t k, size_t i) const { return seqs_[i][k]; }

    // Returns the character at position i (NOT match column i) of sequence k.
    char chr(size_t k, size_t i) const { return Abc::kIntToChar[seqs_[i][k]]; }

    // Returns index of all-alignment column corresponding to match column i.
    size_t col_idx(size_t i) const { return match_idx_[i]; }

    // Returns the number of sequences in the alignment.
    size_t nseqs() const { return seqs_.ncols(); }

    // Returns the total number of alignment columns.
    size_t ncols() const { return seqs_.nrows(); }

    // Returns the number of match columns.
    size_t nmatch() const { return match_idx_.size(); }

    // Returns the number of insert columns.
    int ninsert() const { return ncols() - nmatch(); }

    // Returns the header of sequence k.
    std::string header(size_t k) const { return headers_[k]; }

    // Sets the header of sequence k.
    void set_header(size_t k, const std::string& header) { headers_[k] = header; }

    // Returns the name of the alignment
    std::string name() const { return name_; }

    // Sets the name of the alignment
//    void set_name(const std::string& header) { name_ = name; }

    // Makes all columns with a residue in sequence k match columns.
    void AssignMatchColumnsBySequence(size_t k = 0);

    // Makes all columns with less than X% gaps match columns.
    void AssignMatchColumnsByGapRule(double gap_threshold = 50);

    // Initializes object with an alignment in FASTA format read from given
    // stream.
    void Read(FILE* fin, AlignmentFormat format);

    // Writes the alignment in given format to ouput stream.
    void Write(FILE* fout, AlignmentFormat format, size_t width = 100) const;

    // Returns true if column i is a match column.
    bool is_match(size_t i) const { return is_match_[i]; }

    // Removes all insert columns from the alignment.
    void RemoveInsertColumns();

    // Merges the provided alignment with this alignment considering only sequences
    // that are not already included in this alignment. Warning: Inserts in current
    // alignment are lost!
    void Merge(const Alignment<Abc>& ali);

    // Returns alignment sequence k as Sequence object without gaps.
    Sequence<Abc> GetSequence(size_t k) const;

    // Rearranges the aligned sequences such that they appear in the same order as
    // their unaligned counterparts in the given sequence vector.
    void Rearrange(const std::vector<Sequence<Abc> >& seqs);

    // Prints the Alignment in A2M format for debugging.
    friend std::ostream& operator<< <> (std::ostream& out, const Alignment<Abc>& ali);

};  // Alignment


// Returns the alignment format corresponding to provided filename extension
inline AlignmentFormat AlignmentFormatFromString(const std::string& s);

// Reads all available alignments from the input stream and puts them into vector.
template<class Abc>
static void ReadAll(FILE* fin, AlignmentFormat format, std::vector< Alignment<Abc> >& v);

// Calculates global sequence weights by maximum entropy weighting
// (Henikoff&Henikoff '94).
template<class Abc>
double GlobalWeightsAndDiversity(const Alignment<Abc>& ali, Vector<double>& wg, bool neff_sum_pairs = false);

// Calculates position-dependent sequence weights and number of effective
// sequences on subalignments.
template<class Abc>
Vector<double> PositionSpecificWeightsAndDiversity(const Alignment<Abc>& ali, Matrix<double>& w);

// Converts a character to uppercase and '.' to '-'.
inline char to_match_chr(char c) {
    return isalpha(c) ? toupper(c) : (c == '.' ? '-' : c);
}

// Converts a character to lowercase and '-' to '.'.
inline char to_insert_chr(char c) {
    return isalpha(c) ? tolower(c) : (c == '-' ? '.' : c);
}

// Predicate indicating if character belongs to match column.
inline bool match_chr(char c) {
    return (isalpha(c) && isupper(c)) || c == '-';
}

// Predicate indicating if character belongs to insert column.
inline char insert_chr(char c) {
    return (isalpha(c) && islower(c)) || c == '.';
}

}  // namespace cs

#endif  // CS_ALIGNMENT_H_
