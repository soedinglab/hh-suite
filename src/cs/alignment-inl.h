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

#ifndef CS_ALIGNMENT_INL_H_
#define CS_ALIGNMENT_INL_H_

#include "alignment.h"

#include "blast_hits.h"
#include "sequence-inl.h"

namespace cs {

template<class Abc>
Alignment<Abc>::Alignment(FILE* fin, AlignmentFormat format) {
    Read(fin, format);
}

template<class Abc>
Alignment<Abc>::Alignment(size_t ncols, size_t nseqs) {
    // Fill alignment matrix with gaps and assign empty headers
    Resize(nseqs, ncols);
    for (size_t k = 0; k < nseqs; ++k)
        headers_[k] = "";
    for (size_t i = 0; i < ncols; ++i) {
        col_idx_[i]  = i;
        is_match_[i] = true;
        for (size_t k = 0; k < nseqs; ++k)
            seqs_[i][k]  = Abc::kGap;
    }
    // Initialize index array for match columns
    SetMatchIndices();
}

template<class Abc>
Alignment<Abc>::Alignment(const Sequence<Abc>& seq) {
    std::vector<std::string> headers;
    std::vector<std::string> seqs;
    headers.push_back(seq.header());
    seqs.push_back(seq.ToString());
    Init(headers, seqs);
}

template<class Abc>
Alignment<Abc>::Alignment(const BlastHits& hits, bool best) {
    std::vector<std::string> headers;
    std::vector<std::string> seqs;

    typedef typename BlastHits::ConstHitIter HitIter;
    typedef typename BlastHits::ConstHspIter HspIter;
    for (HitIter hit = hits.begin(); hit != hits.end(); ++hit) {
        for (HspIter hsp = hit->hsps.begin(); hsp != hit->hsps.end(); ++hsp) {
            // Construct query anchored alignment string
            std::string seq(hsp->query_start - 1, '-');
            for (size_t i =  0; i < hsp->length; ++i)
                if (hsp->query_seq[i] != '-')
                    seq += hsp->subject_seq[i];
            seq.append(hits.query_length() - seq.length(), '-');

            headers.push_back(hit->definition);
            seqs.push_back(seq);

            if (best) break;
        }
    }
    Init(headers, seqs);
}

template<class Abc>
void Alignment<Abc>::Init(const std::vector<std::string>& headers,
                          const std::vector<std::string>& seqs) {
    if (headers.size() != seqs.size())
        throw Exception("Bad alignment: unequal number of headers and sequences!");

    const size_t nseqs = seqs.size();
    const size_t ncols = seqs[0].length();
    for (size_t k = 1; k < nseqs; ++k) {
        if (seqs[k].length() != ncols)
            throw Exception("Alignment sequence %i has length %i but should have %i!",
                            k+1, seqs[k].length(), ncols);
    }

    // Validate characters and convert to integer representation
    Resize(seqs.size(), seqs[0].length());
    for (size_t k = 0; k < nseqs; ++k) {
        headers_[k] = headers[k];
        for (size_t i = 0; i < ncols; ++i) {
            const char c = seqs[k][i];
            col_idx_[i] = i;
            is_match_[i] = match_chr(c);
            if (Abc::kValidChar[static_cast<int>(c)] || c == '-' || c == '.')
                seqs_[i][k] = Abc::kCharToInt[static_cast<int>(c)];
            else
                throw Exception("Invalid character %c at position %i of sequence '%s'",
                                c, i, seqs[k].c_str());
        }
    }

    // Replace gap with endgap for all gaps at either end of a sequence
    for (size_t k = 0; k < nseqs; ++k) {
        for (size_t i = 0; i < ncols && seqs_[i][k] == Abc::kGap; ++i)
            seqs_[i][k] = Abc::kEndGap;
        for (int i = ncols - 1; i >= 0 && seqs_[i][k] == Abc::kGap; --i)
            seqs_[i][k] = Abc::kEndGap;
    }

    // Initialize index array for match columns
    SetMatchIndices();
}

template<class Abc>
void Alignment<Abc>::SetMatchIndices() {
    const size_t match_cols = std::count(&is_match_[0], &is_match_[0] + ncols(), true);
    match_idx_.resize(match_cols);
    match_idx_ = col_idx_[is_match_];
}

template<class Abc>
void Alignment<Abc>::Read(FILE* fin, AlignmentFormat format) {
    LOG(DEBUG4) << "Reading alignment from stream ...";

    std::vector<std::string> headers;
    std::vector<std::string> seqs;
    switch (format) {
        case FASTA_ALIGNMENT:
            ReadFasta(fin, headers, seqs);
            break;
        case A2M_ALIGNMENT:
            ReadA2M(fin, headers, seqs);
            break;
        case A3M_ALIGNMENT:
            ReadA3M(fin, headers, seqs);
            break;
        case PSI_ALIGNMENT:
            ReadPsi(fin, headers, seqs);
            break;
        default:
            throw Exception("Unsupported alignment input format %i!", format);
    }
    Init(headers, seqs);

    LOG(DEBUG4) << *this;
}

template<class Abc>
void Alignment<Abc>::ReadFastaFlavors(FILE* fin, std::vector<std::string>& headers, std::vector<std::string>& seqs) {
    headers.clear();
    seqs.clear();

    char buffer[kBufferSize];
    int c = '\0';
    while (!feof(fin)) {
        // Read header
        while (cs::fgetline(buffer, kBufferSize, fin)) {
            if (!strscn(buffer)) continue;
            if (buffer[0] == '#') {
                name_ = std::string(buffer + 1);
            } else if (buffer[0] == '>') {
                if (headers.empty() && 
                    (strstr(buffer, ">ss_") == buffer || strstr(buffer, ">sa_") == buffer)) {
                  while (!feof(fin)) {
                    c = getc(fin);
                    ungetc(c, fin);
                    if (c == '>') break;
                    cs::fgetline(buffer, kBufferSize, fin);
                  }
                  continue;
                }
                headers.push_back(std::string(buffer + 1));
                break;
            } else {
                throw Exception("Header of sequence %i starts with:\n%s",
                                headers.size() + 1, buffer);
            }
        }

        // Read sequence
        seqs.push_back("");
        while (cs::fgetline(buffer, kBufferSize, fin)) {
            seqs.back().append(buffer);

            c = getc(fin);
            if (c == EOF) break;
            ungetc(c, fin);
            if (static_cast<char>(c) == '>') break;
        }
        // Remove whitespace
        seqs.back().erase(remove_if(seqs.back().begin(), seqs.back().end(), isspace), seqs.back().end());

        LOG(DEBUG2) << headers.back();
    }
    //LOG(DEBUG2) << "Number of sequences read: " << headers.size() << std::endl;
    if (headers.empty())
        throw Exception("Bad alignment: no alignment data found in stream!");
}

template<class Abc>
void Alignment<Abc>::ReadPsi(FILE* fin, std::vector<std::string>& headers, std::vector<std::string>& seqs) {
    headers.clear();
    seqs.clear();

    char buffer[kBufferSize];
    const char* ptr;
    size_t block = 1;  // number of block currently read
    size_t n = 0;      // sequence number of first block
    size_t k = 0;      // sequence index to zero for first block

    while (cs::fgetline(buffer, kBufferSize, fin)) {
        if (buffer[0] == '/' && buffer[1] == '/') break;

        // Start of new block
        if (!strscn(buffer)) {
            if (k > 0) {
                if (n > 0 && n != k)
                    throw Exception("Error: different number of sequences in blocks "
                                    "1 and %i.", block);
                ++block;
                n = k;
                k = 0;
            }
            continue;
        }
        // Parse name and sequence characters
        ptr = strchr(buffer, ' ');
        if (!ptr)
            throw Exception("Error: Missing whitespace between identifier and "
                            "sequence characters of sequence %i in block %i.",
                            k+1, block);
        std::string name(buffer, ptr - buffer);
        ptr = strscn(ptr);
        if (!ptr)
            throw Exception("Error: Missing sequence characters in sequence %i of "
                            "block %i.", k+1, block);
        std::string seq(ptr);
        if (block == 1) {
            headers.push_back(name);
            seqs.push_back(seq);
        } else {
            assert(name == headers[k]);
            seqs[k].append(seq);
        }
        ++k;
    }
    if (k > 0 && n > 0 && n != k)
        throw Exception("Error: different number of sequences in blocks "
                        "1 and %i.", block);
}

template<class Abc>
void Alignment<Abc>::ReadFasta(FILE* fin,
                               std::vector<std::string>& headers,
                               std::vector<std::string>& seqs) {
    ReadFastaFlavors(fin, headers, seqs);
    // Convert all characters to match characters
    for (std::vector<std::string>::iterator it = seqs.begin();
         it != seqs.end(); ++it)
        transform(it->begin(), it->end(), it->begin(), to_match_chr);
}

template<class Abc>
void Alignment<Abc>::ReadA2M(FILE* fin, std::vector<std::string>& headers, std::vector<std::string>& seqs) {
    ReadFastaFlavors(fin, headers, seqs);
}

template<class Abc>
void Alignment<Abc>::ReadA3M(FILE* fin, std::vector<std::string>& headers, std::vector<std::string>& seqs) {
    ReadFastaFlavors(fin, headers, seqs);

    // Check number of match states
    const size_t nseqs = seqs.size();
    const size_t nmatch_cols = count_if(seqs[0].begin(), seqs[0].end(), match_chr);
    for (size_t k = 1; k < nseqs; ++k) {
        const size_t nmatch_cols_k = count_if(seqs[k].begin(), seqs[k].end(), match_chr);
        if (nmatch_cols_k != nmatch_cols)
            throw Exception("Sequence %i has %i match columns but should have %i!",
                            k, nmatch_cols_k, nmatch_cols);
        if (count(seqs[k].begin(), seqs[k].end(), '.') > 0)
            throw Exception("Sequence %i in A3M alignment contains gaps!", k);
    }

    // Insert gaps into A3M alignment
    std::vector<std::string> seqs_a2m(seqs.size(), "");
    Matrix<std::string> inserts(seqs.size(), nmatch_cols + 1, "");
    std::vector<size_t> max_insert_len(nmatch_cols + 1, 0);

    // Move inserts before first match state into seqs_a2m and keep track of
    // longest first insert
    size_t max_first_insert_len = 0;
    for (size_t k = 0; k < nseqs; ++k) {
        std::string::iterator i =
            find_if(seqs[k].begin(), seqs[k].end(), match_chr);
        if (i != seqs[k].end()) {
            seqs_a2m[k].append(seqs[k].begin(), i);
            seqs[k].erase(seqs[k].begin(), i);
        }
        max_first_insert_len = MAX(max_first_insert_len,
                                   static_cast<size_t>(i - seqs[k].begin()));
    }

    // Extract all inserts and keep track of longest insert after each match
    // column
    for (size_t k = 0; k < nseqs; ++k) {
        int i = -1;
        std::string::iterator insert_end = seqs[k].begin();
        std::string::iterator insert_start =
            find_if(insert_end, seqs[k].end(), insert_chr);
        while (insert_start != seqs[k].end() && insert_end != seqs[k].end()) {
            i += insert_start - insert_end;
            insert_end = find_if(insert_start, seqs[k].end(), match_chr);
            inserts[k][i] = std::string(insert_start, insert_end);
            max_insert_len[i] = MAX(inserts[k][i].length(), max_insert_len[i]);
            insert_start = find_if(insert_end, seqs[k].end(), insert_chr);
        }
    }

    // Build new A2M alignment
    for (size_t k = 0; k < nseqs; ++k) {
        seqs_a2m[k].append(max_first_insert_len - seqs_a2m[k].length(), '.');
        int i = 0;
        std::string::iterator match = seqs[k].begin();
        while (match != seqs[k].end()) {
            seqs_a2m[k].append(1, *match);
            if (max_insert_len[i] > 0) {
                seqs_a2m[k].append(inserts[k][i]);
                seqs_a2m[k].append(max_insert_len[i] - inserts[k][i].length(), '.');
            }
            match = find_if(match+1, seqs[k].end(), match_chr);
            ++i;
        }
    }

    // Overwrite original A3M alignment with new A2M alignment
    seqs = seqs_a2m;
}

template<class Abc>
void Alignment<Abc>::Write(FILE* fout, AlignmentFormat format, size_t width) const {
    switch (format) {
        case FASTA_ALIGNMENT:
        case A2M_ALIGNMENT:
        case A3M_ALIGNMENT:
            WriteFastaFlavors(fout, format, width);
            break;
        case CLUSTAL_ALIGNMENT:
        case PSI_ALIGNMENT:
            WriteClustalFlavors(fout, format, width);
            break;
        default:
            throw Exception("Unsupported alignment output format %i!", format);
    }
}

template<class Abc>
void Alignment<Abc>::WriteFastaFlavors(FILE* fout, AlignmentFormat format, size_t width) const {
    for (size_t k = 0; k < nseqs(); ++k) {
        fprintf(fout, ">%s\n", headers_[k].c_str());
        size_t j = 0;  // counts printed characters
        for (size_t i = 0; i < ncols(); ++i) {
            switch (format) {
                case FASTA_ALIGNMENT:
                    fputc(to_match_chr(chr(k, i)), fout);
                    ++j;
                    break;
                case A2M_ALIGNMENT:
                    if (is_match_[i])
                        fputc(to_match_chr(chr(k, i)), fout);
                    else
                        fputc(to_insert_chr(chr(k, i)), fout);
                    ++j;
                    break;
                case A3M_ALIGNMENT:
                    if (is_match_[i]) {
                        fputc(to_match_chr(chr(k, i)), fout);
                        ++j;
                    } else if (seq(k, i) != Abc::kGap
                               && seq(k, i) != Abc::kEndGap) {
                        fputc(to_insert_chr(chr(k, i)), fout);
                        ++j;
                    }
                    break;
                default:
                    throw Exception("Unsupported alignment output format %i!", format);
            }
            if (j % width == 0) fputc('\n', fout);
        }
        if (j % width != 0) putc('\n', fout);
    }
}

template<class Abc>
void Alignment<Abc>::WriteClustalFlavors(FILE* fout, AlignmentFormat format, size_t width) const {
    const size_t kHeaderWidth = 15;

    // Convert alignment to character representation
    std::vector<std::string> seqs(nseqs(), "");
    for (size_t k = 0; k < nseqs(); ++k) {
        for (size_t i = 0; i < ncols(); ++i) {
            char c = Abc::kIntToChar[seqs_[i][k]];
            if (c != '-' && !is_match_[i] && format == PSI_ALIGNMENT)
                c = to_insert_chr(c);
            seqs[k].push_back(c);
        }
    }

    // Print alignment in blocks
    if (format == CLUSTAL_ALIGNMENT) fputs("CLUSTAL\n\n", fout);
    while (!seqs.front().empty()) {
        for (size_t k = 0; k < nseqs(); ++k) {
            std::string header(headers_[k].substr(0, headers_[k].find_first_of(' ')));
            if (header.length() <= kHeaderWidth) {
                header += std::string(kHeaderWidth - header.length() + 1, ' ');
                fputs(header.c_str(), fout);
                fputc(' ', fout);  // separator between header and sequence
            } else {
                fputs(header.substr(0, kHeaderWidth).c_str(), fout);
                fputc(' ', fout);  // separator between header and sequence
            }

            size_t len = MIN(width, seqs[k].length());
            fputs(seqs[k].substr(0, len).c_str(), fout);
            fputc('\n', fout);
            seqs[k].erase(0, len);
        }
        fputc('\n', fout);  // blank line after each block
    }
}

template<class Abc>
void Alignment<Abc>::Resize(size_t nseqs, size_t ncols) {
    if (nseqs == 0 || ncols == 0)
        throw Exception("Bad alignment dimensions: nseqs=%i ncols=%i",
                        nseqs, ncols);

    seqs_.Resize(ncols, nseqs);
    col_idx_.resize(ncols);
    match_idx_.resize(ncols);
    is_match_.resize(ncols);
    headers_.resize(nseqs);
}

template<class Abc>
void Alignment<Abc>::AssignMatchColumnsBySequence(size_t k) {
    if (ninsert() == 0) {
        // For alignemnts WITHOUT inserts simply assign columns in which the
        // reference sequence has no residue to be inserts
        for (size_t i = 0; i < ncols(); ++i)
            is_match_[i] = (seqs_[i][k] < Abc::kGap);

    } else {
        // For alignments WITH inserts we remove all columns in which the
        // reference sequence has no residue and remove all inserts in the other
        // sequences.
        const size_t ref_seq_length = GetSequence(k).length();
        Matrix<uint8_t> new_seqs(ref_seq_length, nseqs());
        size_t j = 0;  // column index in new alignment

        for (size_t i = 0; i < ncols(); ++i) {
            if (seqs_[i][k] < Abc::kGap) {
                for (size_t n = 0; n < nseqs(); ++n) {
                    if (n != k && !is_match(i))
                        new_seqs[j][n] = Abc::kGap;
                    else
                        new_seqs[j][n] = seqs_[i][n];
                }
                ++j;
            }
        }
        seqs_.Resize(ref_seq_length, nseqs());
        seqs_ = new_seqs;

        // Update match indices
        col_idx_.resize(ref_seq_length);
        match_idx_.resize(ref_seq_length);
        is_match_.resize(ref_seq_length);
        for (size_t i = 0; i < ref_seq_length; ++i) {
            col_idx_[i] = i;
            is_match_[i]   = true;
        }
    }

    SetMatchIndices();
}

template<class Abc>
void Alignment<Abc>::AssignMatchColumnsByGapRule(double gap_threshold) {
    // FIXME: handle inserts in assignment of match columns by gap rule
    if (ninsert() > 0) return;

    LOG(DEBUG) << "Marking columns with more than " << gap_threshold
               << "% of gaps as insert columns ...";

    // Global weights are sufficient for calculation of gap percentage
    Vector<double> wg;
    GlobalWeightsAndDiversity(*this, wg, gap_threshold <= 1.0);
    for (size_t i = 0; i < ncols(); ++i) {
        double gap = 0.0f;
        double res = 0.0f;

        for (size_t k = 0; k < nseqs(); ++k)
            if (seqs_[i][k] < Abc::kAny)
                res += wg[k];
            else if (seqs_[i][k] != Abc::kEndGap)  // ENDGAPs are ignored
                gap += wg[k];

        if (gap_threshold > 1.0) { // interpret as number between 1 and 100
            double percent_gaps = 100.0 * gap / (res + gap);
            is_match_[i] = (percent_gaps <= gap_threshold);
        } else { // interpret as decimal number
            double frac_res = res / (res + gap);
            is_match_[i] = (frac_res > gap_threshold);
        }
    }
    SetMatchIndices();
}

template<class Abc>
void Alignment<Abc>::RemoveInsertColumns() {
    // Create new sequence matrix
    const size_t match_cols = nmatch();
    Matrix<uint8_t> new_seqs(match_cols, nseqs());
    for (size_t i = 0; i < match_cols; ++i) {
        for (size_t k = 0; k < nseqs(); ++k) {
            new_seqs[i][k] = seqs_[match_idx_[i]][k];
        }
    }
    seqs_.Resize(match_cols, nseqs());
    seqs_ = new_seqs;

    // Update match indices
    col_idx_.resize(match_cols);
    match_idx_.resize(match_cols);
    is_match_.resize(match_cols);
    for (size_t i = 0; i < match_cols; ++i) {
        col_idx_[i] = i;
        is_match_[i]   = true;
    }
    SetMatchIndices();
}

template<class Abc>
void Alignment<Abc>::Merge(const Alignment<Abc>& ali) {
    if (nmatch() != ali.nmatch()) return;
    // FIXME: Keep insert columns when merging two alignments
    RemoveInsertColumns();

    // Copy and keep track of headers that are not already contained in master
    // alignment.
    Vector<bool> include_seq(ali.nseqs(), false);
    std::vector<std::string> headers_merged(headers_);
    for (size_t k = 0; k < ali.nseqs(); ++k) {
        if (find(headers_.begin(), headers_.end(), ali.header(k)) == headers_.end()) {
            include_seq[k] = true;
            headers_merged.push_back(ali.header(k));
        }
    }

    // Append new sequences to alignment
    Matrix<uint8_t> seqs_merged(nmatch(), headers_merged.size());
    for (size_t i = 0; i < nmatch(); ++i) {
        for (size_t k = 0; k < nseqs(); ++k) {
            seqs_merged[i][k] = seqs_[match_idx_[i]][k];
        }
        size_t l = nseqs(); // index of sequence k in merged alignment
        for (size_t k = 0; k < ali.nseqs(); ++k) {
            if (include_seq[k]) {
                seqs_merged[i][l] = ali.seqs_[ali.match_idx_[i]][k];
                ++l;
            }
        }
    }

    // Apply the merging
    headers_ = headers_merged;
    seqs_    = seqs_merged;
}

template<class Abc>
Sequence<Abc> Alignment<Abc>::GetSequence(size_t k) const {
    size_t nres = 0;
    for (size_t i = 0; i < ncols(); ++i)
        if (seq(k,i) < Abc::kGap) nres++;
    Sequence<Abc> sequence(nres);
    size_t j = 0;
    for (size_t i = 0; i < ncols(); ++i)
        if (seq(k,i) < Abc::kGap) sequence[j++] = seq(k,i);
    sequence.set_header(header(k));
    return sequence;
}

template<class Abc>
void Alignment<Abc>::Rearrange(const std::vector<Sequence<Abc> >& seqs) {
    assert_eq(nseqs(), seqs.size());

    // Determine mapping of alignment sequecnes to their position in 'seqs'
    Vector<int> mapping(nseqs());
    for (size_t k = 0; k < nseqs(); ++k) {
        for (size_t l = 0; l < seqs.size(); ++l) {
            if (seqs[l].header() == headers_[k]) {
                mapping[l] = k;
                break;
            }
        }
    }

    // Sort headers and sequence matrix
    std::vector<std::string> new_headers(nseqs());
    Matrix<uint8_t> new_seqs(ncols(), nseqs());
    for (size_t k = 0; k < nseqs(); ++k) {
        new_headers[k] = headers_[mapping[k]];
        for (size_t i = 0; i < ncols(); ++i)
            new_seqs[i][k] = seqs_[i][mapping[k]];
    }
    headers_ = new_headers;
    seqs_ = new_seqs;
}


template<class Abc>
std::ostream& operator<< (std::ostream& out, const Alignment<Abc>& ali) {
    const size_t kHeaderWidth = 18;
    const size_t kWrap     = 100;

    // Convert alignment to character representation
    Vector<std::string> seqs(ali.nseqs(), "");
    for (size_t k = 0; k < ali.nseqs(); ++k) {
        for (size_t i = 0; i < ali.ncols(); ++i) {
            char c = Abc::kIntToChar[ali.seqs_[i][k]];
            seqs[k].push_back(c);
        }
    }

    // Print alignment in blocks
    out << "CLUSTAL\n\n";
    while (!seqs[0].empty()) {
        for (size_t k = 0; k < ali.nseqs(); ++k) {
            std::string header =
                ali.headers_[k].substr(0, ali.headers_[k].find_first_of(' '));
            if (header.length() <= kHeaderWidth)
                out << header + std::string(kHeaderWidth - header.length(), ' ') << ' ';
            else
                out << header.substr(0, kHeaderWidth) << ' ';

            out << seqs[k].substr(0, MIN(kWrap, seqs[k].length()))
                << std::endl;
            seqs[k].erase(0, MIN(kWrap, seqs[k].length()));
        }
        out << std::endl;  // blank line after each block
    }
    return out;
}


template<class Abc>
void ReadAll(FILE* fin, AlignmentFormat format, std::vector< Alignment<Abc> >& v) {
    while (!feof(fin)) {
        v.push_back(Alignment<Abc>(fin, format));
        uint8_t c = fgetc(fin);
        if (c == EOF) break;
        ungetc(c, fin);
    }
}


// Returns the alignment format corresponding to provided filename extension
inline AlignmentFormat AlignmentFormatFromString(const std::string& s) {
    if (s == "fas" || s == "mfa")
        return FASTA_ALIGNMENT;
    else if (s == "a2m")
        return A2M_ALIGNMENT;
    else if (s == "a3m" || s == "ca3m")
        return A3M_ALIGNMENT;
    else if (s == "clu")
        return CLUSTAL_ALIGNMENT;
    else if (s == "psi")
        return PSI_ALIGNMENT;
    else
        throw Exception("Unknown alignment format extension '%s'!", s.c_str());
}

template<class Abc>
double GlobalWeightsAndDiversity(const Alignment<Abc>& ali, Vector<double>& wg, bool neff_sum_pairs) {
    const double kZero = 1E-10;  // for calculation of entropy
    const size_t nseqs = ali.nseqs();
    const size_t ncols = ali.nmatch();
    const size_t alphabet_size = Abc::kSize;
    const uint8_t any = Abc::kAny;

    // Return values
    wg.Assign(nseqs, 0.0f);  // global sequence weights
    double neff = 0.0f;      // diversity of alignment

    Vector<int> n(nseqs, 0);                      // number of residues in sequence i
    Vector<double> fj(alphabet_size, 0.0f);       // to calculate entropy
    Vector<int> adiff(ncols, 0);                  // different letters in each column
    Matrix<int> counts(ncols, alphabet_size, 0);  // column counts (excl. ANY)

    LOG(INFO) << "Calculation of global weights and alignment diversity ...";

    // Count number of residues in each column
    for (size_t i = 0; i < ncols; ++i) {
        for (size_t k = 0; k < nseqs; ++k) {
            if (ali[i][k] < any) {
                ++counts[i][ali[i][k]];
                ++n[k];
            }
        }
    }
    // Count number of different residues in each column
    for (size_t i = 0; i < ncols; ++i) {
        for (size_t a = 0; a < alphabet_size; ++a) {
            if (counts[i][a]) ++adiff[i];
        }
        if (adiff[i] == 0) adiff[i] = 1;  // col consists of only gaps and ANYs
    }
    // Calculate weights
    for (size_t i = 0; i < ncols; ++i) {
        for (size_t k = 0; k < nseqs; ++k) {
            if (adiff[i] > 0 && ali[i][k] < any)
                wg[k] += 1.0 / (adiff[i] * counts[i][ali[i][k]] * n[k]);
            //wg[k] += (adiff[i] - 1.0) / (counts[i][ali[i][k]] * n[k]);
        }
    }
    Normalize(&wg[0], nseqs);
    // Calculate number of effective sequences
    if (!neff_sum_pairs) {
        for (size_t i = 0; i < ncols; ++i) {
            Reset(&fj[0], alphabet_size);
            for (size_t k = 0; k < nseqs; ++k)
                if (ali[i][k] < any) fj[ali[i][k]] += wg[k];
            Normalize(&fj[0], alphabet_size);
            for (size_t a = 0; a < alphabet_size; ++a)
                if (fj[a] > kZero) neff -= fj[a] * log2(fj[a]);
        }
        neff = pow(2.0, neff / ncols);
    } else {
        for (size_t i = 0; i < ncols; ++i) {
            Reset(&fj[0], alphabet_size);
            for (size_t k = 0; k < nseqs; ++k)
                if (ali[i][k] < any) fj[ali[i][k]] += wg[k];
            Normalize(&fj[0], alphabet_size);
            double ni = nseqs + 1.0;
            for (size_t a = 0; a < alphabet_size; ++a)
                ni -= SQR(fj[a]) * nseqs;
            neff += ni;
        }
        neff /= ncols;
    }

    LOG(DEBUG) << "neff=" << neff;

    return neff;
}

template<class Abc>
Vector<double> PositionSpecificWeightsAndDiversity(const Alignment<Abc>& ali, Matrix<double>& w) {
    // Maximal fraction of sequences with an endgap
    const double kMaxEndgapFraction = 0.1;
    // Minimum number of columns in subalignments
    const size_t kMinCols = 10;
    // Zero value for calculation of entropy
    const double kZero = 1E-10;
    const size_t nseqs  = ali.nseqs();
    const size_t ncols  = ali.nmatch();
    const size_t alphabet_size = Abc::kSize;
    const uint8_t any    = Abc::kAny;
    const uint8_t endgap = Abc::kEndGap;

    // Return values
    Vector<double> neff(ncols, 0.0f);   // diversity of subalignment i
    w.Assign(ncols, nseqs, 0.0f);      // weight of seq k in column i
    // Helper variables
    size_t ncoli = 0;        // number of columns j that contribute to neff[i]
    size_t nseqi = 0;        // number of sequences in subalignment i
    size_t ndiff = 0;        // number of different alphabet letters
    bool change = false;  // has the set of sequences in subalignment changed?
    // Number of seqs with some residue in column i AND a at position j
    Matrix<int> n(ncols, endgap + 1, 0);
    // To calculate entropy
    Vector<double> fj(alphabet_size, 0.0f);
    // Weight of sequence k in column i, calculated from subalignment i
    Vector<double> wi(nseqs, 0.0f);
    Vector<double> wg;                   // global weights
    Vector<int> nseqi_debug(ncols, 0);   // debugging
    Vector<int> ncoli_debug(ncols, 0);   // debugging

    // Calculate global weights for fallback
    GlobalWeightsAndDiversity(ali, wg);

    for (size_t i = 0; i < ncols; ++i) {
        change = false;
        for (size_t k = 0; k < nseqs; ++k) {
            if ((i == 0 || ali[i-1][k] >= any) && ali[i][k] < any) {
                change = true;
                ++nseqi;
                for (size_t j = 0; j < ncols; ++j)
                    ++n[j][ali[j][k]];
            } else if (i > 0 && ali[i-1][k] < any && ali[i][k] >= any) {
                change = true;
                --nseqi;
                for (size_t j = 0; j < ncols; ++j)
                    --n[j][ali[j][k]];
            }
        }  // for k over nseqs
        nseqi_debug[i] = nseqi;

        if (change) {  // set of sequences in subalignment has changed
            ncoli = 0;
            Reset(&wi[0], nseqs);

            for (size_t j = 0; j < ncols; ++j) {
                if (n[j][endgap] > kMaxEndgapFraction * nseqi) continue;
                ndiff = 0;
                for (size_t a = 0; a < alphabet_size; ++a) if (n[j][a]) ++ndiff;
                if (ndiff == 0) continue;
                ++ncoli;
                for (size_t k = 0; k < nseqs; ++k) {
                    if (ali[i][k] < any && ali[j][k] < any) {
                        assert_ne(0, n[j][ali[j][k]]);
                        wi[k] += 1.0 / static_cast<double>((n[j][ali[j][k]] * ndiff));
                        //wi[k] += (ndiff - 1.0) / static_cast<double>(n[j][ali[j][k]]);
                    }
                }
            }  // for j over ncols
            Normalize(&wi[0], nseqs);

            if (ncoli < kMinCols)  // number of columns in subalignment insufficient?
                for (size_t k = 0; k < nseqs; ++k)
                    if (ali[i][k] < any)
                        wi[k] = wg[k];
                    else
                        wi[k] = 0.0f;

            neff[i] = 0.0f;
            for (size_t j = 0; j < ncols; ++j) {
                if (n[j][endgap] > kMaxEndgapFraction * nseqi) continue;
                Reset(&fj[0], alphabet_size);

                for (size_t k = 0; k < nseqs; ++k)
                    if (ali[i][k] < any && ali[j][k] < any)
                        fj[ali[j][k]] += wi[k];
                Normalize(&fj[0], alphabet_size);
                for (size_t a = 0; a < alphabet_size; ++a)
                    if (fj[a] > kZero) neff[i] -= fj[a] * log2(fj[a]);
            }  // for j over ncols

            neff[i] = (ncoli > 0 ? pow(2.0, neff[i] / ncoli) : 1.0f);

        } else {  // set of sequences in subalignment has NOT changed
            neff[i] = (i == 0 ? 0.0 : neff[i-1]);
        }

        for (size_t k = 0; k < nseqs; ++k) w[i][k] = wi[k];
        ncoli_debug[i] = ncoli;
    }  // for i over ncols

    return neff;
}

}  // namespace cs

#endif  // CS_ALIGNMENT_INL_H_
