/*
  Copyright 2009-2012 Andreas Biegert, Christof Angermueller

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

#ifndef CS_CONTEXT_WEIGHT_STATE_INL_H_
#define CS_CONTEXT_WEIGHT_STATE_INL_H_

#include "crf_state.h"

namespace cs {

template<class Abc>
void CrfState<Abc>::Read(FILE* fin) {
    // Parse and check header information
    if (!StreamStartsWith(fin, "CrfState"))
        throw Exception("Stream does not start with class id 'CrfState'!");

    char buffer[KB];
    cs::fgetline(buffer, KB, fin);
    if (strstr(buffer, "NAME")) {
        name = ReadString(buffer, "NAME", "Unable to parse CRF state 'NAME'!");
        cs::fgetline(buffer, KB, fin);
    }
    bias_weight = ReadDouble(buffer, "BIAS", "Unable to parse CRF state 'BIAS'!");
    cs::fgetline(buffer, KB, fin);
    size_t len = ReadInt(buffer, "LENG", "Unable to parse CRF state 'LENG'!");
    cs::fgetline(buffer, KB, fin);
    size_t nalph = ReadInt(buffer, "ALPH", "Unable to parse CRF state 'ALPH'!");
    assert(len & 1);
    if (nalph != Abc::kSize)
        throw Exception("Alphabet size of serialized CRF state should be %d"
                        "but is acutally %d!", Abc::kSize, nalph);

    // If everything went fine we can resize our data memmbers
    context_weights.Resize(len);

    // Read context weights and pseudocount weights
    const char* ptr = buffer;
    size_t i = 0;
    cs::fgetline(buffer, KB, fin);  // skip alphabet description line
    while (cs::fgetline(buffer, KB, fin) && buffer[0] != '/' && buffer[1] != '/') {
        ptr = buffer;
        if (buffer[0] != 'P' && buffer[1] != 'C') {
            i = strtoi(ptr) - 1;
            assert(i < len);
            // TODO: include ANY char in serialization
            for (size_t a = 0; a < Abc::kSize; ++a)
                context_weights[i][a] = static_cast<double>(strastoi(ptr)) / kScale;
            context_weights[i][Abc::kAny] = 0.0;
        } else {
            for (size_t a = 0; a < Abc::kSize; ++a)
                pc_weights[a] = static_cast<double>(strastoi(ptr)) / kScale;
        }
        context_weights[i][Abc::kAny] = 0.0;
        pc_weights[Abc::kAny] = 0.0;
    }
    if (i != len - 1)
        throw Exception("CRF state should have %i columns but actually has %i!",
                        len, i+1);
    UpdatePseudocounts(*this);
}

template<class Abc>
void CrfState<Abc>::Write(FILE* fout) const {
    // Print header section
    fputs("CrfState\n", fout);
    if (!name.empty()) fprintf(fout, "NAME\t%s\n", name.c_str());
    fprintf(fout, "BIAS\t%-10.8g\n", bias_weight);
    fprintf(fout, "LENG\t%d\n", static_cast<int>(context_weights.length()));
    fprintf(fout, "ALPH\t%d\n", static_cast<int>(Abc::kSize));

    // Print alphabet description line
    fputs("WEIGHTS", fout);
    for (size_t a = 0; a < Abc::kSize; ++a)
        fprintf(fout, "\t%c", Abc::kIntToChar[a]);
    fputs("\n", fout);
    // Print context weights scaled by 'kScale'
    for (size_t i = 0; i < context_weights.length(); ++i) {
        fprintf(fout, "%i", static_cast<int>(i+1));
        // TODO: include ANY char in serialization
        for (size_t a = 0; a < Abc::kSize; ++a) {
            if (context_weights[i][a] == -INFINITY) fputs("\t*", fout);
            else fprintf(fout, "\t%i", iround(context_weights[i][a] * kScale));
        }
        fputs("\n", fout);
    }
    // Print pseudocount weights
    fputs("PC", fout);
    for (size_t a = 0; a < Abc::kSize; ++a)
        fprintf(fout, "\t%i", iround(pc_weights[a] * kScale));
    fputs("\n//\n", fout);
}

// Prints CRF state weights probabilities in human-readable format for debugging.
template<class Abc>
std::ostream& operator<< (std::ostream& out, const CrfState<Abc>& state) {
    out << "CrfState" << std::endl;
    out << "name:\t" << state.name << std::endl;
    out << "bias:\t" << strprintf("%-10.8g", state.bias_weight) << std::endl;

    const int c = (state.context_weights.length() - 1) / 2;
    for (size_t i = 0; i < state.context_weights.length(); ++i)
        out << "\t   " << abs(static_cast<int>(i) - c);
    out << "\t  PCW\t  PC" << std::endl;

    for (size_t a = 0; a < Abc::kSizeAny; ++a) {
        out << Abc::kIntToChar[a];
        for (size_t i = 0; i < state.context_weights.length(); ++i) {
            out << strprintf("\t%+6.2f", state.context_weights[i][a]);
        }
        out << strprintf("\t%+6.2f\t%6.4f", state.pc_weights[a], state.pc[a])
            << std::endl;
    }

    return out;
}

// Updates pseudocount emission probs in given CRF state based on 'pc_weights' and
// rescales 'pc_weights'
template<class Abc>
inline void UpdatePseudocounts(CrfState<Abc>& state) {
    // Calculate maximum of pseudocount weights
    double max = -DBL_MAX;
    double mean = 0.0;
    for (size_t a = 0; a < Abc::kSize; ++a) {
        mean += state.pc_weights[a];
        if (state.pc_weights[a] > max) max = state.pc_weights[a];
    }
    mean /= Abc::kSize;

    // Rescale pseudocount weights and calculate their sum in lin-space
    long double sum = 0.0;
    for (size_t a = 0; a < Abc::kSize; ++a)
        sum += exp(state.pc_weights[a] - max);

    // Update emission pseudocount vector
    double tmp = max + log(sum);
    for (size_t a = 0; a < Abc::kSize; ++a) {
        state.pc[a] = DBL_MIN + exp(state.pc_weights[a] - tmp);        
        // state.pc_weights[a] -= mean; // Not necessary if pc_weights are centered on central context weights
    }
}

// Calculates context score between a CRF state and a sequence window
template<class Abc>
inline double ContextScore(const Profile<Abc>& context_weights,
                           const Sequence<Abc>& seq,
                           size_t idx,
                           size_t center) {
    assert(context_weights.length() & 1);
    const size_t beg = MAX(0, static_cast<int>(idx - center));
    const size_t end = MIN(seq.length(), idx + center + 1);
    double score = 0.0;
    for(size_t i = beg, j = beg - idx + center; i < end; ++i, ++j)
        score += context_weights[j][seq[i]];
    return score;
}

// Calculates context score between a CRF state and a count profile window
template<class Abc>
inline double ContextScore(const Profile<Abc>& context_weights,
                           const CountProfile<Abc>& cp,
                           size_t idx,
                           size_t center) {
    assert(context_weights.length() & 1);
    const size_t beg = MAX(0, static_cast<int>(idx - center));
    const size_t end = MIN(cp.counts.length(), idx + center + 1);
    double score = 0.0;
    for(size_t i = beg, j = beg - idx + center; i < end; ++i, ++j) {
        for (size_t a = 0; a < Abc::kSize; ++a)
            score += context_weights[j][a] * cp.counts[i][a];
    }
    return score;
}

}  // namespace cs

#endif  // CS_CRF_STATE_INL_H_
