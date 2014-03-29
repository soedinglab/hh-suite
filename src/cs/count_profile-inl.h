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

#ifndef CS_COUNT_PROFILE_INL_H_
#define CS_COUNT_PROFILE_INL_H_

#include "count_profile.h"

#include "alignment-inl.h"
#include "profile-inl.h"
#include "sequence-inl.h"

namespace cs {

template<class Abc>
CountProfile<Abc>::CountProfile(const Alignment<Abc>& ali, bool pos_weights, bool neff_sum_pairs)
  : counts(ali.nmatch(), 0.0f),
    neff(ali.nmatch()) {
  // Add counts and neff from alignment to count profile
  if (pos_weights) {  // use position-specific sequence weights
    Matrix<double> w;
    Vector<double> neff_tmp(PositionSpecificWeightsAndDiversity(ali, w));
    assert(neff.size());

    for (size_t i = 0; i < counts.length(); ++i) {
      neff[i] = neff_tmp[i];
      for (size_t k = 0; k < ali.nseqs(); ++k)
        if (ali[i][k] < Abc::kAny)
          counts[i][ali[i][k]] += w[i][k];
    }
  } else {  // use faster global sequence weights
    Vector<double> wg;
    double neff_glob = GlobalWeightsAndDiversity(ali, wg, neff_sum_pairs);
    for (size_t i = 0; i < counts.length(); ++i) {
      neff[i] = neff_glob;
      for (size_t k = 0; k < ali.nseqs(); ++k)
        if (ali[i][k] < Abc::kAny)
          counts[i][ali[i][k]] += wg[k];
    }
  }
  // Normalize counts to effective number of sequences
  Normalize(counts, neff);
}

template<class Abc>
void CountProfile<Abc>::Read(FILE* fin) {
  // Parse and check header information
  if (!StreamStartsWith(fin, "CountProfile"))
      throw Exception("Stream does not start with class id 'CountProfile'!");

  char buffer[KB];
  cs::fgetline(buffer, KB, fin);
  if (strstr(buffer, "NAME")) {
    name = ReadString(buffer, "NAME", "Unable to parse count profile 'NAME'!");
    cs::fgetline(buffer, KB, fin);
  }
  size_t len = ReadInt(buffer, "LENG", "Unable to parse count profile 'LENG'!");
  cs::fgetline(buffer, KB, fin);
  size_t nalph = ReadInt(buffer, "ALPH", "Unable to parse count profile 'ALPH'!");
  if (nalph != Abc::kSize)
    throw Exception("Alphabet size of serialized count profile should be %d "
                    "but is actually %d!", Abc::kSize, nalph);

  // If everything went fine we can resize our data memmbers
  counts.Resize(len);
  neff.Resize(len);

  // Read counts and effective number of sequences for each column
  size_t i = 0;
  const char* ptr = buffer;
  cs::fgetline(buffer, KB, fin);  // skip alphabet description line
  while (cs::fgetline(buffer, KB, fin) && buffer[0] != '/' && buffer[1] != '/') {
    ptr = buffer;
    i = strtoi(ptr) - 1;
    assert(i < len);
    // TODO: include ANY char in seialization
    for (size_t a = 0; a < Abc::kSize; ++a)
      // TODO: save counts as log base e instead of log base 2
      counts[i][a] = pow(2, static_cast<double>(-strastoi(ptr)) / kScale);
    counts[i][Abc::kAny] = 0.0f;
    neff[i] = static_cast<double>(strtoi(ptr)) / kScale;
  }
  Normalize(counts, neff);  // normalize probs to counts
  if (i != len - 1)
    throw Exception("Count profile should have %i columns but actually has %i!",
                    len, i+1);
}

template<class Abc>
void CountProfile<Abc>::Write(FILE* fout) const {
  // Print header section
  fputs("CountProfile\n", fout);
  if (!name.empty()) fprintf(fout, "NAME\t%s\n", name.c_str());
  fprintf(fout, "LENG\t%zu\n", counts.length());
  fprintf(fout, "ALPH\t%zu\n", Abc::kSize);

  // Print alphabet description line
  fputs("COUNTS", fout);
  for (size_t a = 0; a < Abc::kSize; ++a)
    fprintf(fout, "\t%c", Abc::kIntToChar[a]);
  fputs("\tNEFF\n", fout);

  // Print counts matrix and neff vector as negative logs scaled by 'kScale'
  for (size_t i = 0; i < counts.length(); ++i) {
    fprintf(fout, "%zu", i+1);
    // TODO: include ANY char in seialization
    for (size_t a = 0; a < Abc::kSize; ++a) {
      // TODO: save counts as log base e instead of log base 2
      if (counts[i][a] == 0.0) fputs("\t*", fout);
      else fprintf(fout, "\t%d", -iround(log2(counts[i][a] / neff[i]) * kScale));
    }
    fprintf(fout, "\t%d\n", iround(neff[i] * kScale));
  }
  fputs("//\n", fout);
}


// Returns the average Neff in given count profile.
template<class Abc>
inline double Neff(const CountProfile<Abc>& cp) {
  double neff_sum = 0.0;
  for (size_t i = 0; i < cp.neff.size(); ++i)
    neff_sum += cp.neff[i];
  return neff_sum / cp.neff.size();
}

// Builds and returns a consensus string of the given count profile by
// calculating at each position the alphabet character that deviates most strongly
// from its background probability.
template<class Abc>
std::string ConsensusSequence(const CountProfile<Abc>& cp,
                              const SubstitutionMatrix<Abc>& sm) {
  Profile<Abc> prof(cp.counts);
  Normalize(prof, 1.0);
  std::string cons(prof.length(), ' ');

  for (size_t i = 0; i < prof.length(); ++i) {
    double maxw = 0.0;
    size_t maxa = 0;
    for (size_t a = 0; a < Abc::kSize; ++a) {
      if (prof[i][a] - sm.p(a) > maxw) {
        maxw = prof[i][a] - sm.p(a);
        maxa = a;
      }
    }
    cons[i] = Abc::kIntToChar[maxa];
  }
  return cons;
}

// Builds and returns a conservation string for given count profile that
// indicates conservation of residues by uppercase, lowercase, and '~'
template<class Abc>
std::string ConservationSequence(const CountProfile<Abc>& cp,
                                 const SubstitutionMatrix<Abc>& sm) {
  static const char kMixedColChar = '~';

  // Precompute similarity matrix for amino acid pairs
  Matrix<double> sim(Abc::kSize, Abc::kSize);
  for (size_t a = 0; a < Abc::kSize; ++a)
    for (size_t b = 0; b < Abc::kSize; ++b)
      sim[a][b] = sm.q(a,b) * sm.q(a,b) / sm.q(a,a) / sm.q(b,b);

  // Normalize count profile
  Profile<Abc> prof(cp.counts);
  Normalize(prof, 1.0);
  std::string cons(prof.length(), ' ');

  // Now compute conservation string
  for (size_t i = 0; i < prof.length(); ++i) {
    double maxw = 0.0;
    size_t maxa = 0;
    for (size_t a = 0; a < Abc::kSize; ++a) {
      if (prof[i][a] - sm.p(a) > maxw) {
        maxw = prof[i][a] - sm.p(a);
        maxa = a;
      }
    }

    maxw = 0.0;
    for (size_t b = 0; b < Abc::kSize; ++b)
      maxw += prof[i][b] * sim[maxa][b] * sim[maxa][b];

    if (maxw > 0.6)
      cons[i] = toupper(Abc::kIntToChar[maxa]);
    else if (maxw > 0.4)
      cons[i] = tolower(Abc::kIntToChar[maxa]);
    else
      cons[i] = kMixedColChar;
  }

  return cons;
}

// Prints counts and neff in human-readable format for debugging.
template<class Abc>
std::ostream& operator<< (std::ostream& out, const CountProfile<Abc>& cp) {
  out << "CountProfile" << std::endl;
  out << "name:\t" << cp.name << std::endl;
  for (size_t a = 0; a < Abc::kSizeAny; ++a)
    out << "\t" << Abc::kIntToChar[a];
  out << "\tNeff" << std::endl;
  for (size_t i = 0; i < cp.counts.length(); ++i) {
    out << i+1;
    for (size_t a = 0; a < Abc::kSizeAny; ++a)
      out << strprintf("\t%6.4f", cp.counts[i][a]);
    out << strprintf("\t%-5.2f\n", cp.neff[i]);
  }
  return out;
}

}  // namespace cs

#endif  // CS_COUNT_PROFILE_INL_H_
