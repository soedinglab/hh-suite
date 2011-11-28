// Copyright 2009, Andreas Biegert

#ifndef CS_PSSM_H_
#define CS_PSSM_H_

#include "profile-inl.h"
#include "sequence-inl.h"

namespace cs {

// Container class for a PSI-BLAST position specific scoring matrix (PSSM) and its
// associated query sequence.
struct Pssm {
  // Constructor to create a PSSM from query string and a sequence profile.
  Pssm(const Sequence<AA>& seq, const Profile<AA>& prof)
      : query(seq), profile(prof) {}

  // Constructor to create a PSSM from a PSI-BLAST checkpoint file.
  Pssm(FILE* fin);

  // Overwrites existing PSSM with PSSM from PSI-BLAST checkpoint.
  void Read(FILE* fin) {
    size_t count = 0;
    int query_length = 0;

    count += fread(reinterpret_cast<char*>(&query_length), kIntSize, 1, fin);
    assert(count == 1);
    assert(query_length > 0);

    query.Resize(query_length);
    profile.Resize(query_length);

    for (int i = 0; i < query_length; ++i) {
      char c;
      count += fread(&c, kCharSize, 1, fin);
      assert(AA::kValidChar[c]);
      query[i] = AA::kCharToInt[static_cast<int>(c)];
    }

    for(int i = 0; i < query_length; ++i) {
      for(size_t a = 0; a < AA::kSize; ++a) {
        double p;
        count += fread(reinterpret_cast<char*>(&p), kDoubleSize, 1, fin);
        profile[i][a] = p;
      }
    }
    assert_eq(1 + query_length + query_length * AA::kSize, count);
  }

  // Writes PSSM in PSI-BLAST binary checkpoint format to stream.
  void Write(FILE* fout) const {
    LOG(INFO) << "Writing PSI-BLAST checkpoint ...";
    LOG(INFO) << query;
    LOG(INFO) << profile;

    int query_length = query.length();
    size_t count = 0;

    count += fwrite(reinterpret_cast<char*>(&query_length), kIntSize, 1, fout);
    for(int i = 0; i < query_length; ++i) {
      char c = query.chr(i);
      count += fwrite(&c, kCharSize, 1, fout);
    }

    for(int i = 0; i < query_length; ++i) {
      for(size_t a = 0; a < AA::kSize; ++a) {
        double p = profile[i][a];
        count += fwrite(reinterpret_cast<char*>(&p), kDoubleSize, 1, fout);
      }
    }
    assert_eq(1 + query_length + query_length * AA::kSize, count);
  }

  // Query sequence with which search was started
  Sequence<AA> query;
  // Evolving sequence profile
  Profile<AA> profile;
};  // class Pssm

}  // namespace cs

#endif  // CS_PSSM_H_
