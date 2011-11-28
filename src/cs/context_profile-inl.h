// Copyright 2009, Andreas Biegert

#ifndef CS_CONTEXT_PROFILE_INL_H_
#define CS_CONTEXT_PROFILE_INL_H_

#include "context_profile.h"

namespace cs {

template<class Abc>
void ContextProfile<Abc>::Read(FILE* fin) {
    // Parse and check header information
    if (!StreamStartsWith(fin, "ContextProfile"))
        throw Exception("Stream does not start with class id 'ContextProfile'!");

    char buffer[KB];
    cs::fgetline(buffer, KB, fin);
    if (strstr(buffer, "NAME")) {
        name = ReadString(buffer, "NAME", "Unable to parse context profile 'NAME'!");
        cs::fgetline(buffer, KB, fin);
    }
    prior = ReadDouble(buffer, "PRIOR", "Unable to parse context profile 'PRIOR'!");
    cs::fgetline(buffer, KB, fin);
    if (strstr(buffer, "COLOR")) {
        std::string coldef;
        coldef = ReadString(buffer, "COLOR", "Unable to parse context profile 'COLOR'!");
        color = Color(coldef);
        cs::fgetline(buffer, KB, fin);
    }
    is_log = ReadBool(buffer, "ISLOG", "Unable to parse context profile 'ISLOG'!");
    cs::fgetline(buffer, KB, fin);
    size_t  len = ReadInt(buffer, "LENG", "Unable to parse context profile 'LENG'!");
    cs::fgetline(buffer, KB, fin);
    size_t nalph = ReadInt(buffer, "ALPH", "Unable to parse context profile 'ALPH'!");
    if (is_log) prior = log(prior);
    assert(len & 1);
    if (nalph != Abc::kSize)
        throw Exception("Alphabet size of serialized context profile should be %d"
                        "but is acutally %d!", Abc::kSize, nalph);

    // If everything went fine we can resize our data memmbers
    probs.Resize(len);

    // Read counts and effective number of sequences for each column
    const char* ptr = buffer;
    const size_t center = (len - 1) / 2;
    size_t i = 0;
    cs::fgetline(buffer, KB, fin);  // skip alphabet description line
    while (cs::fgetline(buffer, KB, fin) && buffer[0] != '/' && buffer[1] != '/') {
        ptr = buffer;
        i = strtoi(ptr) - 1;
        assert(i < len);
        // TODO: include ANY char in seialization
        for (size_t a = 0; a < Abc::kSize; ++a) {
            // TODO: save probs as log base e instead of log base 2
            if (is_log)
                probs[i][a] = log(pow(2 ,static_cast<double>(-strastoi(ptr)) / kScale));
            // probs[i][a] = static_cast<double>(-strastoi(ptr)) / kScale;
            else
                probs[i][a] = pow(2, static_cast<double>(-strastoi(ptr)) / kScale);
            // probs[i][a] = exp(static_cast<double>(-strastoi(ptr)) / kScale);
        }
        probs[i][Abc::kAny] = is_log ? 0.0 : 1.0;
        // Set pseudocounts of central column
        if (i == center) {
            for (size_t a = 0; a < Abc::kSize; ++a)
                pc[a] = is_log ? exp(probs[i][a]) : probs[i][a];
            pc[Abc::kAny] = 1.0;
        }
    }
    if (i != len - 1)
        throw Exception("Context profile should have %i columns but actually has %i!",
                        len, i+1);
}

template<class Abc>
void ContextProfile<Abc>::Write(FILE* fout) const {
    // Print header section
    fputs("ContextProfile\n", fout);
    if (!name.empty()) fprintf(fout, "NAME\t%s\n", name.c_str());
    fprintf(fout, "PRIOR\t%-10.8g\n", is_log ? exp(prior) : prior);
    fprintf(fout, "COLOR\t%s\n", color.ToString().c_str());
    fprintf(fout, "ISLOG\t%c\n", is_log ? 'T' : 'F');
    fprintf(fout, "LENG\t%d\n",  static_cast<int>(probs.length()));
    fprintf(fout, "ALPH\t%d\n",  static_cast<int>(Abc::kSize));

    // Print alphabet description line
    fputs("PROBS", fout);
    for (size_t a = 0; a < Abc::kSize; ++a)
        fprintf(fout, "\t%c", Abc::kIntToChar[a]);
    fputs("\n", fout);
    // Print probs as negative logs scaled by 'kScale'
    for (size_t i = 0; i < probs.length(); ++i) {
        fprintf(fout, "%zu", i+1);
        // TODO: include ANY char in serialization
        for (size_t a = 0; a < Abc::kSize; ++a) {
            // TODO: save probs as log base e instead of log base 2a
            // double log_val = is_log ? probs[i][a] : log(probs[i][a]);
            double log_val = is_log ? log2(exp(probs[i][a])) : log2(probs[i][a]);
            if (log_val == -INFINITY) fputs("\t*", fout);
            else fprintf(fout, "\t%i", -iround(log_val * kScale));
        }
        fputs("\n", fout);
    }
    fputs("//\n", fout);
}

// Prints context profile probabilities in human-readable format for debugging.
template<class Abc>
std::ostream& operator<< (std::ostream& out, const ContextProfile<Abc>& cp) {
    out << "ContextProfile" << std::endl;
    out << "name:\t" << cp.name << std::endl;
    out << "prior:\t" << strprintf("%-10.8g\n", cp.is_log ? exp(cp.prior) : cp.prior);
    out << "color:\t" << cp.color << std::endl;
    out << "is_log:\t" << (cp.is_log ? "true" : "false") << std::endl;

    const int c = (cp.probs.length() - 1) / 2;
    for (size_t i = 0; i < cp.probs.length(); ++i)
        out << "\t   " << abs(static_cast<int>(i) - c);
    out << std::endl;

    for (size_t a = 0; a < Abc::kSizeAny; ++a) {
        out << Abc::kIntToChar[a];
        for (size_t i = 0; i < cp.probs.length(); ++i) {
            double val = cp.is_log ? exp(cp.probs[i][a]) : cp.probs[i][a];
            out << strprintf("\t%6.4f", val);
        }
        out << std::endl;
    }

    return out;
}

// Transforms probabilites in context profile to log-space and sets 'is_log' flag.
template<class Abc>
inline void TransformToLog(ContextProfile<Abc>& cp) {
    if (!cp.is_log) {
        cp.prior = log(cp.prior);
        for (size_t i = 0; i < cp.probs.length(); ++i)
            for (size_t a = 0; a < Abc::kSizeAny; ++a)
                cp.probs[i][a] = log(cp.probs[i][a]);
        cp.is_log = true;
    }
}

// Transforms probabilites in context profile to lin-space and sets 'is_log' flag.
template<class Abc>
inline void TransformToLin(ContextProfile<Abc>& cp) {
    if (cp.is_log) {
        cp.prior = exp(cp.prior);
        for (size_t i = 0; i < cp.probs.length(); ++i)
            for (size_t a = 0; a < Abc::kSizeAny; ++a)
                cp.probs[i][a] = exp(cp.probs[i][a]);
        cp.is_log = false;
    }
}

}  // namespace cs

#endif  // CS_CONTEXT_PROFILE_INL_H_
