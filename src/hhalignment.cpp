// hhalignment.C

/////////////////////////////////////////////////////////////////////////////////////
//// Class Alignment
/////////////////////////////////////////////////////////////////////////////////////

// hhalignment.C
#include "hhalignment.h"
#include "util.h"
#include <vector>
#include <map>

using std::cout;
using std::cerr;
using std::endl;
using std::ios;
using std::ifstream;
using std::ofstream;

/////////////////////////////////////////////////////////////////////////////////////
// Class Alignment
/////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////
// Object constructor
/////////////////////////////////////////////////////////////////////////////////////
Alignment::Alignment(int maxseq, int maxres) {
  longname = new char[DESCLEN];
  sname = new char*[maxseq + 2];
  seq = new char*[maxseq + 2];
  l = new int[maxres];
  X = new char*[maxseq + 2];
  I = new short unsigned int*[maxseq + 2];
  keep = new char[maxseq + 2];
  display = new char[maxseq + 2];
  wg = new float[maxseq + 2];
  this->maxseq = maxseq;
  nseqs = new int[maxres + 2];
  N_in = L = 0;
  nres = NULL;           // number of residues per sequence k
  first = NULL;          // first residue in sequence k
  last = NULL;           // last  residue in sequence k
  ksort = NULL;          // sequence indices sorted by descending nres[k]
  name[0] = '\0';        // no name defined yet
  longname[0] = '\0';    // no name defined yet
  fam[0] = '\0';         // no name defined yet
  file[0] = '\0';        // no name defined yet
  readCommentLine = '0';
}

/////////////////////////////////////////////////////////////////////////////////////
// Object destructor
/////////////////////////////////////////////////////////////////////////////////////
Alignment::~Alignment() {
  delete[] longname;
  for (int k = 0; k < N_in; ++k) {
    delete[] sname[k];
    delete[] seq[k];
    free(X[k]);
    delete[] I[k];
  }
  delete[] sname;
  delete[] seq;
  delete[] X;
  delete[] I;
  delete[] l;
  delete[] keep;
  delete[] display;
  delete[] wg;
  delete[] nseqs;
  delete[] nres;
  delete[] first;
  delete[] last;
  delete[] ksort;
}

char * Alignment::initX(int len) {
  int seqSimdLength = (len) / (VECSIZE_INT * 4) + 2;
  seqSimdLength *= (VECSIZE_INT * 4);
  char * ptr = (char *) malloc_simd_int(seqSimdLength);
  std::fill(ptr, ptr + seqSimdLength, GAP);
  return ptr;
}

/////////////////////////////////////////////////////////////////////////////////////
// Deep-copy constructor
/////////////////////////////////////////////////////////////////////////////////////
Alignment& Alignment::operator=(Alignment& ali) {

  // First delete all arrays
  for (int k = 0; k < N_in; ++k) {
    delete[] sname[k];
    delete[] seq[k];
    delete[] X[k];
    delete[] I[k];
  }

  L = ali.L;
  N_in = ali.N_in;
  N_filtered = ali.N_filtered;
  N_ss = ali.N_ss;
  n_display = ali.n_display;

  if (ksort) {
    delete[] ksort;
    ksort = NULL;
  }

  if (first) {
    delete[] first;
    first = NULL;
  }

  if (last) {
    delete[] last;
    last = NULL;
  }

  if (nres) {
    delete[] nres;
    nres = NULL;
  }

  // Then allocate new space and copy stuff from source alignment
  if (ali.nres) {
    nres = new int[N_in];
    for (int k = 0; k < N_in; ++k)
      nres[k] = ali.nres[k];
  }
  for (int k = 0; k < N_in; ++k) {
    sname[k] = new char[strlen(ali.sname[k]) + 1];
    if (!sname[k])
      MemoryError("array of names for sequences to display", __FILE__,
      __LINE__,
                  __func__);
    strcpy(sname[k], ali.sname[k]);
  }
  for (int k = 0; k < N_in; ++k) {
    seq[k] = new char[strlen(ali.seq[k]) + 1];
    if (!seq[k])
      MemoryError("array of names for sequences to display", __FILE__,
      __LINE__,
                  __func__);
    strcpy(seq[k], ali.seq[k]);
  }
  for (int k = 0; k < N_in; ++k) {
    X[k] = initX(strlen(ali.seq[k]) + 2);
    if (!X[k])
      MemoryError("array for input sequences", __FILE__, __LINE__, __func__);
  }
  for (int k = 0; k < N_in; ++k) {
    I[k] = new short unsigned int[strlen(ali.seq[k]) + 2];
    if (!I[k])
      MemoryError("array for input sequences", __FILE__, __LINE__, __func__);
  }
  for (int k = 0; k < N_in; ++k) {
    for (int i = 0; i <= L; ++i) {
      X[k][i] = ali.X[k][i];
      I[k][i] = ali.I[k][i];
    }
  }
  for (int k = 0; k < N_in; ++k) {
    keep[k] = ali.keep[k];
    display[k] = ali.display[k];
    wg[k] = ali.wg[k];
  }

  kss_dssp = ali.kss_dssp;
  ksa_dssp = ali.ksa_dssp;
  kss_pred = ali.kss_pred;
  kss_conf = ali.kss_conf;
  kfirst = ali.kfirst;

  strcpy(longname, ali.longname);
  strmcpy(name, ali.name, NAMELEN - 1);
  strmcpy(fam, ali.fam, NAMELEN - 1);
  strmcpy(file, ali.file, NAMELEN - 1);

  for (int i = 1; i <= L; ++i)
    l[i] = ali.l[i];

  return (Alignment&) (*this);
}

/////////////////////////////////////////////////////////////////////////////////////
// Reads in an alignment from file into matrix seq[k][l] as ASCII
/////////////////////////////////////////////////////////////////////////////////////
void Alignment::Read(FILE* inf, char infile[], const char mark,
                     const int maxcol, const int nseqdis, char* firstline) {
  int l;                  // Postion in alignment incl. gaps (first=1)
  int h;                  // Position in input line (first=0)
  int k;                  // Index of sequence being read currently (first=0)
  char line[LINELEN] = "";  // input line
  char cur_seq[maxcol];   // Sequence currently read in
  char* cur_name;         // Sequence currently read in
  int linenr = 0;           // current line number in input file
  char skip_sequence = 0;
  char no_name[] = "no_name_1234567890";
//	std::cout << infile << std::endl;
  RemoveExtension(file, infile);  //copy rootname (w/o path) of infile into file variable of class object
//	std::cout << file << std::endl;

  kss_dssp = ksa_dssp = kss_pred = kss_conf = kfirst = -1;
  n_display = 0;
  N_in = 0;
  N_filtered = 0;
  N_ss = 0;
  cur_seq[0] = ' ';  // overwrite '\0' character at beginning to be able to do strcpy(*,cur_seq)
  l = 1;
  k = -1;

  // Does firstline already contain first line of file?
  if (firstline != NULL)
    strcpy(line, firstline);

  /////////////////////////////////////////////////////////////////////////
  // Read infile line by line
  while (firstline || (fgetline(line, LINELEN, inf) && k < MAXSEQ)) {
    linenr++;
    firstline = NULL;
    if (line[0] == '>')             //line contains sequence name
        {
      if (k >= MAXSEQ - 1) {
        if (k >= MAXSEQ) {
          HH_LOG(WARNING) << "Maximum number " << MAXSEQ
                                    << " of sequences exceeded in file "
                                    << infile << std::endl;
        }
        break;
      }

      if (k >= 0)  //if this is at least the second name line
          {
        // 1, because the sequence in cur_seq starts at position 1 => no residues = length 1
        if (strlen(cur_seq) <= 1) {
          HH_LOG(ERROR) << "In " << __FILE__ << ":" << __LINE__ << ": " << __func__ << ":" << std::endl;
          HH_LOG(ERROR) << "\tsequence " << sname[k] << " contains no residues." << std::endl;
          exit(1);
        }

        // Create space for residues and paste new sequence in
        seq[k] = new char[strlen(cur_seq) + 2];
        if (!seq[k])
          MemoryError("array for input sequences", __FILE__, __LINE__,
                      __func__);
        X[k] = initX(strlen(cur_seq) + 2);
        if (!X[k])
          MemoryError("array for input sequences", __FILE__, __LINE__,
                      __func__);
        I[k] = new short unsigned int[strlen(cur_seq) + 2];
        if (!I[k])
          MemoryError("array for input sequences", __FILE__, __LINE__,
                      __func__);
        strcpy(seq[k], cur_seq);
      }
      skip_sequence = 0;

      // Set name of new sequence
      ++k;
      cur_name = strscn(line + 1);  // go to beginning of current sequence name (skip white space)
      if (!cur_name)  // if no name is given, name sequence >no_name_123
      {
        cur_name = no_name;
        sprintf(cur_name, "no_name_%i", k);
      } else if (cur_name[0] == '@')
        cur_name = strscn(cur_name + 1);  //skip @-character in name

      l = 1;  //position in current sequence (first=1)
      cur_seq[l] = '\0';  // avoids taking wrong sequence in case the input alignment is corrupted (two header lines with no sequence line between)

      // display[k]= 0: do not show in Q-T alignments  1: show if not filtered out later     2: show in any case    (do not filter out)
      // keep[k]   = 0: do not include in profile      1: include if not filtered out later  2: include in any case (do not filter out)
      if (!strncmp(line, ">ss_dssp", 8)) {
        if (kss_dssp < 0) {
          display[k] = 2;
          n_display++;
          keep[k] = 0;
          kss_dssp = k;
          N_ss++;
        } else {
          skip_sequence = 1;
          k--;
          continue;
        }
      } else if (!strncmp(line, ">sa_dssp", 8)) {
        if (ksa_dssp < 0) {
          display[k] = 2;
          n_display++;
          keep[k] = 0;
          ksa_dssp = k;
          N_ss++;
        } else {
          skip_sequence = 1;
          k--;
          continue;
        }
      } else if (!strncmp(line, ">ss_pred", 8)) {
        if (kss_pred < 0) {
          display[k] = 2;
          n_display++;
          keep[k] = 0;
          kss_pred = k;
          N_ss++;
        } else {
          skip_sequence = 1;
          k--;
          continue;
        }
      } else if (!strncmp(line, ">ss_conf", 8)) {
        if (kss_conf < 0) {
          display[k] = 2;
          n_display++;
          keep[k] = 0;
          kss_conf = k;
          N_ss++;
        } else {
          skip_sequence = 1;
          k--;
          continue;
        }
      } else if (!strncmp(line, ">ss_", 4) || !strncmp(line, ">sa_", 4)) {
        display[k] = 2;
        n_display++;
        keep[k] = 0;
        N_ss++;
      } else if (!strncmp(line, ">aa_", 4)) {  // ignore sequences beginning with ">aa_"
        skip_sequence = 1;
        k--;
        continue;
      }
      //store first real seq
      else if (kfirst < 0) {
        char word[NAMELEN];
        strwrd(word, line, NAMELEN);  // Copies first word in ptr to str
        if (strstr(word, "_consensus")) {
          display[k] = 2;
          keep[k] = 0;
          n_display++;
          kfirst = k;
        } else {
          display[k] = keep[k] = 2;
          n_display++;
          kfirst = k;
        }
      }
      //store all sequences
      else if (mark == 0) {
        display[k] = keep[k] = 1;
        n_display++;
      }
      //store sequences up to nseqdis
      else if (line[1] == '@' && n_display - N_ss < nseqdis) {
        display[k] = keep[k] = 2;
        n_display++;
      } else if (mark == 1) {
        display[k] = keep[k] = 1;
        n_display++;
      }
      //store sequences up to nseqdis
      else {
        display[k] = 0;
        keep[k] = 1;
      }

      // store sequence name
      HH_LOG(DEBUG1) << "Reading seq " << cur_name << " k=" << k
                               << "  n_displ=" << n_display << "  display[k]="
                               << display[k] << " keep[k]=" << keep[k]
                               << std::endl;
      sname[k] = new char[strlen(cur_name) + 1];
      if (!sname[k]) {
        MemoryError("array for sequence names", __FILE__, __LINE__, __func__);
      }
      strcpy(sname[k], cur_name);
    }  // end if(line contains sequence name)

    else if (line[0] == '#')  // Commentary line?
        {
      // #PF01367.9 5_3_exonuc: 5'-3' exonuclease, C-terminal SAM fold; PDB 1taq, 1bgx (T:271-174), 1taq (271-174)
      if (name[0])
        continue;  // if already name defined: skip commentary line
      char *ptr1;
      ptr1 = strscn_(line + 1);  // set ptr1 to first non-whitespace character after '#' -> AC number
      strncpy(longname, ptr1, DESCLEN - 1);  // copy whole commentary line after '# ' into longname
      longname[DESCLEN - 1] = '\0';
      strtr(longname, "", " ");
      strcut_(ptr1);  // cut after AC number and set ptr2 to first non-whitespace character after AC number
//        char* ptr2=strcut_(ptr1); // **instead** of previous line
//        strcpy(fam,ptr1);    // copy AC number to fam
//        if (!strncmp(fam,"PF",2)) strcut_(fam,'.'); // if PFAM identifier contains '.' cut it off
//        strcut_(ptr2);       // cut after first word ...
      strmcpy(name, ptr1, NAMELEN - 1);  // ... and copy first word into name
      readCommentLine = '1';
    }

    //line contains sequence residues or SS information and does not belong to a >aa_ sequence
    else if (!skip_sequence) {
      if (k == -1) {
        HH_LOG(WARNING)
            << "No sequence name preceding following line in "
            << infile << ":\n\'" << line << "\'\n";
        continue;
      }

      h = 0;  //counts characters in current line

      // Check whether all characters are correct; store into cur_seq
      if (keep[k] || k == kfirst)  // normal line containing residues
          {
        while (h < LINELEN && line[h] > '\0' && l < maxcol - 1) {
          if (aa2i(line[h]) >= 0)  // ignore white-space characters ' ', \t and \n (aa2i()==-1)
              {
            cur_seq[l] = line[h];
            l++;
          } else if (aa2i(line[h]) == -2) {
            HH_LOG(WARNING) << "Ignoring invalid symbol \'"
                                      << line[h] << "\' at pos. " << h
                                      << " in line " << linenr << " of "
                                      << infile << "\n";
          }
          h++;
        }
      } else if (k == kss_dssp)  // lines with dssp secondary structure states (. - H E C S T G B)
          {
        while (h < LINELEN && line[h] > '\0' && l < maxcol - 1) {
          if (ss2i(line[h]) >= 0 && ss2i(line[h]) <= 7) {
            cur_seq[l] = ss2ss(line[h]);
            l++;
          } else
          HH_LOG(WARNING) << "Ignoring invalid symbol \'"
                                    << line[h] << "\' at pos. " << h
                                    << " in line " << linenr << " of " << infile
                                    << "\n";
          h++;
        }
      } else if (k == ksa_dssp)  // lines with dssp solvent accessibility states (. - ???)
          {
        while (h < LINELEN && line[h] > '\0' && l < maxcol - 1) {
          if (sa2i(line[h]) >= 0)
            cur_seq[l++] = line[h];
          else {
            HH_LOG(WARNING) << "Ignoring invalid symbol \'"
                                      << line[h] << "\' at pos. " << h
                                      << " in line " << linenr << " of "
                                      << infile << "\n";
          }
          h++;
        }
      } else if (k == kss_pred)  // lines with predicted secondary structure (. - H E C)
          {
        while (h < LINELEN && line[h] > '\0' && l < maxcol - 1) {
          if (ss2i(line[h]) >= 0 && ss2i(line[h]) <= 3) {
            cur_seq[l] = ss2ss(line[h]);
            l++;
          } else
          HH_LOG(WARNING) << "ignoring invalid symbol \'" << line[h]
                                    << "\' at pos. " << h << " in line "
                                    << linenr << " of " << infile << "\n";
          h++;
        }
      } else if (k == kss_conf)  // lines with confidence values should contain only 0-9, '-', or '.'
          {
        while (h < LINELEN && line[h] > '\0' && l < maxcol - 1) {
          if (line[h] == '-' || line[h] == '.'
              || (line[h] >= '0' && line[h] <= '9')) {
            cur_seq[l] = line[h];
            l++;
          } else
          HH_LOG(WARNING) << "ignoring invalid symbol \'" << line[h]
                                    << "\' at pos. " << l << " in line "
                                    << linenr << " of " << infile << "\n";
          h++;
        }
      } else if (display[k])  // other lines such as >sa_pred etc
      {
        while (h < LINELEN && line[h] > '\0' && l < maxcol - 1) {
          if (line[h] == '-' || line[h] == '.'
              || (line[h] >= '0' && line[h] <= '9')
              || (line[h] >= 'A' && line[h] <= 'B')) {
            cur_seq[l] = line[h];
            l++;
          } else
          HH_LOG(WARNING) << "ignoring invalid symbol \'" << line[h]
                                    << "\' at pos. " << l << " in line "
                                    << linenr << " of " << infile << "\n";
          h++;
        }
      }

      if (l >= maxcol - 1) {
        HH_LOG(WARNING) << "maximum number of residues " << maxcol - 2
                                  << " exceeded in sequence " << sname[k]
                                  << "\n";
        l = maxcol - 1; // Bug removed 29.10.14 by MM & JS by inserting this line
        skip_sequence = 1;
      }
      cur_seq[l] = '\0';  //Ensure that cur_seq ends with a '\0' character
    }  //end else
  }
  /////////////////////////////////////////////////////////////////////////

  if (k >= 0)  //if  at least one sequence was read in
      {
    seq[k] = new char[strlen(cur_seq) + 2];
    if (!seq[k])
      MemoryError("array for input sequences", __FILE__, __LINE__, __func__);
    X[k] = initX(strlen(cur_seq) + 2);
    if (!X[k])
      MemoryError("array for input sequences", __FILE__, __LINE__, __func__);
    I[k] = new short unsigned int[strlen(cur_seq) + 2];
    if (!I[k])
      MemoryError("array for input sequences", __FILE__, __LINE__, __func__);
    strcpy(seq[k], cur_seq);
  } else {
    HH_LOG(ERROR) << "In " << __FILE__ << ":" << __LINE__ << ": " << __func__ << ":" << std::endl;
    HH_LOG(ERROR) << "\tNo sequences found in file " << infile << std::endl;
    exit(1);
  }

  N_in = k + 1;

  // Warn if there are only special sequences but no master sequence (consensus seq given if keep[kfirst]==0)
  if (kfirst < 0 || (N_in - N_ss - (keep[kfirst] == 0 ? 1 : 0)) == 0) {
    HH_LOG(ERROR) << "In " << __FILE__ << ":" << __LINE__ << ": " << __func__ << ":" << std::endl;
    HH_LOG(ERROR) << "\tMSA file " << infile << " contains no master sequence!" << std::endl;
    exit(1);
  }

  // Set name, longname, fam
  if (!*name)  // longname, name and family were not set by '#...' line yet -> extract from first sequence
  {
    char* ptr;
//    strtr(sname[kfirst],"~"," ");              // 'transpose': replaces the tilde with a blanc everywhere in sname[kfirst]
    strncpy(longname, sname[kfirst], DESCLEN - 1);  // longname is name of first sequence
    longname[DESCLEN - 1] = '\0';
    strncpy(name, sname[kfirst], NAMELEN - 1);  // Shortname is first word of longname...
    name[NAMELEN - 1] = '\0';
    ptr = strcut(name);              // ...until first white-space character
    if (ptr && islower(ptr[0]) && ptr[1] == '.' && isdigit(ptr[2]))  //Scop family code present as second word?
        {
      lwrstr(name);                  // Transform upper case to lower case
      strcut(ptr);  // Non-white-space characters until next white-space character..
      strcpy(fam, ptr);                    // ...are the SCOP familiy code
    } else if (name[0] == 'P' && name[1] == 'F' && isdigit(name[2])
        && isdigit(name[3]))  //Pfam code
            {
      strcpy(fam, name);                    // set family name = Pfam code
    }
  }

  HH_LOG(DEBUG) << "Read " << infile << " with " << N_in
                          << " sequences\n";
  HH_LOG(DEBUG1) << "Query sequence for alignment has number "
                           << kfirst << " (0 is first)\n";
}

void Alignment::ReadCompressed(ffindex_entry_t* entry, char* data,
                               ffindex_index_t* ffindex_sequence_database_index,
                               char* ffindex_sequence_database_data,
                               ffindex_index_t* ffindex_header_database_index,
                               char* ffindex_header_database_data,
                               const char mark, const int maxcol) {

  char cur_seq[maxcol];
  cur_seq[0] = ' ';
  char cur_header[NAMELEN];

  RemoveExtension(file, entry->name);
  size_t data_size = entry->length - 1;
  size_t index = 0;

  kss_dssp = ksa_dssp = kss_pred = kss_conf = kfirst = -1;
  n_display = 0;
  N_in = 0;
  N_filtered = 0;
  N_ss = 0;

  // Commentary line?
  if ((*data) == '#' && !name[0]) {
    //skip #
    index++;
    data++;

    //skip first white spaces
    while (isspace(*data) && index < data_size) {
      index++;
      data++;
    }

    //copy name in name and complete header in longname
    size_t longname_index = 0;
    size_t name_index = 0;
    size_t count_ws = 0;
    while ((*data) != '\n' && index < data_size) {
      if (isspace(*data)) {
        count_ws++;
        if (count_ws == 1) {
          name[name_index++] = '\0';
        }
      }

      if (count_ws == 0) {
        name[name_index++] = (*data);
      }

      longname[longname_index++] = (*data);

      data++;
      index++;
    }
    longname[longname_index] = '\0';

    readCommentLine = '1';

    //skip new line after commentary line
    data++;
    index++;
  }

  char last_char = '\0';
  char inConsensus = 0;
  size_t consensus_length = 0;
  size_t name_length = 0;

  while (!(last_char == '\n' && (*data) == ';') && index < data_size) {
    if ((*data) == '\n') {
      inConsensus++;
    } else {
      if (inConsensus == 0) {
        cur_header[name_length++] = (*data);
      } else if (inConsensus == 1) {
        cur_seq[++consensus_length] = (*data);
      }
    }

    last_char = (*data);
    data++;
    index++;
  }
  cur_seq[consensus_length + 1] = '\0';
  cur_header[name_length] = '\0';

  //get past ';'
  data++;
  index++;

  // Index of sequence being read currently (first=0)
  int k = 0;
  //process consensus
  display[k] = 2;
  keep[k] = 0;
  n_display++;
  kfirst = k;

  X[k] = initX(consensus_length + 2);
  I[k] = new short unsigned int[consensus_length + 2];

  seq[k] = new char[consensus_length + 2];
  seq[k][0] = ' ';
  int copy_pos = 1;
  for (size_t string_pos = 1; string_pos <= consensus_length; string_pos++) {
    if (aa2i(cur_seq[string_pos]) >= 0) {
      seq[k][copy_pos++] = cur_seq[string_pos];
    }
  }
  seq[k][copy_pos] = '\0';  //Ensure that cur_seq ends with a '\0' character

  char* cur_name = strscn(cur_header + 1);
  sname[k] = new char[strlen(cur_name) + 1];
  strcpy(sname[k], cur_name);

  if (copy_pos >= maxcol - 1) {
    HH_LOG(WARNING) << "Maximum number of residues "
                              << maxcol - 2 << " exceeded in sequence "
                              << sname[k] << std::endl;
  }

  k++;

  while (index < data_size) {
    unsigned int entry_index;
    unsigned short int nr_blocks;
    unsigned short int start_pos;

    readU32(&data, entry_index);
    index += 4;

    ffindex_entry_t* sequence_entry = ffindex_get_entry_by_index(
        ffindex_sequence_database_index, entry_index);
    if (sequence_entry == NULL) {
      HH_LOG(ERROR) << "Could not fetch sequence entry: " << entry_index << " for alignment " << entry->name << std::endl;
      exit(1);
    }

    char* sequence_data = ffindex_get_data_by_entry(
        ffindex_sequence_database_data, sequence_entry);
    if (sequence_data == NULL) {
      HH_LOG(ERROR) << "Could not fetch sequence data: " << entry_index << " for alignment " << entry->name << std::endl;
      exit(1);
    }

    ffindex_entry_t* header_entry = ffindex_get_entry_by_index(
        ffindex_header_database_index, entry_index);
    if (header_entry == NULL) {
      HH_LOG(ERROR) << "Could not fetch header entry: " << entry_index << " for alignment " << entry->name << std::endl;
      exit(1);
    }

    char* header_data = ffindex_get_data_by_entry(ffindex_header_database_data,
                                                  header_entry);
    if (header_data == NULL) {
      HH_LOG(ERROR) << "Could not fetch header data: " << entry_index << " for alignment " << entry->name << std::endl;
      exit(1);
    }

    readU16(&data, start_pos);
    index += 2;

    readU16(&data, nr_blocks);
    index += 2;

    size_t actual_pos = start_pos;
    size_t alignment_length = 0;
    size_t alignment_index = 1;
    for (unsigned short int block_index = 0; block_index < nr_blocks;
        block_index++) {
      unsigned char nr_matches = (unsigned char) (*data);
      data++;
      index++;

      for (int i = 0; i < nr_matches; i++) {
        cur_seq[alignment_index++] = sequence_data[actual_pos - 1];
        actual_pos++;
        alignment_length++;
      }

      char nr_insertions_deletions = (*data);
      data++;
      index++;

      if (nr_insertions_deletions > 0) {
        for (int i = 0; i < nr_insertions_deletions; i++) {
          cur_seq[alignment_index++] = tolower(sequence_data[actual_pos - 1]);
          actual_pos++;
        }
      } else {
        for (int i = 0; i < -nr_insertions_deletions; i++) {
          cur_seq[alignment_index++] = '-';
          alignment_length++;
        }
      }
    }

    while (alignment_length < consensus_length) {
      cur_seq[alignment_index++] = '-';
      alignment_length++;
    }

    cur_seq[alignment_index] = '\0';

    //process sequence with header
    if (mark == 0) {
      display[k] = keep[k] = 1;
      n_display++;
    } else if (mark == 1) {
      display[k] = keep[k] = 1;
      n_display++;
    } else {
      display[k] = 0;
      keep[k] = 1;
    }

    X[k] = initX(alignment_index + 1);
    I[k] = new short unsigned int[alignment_index + 1];

    seq[k] = new char[alignment_index + 1];
    seq[k][0] = ' ';
    size_t copy_pos = 1;
    for (size_t string_pos = 1; string_pos < alignment_index; string_pos++) {
      if (aa2i(cur_seq[string_pos]) >= 0) {
        seq[k][copy_pos++] = cur_seq[string_pos];
      }
    }
    seq[k][copy_pos] = '\0';  //Ensure that cur_seq ends with a '\0' character

    char* cur_name = strscn(header_data + 1);
    sname[k] = new char[strlen(cur_name) + 1];
    strcpy(sname[k], cur_name);

    k++;
  }

  N_in = k;

  // Warn if there are only special sequences but no master sequence (consensus seq given if keep[kfirst]==0)
  if (kfirst < 0 || (N_in - N_ss - (keep[kfirst] == 0 ? 1 : 0)) == 0) {
    HH_LOG(ERROR) << "In " << __FILE__ << ":" << __LINE__ << ": " << __func__ << ":" << std::endl;
    HH_LOG(ERROR) << "\tMSA file " << entry->name << " contains no master sequence!" << std::endl;
    exit(1);
  }

  // Set name, longname, fam
  // longname, name and family were not set by '#...' line yet -> extract from first sequence
  if (!*name) {
    char* ptr;
    strncpy(longname, sname[kfirst], DESCLEN - 1);  // longname is name of first sequence
    longname[DESCLEN - 1] = '\0';
    strncpy(name, sname[kfirst], NAMELEN - 1);  // Shortname is first word of longname...
    name[NAMELEN - 1] = '\0';
    ptr = strcut(name);              // ...until first white-space character
    if (ptr && islower(ptr[0]) && ptr[1] == '.' && isdigit(ptr[2]))  //Scop family code present as second word?
        {
      lwrstr(name);                  // Transform upper case to lower case
      strcut(ptr);  // Non-white-space characters until next white-space character..
      strcpy(fam, ptr);                    // ...are the SCOP familiy code
    } else if (name[0] == 'P' && name[1] == 'F' && isdigit(name[2])
        && isdigit(name[3]))  //Pfam code
            {
      strcpy(fam, name);                    // set family name = Pfam code
    }
  }

  HH_LOG(DEBUG) << "Read " << entry->name << " with " << N_in
                          << " sequences\n";
  HH_LOG(DEBUG1) << "Query sequence for alignment has number "
                           << kfirst << " (0 is first)\n";
}

/////////////////////////////////////////////////////////////////////////////////////
// Convert ASCII in seq[k][l] to int (0-20) in X[k][i], throw out all insert states, record their number in I[k][i]
// and store sequences to be displayed in seq[k]
/////////////////////////////////////////////////////////////////////////////////////
void Alignment::Compress(const char infile[], const char cons, const int maxres,
                         const int maxcol, const int par_M, const int Mgaps) {
  int i;                  //Index for match state (first=1)
  int l;                  //Postion in alignment incl. gaps (first=1)
  int k;                  //Index for sequences (first=0)
  int a;                  //amino acid index
  char c;
  int unequal_lengths = 0;  //k: seq k doesn't have same number of match states as seq 0 => WARNING
  short unsigned int h[MAXSEQ];  //points to next character in seq[k] to be written

  int M = par_M;

  // for (k=0;k<N_in; ++k) printf("k=%i >%s\n%s\n",k,sname[k],seq[k]); // DEBUG

  // Initialize
  for (k = 0; k < N_in; ++k) {
    I[k][0] = 0;
    X[k][0] = ANY;
  }

  if (M == 1) {
    HH_LOG(DEBUG1)
        << "Using match state assignment by capital letters (a2m/a3m format)\n";
  } else if (M == 2) {
    HH_LOG(DEBUG1)
        << "Using percentage-rule match state assignment\n";
  } else if (M == 3) {
    HH_LOG(DEBUG1)
        << "Using residues of first sequence as match states\n";
  }

  // Warn, if there are gaps in a single sequence
  if (N_in - N_ss == 1 && M != 2 && strchr(seq[kfirst] + 1, '-') != NULL) {
    HH_LOG(DEBUG)
        << "File "
        << infile
        << " has a single sequence containing gaps, which will be ignored.\nIf you want to treat the gaps as match states, use the '-M 100' option."
        << std::endl;
  }

  // Too few match states?
  if (M == 1) {
    int match_states = strcount(seq[kfirst] + 1, 'A', 'Z')
        + strcount(seq[kfirst] + 1, '-', '-');
    if (match_states < 6) {
      if (N_in - N_ss <= 1) {
        M = 3;  // if only single sequence in input file, use M=3 (match states by first seq)
        HH_LOG(DEBUG)
            << "Single sequence in file " << infile
            << " contains only " << match_states
            << " match_states! Switching to option -M first\n seq="
            << seq[kfirst] << std::endl;
      }
    } else {
      HH_LOG(DEBUG)
          << "Single sequence in file " << infile << " contains only "
          << match_states
          << " match_states! Switching to option -M first\n seq=" << seq[kfirst]
          << std::endl;
    }
  }

  // Create matrices X and I with amino acids represented by integer numbers
  switch (M) {

    /////////////////////////////////////////////////////////////////////////////////////
    // a2m/a3m format: match states capital case, inserts lower case, delete states '-', inserted gaps '.'
    // The confidence values for ss prediction are interpreted as follows: 0-9:match states(!)  '-' :match state  '.':insert
    case 1:
    default:

      // Warn if alignment is meant to be -M first or -M <%> instead of A2M/A3M
      // Seed/query sequence contains a gap ...
      if (strchr(seq[kfirst], '-')) {
        unsigned int len = strlen(seq[kfirst]) - 1;
        for (k = 1; k < N_in; ++k) {
          if (keep[k] && strcount(seq[k], 'a', 'z'))
            break;
          if (strlen(seq[k]) != len)
            k = N_in;
        }
        // ... but alignment contains no lower case residue
        if (k >= N_in) {
          HH_LOG(WARNING)
              << "Input alignment "
              << infile
              << " looks like aligned FASTA instead of A2M/A3M format. Consider using '-M first' or '-M 50'"
              << std::endl;
        }
      }

      // Remove '.' characters from seq[k]
      for (k = 0; k < N_in; ++k) {
        char* ptrS = seq[k];       // pointer to source: character in seq[k]
        char* ptrD = seq[k];         // pointer to destination: seq[k]
        while (1)                   // omit '.' symbols
        {
          if (*ptrS != '.') {
            *ptrD = *ptrS;
            ptrD++;
          }  //leave out '.' symbols
          if (!*ptrS)
            break;
          ptrS++;
        }
      }
      L = maxres - 2;  // needed because L=imin(L,i)
      for (k = 0; k < N_in; ++k) {
        i = 1;
        l = 1;  // start at i=1, not i=0!
        if (keep[k])  //skip >ss_dssp, >ss_pred, >ss_conf, >aa_... sequences
        {
          while ((c = seq[k][l++]))  // assign residue to c at same time
          {
            if (c >= 'a' && c <= 'z')
              I[k][i - 1]++;  //insert state = lower case character
            else if (c != '.')  //match state = upper case character
                {
              X[k][i] = aa2i(c);
              I[k][i] = 0;
              ++i;
            }
          }
        } else if (k == kss_dssp || k == kss_pred)  // does alignment contain sequence of secondary structure states?
            {
          while ((c = seq[k][l++]))  // assign residue to c at same time
            if (c != '.' && !(c >= 'a' && c <= 'z'))
              X[k][i++] = ss2i(c);  //match state = upper case character
        } else if (k == ksa_dssp)  // does alignment contain sequence of prediction confidence values?
            {
          while ((c = seq[k][l++]))  // assign residue to c at same time
            if (c != '.' && !(c >= 'a' && c <= 'z'))
              X[k][i++] = sa2i(c);  //match state = upper case character
        } else if (k == kss_conf)  // does alignment contain sequence of prediction confidence values?
            {
          while ((c = seq[k][l++]))  // assign residue to c at same time
            if (c != '.')
              X[k][i++] = cf2i(c);  //match state = 0-9 or '-'
        } else if (k == kfirst)  // does alignment contain sequence of prediction confidence values?
            {
          while ((c = seq[k][l++]))  // assign residue to c at same time
            if (c != '.') {
              X[k][i] = aa2i(c);
              I[k][i] = 0;
              ++i;
            }
        } else
          continue;
        i--;
        if (L != i && L != maxres - 2 && !unequal_lengths)
          unequal_lengths = k;   //sequences have different lengths

        L = imin(L, i);
      }
      if (unequal_lengths)
        break;

      //Replace GAP with ENDGAP for all end gaps
      for (k = 0; k < N_in; ++k) {
        if (!keep[k])
          continue;
        for (i = 1; i <= L && X[k][i] == GAP; ++i)
          X[k][i] = ENDGAP;
        for (i = L; i >= 1 && X[k][i] == GAP; i--)
          X[k][i] = ENDGAP;
      }

      for (i = 1; i <= L; ++i)
        this->l[i] = i;  //assign column indices to match states
      if (L <= 0) {
        HH_LOG(ERROR) << "In " << __FILE__ << ":" << __LINE__ << ": " << __func__ << ":" << std::endl;
        HH_LOG(ERROR) << "\tAlignment in " << infile << " contains no match states. Consider using -M first or -M <int> option" << std::endl;
        exit(1);
      }

      if (L == maxres - 2) {
        HH_LOG(WARNING)
            << "Number of match columns too large. Only first " << L
            << " match columns will be kept!\n";
        break;
      }

      HH_LOG(DEBUG) << "Alignment in " << infile << " contains " << L
                              << " match states\n";
      break;

      /////////////////////////////////////////////////////////////////////////////////////
      // gap-rule assignment of match states
    case 2: {
      int nl[NAA + 2];  //nl[a] = number of seq's with amino acid a at position l
      float* percent_gaps = new float[maxcol];  //percentage of gaps in column k (with weighted sequences)

      //determine number of columns L in alignment
      L = strlen(seq[kfirst]) - 1;

      // Conversion to integer representation, checking for unequal lengths and initialization
      if (nres == NULL)
        nres = new int[N_in];
      for (k = 0; k < N_in; ++k) {
        if (!keep[k])
          continue;
        int nr = 0;
        wg[k] = 0;
        nres[k] = 0;
        for (l = 1; l <= L; l++) {
          X[k][l] = aa2i(seq[k][l]);
          if (X[k][l] < NAA)
            nr++;
        }
        nres[k] = nr;
        if (seq[k][L + 1] != '\0' && !unequal_lengths)
          unequal_lengths = k;
      }
      if (unequal_lengths)
        break;

      // Quick and dirty calculation of the weight per sequence wg[k]
      for (l = 1; l <= L; l++)  // for all positions l in alignment
          {
        int naa = 0;            //number of different amino acids
        for (a = 0; a < 20; ++a)
          nl[a] = 0;
        for (k = 0; k < N_in; ++k)
          if (keep[k])
            nl[(int) X[k][l]]++;
        for (a = 0; a < 20; ++a)
          if (nl[a])
            ++naa;
        if (!naa)
          naa = 1;  //naa=0 when column consists of only gaps and Xs (=ANY)
        for (k = 0; k < N_in; ++k)
          if (keep[k] && X[k][l] < 20)
            wg[k] += 1.0 / float(nl[(int) X[k][l]] * naa * (nres[k] + 30.0));
        // wg[k] += 1.0/float(nl[ (int)X[k][l]]*(nres[k]+30.0));
        // wg[k] += (naa-1.0)/float(nl[ (int)X[k][l]]*(nres[k]+30.0));
      }

      //Replace GAP with ENDGAP for all end gaps
      for (k = 0; k < N_in; ++k) {
        if (!keep[k])
          continue;
        for (i = 1; i <= L && X[k][i] == GAP; ++i)
          X[k][i] = ENDGAP;
        for (i = L; i >= 1 && X[k][i] == GAP; i--)
          X[k][i] = ENDGAP;
      }

      // Add up percentage of gaps
      for (l = 1; l <= L; l++) {
        float res = 0;
        float gap = 0;
        for (k = 0; k < N_in; ++k) {
          if (keep[k]) {
            if (X[k][l] < GAP)
              res += wg[k];   // AA or ANY; Changed from <ANY
            else if (X[k][l] != ENDGAP)
              gap += wg[k];  // else: GAP. ENDGAPs are ignored for counting percentage (multi-domain proteins)
          }
        }
        percent_gaps[l] = 100. * gap / (res + gap);
        HH_LOG(DEBUG1) << "percent gaps[" << l << "]="
                                 << percent_gaps[l] << " first seq:"
                                 << seq[0][l] << "\n";
      }
      // Throw out insert states and keep only match states
      i = 0;
      for (k = 0; k < N_in; ++k) {
        h[k] = 1;
        seq[k][0] = '-';
      }

      for (l = 1; l <= L; l++) {
        // for first and last columns look for slightly shifted columns for nr gaps
        if (percent_gaps[l] <= float(Mgaps)) {
          if (i >= maxres - 2) {
            HH_LOG(WARNING)
                << "Number of match columns too large. Only first "
                << i << " match columns will be kept!\n";
            break;
          }
          ++i;
          this->l[i] = l;
          for (k = 0; k < N_in; ++k) {
            if (keep[k]) {
              seq[k][h[k]++] = MatchChr(seq[k][l]);
              X[k][i] = X[k][l];
              I[k][i] = 0;

              //kfirst might not be copied in the case of consensus sequences
              //so it will be deleted and kfirst will be set to 0
              //kfirst needs to be set to the next following sequence
              if (kfirst == -1) {
                kfirst = k;
              }
            } else if (k == kss_dssp || k == kss_pred) {
              seq[k][h[k]++] = MatchChr(seq[k][l]);
              X[k][i] = ss2i(seq[k][l]);
            } else if (k == ksa_dssp) {
              seq[k][h[k]++] = MatchChr(seq[k][l]);
              X[k][i] = sa2i(seq[k][l]);
            } else if (k == kss_conf) {
              seq[k][h[k]++] = seq[k][l];
              X[k][i] = cf2i(seq[k][l]);
            }
            //consensus sequence, keep[kfirst] == 0
            else if (k == kfirst) {
              kfirst = -1;
            }
          }
        } else {
          for (k = 0; k < N_in; ++k)
            if (keep[k] && X[k][l] < GAP) {
              I[k][i]++;
              seq[k][h[k]++] = InsertChr(seq[k][l]);
            }
        }
      }
      for (k = 0; k < N_in; ++k)
        seq[k][h[k]] = '\0';

      HH_LOG(DEBUG) << "Alignment in " << infile << " contains " << L
                              << " columns and " << i << " match states\n";
      L = i;        //Number of match states
      delete[] percent_gaps;
      break;
    }

      /////////////////////////////////////////////////////////////////////////////////////
      // Using residues of first sequence as match states
    case 3: {
      char* match_state = new char[maxcol];  //1: column assigned to match state 0: insert state

      // Determine number of columns L in alignment
      L = strlen(seq[0] + 1);
      HH_LOG(DEBUG) << "Length of first seq = " << L << std::endl;
      // Check for sequences with unequal lengths
      for (k = 1; k < N_in; ++k)
        if (int(strlen(seq[k] + 1)) != L) {
          unequal_lengths = k;
          break;
        }
      if (unequal_lengths)
        break;

      // Determine match states: seq kfirst has residue at pos l -> match state
      for (l = 1; l <= L; l++)
        if (isalpha(seq[kfirst][l]))
          match_state[l] = 1;
        else
          match_state[l] = 0;
      // Throw out insert states and keep only match states
      for (k = 0; k < N_in; ++k) {
        h[k] = 1;
        seq[k][0] = '-';
      }
      i = 0;
      for (l = 1; l <= L; l++) {
        if (match_state[l])  // does sequence 0 have residue at position l?
        {
          if (i >= maxres - 2) {
            HH_LOG(WARNING)
                << "Number of match columns too large. Only first "
                << i << " match columns will be kept!\n";
            break;
          }
          ++i;
          this->l[i] = l;
          for (k = 0; k < N_in; ++k) {
            if (keep[k]) {
              seq[k][h[k]++] = MatchChr(seq[k][l]);
              X[k][i] = aa2i(seq[k][l]);
              I[k][i] = 0;
            } else if (k == kss_dssp || k == kss_pred) {
              seq[k][h[k]++] = MatchChr(seq[k][l]);
              X[k][i] = ss2i(seq[k][l]);
            } else if (k == ksa_dssp) {
              seq[k][h[k]++] = MatchChr(seq[k][l]);
              X[k][i] = sa2i(seq[k][l]);
            } else if (k == kss_conf) {
              seq[k][h[k]++] = seq[k][l];
              X[k][i] = cf2i(seq[k][l]);
            }
          }
        } else {
          for (k = 0; k < N_in; ++k)
            if (keep[k] && aa2i(seq[k][l]) < GAP) {
              I[k][i]++;
              seq[k][h[k]++] = InsertChr(seq[k][l]);
            }
        }
      }
      for (k = 0; k < N_in; ++k)
        seq[k][h[k]] = '\0';

      HH_LOG(DEBUG) << "Alignment in " << infile << " contains " << L
                              << " columns and " << i << " match states\n";
      L = i;        //Number of match states

      //Replace GAP with ENDGAP for all end gaps
      for (k = 0; k < N_in; ++k) {
        if (!keep[k])
          continue;
        for (i = 1; i <= L && X[k][i] == GAP; ++i)
          X[k][i] = ENDGAP;
        for (i = L; i >= 1 && X[k][i] == GAP; i--)
          X[k][i] = ENDGAP;
      }

      delete[] match_state;

      break;
    }

  }  //end switch()
     /////////////////////////////////////////////////////////////////////////////////////

  // Error
  if (unequal_lengths) {
    HH_LOG(ERROR) << strcut(sname[unequal_lengths]);
    HH_LOG(ERROR) << "Error in " << __FILE__ << ":" << __LINE__
                            << ": " << __func__ << ":" << std::endl;
    HH_LOG(ERROR) << "\tsequences in " << infile
                            << " do not all have the same number of columns, "
                            << std::endl;
    HH_LOG(ERROR) << "\t\ne.g. first sequence and sequence "
                            << sname[unequal_lengths] << ".\n";

    if (M == 1) {
      HH_LOG(ERROR)
          << "Check input format for '-M a2m' option and consider using '-M first' or '-M 50'\n";
    }

    if (!strncmp(infile, "merged A3M", 10)) {
      HH_LOG(ERROR) << "Merged MSA:\n";
      for (k = 0; k <= unequal_lengths; ++k)
        HH_LOG(ERROR) << k << "\n" << sname[k] << "\n" << seq[k]
                                << std::endl;
    }
    exit(1);
  }

  if (L == 0) {
    HH_LOG(ERROR) << "Error in " << __FILE__ << ":" << __LINE__
                            << ": " << __func__ << ":" << std::endl;
    HH_LOG(ERROR) << "\tno match states found in " << infile << "!\n";
    HH_LOG(ERROR)
        << "Better use '-M first' option or reduce gap percentage threshold for match states\n";
    exit(1);
  }

  // Avert user about -add_cons option?
  if (!cons) {
    for (i = 1; i <= L; ++i)
      if (X[kfirst][i] == GAP) {
        HH_LOG(INFO)
            << "NOTE: Use the '-add_cons' option to calculate a consensus sequence as first sequence of the alignment with hhconsensus or hhmake.\n";
        break;
      }
  }

  // DEBUG
  for (k = 0; k < N_in; ++k) {
    if (!display[k])
      continue;
    HH_LOG(DEBUG1) << "k=" << k << " >" << sname[k] << "\n";
    if (k == kss_dssp || k == kss_pred) {
      for (i = 1; i <= L; ++i)
        HH_LOG(DEBUG1) << char(i2ss(X[k][i]));
    } else if (k == kss_conf) {
      for (i = 1; i <= L; ++i)
        HH_LOG(DEBUG1) << char(i2cf(X[k][i]));
    } else if (k == ksa_dssp) {
      for (i = 1; i <= L; ++i)
        HH_LOG(DEBUG1) << char(i2sa(X[k][i]));
    } else {
      for (i = 1; i <= L; ++i)
        HH_LOG(DEBUG1) << char(i2aa(X[k][i]));
      HH_LOG(DEBUG1) << "\n";
      for (i = 1; i <= L; ++i) {
        if (I[k][i] == 0) {
          HH_LOG(DEBUG1) << "-";
        }
        else if (I[k][i] > 9) {
          HH_LOG(DEBUG1) << "X";
        }
        else {
          HH_LOG(DEBUG1) << I[k][i];
        }
      }
    }
    HH_LOG(DEBUG1) << "\n";
  }
}

/////////////////////////////////////////////////////////////////////////////////////
// Copy sequences from HMM into matrix seq[k][l] as ASCII
/////////////////////////////////////////////////////////////////////////////////////
void Alignment::GetSeqsFromHMM(HMM* q) {
  int qk;                 // Index of sequence in HHM
  int k;                  // Index of sequence being read currently (first=0)

  kss_dssp = ksa_dssp = kss_pred = kss_conf = kfirst = -1;
  n_display = 0;
  N_in = 0;
  N_filtered = 0;
  N_ss = 0;
  k = 0;

  for (qk = 0; qk < q->n_seqs; ++qk) {
    if (qk == q->ncons) {
      continue;
    }

    // Create space for residues and paste new sequence in
    seq[k] = new char[strlen(q->seq[qk]) + 1];
    if (!seq[k])
      MemoryError("array for input sequences", __FILE__, __LINE__, __func__);
    strcpy(seq[k], q->seq[qk]);
    X[k] = initX(strlen(q->seq[qk]) + 1);
    if (!X[k])
      MemoryError("array for input sequences", __FILE__, __LINE__, __func__);
    I[k] = new short unsigned int[strlen(q->seq[qk]) + 1];
    if (!I[k])
      MemoryError("array for input sequences", __FILE__, __LINE__, __func__);

    if (qk == q->nss_dssp) {
      display[k] = 2;
      n_display++;
      keep[k] = 0;
      kss_dssp = k;
      N_ss++;
    } else if (qk == q->nsa_dssp) {
      display[k] = 2;
      n_display++;
      keep[k] = 0;
      ksa_dssp = k;
      N_ss++;
    } else if (qk == q->nss_pred) {
      display[k] = 2;
      n_display++;
      keep[k] = 0;
      kss_pred = k;
      N_ss++;
    } else if (qk == q->nss_conf) {
      display[k] = 2;
      n_display++;
      keep[k] = 0;
      kss_conf = k;
      N_ss++;
    } else {
      if (kfirst < 0) {
        display[k] = keep[k] = 2;
        n_display++;
        kfirst = k;
      } else {
        display[k] = keep[k] = 1;
        n_display++;
      }
    }

    // store sequence name
    HH_LOG(DEBUG) << "Reading seq " << q->sname[qk] << " k=" << k
                            << "  n_displ=" << n_display << "  display[k]="
                            << display[k] << " keep[k]=" << keep[k] << "\n";
    sname[k] = new char[strlen(q->sname[qk]) + 1];
    if (!sname[k]) {
      MemoryError("array for sequence names", __FILE__, __LINE__, __func__);
    }
    strcpy(sname[k], q->sname[qk]);

    ++k;
  }

  N_in = k;

  strcpy(longname, q->longname);
  strmcpy(name, q->name, NAMELEN - 1);
  strmcpy(fam, q->fam, NAMELEN - 1);
}

/////////////////////////////////////////////////////////////////////////////////////
// Remove sequences with seq. identity larger than seqid percent (remove the shorter of two) or coverage<cov_thr
/////////////////////////////////////////////////////////////////////////////////////
int Alignment::FilterForDisplay(int max_seqid, const char mark,
                                const float S[20][20], int coverage, int qid,
                                float qsc, int N) {
  if (mark)
    return n_display;
  char *dummy = new char[N_in + 1];
  int seqid;
  n_display = 0;
  if (kss_dssp >= 0)
    display[kss_dssp] = 0;
  if (ksa_dssp >= 0)
    display[ksa_dssp] = 0;
  if (kss_pred >= 0)
    display[kss_pred] = 0;
  if (kss_conf >= 0)
    display[kss_conf] = 0;
  for (seqid = imin(10, max_seqid); n_display < N && seqid <= max_seqid;
      seqid++) {
    for (int k = 0; k < N_in; ++k)
      dummy[k] = display[k];
    n_display = Filter2(dummy, coverage, qid, qsc, 20, seqid, 0, S);
    //printf("Seqid=%3i  n_display=%4i\n",seqid,n_display);
  }
  if (n_display > N) {
    for (int k = 0; k < N_in; ++k)
      dummy[k] = display[k];
    n_display = Filter2(dummy, coverage, qid, qsc, 20, --(--seqid), 0, S);
  }
  for (int k = 0; k < N_in; ++k)
    display[k] = dummy[k];
  if (kss_dssp >= 0) {
    display[kss_dssp] = 1;
    n_display++;
  }
  if (ksa_dssp >= 0) {
    display[ksa_dssp] = 1;
    n_display++;
  }
  if (kss_pred >= 0) {
    display[kss_pred] = 1;
    n_display++;
  }
  if (kss_conf >= 0) {
    display[kss_conf] = 1;
    n_display++;
  }
  delete[] dummy;

  return n_display;
}

/////////////////////////////////////////////////////////////////////////////////////
// Remove sequences with seq. identity larger than seqid percent (remove the shorter of two) or coverage<cov_thr
/////////////////////////////////////////////////////////////////////////////////////
int Alignment::Filter(int max_seqid, const float S[20][20], int coverage,
                      int qid, float qsc, int N) {
  return Filter2(keep, coverage, qid, qsc, 20, max_seqid, N, S);
}

void Alignment::Shrink() {
  char** new_X = new char*[maxseq + 2];
  short unsigned int** new_I = new short unsigned int*[maxseq + 2];;
  char** new_sname = new char*[maxseq + 2];
  char** new_seq = new char*[maxseq + 2];

  char* new_keep = new char[maxseq + 2];
  char* new_display = new char[maxseq + 2];
  float* new_wg = new float[maxseq + 2];

  int new_kss_dssp = -1;
  int new_ksa_dssp = -1;
  int new_kss_pred = -1;
  int new_kss_conf = -1;
  int new_kfirst = -1;

  int new_N_in = N_in;
  int new_k = 0;
  for(int k = 0; k < N_in; k++) {
	if(keep[k] == 0 && k != kss_dssp && k != ksa_dssp && k != kss_pred && k != kss_conf && k != kfirst) {
	  free(X[k]);
	  delete[] I[k];
	  delete[] sname[k];
	  delete[] seq[k];
	  new_N_in--;
	}
	else {
	  new_X[new_k] = X[k];
	  new_I[new_k] = I[k];
	  new_sname[new_k] = sname[k];
	  new_seq[new_k] = seq[k];

	  new_keep[new_k] = keep[k];
	  new_display[new_k] = display[k];
	  new_wg[new_k] = wg[k];

	  if (k == kss_dssp) {
		  new_kss_dssp = new_k;
	  }
	  else if (k == ksa_dssp) {
		  new_ksa_dssp = new_k;
	  }
	  else if (k == kss_pred) {
		  new_kss_pred = new_k;
	  }
	  else if (k == kss_conf) {
		  new_kss_conf = new_k;
	  }
	  else if (k == kfirst) {
		  new_kfirst = new_k;
	  }

	  new_k++;
	}
  }

  delete[] X;
  X= new_X;

  delete[] I;
  I = new_I;

  delete[] sname;
  sname = new_sname;

  delete[] seq;
  seq = new_seq;

  delete[] keep;
  keep = new_keep;

  delete[] display;
  display = new_display;

  delete[] wg;
  wg = new_wg;

  kss_dssp = new_kss_dssp;
  ksa_dssp = new_ksa_dssp;
  kss_pred = new_kss_pred;
  kss_conf = new_kss_conf;
  kfirst = new_kfirst;

  N_in = new_N_in;

  if(ksort != NULL) {
    delete[] ksort;
    ksort = NULL;
  }
  if(first != NULL) {
    delete[] first;
	first = NULL;
  }
  if(last != NULL) {
    delete[] last;
    last = NULL;
  }
  if(nres != NULL) {
	delete[] nres;
	nres = NULL;
  }
}



/////////////////////////////////////////////////////////////////////////////////////
// Select set of representative sequences in the multiple sequence alignment
// Filter criteria:
//   * Remove sequences with coverage of query less than "coverage" percent
//   * Remove sequences with sequence identity to query of less than "qid" percent
//   * If Ndiff==0, remove sequences with seq. identity larger than seqid2(=max_seqid) percent
//   * If Ndiff>0, remove sequences with minimum-sequence-identity filter of between seqid1
//     and seqid2 (%), where the minimum seqid threshold is determined such that,
//     in all column blocks of at least WMIN=25 residues, at least Ndiff sequences are left.
//     This ensures that in multi-domain proteins sequences covering one domain are not
//     removed completely because sequences covering other domains are more diverse.
//
// Allways the shorter of two compared sequences is removed (=> sort sequences by length first).
// Please note: sequence identity of sequence x with y when filtering x is calculated as
// number of residues in sequence x that are identical to an aligned residue in y / number of residues in x
// Example: two sequences x and y are 100% identical in their overlapping region but one overlaps by 10% of its
// length on the left and the other by 20% on the right. Then x has 10% seq.id with y and y has 20% seq.id. with x.
/////////////////////////////////////////////////////////////////////////////////////
int Alignment::Filter2(char keep[], int coverage, int qid, float qsc,
                       int seqid1, int seqid2, int Ndiff,
                       const float S[20][20]) {

  // In the beginnning, keep[k] is 1 for all regular amino acid sequences and 0 for all others (ss_conf, ss_pred,...)
  // In the end, keep[k] will be 1 for all regular representative sequences kept in the alignment, 0 for all others
  // Sequences with keep[k] = 2 will cannot be filtered out and will remain in the alignment.
  // If a consensus sequence exists it has k = kfirst and keep[k] = 0, since it should not enter into the profile calculation.
  char* in = new char[N_in + 1];  // in[k]=1: seq k has been accepted; in[k]=0: seq k has not yet been accepted at current seqid
  char* inkk = new char[N_in + 1];  // inkk[k]=1 iff in[ksort[k]]=1 else 0;
  int* Nmax = new int[L + 2];  // position-dependent maximum-sequence-identity threshold for filtering? (variable used in former version was idmax)
  int* idmaxwin = new int[L + 2];     // minimum value of idmax[i-WFIL,i+WFIL]
  int* seqid_prev = new int[N_in + 1];  // maximum-sequence-identity threshold used in previous round of filtering (with lower seqid)
  int* N = new int[L + 2];  // N[i] number of already accepted sequences at position i
  const int WFIL = 25;            // see previous line

  int diffNmax = Ndiff;    // current  maximum difference of Nmax[i] and Ndiff
  int diffNmax_prev = 0;   // previous maximum difference of Nmax[i] and Ndiff
  int seqid;  // current  maximum value for the position-dependent maximum-sequence-identity thresholds in idmax[]
  int seqid_step = 0;         // previous increment of seqid

  float diff_min_frac;  // minimum fraction of differing positions between sequence j and k needed to accept sequence k
  float qdiff_max_frac = 0.9999 - 0.01 * qid;  // maximum allowable number of residues different from query sequence
  int diff;  // number of differing positions between sequences j and k (counted so far)
  int diff_suff;  // number of differing positions between sequences j and k that would be sufficient
  int qdiff_max;  // maximum number of residues required to be different from query
  int cov_kj;  // upper limit of number of positions where both sequence k and j have a residue
  int first_kj;             // first non-gap position in sequence j AND k
  int last_kj;              // last  non-gap position in sequence j AND k
  int kk, jj;               // indices for sequence from 1 to N_in
  int k, j;                 // kk=ksort[k], jj=ksort[j]
  int i;                    // counts residues
  int n;                    // number of sequences accepted so far

  // Initialize in[k]
  for (n = k = 0; k < N_in; ++k) {
    if (keep[k] == 2) {
      in[k] = 2;
      n++;
    } else
      in[k] = 0;
  }
  // Determine first[k], last[k]?
  if (first == NULL) {
    first = new int[N_in];         // first non-gap position in sequence k
    last = new int[N_in];          // last  non-gap position in sequence k
    for (k = 0; k < N_in; ++k)  // do this for ALL sequences, not only those with in[k]==1 (since in[k] may be display[k])
    {
      for (i = 1; i <= L; ++i)
        if (X[k][i] < NAA)
          break;
      first[k] = i;
      for (i = L; i >= 1; i--)
        if (X[k][i] < NAA)
          break;
      last[k] = i;
    }
  }

  // Determine number of residues nres[k]?
  if (nres == NULL || sizeof(nres) < N_in * sizeof(int)) {
    if (nres)
      delete[] nres;
    nres = new int[N_in];
    for (k = 0; k < N_in; ++k)  // do this for ALL sequences, not only those with in[k]==1 (since in[k] may be display[k])
    {
      int nr = 0;
      for (i = first[k]; i <= last[k]; ++i)
        if (X[k][i] < NAA)
          nr++;
      nres[k] = nr;
      //printf("%20.20s nres=%3i  first=%3i  last=%3i\n",sname[k],nr,first[k],last[k]);
      if (nr == 0)
        keep[k] = 0;
    }
  }

  // Sort sequences according to length; afterwards, nres[ksort[kk]] is sorted by size
  if (ksort == NULL) {
    ksort = new int[N_in];  // never reuse alignment object for new alignment with more sequences
    for (k = 0; k < N_in; ++k)
      ksort[k] = k;
    QSortInt(nres, ksort, kfirst + 1, N_in - 1, -1);  //Sort sequences after kfirst (query) in descending order
  }
  for (kk = 0; kk < N_in; ++kk) {
    inkk[kk] = in[ksort[kk]];
  }

  // Initialize N[i], idmax[i], idprev[i]
  for (i = 1; i < first[kfirst]; ++i)
    N[i] = 0;
  for (i = first[kfirst]; i <= last[kfirst]; ++i)
    N[i] = 1;
  for (i = last[kfirst] + 1; i <= L; ++i)
    N[i] = 0;
  for (i = 1; i <= L; ++i) {
    Nmax[i] = 0;
    idmaxwin[i] = -1;
  }
  for (k = 0; k < N_in; ++k)
    seqid_prev[k] = -1;
  if (Ndiff <= 0 || Ndiff >= N_in) {
    seqid1 = seqid2;
    Ndiff = N_in;
    diffNmax = Ndiff;
  }

  // Check coverage and sim-to-query criteria for each sequence k
  for (k = 0; k < N_in; ++k) {
    if (keep[k] == 0 || keep[k] == 2)
      continue;  // seq k not regular sequence OR is marked sequence
    if (100 * nres[k] < coverage * L) {
      keep[k] = 0;
      continue;
    }  // coverage too low? => reject once and for all

    float qsc_sum = 0.0;

    // Check if score-per-column with query is at least qsc
    if (qsc > -10) {
      float qsc_min = qsc * nres[k];  // minimum total score of seq k with query

      int gapq = 0, gapk = 0;  // number of consecutive gaps in query or k'th sequence at position i
      for (int i = first[k]; i <= last[k]; ++i) {
        if (X[k][i] < 20) {
          gapk = 0;
          if (X[kfirst][i] < 20) {
            gapq = 0;
            qsc_sum += S[(int) X[kfirst][i]][(int) X[k][i]];
          } else if (X[kfirst][i] == ANY)
            // Treat score of X with other amino acid as 0.0
            continue;
          else if (gapq++)
            qsc_sum -= PLTY_GAPEXTD;
          else
            qsc_sum -= PLTY_GAPOPEN;
        } else if (X[k][i] == ANY)
          // Treat score of X with other amino acid as 0.0
          continue;
        else if (X[kfirst][i] < 20) {
          gapq = 0;
          if (gapk++)
            qsc_sum -= PLTY_GAPEXTD;
          else
            qsc_sum -= PLTY_GAPOPEN;
        }
      }
//        printf("k=%3i qsc=%6.2f\n",k,qsc_sum);
      if (qsc_sum < qsc_min) {
        keep[k] = 0;
        continue;
      }  // too different from query? => reject once and for all
    }

    //Check if sequence similarity with query at least qid?
    if (qdiff_max_frac < 0.999) {
      qdiff_max = int(qdiff_max_frac * nres[k] + 0.9999);
      //      printf("k=%-4i  nres=%-4i  qdiff_max=%-4i first=%-4i last=%-4i",k,nres[k],qdiff_max,first[k],last[k]);
      diff = 0;
      for (int i = first[k]; i <= last[k]; ++i)
        // enough different residues to reject based on minimum qid with query? => break
        if (X[k][i] < 20 && X[k][i] != X[kfirst][i] && ++diff >= qdiff_max)
          break;
      //      printf("  diff=%4i\n",diff);
      if (diff >= qdiff_max) {
        keep[k] = 0;
        continue;
      }  // too different from query? => reject once and for all
    }
//      printf("  qsc=%6.2f     qid=%6.2f  \n",qsc_sum/nres[k],100.0*(1.0-(float)(diff)/nres[k]));
  }

  // If no sequence left, issue warning and put back first real sequence into alignment
  int nn=0;
  for (k = 0; k < N_in; ++k)
    if (keep[k] > 0)
      nn++;

  if (nn == 0) {
    for (k = 0; k < N_in; k++) {
      if (display[k] != 2) {
        keep[k] = 1;
        break;
      }
    }
    if (keep[k] == 1) {
      HH_LOG(WARNING) << "Warning in " << __FILE__ << ":" << __LINE__
                                << ": " << __func__ << ":" << std::endl;
      HH_LOG(WARNING)
          << "\tFiltering removed all sequences in alignment " << name
          << ". Inserting back first sequence.\n";
    } else if (display[kfirst] == 2) {  // the only sequence in the alignment is the consensus sequence :-(
      HH_LOG(WARNING) << "Warning in " << __FILE__ << ":" << __LINE__
                                << ": " << __func__ << ":" << std::endl;
      HH_LOG(WARNING)
          << "\tAlignment "
          << name
          << " contains no sequence except consensus sequence. Using consensus sequence for searching.\n";
    } else {
      const char unknown[] = "'unknown'";
      char details[100];
      sprintf(details, "The alingment %s does not contain any sequences.",
              name);
      FormatError(unknown, __FILE__, __LINE__, __func__, details);
    }
  }

  // If min required seqid larger than max required seqid, return here without doing pairwise seqid filtering
  if (seqid1 > seqid2)
    return nn;

  // Successively increment idmax[i] at positons where N[i]<Ndiff
  seqid = seqid1;
  while (seqid <= seqid2) {
    char stop = 1;
    // Update Nmax[i]
    diffNmax_prev = diffNmax;
    diffNmax = 0;
    for (i = 1; i <= L; ++i) {
      int max = 0;
      for (j = imax(1, imin(L - 2 * WFIL + 1, i - WFIL));
          j <= imin(L, imax(2 * WFIL, i + WFIL)); ++j)
        if (N[j] > max)
          max = N[j];
      if (Nmax[i] < max)
        Nmax[i] = max;
      if (Nmax[i] < Ndiff) {
        stop = 0;
        idmaxwin[i] = seqid;
        if (diffNmax < Ndiff - Nmax[i])
          diffNmax = Ndiff - Nmax[i];
      }
    }

    //printf("seqid=%3i  diffNmax_prev= %-4i   diffNmax= %-4i   n=%-5i  N_in-N_ss=%-5i\n",seqid,diffNmax_prev,diffNmax,n,N_in-N_ss);
    if (stop)
      break;

    // Loop over all candidate sequences kk (-> k)
    for (kk = 0; kk < N_in; ++kk) {
      if (inkk[kk])
        continue;   // seq k already accepted
      k = ksort[kk];
      if (!keep[k])
        continue;  // seq k is not regular aa sequence or already suppressed by coverage or qid criterion
      if (keep[k] == 2) {
        inkk[kk] = 2;
        continue;
      }  // accept all marked sequences (no n++, since this has been done already)

      // Calculate max-seq-id threshold seqidk for sequence k (as maximum over idmaxwin[i])
      if (seqid >= 100) {
        in[k] = inkk[kk] = 1;
        n++;
        continue;
      }

      float seqidk = seqid1;
      for (i = first[k]; i <= last[k]; ++i)
        if (idmaxwin[i] > seqidk)
          seqidk = idmaxwin[i];

      if (seqid == seqid_prev[k])
        continue;  // sequence has already been rejected at this seqid threshold => reject this time

      seqid_prev[k] = seqid;
      diff_min_frac = 0.9999 - 0.01 * seqidk;  // min fraction of differing positions between sequence j and k needed to accept sequence k

      // Loop over already accepted sequences
      for (jj = 0; jj < kk; ++jj) {
        if (!inkk[jj])
          continue;
        j = ksort[jj];

        first_kj = imax(first[k], first[j]);
        last_kj = imin(last[k], last[j]);
        cov_kj = last_kj - first_kj + 1;
        diff_suff = int(diff_min_frac * imin(nres[k], cov_kj) + 0.999);  // nres[j]>nres[k] anyway because of sorting
        diff = 0;
        const simd_int * XK = (simd_int *) X[k];
        const simd_int * XJ = (simd_int *) X[j];
        const int first_kj_simd = first_kj / (VECSIZE_INT * 4);
        const int last_kj_simd = last_kj / (VECSIZE_INT * 4) + 1;
        // coverage correction for simd
        // because we do not always hit the right start with simd.
        // This works because all sequence vector are initialized with GAPs so the sequnces is surrounded by GAPs
        const int first_diff_simd_scalar = std::abs(
            first_kj_simd * (VECSIZE_INT * 4) - first_kj);
        const int last_diff_simd_scalar = std::abs(
            last_kj_simd * (VECSIZE_INT * 4) - (last_kj + 1));

        cov_kj += (first_diff_simd_scalar + last_diff_simd_scalar);

        // _mm_set1_epi8 pseudo-instruction is slow!
        const simd_int NAAx16 = simdi8_set(NAA - 1);
        for (int i = first_kj_simd; i < last_kj_simd && diff < diff_suff; ++i) {
          // None SIMD function
          // enough different residues to accept? => break
          // if (X[k][i] >= NAA || X[j][i] >= NAA)
          //    cov_kj--;
          // else if (X[k][i] != X[j][i] && ++diff >= diff_suff)
          //    break; // accept (k,j)

          const simd_int NO_AA_K = simdi8_gt(XK[i], NAAx16);  // pos without amino acid in seq k
          const simd_int NO_AA_J = simdi8_gt(XJ[i], NAAx16);  // pos without amino acid in seq j

          // Compute 16 bits indicating positions with GAP, ANY or ENDGAP in seq k or j
          // int _mm_movemask_epi8(__m128i a) creates 16-bit mask from most significant bits of
          // the 16 signed or unsigned 8-bit integers in a and zero-extends the upper bits.
          int res = simdi8_movemask(simdi_or(NO_AA_K, NO_AA_J));

          cov_kj -= NumberOfSetBits(res);  // subtract positions that should not contribute to coverage

          // Compute 16 bit mask that indicates positions where k and j have identical residues
          int c = simdi8_movemask(simdi8_eq(XK[i], XJ[i]));

          // Count positions where  k and j have different amino acids, which is equal to 16 minus the
          //  number of positions for which either j and k are equal or which contain ANY, GAP, or ENDGAP
          diff += (VECSIZE_INT * 4) - NumberOfSetBits(c | res);

        }

        //dissimilarity < acceptace threshold? Reject!
        if (diff < diff_suff && float(diff) <= diff_min_frac * cov_kj)
          break;
      }

      // did loop reach end? => accept k. Otherwise reject k (the shorter of the two)
      if (jj >= kk) {
        in[k] = inkk[kk] = 1;
        n++;
        for (i = first[k]; i <= last[k]; ++i)
          N[i]++;  // update number of sequences at position i
      }
    }  // End Loop over all candidate sequences kk

    // Increment seqid
    seqid_step = imax(1, imin(5, diffNmax / (diffNmax_prev - diffNmax + 1) * seqid_step / 2));
    seqid += seqid_step;
  }  // End Loop over seqid

  HH_LOG(DEBUG) << n << " out of " << N_in - N_ss
                          << " sequences passed filter (";
  if (coverage) {
    HH_LOG(DEBUG) << coverage << "% min coverage, ";
  }
  if (qid) {
    HH_LOG(DEBUG) << qid << "% min sequence identity to query, ";
  }
  if (qsc > -10) {
    HH_LOG(DEBUG) << qsc << " bits min score per column to query, ";
  }
  if (Ndiff < N_in && Ndiff > 0) {
    HH_LOG(DEBUG)
        << "up to " << seqid
        << "% position-dependent max pairwise sequence identity)\n";
  } else {
    HH_LOG(DEBUG) << seqid1 << "% max pairwise sequence identity)\n";
  }

  for (k = 0; k < N_in; ++k)
    keep[k] = in[k];
  delete[] in;
  delete[] inkk;
  //  delete[] idmax;
  delete[] Nmax;
  delete[] idmaxwin;
  delete[] seqid_prev;
  delete[] N;
  return n;
}

/////////////////////////////////////////////////////////////////////////////////////
// Filter out all sequences below a minimum score per column with profile qcore
/////////////////////////////////////////////////////////////////////////////////////
int Alignment::FilterWithCoreHMM(char in[], float coresc, HMM* qcore,
                                 const float* pb) {
  int k;     // count sequences in alignment
  int i;     // column in query alignment
  int a;     // amino acid (0..19)
  int n = 1;   // number of sequences that passed filter
  float** logodds = new float*[L + 1];  // log-odds ratios for HMM qcore
  char gap;  // 1: previous state in seq k was a gap  0: previous state in seq k was an amino acid
  float score;  // score of sequence k aligned with qcore

  for (i = 1; i <= L; ++i)
    logodds[i] = new float[21];

  // Determine first[k], last[k]?
  if (first == NULL) {
    first = new int[N_in];  // first non-gap position in sequence k
    last = new int[N_in];  // last  non-gap position in sequence k
    for (k = 0; k < N_in; ++k)  // do this for ALL sequences, not only those with in[k]==1 (since in[k] may be display[k])
        {
      for (i = 1; i <= L; ++i)
        if (X[k][i] < NAA)
          break;
      first[k] = i;
      for (i = L; i >= 1; i--)
        if (X[k][i] < NAA)
          break;
      last[k] = i;
    }
  }

  // Determine number of residues nres[k]?
  if (nres == NULL) {
    nres = new int[N_in];
    for (k = 0; k < N_in; ++k)  // do this for ALL sequences, not only those with in[k]==1 (since in[k] may be display[k])
        {
      int nr = 0;
      for (i = first[k]; i <= last[k]; ++i)
        if (X[k][i] < NAA)
          nr++;
      nres[k] = nr;
      //    printf("%20.20s nres=%3i  first=%3i  last=%3i\n",sname[k],nr,f,l);
    }
  }

  // Precalculate the log-odds for qcore
  for (i = 1; i <= L; ++i) {
    for (a = 0; a < NAA; ++a)
      logodds[i][a] = fast_log2(qcore->p[i][a] / pb[a]);
    logodds[i][ANY] = -0.5;  // half a bit penalty for X

//       printf("         A     R     N     D     C     Q     E     G     H     I     L     K     M     F     P     S     T     W     Y     V\n");
//       printf("%6i ",i);
//       for (a=0; a<20; ++a) fprintf(stdout,"%5.1f ",100 * qcore->f[i][a]);
//       printf("\n");
//       printf("       ");
//       for (a=0; a<20; ++a) fprintf(stdout,"%5.1f ",100 * qcore->g[i][a]);
//       printf("\n");
//       printf("       ");
//       for (a=0; a<20; ++a) fprintf(stdout,"%5.1f ",100 * qcore->p[i][a]);
//       printf("\n");
//       printf("       ");
//       for (a=0; a<20; ++a) fprintf(stdout,"%5.1f ",100*pb[a]);
//       printf("\n");
//       printf("       ");
//       for (a=0; a<20; ++a) fprintf(stdout,"%5.2f ",fast_log2(qcore->p[i][a]/pb[a]));
//       printf("\n");
  }

  // Main loop: test all sequences k
  for (k = kfirst + 1; k < N_in; ++k) {
    if (!in[k])
      continue;  // if in[k]==0 sequence k will be suppressed directly

    // float score_M=0.0;
    // float score_prev=0.0;

    // Calculate score of sequence k with core HMM
    score = 0;
    gap = 0;
    for (i = first[k]; i <= last[k]; ++i) {
      // score_M=0.0;
      if (X[k][i] <= ANY)       // current state is Match
          {
        // score_M=logodds[i][ (int)X[k][i]];
        score += logodds[i][(int) X[k][i]];
        if (gap)
          score += qcore->tr[i][D2M];
        else
          score += qcore->tr[i][M2M];
        gap = 0;
      } else if (X[k][i] == GAP)  // current state is Delete (ignore ENDGAPs)
          {
        if (gap)
          score += qcore->tr[i][D2D];
        else
          score += qcore->tr[i][M2D];
        gap = 1;
      }
      if (I[k][i])
        score += qcore->tr[i][M2I] + (I[k][i] - 1) * qcore->tr[i][I2I]
            + qcore->tr[i][I2M];
//        if (k==2) printf("i=%3i %c:%c   score_M=%6.2f   score=%6.2f  score_sum=%6.2f \n",i,i2aa(X[kfirst][i]),i2aa(X[k][i]),score_M,score-score_prev,score);
      // score_prev=score;
    }

    printf("k=%3i score=%6.2f\n", k, score);
    if (score < nres[k] * coresc)
      in[k] = 0;
    else
      n++;          // reject sequence k?
  }
  for (i = 1; i <= L; ++i)
    delete[] logodds[i];
  delete[] logodds;
  return n;
}

/////////////////////////////////////////////////////////////////////////////////////
// Filter alignment to given diversity/Neff
/////////////////////////////////////////////////////////////////////////////////////
void Alignment::FilterNeff(char use_global_weights, const char mark,
                           const char cons, const char showcons,
                           const int maxres, const int max_seqid,
                           const int coverage, const float Neff,
                           const float* pb, const float S[20][20],
                           const float Sim[20][20]) {
  const float TOLX = 0.01;
  const float TOLY = 0.02;
  char keep_orig[N_in + 1];
  for (int k = 0; k < N_in; ++k)
    keep_orig[k] = keep[k];
  float x0 = -1.0;
  float x1 = +1.0;
  float x = 0.0;
  float y0;
  float y1 = 1.0;
  float y;

  HMM q;
  FrequenciesAndTransitions(&q, use_global_weights, mark, cons, showcons,
                            maxres, pb, Sim);
  y = y0 = q.Neff_HMM;
  if (fabs(Neff - y0) < TOLY) {
    return;
  }
  if (y0 < Neff) {
    HH_LOG(DEBUG) << "Diversity of alignment before Neff filter "
                            << y0 << " is below target diversity Neff = "
                            << Neff << std::endl;
    return;
  }

  y1 = filter_by_qsc(x1, use_global_weights, mark, cons, showcons, maxres,
                     max_seqid, coverage, keep_orig, pb, S, Sim);
  if (fabs(Neff - y1) < TOLY) {
    return;
  }

  // Contract interval (x0,x1) for par.sc around target value: Neff(x0) > Neff_target > Neff(x1)
  do {
    const float w = 0.5;  // mixture coefficient
    x = w * (0.5 * (x0 + x1))
        + (1 - w) * (x0 + (Neff - y0) * (x1 - x0) / (y1 - y0));  // mixture of bisection and linear interpolation
    y = filter_by_qsc(x, use_global_weights, mark, cons, showcons, maxres,
                      max_seqid, coverage, keep_orig, pb, S, Sim);
    if (y > Neff) {
      x0 = x;
      y0 = y;
    } else {
      x1 = x;
      y1 = y;
    }
  } while (fabs(Neff - y) > TOLY && x1 - x0 > TOLX);

  // Write filtered alignment WITH insert states (lower case) to alignment file
  HH_LOG(DEBUG) << "Found Neff=" << y << " at filter threshold qsc="
                          << x << std::endl;
}

float Alignment::filter_by_qsc(float qsc, char use_global_weights,
                               const char mark, const char cons,
                               const char showcons, const int maxres,
                               const int max_seqid, const int coverage,
                               char* keep_orig, const float* pb,
                               const float S[20][20], const float Sim[20][20]) {
  HMM q;
  for (int k = 0; k < N_in; ++k)
    keep[k] = keep_orig[k];
  Filter2(keep, coverage, 0, qsc, max_seqid + 1, max_seqid, 0, S);
  FrequenciesAndTransitions(&q, use_global_weights, mark, cons, showcons,
                            maxres, pb, Sim);  // Might be sped up by calculating wg and calling only Amino_acid_frequencies_and_transitions_from_M_state(q,in);
  return q.Neff_HMM;
}

/////////////////////////////////////////////////////////////////////////////////////
// Calculate AA frequencies q->p[i][a] and transition probabilities q->tr[i][a] from alignment
/////////////////////////////////////////////////////////////////////////////////////
void Alignment::FrequenciesAndTransitions(HMM* q, char use_global_weights,
                                          const char mark, const char cons,
                                          const char showcons, const int maxres,
                                          const float* pb,
                                          const float Sim[20][20], char* in,
                                          bool time) {
  int k;                // index of sequence
  int i;                // position in alignment
  int a;                // amino acid (0..19)
  int ni[NAA + 3];        // number of times amino acid a occurs at position i
  int naa;              // number of different amino acids

  //if (time) { ElapsedTimeSinceLastCall("begin freq and trans"); }

  //Delete name and seq matrices of old HMM q
  if (!q->dont_delete_seqs)  // don't delete sname and seq if flat copy to hit object has been made
  {
    for (k = 0; k < q->n_seqs; k++)
      delete[] q->sname[k];
    for (k = 0; k < q->n_seqs; k++)
      delete[] q->seq[k];
  } else  // Delete all not shown sequences (lost otherwise)
  {
    if (q->n_seqs > q->n_display) {
      for (k = q->n_display; k < q->n_seqs; k++)
        delete[] q->sname[k];
      for (k = q->n_display; k < q->n_seqs; k++)
        delete[] q->seq[k];
    }
  }

  HH_LOG(DEBUG)
      << "Calculating position-dependent weights on subalignments\n";

  if (in == NULL)
    in = keep;  //why not in declaration?

  if (N_filtered > 1) {
    for (k = 0; k < N_in; ++k)
      wg[k] = 1e-6;  // initialized wg[k] with tiny pseudocount
    // Calculate global weights
    for (i = 1; i <= L; ++i)  // for all positions i in alignment
        {
      for (a = 0; a < 20; ++a)
        ni[a] = 0;
      for (k = 0; k < N_in; ++k)
        if (in[k])
          ni[(int) X[k][i]]++;
      naa = 0;
      for (a = 0; a < 20; ++a)
        if (ni[a])
          naa++;
      if (!naa)
        naa = 1;  //naa=0 when column consists of only gaps and Xs (=ANY)
      for (k = 0; k < N_in; ++k)
        if (in[k] && X[k][i] < 20)
          wg[k] += 1.0 / float(ni[(int) X[k][i]] * naa * (nres[k] + 30.0));
      // wg[k] += 1.0/float(ni[ (int)X[k][i]]*(nres[k]+30.0));
      // wg[k] += (naa-1.0)/float(ni[ (int)X[k][i]]*(nres[k]+30.0));
      // ensure that each residue of a short sequence contributes as much as a residue of a long sequence:
      // contribution is proportional to one over sequence length nres[k] plus 30.
    }
    NormalizeTo1(wg, N_in);
    //if (time) { ElapsedTimeSinceLastCall("Calc global weights"); }

    // Do pos-specific sequence weighting and calculate amino acid frequencies and transitions
    for (k = 0; k < N_in; ++k)
      X[k][0] = ENDGAP;  // make sure that sequences ENTER subalignment j for j=1
    for (k = 0; k < N_in; ++k)
      X[k][L + 1] = ENDGAP;  // does it have an influence?

    Amino_acid_frequencies_and_transitions_from_M_state(q, use_global_weights,
                                                        in, maxres, pb);  // use subalignments of seqs with residue in i
    Transitions_from_I_state(q, in, maxres);  // use subalignments of seqs with insert in i
    Transitions_from_D_state(q, in, maxres);  // use subalignments of seqs with delete in i. Must be last of these three calls if par.wg==1!
    //if (time) { ElapsedTimeSinceLastCall("Do pos-specific sequence weighting and calculate amino acid frequencies and transitions"); }
  } else  // N_filtered==1
  {
    //use first useful sequence (MM,JS 24.11.2013: removed bug that used consensus seq instead of first real seq)
    for (k = 0; k < N_in; k++) {
      if (in[k]) {
        break;
      }
    }

    X[k][0] = X[k][L + 1] = ANY;  // (to avoid unallowed access within loop)
    q->Neff_HMM = 1.0f;

    // for all positions i in alignment
    for (i = 0; i <= L + 1; ++i) {
      q->Neff_M[i] = 1.0f;
      q->Neff_I[i] = q->Neff_D[i] = 0.0f;
      for (a = 0; a < 20; ++a)
        q->f[i][a] = 0.0;
      if (X[k][i] < ANY)
        q->f[i][(unsigned int) X[k][i]] = 1.0;
      else
        for (a = 0; a < 20; ++a)
          q->f[i][a] = pb[a];
      q->tr[i][M2M] = 0;
      q->tr[i][M2I] = -100000.0;
      q->tr[i][M2D] = -100000.0;
      q->tr[i][I2M] = -100000.0;
      q->tr[i][I2I] = -100000.0;
      q->tr[i][D2M] = -100000.0;
      q->tr[i][D2D] = -100000.0;
    }
    q->tr[0][I2M] = 0;
    q->tr[L][I2M] = 0;
    q->tr[0][D2M] = 0;
    q->Neff_M[0] = q->Neff_I[0] = q->Neff_D[0] = 99.999;  // Neff_av[0] is used for calculation of transition pseudocounts for the start state
  }

  if (Log::reporting_level() >= DEBUG) {
    HH_LOG(DEBUG) << "Matches:\n";
    HH_LOG(DEBUG) << "col\tNeff\tnseqs\n";
    for (i = 1; i <= imin(L, 100); ++i)
      HH_LOG(DEBUG) << i << "\t" << q->Neff_M[i] << "\t" << nseqs[i]
                              << std::endl;

    HH_LOG(DEBUG) << "Inserts:";
    HH_LOG(DEBUG) << "col\tNeff\tnseqs" << std::endl;
    for (i = 1; i <= imin(L, 100); ++i)
      HH_LOG(DEBUG) << i << "\t" << q->Neff_I[i] << "\t" << nseqs[i]
                              << std::endl;

    HH_LOG(DEBUG) << "Deletes:\n";
    HH_LOG(DEBUG) << "col\tNeff\tnseqs\n";
    for (i = 1; i <= imin(L, 100); ++i)
      HH_LOG(DEBUG) << i << "\t" << q->Neff_D[i] << "\t" << nseqs[i]
                              << std::endl;
  }

  // Copy column information into HMM q
  q->L = L;
  q->N_in = N_in;
  q->N_filtered = N_filtered;
  for (i = 1; i <= L; ++i)
    q->l[i] = l[i];

  // Set names in HMM q
  if (strlen(q->name) == 0)
    strcpy(q->name, name);
  if (strlen(q->longname) == 0)
    strcpy(q->longname, longname);
  if (strlen(q->fam) == 0)
    strcpy(q->fam, fam);
  ScopID(q->cl, q->fold, q->sfam, q->fam);  // derive superfamily, fold and class code from family name
  strcpy(q->file, file);   // Store basename of alignment file name in q->file

  // Copy sequences to be displayed into HMM
  q->nss_dssp = q->nsa_dssp = q->nss_pred = q->nss_conf = q->nfirst = -1;
  int n = 0;
  if (kss_dssp >= 0)
    q->nss_dssp = n++;  // copy dssp sequence?
  if (ksa_dssp >= 0)
    q->nsa_dssp = n++;  // copy dssp sequence?
  if (kss_pred >= 0)
    q->nss_pred = n++;  // copy psipred sequence?
  if (kss_conf >= 0)
    q->nss_conf = n++;  // copy confidence value sequence?

  //if (time) { ElapsedTimeSinceLastCall("Copy to HMM"); }

  // Calculate consensus sequence?
  if (showcons || cons) {
    float maxw;
    int maxa;
    if (showcons) {
      // Reserve space for consensus/conservation sequence as Q-T alignment mark-up
      q->ncons = n++;
      q->sname[q->ncons] = new char[10];
      if (!q->sname[q->ncons]) {
        MemoryError("array of names for displayed sequences", __FILE__,
        __LINE__,
                    __func__);
      }
      strcpy(q->sname[q->ncons], "Consensus");
      q->seq[q->ncons] = new char[L + 2];
      if (!q->seq[q->ncons]) {
        MemoryError("array of names for displayed sequences", __FILE__,
        __LINE__,
                    __func__);
      }
    }
    if (cons) {
      // Reserve space for consensus sequence as first sequence in alignment
      q->nfirst = n++;
      kfirst = -1;
      q->sname[q->nfirst] = new char[strlen(name) + 11];
      if (!q->sname[q->nfirst]) {
        MemoryError("array of names for displayed sequences", __FILE__,
        __LINE__,
                    __func__);
      }
      strcpy(q->sname[q->nfirst], name);
      strcat(q->sname[q->nfirst], "_consensus");
      q->seq[q->nfirst] = new char[L + 2];
      if (!q->seq[q->nfirst]) {
        MemoryError("array of names for displayed sequences", __FILE__,
        __LINE__,
                    __func__);
      }
    }
    // Calculate consensus amino acids using similarity matrix
    for (i = 1; i <= L; ++i) {
      maxw = 0;
      maxa = ANY;
      for (a = 0; a < 20; ++a)
        if (q->f[i][a] - pb[a] > maxw) {
          maxw = q->f[i][a] - pb[a];
          maxa = a;
        }

      if (showcons) {
        maxw = 0.0;
        for (int b = 0; b < 20; b++)
          maxw += q->f[i][b] * Sim[maxa][b] * Sim[maxa][b];
        maxw *= q->Neff_M[i] / (q->Neff_HMM + 1);  // columns with many gaps don't get consensus symbol
        if (maxw > 0.6)
          q->seq[q->ncons][i] = uprchr(i2aa(maxa));
        else if (maxw > 0.4)
          q->seq[q->ncons][i] = lwrchr(i2aa(maxa));
        else
          q->seq[q->ncons][i] = 'x';
      }
      if (cons)
        q->seq[q->nfirst][i] = uprchr(i2aa(maxa));
    }
    if (showcons) {
      q->seq[q->ncons][0] = '-';
      q->seq[q->ncons][L + 1] = '\0';
    }
    if (cons) {
      q->seq[q->nfirst][0] = '-';
      q->seq[q->nfirst][L + 1] = '\0';
    }
  }

  //if (time) { ElapsedTimeSinceLastCall("Calc consensus sequence"); }

  // Copy sequences to be displayed from alignment to HMM
  for (k = 0; k < N_in; ++k) {
    int nn;
    if (display[k]) {
      if (n >= q->maxseqdis) { //MAXSEQDIS) {
        if (mark)
          HH_LOG(WARNING) << "Maximum number " << q->maxseqdis << " of sequences for display of alignment exceeded" << std::endl;
        break;
      }
      if (k == kss_dssp)
        nn = q->nss_dssp;  // copy dssp sequence to nss_dssp
      else if (k == ksa_dssp)
        nn = q->nsa_dssp;
      else if (k == kss_pred)
        nn = q->nss_pred;
      else if (k == kss_conf)
        nn = q->nss_conf;
      else if (k == kfirst)
        nn = q->nfirst = n++;
      else
        nn = n++;
//        strcut(sname[k],"  "); // delete rest of name line beginning with two spaces "  " // Why this?? Problem for pdb seqs without chain
      q->sname[nn] = new char[strlen(sname[k]) + 1];
      if (!q->sname[nn]) {
        MemoryError("array of names for displayed sequences", __FILE__,
        __LINE__,
                    __func__);
      }
      strcpy(q->sname[nn], sname[k]);
      q->seq[nn] = new char[strlen(seq[k]) + 1];
      if (!q->seq[nn]) {
        MemoryError("array of names for displayed sequences", __FILE__,
        __LINE__,
                    __func__);
      }
      strcpy(q->seq[nn], seq[k]);
    }
  }
  q->n_display = n;  // how many sequences to be displayed in alignments?
  q->n_seqs = n;

  // Copy secondary structure information into HMM
  if (kss_dssp >= 0)
    for (i = 1; i <= L; ++i)
      q->ss_dssp[i] = X[kss_dssp][i];
  if (ksa_dssp >= 0)
    for (i = 1; i <= L; ++i)
      q->sa_dssp[i] = X[ksa_dssp][i];
  if (kss_pred >= 0) {
    for (i = 1; i <= L; ++i)
      q->ss_pred[i] = X[kss_pred][i];
    if (kss_conf >= 0)
      for (i = 1; i <= L; ++i)
        q->ss_conf[i] = X[kss_conf][i];
    else
      for (i = 1; i <= L; ++i)
        q->ss_conf[i] = 5;
  }

  q->lamda = 0.0;
  q->mu = 0.0;

  q->trans_lin = 0;  // transition probs in log space
  q->has_pseudocounts = false;
  q->dont_delete_seqs = false;
  q->divided_by_local_bg_freqs = false;

  //if (time) { ElapsedTimeSinceLastCall("Copy sequences and SS"); }

  // Debug: print occurrence of amino acids for each position i
  if (Log::reporting_level() >= DEBUG) {
    HH_LOG(DEBUG) << "Effective number of sequences exp(entropy) = "
                            << q->Neff_HMM << std::endl;

    HH_LOG(DEBUG) << "Matr: ";
    for (a = 0; a < 20; ++a)
      HH_LOG(DEBUG) << 100 * pb[a] << " ";
    HH_LOG(DEBUG)
        << "\nAmino acid frequencies without pseudocounts:\n";
    HH_LOG(DEBUG)
        << "         A    R    N    D    C    Q    E    G    H    I    L    K    M    F    P    S    T    W    Y    V\n";
    for (i = 1; i <= L; ++i) {
      HH_LOG(DEBUG) << i << " ";
      for (a = 0; a < 20; ++a)
        HH_LOG(DEBUG) << 100 * q->f[i][a] << " ";
      HH_LOG(DEBUG) << std::endl;
    }
    HH_LOG(DEBUG) << std::endl;

    HH_LOG(DEBUG)
        << "Listing transition probabilities without pseudocounts:"
        << std::endl;
    HH_LOG(DEBUG)
        << "   i    M->M   M->I   M->D   I->M   I->I   D->M   D->D  Neff_M Neff_I Neff_D"
        << std::endl;
    for (i = 0; i <= L; ++i) {
      HH_LOG(DEBUG) << i << "  " << fpow2(q->tr[i][M2M]) << " "
                              << fpow2(q->tr[i][M2I]) << " "
                              << fpow2(q->tr[i][M2D]) << " ";
      HH_LOG(DEBUG) << fpow2(q->tr[i][I2M]) << " "
                              << fpow2(q->tr[i][I2I]) << " ";
      HH_LOG(DEBUG) << fpow2(q->tr[i][D2M]) << " "
                              << fpow2(q->tr[i][D2D]) << "  ";
      HH_LOG(DEBUG) << q->Neff_M[i] << " " << q->Neff_I[i] << " "
                              << q->Neff_D[i] << std::endl;
    }
  }
}

/////////////////////////////////////////////////////////////////////////////////////
// Calculate freqs q->f[i][a] and transitions q->tr[i][a] (a=MM,MI,MD) with pos-specific subalignments
// Pos-specific weights are calculated like in "GetPositionSpecificWeights()"
/////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////
// Calculate freqs q->f[i][a] and transitions q->tr[i][a] (a=MM,MI,MD) with pos-specific subalignments
// Pos-specific weights are calculated like in "GetPositionSpecificWeights()"
/////////////////////////////////////////////////////////////////////////////////////
void Alignment::Amino_acid_frequencies_and_transitions_from_M_state(
    HMM* q, char use_global_weights, char* in, const int maxres,
    const float* pb) {
  // Calculate position-dependent weights wi[k] for each i.
  // For calculation of weights in column i use sub-alignment
  // over sequences which have a *residue* in column i (no gap, no end gap)
  // and over columns where none of these sequences has an end gap.
  // This is done by updating the arrays n[j][a] at each step i-1->i while letting i run from 1 to L.
  // n[j][a] = number of occurences of amino acid a at column j of the subalignment,
  //        => only columns with n[j][ENDGAP]=0 are contained in the subalignment!
  // If no sequences enter or leave the subalignment at the step i-1 -> i (i.e. change=0)
  // then the old values wi[k], Neff[i-1], and ncol are used for the new position i.
  // Index a can be an amino acid (0-19), ANY=20, GAP=21, or ENDGAP=22

  int k;                      // index of sequence
  int i, j;                    // position in alignment
  int a;                      // amino acid (0..19)
  int** n = NULL;  // n[j][a] = number of seq's with some non-gap amino acid at column i AND residue a at position j
  float** w_contrib = NULL;  // weight contribution of amino acid a at pos. j to weight of a sequence
  float wi[MAXSEQ];  // weight of sequence k in column i, calculated from subalignment i
  float Neff[maxres];         // diversity of subalignment i
  int nseqi = 0;                // number of sequences in subalignment i
  int ncol = 0;                // number of columns j that contribute to Neff[i]
  char change;  // has the set of sequences in subalignment changed? 0:no  1:yes
  float sum;

  int* naa = new int[L + 1];   // number of different amino acids

  // Allocate memory for f[j]
  float** f = malloc_matrix<float>(L+1, NAA+3);


  // Initialization
  if (use_global_weights == 1)  // If global weights
    for (k = 0; k < N_in; ++k)
      wi[k] = wg[k];
  else {
    unsigned int NAA_VECSIZE = ((NAA+ 3 + VECSIZE_INT - 1) / VECSIZE_INT) * VECSIZE_INT; // round NAA+3 up to next multiple of VECSIZE_INT
    n = new int*[L + 2];
    for (j = 1; j <= L; ++j)
      n[j] = (int *) malloc_simd_int(NAA_VECSIZE * sizeof(int));
    for (j = 1; j <= L; ++j)
      for (a = 0; a < NAA + 3; ++a)
        n[j][a] = 0;
    w_contrib = new float*[L + 2];
    for (j = 1; j <= L; ++j){
      w_contrib[j] = (float *) malloc_simd_int(NAA_VECSIZE * sizeof(float));
      memset(w_contrib[j], 0,NAA_VECSIZE * sizeof(int));
    }
  }
  q->Neff_HMM = 0.0f;
  Neff[0] = 0.0;  // if the first column has no residues (i.e. change==0), Neff[i]=Neff[i-1]=Neff[0]

  //////////////////////////////////////////////////////////////////////////////////////////////
  // Main loop through alignment columns
  for (i = 1; i <= L; ++i)  // Calculate wi[k] at position i as well as Neff[i]
      {

    if (use_global_weights == 0) {

      change = 0;
      // Check all sequences k and update n[j][a] and ri[j] if necessary
      for (k = 0; k < N_in; ++k) {
        if (!in[k])
          continue;

        // Update amino acid and GAP / ENDGAP counts for sequences with AA in i-1 and GAP/ENDGAP in i or vice versa
        if (X[k][i - 1] >= ANY && X[k][i] < ANY) {  // ... if sequence k was NOT included in i-1 and has to be included for column i
          change = 1;
          nseqi++;
          for (int j = 1; j <= L; ++j)
            n[j][(int) X[k][j]]++;
        } else if (X[k][i - 1] < ANY && X[k][i] >= ANY) {  // ... if sequence k WAS included in i-1 and has to be thrown out for column i
          change = 1;
          nseqi--;
          for (int j = 1; j <= L; ++j)
            n[j][(int) X[k][j]]--;
        }
      }  //end for (k)
      nseqs[i] = nseqi;

      // Only if subalignment changed we need to update weights wi[k] and Neff[i]
      if (change) {
        // We gained a factor ~8.0 for the following computation of weights
        // and profile by exchanging the inner two loops (j, k => k, j)
        // and precomputing the weight contributions w_contrib[j][a].
        // M. Steinegger and J. Soeding (29 July 2014)

        // Initialize weights and numbers of residues for subalignment i
        ncol = 0;
        for (k = 0; k < N_in; ++k)
          wi[k] = 1E-8;  // for pathological alignments all wi[k] can get 0;

        // Find min and max borders between which > fraction MAXENDGAPFRAC of sequences in subalignment contain an aa
        int jmin;
        int jmax;
        for (jmin = 1; jmin <= L && n[jmin][ENDGAP] > MAXENDGAPFRAC * nseqi;
            ++jmin) {
        };
        for (jmax = L; jmax >= 1 && n[jmax][ENDGAP] > MAXENDGAPFRAC * nseqi;
            --jmax) {
        };
        ncol = jmax - jmin + 1;

        // Check whether number of columns in subalignment is sufficient
        if (ncol < NCOLMIN) {
          // Take global weights
          for (k = 0; k < N_in; ++k)
            wi[k] = (in[k] && X[k][i] < ANY)? wg[k] : 0.0f;
        } else {
          // Count number of different amino acids in column j
          for (j = jmin; j <= jmax; ++j)
            for (naa[j] = a = 0; a < ANY; ++a)
              naa[j] += (n[j][a] ? 1 : 0);

            // Compute the contribution of amino acid a to the weight
            //for (a = 0; a < ANY; ++a)
            //      w_contrib[j][a] = (n[j][a] > 0) ? 1.0/ float(naa[j]*n[j][a]): 0.0f;
            for (j = jmin; j <= jmax; ++j) {
              simd_float naa_j = simdi32_i2f(simdi32_set(naa[j]));
              const simd_int *nj = (const simd_int *) n[j];
              const int aa_size = (ANY + VECSIZE_INT - 1) / VECSIZE_INT;
              for (a = 0; a < aa_size; ++a) {
                simd_float nja = simdi32_i2f(simdi_load(nj + a));
                simd_float res = simdf32_mul(nja, naa_j);
                simdf32_store(w_contrib[j] + (a * VECSIZE_INT), simdf32_rcp(res));
              }
              for (a = ANY; a < NAA + 3; ++a)
                w_contrib[j][a] = 0.0f;  // set non-amino acid values to 0 to avoid checking in next loop for X[k][j]<ANY
            }

            // Compute pos-specific weights wi[k]
            for (k = 0; k < N_in; ++k) {
              if (!in[k] || X[k][i] >= ANY)
                continue;
              for (j = jmin; j <= jmax; ++j)  // innermost, time-critical loop; O(L*N_in*L)
                wi[k] += w_contrib[j][(int) X[k][j]];
            }
        }

        // Calculate Neff[i]
        Neff[i] = 0.0;

        // Allocate and reset amino acid frequencies
        for (j = jmin; j <= jmax; ++j)
            memset(f[j], 0, ANY * sizeof(float));

        // Update f[j][a]
        for (k = 0; k < N_in; ++k) {
          if (!in[k] || X[k][i] >= ANY)
            continue;
          for (j = jmin; j <= jmax; ++j)  // innermost loop; O(L*N_in*L)
            f[j][(int) X[k][j]] += wi[k];
        }

        // Add contributions to Neff[i]
        for (j = jmin; j <= jmax; ++j) {
          NormalizeTo1(f[j], NAA);
          for (a = 0; a < 20; ++a)
            if (f[j][a] > 1E-10)
              Neff[i] -= f[j][a] * fast_log2(f[j][a]);
        }

        if (ncol > 0)
          Neff[i] = fpow2(Neff[i] / ncol);
        else
          Neff[i] = 1.0;

      }

      else  //no update was necessary; copy values for i-1
      {
        Neff[i] = Neff[i - 1];
      }
    }

    // Calculate amino acid frequencies q->f[i][a] from weights wi[k]
    for (a = 0; a < 20; ++a)
      q->f[i][a] = 0;
    for (k = 0; k < N_in; ++k)
      if (in[k])
        q->f[i][(int) X[k][i]] += wi[k];
    NormalizeTo1(q->f[i], NAA, pb);

    // Calculate transition probabilities from M state
    q->tr[i][M2M] = q->tr[i][M2D] = q->tr[i][M2I] = 0.0;
    for (k = 0; k < N_in; ++k)  //for all sequences
        {
      if (!in[k])
        continue;
      //if input alignment is local ignore transitions from and to end gaps
      if (X[k][i] < ANY)            //current state is M
          {
        if (I[k][i])             //next state is I
          q->tr[i][M2I] += wi[k];
        else if (X[k][i + 1] <= ANY)  //next state is M
          q->tr[i][M2M] += wi[k];
        else if (X[k][i + 1] == GAP)  //next state is D
          q->tr[i][M2D] += wi[k];
      }
    }  // end for(k)
    // Normalize and take log
    sum = q->tr[i][M2M] + q->tr[i][M2I] + q->tr[i][M2D] + FLT_MIN;
    q->tr[i][M2M] = flog2(q->tr[i][M2M] / sum);
    q->tr[i][M2I] = flog2(q->tr[i][M2I] / sum);
    q->tr[i][M2D] = flog2(q->tr[i][M2D] / sum);

//       for (k=0; k<N_in; ++k) if (in[k]) w[k][i]=wi[k];
  }
  // end loop through alignment columns i
  //////////////////////////////////////////////////////////////////////////////////////////////

  if (use_global_weights == 0) {
    // Delete n[][]
    for (j = 1; j <= L; ++j)
      free(n[j]);
    delete[] (n);
    // Delete w_contib[][]
    for (j = 1; j <= L; ++j)
      free(w_contrib[j]);
    delete[] (w_contrib);
  }

  // Delete f[j]
  free(f);
  delete[] naa;

  q->tr[0][M2M] = 0;
  q->tr[0][M2I] = -100000;
  q->tr[0][M2D] = -100000;
  q->tr[L][M2M] = 0;
  q->tr[L][M2I] = -100000;
  q->tr[L][M2D] = -100000;
  q->Neff_M[0] = 99.999;  // Neff_av[0] is used for calculation of transition pseudocounts for the start state

  // Set emission probabilities of zero'th (begin) state and L+1st (end) state to background probabilities
  for (a = 0; a < 20; ++a)
    q->f[0][a] = q->f[L + 1][a] = pb[a];

  // Assign Neff_M[i] and calculate average over alignment, Neff_M[0]
  if (use_global_weights == 1) {
    for (i = 1; i <= L; ++i) {
      float sum = 0.0f;
      for (a = 0; a < 20; ++a)
        if (q->f[i][a] > 1E-10)
          sum -= q->f[i][a] * fast_log2(q->f[i][a]);
      q->Neff_HMM += fpow2(sum);
    }
    q->Neff_HMM /= L;
    float Nlim = fmax(10.0, q->Neff_HMM + 1.0);    // limiting Neff
    float scale = flog2((Nlim - q->Neff_HMM) / (Nlim - 1.0));  // for calculating Neff for those seqs with inserts at specific pos
    for (i = 1; i <= L; ++i) {
      float w_M = -1.0 / N_filtered;
      for (k = 0; k < N_in; ++k)
        if (in[k] && X[k][i] <= ANY)
          w_M += wg[k];
      if (w_M < 0)
        q->Neff_M[i] = 1.0;
      else
        q->Neff_M[i] = Nlim - (Nlim - 1.0) * fpow2(scale * w_M);
//        fprintf(stderr,"M  i=%3i  ncol=---  Neff_M=%5.2f  Nlim=%5.2f  w_M=%5.3f  Neff_M=%5.2f\n",i,q->Neff_HMM,Nlim,w_M,q->Neff_M[i]);
    }
  } else {
    for (i = 1; i <= L; ++i) {
      q->Neff_HMM += Neff[i];
      q->Neff_M[i] = Neff[i];
      if (q->Neff_M[i] == 0) {
        q->Neff_M[i] = 1;
      }
    }
    q->Neff_HMM /= L;
  }

  // printf("Neff = %g\n",q->Neff_HMM);

  return;
}

//void Alignment::Amino_acid_frequencies_and_transitions_from_M_state(
//    HMM* q, char use_global_weights, char* in, const int maxres,
//    const float* pb) {
//  // Calculate position-dependent weights wi[k] for each i.
//  // For calculation of weights in column i use sub-alignment
//  // over sequences which have a *residue* in column i (no gap, no end gap)
//  // and over columns where none of these sequences has an end gap.
//  // This is done by updating the arrays n[j][a] at each step i-1->i while letting i run from 1 to L.
//  // n[j][a] = number of occurences of amino acid a at column j of the subalignment,
//  //        => only columns with n[j][ENDGAP]=0 are contained in the subalignment!
//  // If no sequences enter or leave the subalignment at the step i-1 -> i (i.e. change=0)
//  // then the old values wi[k], Neff[i-1], and ncol are used for the new position i.
//  // Index a can be an amino acid (0-19), ANY=20, GAP=21, or ENDGAP=22
//
//  int k;                      // index of sequence
//  int i, j;                    // position in alignment
//  int a;                      // amino acid (0..19)
//  int naa;                    // number of different amino acids
//  int** n;  // n[j][a] = number of seq's with some residue at column i AND a at position j
//  float wi[MAXSEQ];  // weight of sequence k in column i, calculated from subalignment i
//  float Neff[maxres];         // diversity of subalignment i
//  int nseqi = 0;                // number of sequences in subalignment i
//  int ncol = 0;              // number of columns j that contribute to Neff[i]
//  char change;  // has the set of sequences in subalignment changed? 0:no  1:yes
//  float sum;
//
//  // Global weights?
//  if (use_global_weights == 1)
//    for (k = 0; k < N_in; ++k)
//      wi[k] = wg[k];
//
//  // Initialization
//  q->Neff_HMM = 0.0f;
//  Neff[0] = 0.0;  // if the first column has no residues (i.e. change==0), Neff[i]=Neff[i-1]=Neff[0]
//
//  n = new int*[L + 2];
//  for (j = 1; j <= L; ++j)
//    n[j] = new int[NAA + 3];
//  for (j = 1; j <= L; ++j)
//    for (a = 0; a < NAA + 3; ++a)
//      n[j][a] = 0;
//
//  //////////////////////////////////////////////////////////////////////////////////////////////
//  // Main loop through alignment columns
//  for (i = 1; i <= L; ++i)  // Calculate wi[k] at position i as well as Neff[i]
//      {
//
//    if (use_global_weights == 0) {
//
//      change = 0;
//      // Check all sequences k and update n[j][a] and ri[j] if necessary
//      for (k = 0; k < N_in; ++k) {
//        if (!in[k])
//          continue;
//        if (X[k][i - 1] >= ANY && X[k][i] < ANY) {  // ... if sequence k was NOT included in i-1 and has to be included for column i
//          change = 1;
//          nseqi++;
//          for (int j = 1; j <= L; ++j)
//            n[j][(int) X[k][j]]++;
//        } else if (X[k][i - 1] < ANY && X[k][i] >= ANY) {  // ... if sequence k WAS included in i-1 and has to be thrown out for column i
//          change = 1;
//          nseqi--;
//          for (int j = 1; j <= L; ++j)
//            n[j][(int) X[k][j]]--;
//        }
//      }  //end for (k)
//      nseqs[i] = nseqi;
//
//      // If subalignment changed: update weights wi[k] and Neff[i]
//      if (change) {
//        // Initialize weights and numbers of residues for subalignment i
//        ncol = 0;
//        for (k = 0; k < N_in; ++k)
//          wi[k] = 1E-8;  // for pathological alignments all wi[k] can get 0;
//
//        // sum wi[k] over all columns j and sequences k of subalignment
//        for (j = 1; j <= L; ++j) {
//          //  do at least a fraction MAXENDGAPFRAC of sequences in subalignment contain an end gap in j?
//          if (n[j][ENDGAP] > MAXENDGAPFRAC * nseqi)
//            continue;
//          naa = 0;
//          for (a = 0; a < 20; ++a)
//            if (n[j][a])
//              naa++;
//          if (naa == 0)
//            continue;
//          ncol++;
//          for (k = 0; k < N_in; ++k) {
//            if (in[k] && X[k][i] < ANY && X[k][j] < ANY) {
//              //if (!n[j][ (int)X[k][j]]) {fprintf(stderr,"Error in "<<par.argv[0]<<": Mi=%i: n[%i][X[%i]]=0! (X[%i]=%i)\n",i,j,k,k,X[k][j]);}
//              wi[k] += 1.0 / float(n[j][(int) X[k][j]] * naa);
//            }
//          }
//        }
//
//        // Check whether number of columns in subalignment is sufficient
//        if (ncol < NCOLMIN)
//          // Take global weights
//          for (k = 0; k < N_in; ++k)
//            if (in[k] && X[k][i] < ANY)
//              wi[k] = wg[k];
//            else
//              wi[k] = 0.0;
//
//        // Calculate Neff[i]
//        Neff[i] = 0.0;
//
//        // New Code -- 07.08.2014
//        // Find min and max borders of MSA between which > fraction MAXENDGAPFRAC of sequences in subalignment contain an aa
//        int jmin = 0;
//        int jmax = 0;
//        for (jmin = 1; jmin <= L && n[jmin][ENDGAP] > MAXENDGAPFRAC * nseqi;
//            ++jmin) {
//        };
//        for (jmax = L; jmax >= 1 && n[jmax][ENDGAP] > MAXENDGAPFRAC * nseqi;
//            --jmax) {
//        };
//
//        // Allocate and reset amino acid frequencies
//        float* f[L + 1];
//        for (j = 0; j <= L; ++j)
//          f[j] = new float[NAA + 3];
//
//        for (j = jmin; j <= jmax; ++j)
//          for (a = 0; a < ANY; ++a)
//            f[j][a] = 0.0;
//
//        // Update f[j][a]
//        for (k = 0; k < N_in; ++k) {
//          if (!in[k] || X[k][i] >= ANY)
//            continue;
//          for (j = jmin; j <= jmax; ++j)
//            f[j][(int) X[k][j]] += wi[k];
//        }
//
//        // Add contributions to Neff[i]
//        for (j = jmin; j <= jmax; ++j) {
//          NormalizeTo1(f[j], NAA);
//          for (a = 0; a < 20; ++a)
//            if (f[j][a] > 1E-10)
//              Neff[i] -= f[j][a] * fast_log2(f[j][a]);
//        }
//
//        for (int j = 0; j <= L; ++j) {
//          delete [] f[j];
//        }
////        delete [] f;
//
//
////        old code fragment
////        float fj [NAA + 3];
////				for (j = 1; j <= L; ++j) {
////					//  do at least a fraction MAXENDGAPFRA of sequences in subalignment contain an end gap in j?
////					if (n[j][ENDGAP] > MAXENDGAPFRAC * nseqi)
////						continue;
////					for (a = 0; a < 20; ++a)
////						fj[a] = 0;
////					for (k = 0; k < N_in; ++k)
////						if (in[k] && X[k][i] < ANY && X[k][j] < ANY)
////							fj[(int) X[k][j]] += wi[k];
////					NormalizeTo1(fj, NAA);
////					for (a = 0; a < 20; ++a)
////						if (fj[a] > 1E-10)
////							Neff[i] -= fj[a] * fast_log2(fj[a]);
////				}
//
//				if (ncol > 0)
//					Neff[i] = pow(2.0, Neff[i] / ncol);
//				else
//					Neff[i] = 1.0;
//      }
//
//      else  //no update was necessary; copy values for i-1
//      {
//        Neff[i] = Neff[i - 1];
//      }
//    }
//
//    // Calculate amino acid frequencies q->f[i][a] from weights wi[k]
//    for (a = 0; a < 20; ++a)
//      q->f[i][a] = 0;
//    for (k = 0; k < N_in; ++k)
//      if (in[k])
//        q->f[i][(int) X[k][i]] += wi[k];
//    NormalizeTo1(q->f[i], NAA, pb);
//
//    // Calculate transition probabilities from M state
//    q->tr[i][M2M] = q->tr[i][M2D] = q->tr[i][M2I] = 0.0;
//    for (k = 0; k < N_in; ++k)  //for all sequences
//        {
//      if (!in[k])
//        continue;
//      //if input alignment is local ignore transitions from and to end gaps
//      if (X[k][i] < ANY)            //current state is M
//          {
//        if (I[k][i])             //next state is I
//          q->tr[i][M2I] += wi[k];
//        else if (X[k][i + 1] <= ANY)  //next state is M
//          q->tr[i][M2M] += wi[k];
//        else if (X[k][i + 1] == GAP)  //next state is D
//          q->tr[i][M2D] += wi[k];
//      }
//    }  // end for(k)
//       // Normalize and take log
//    sum = q->tr[i][M2M] + q->tr[i][M2I] + q->tr[i][M2D] + FLT_MIN;
//    q->tr[i][M2M] = log2(q->tr[i][M2M] / sum);
//    q->tr[i][M2I] = log2(q->tr[i][M2I] / sum);
//    q->tr[i][M2D] = log2(q->tr[i][M2D] / sum);
//
////       for (k=0; k<N_in; ++k) if (in[k]) w[k][i]=wi[k];
//  }
//
//  // end loop through alignment columns i
//  //////////////////////////////////////////////////////////////////////////////////////////////
//
//  // delete n[][]
//  for (j = 1; j <= L; ++j)
//    delete[] (n[j]);
//  delete[] (n);
//
//  q->tr[0][M2M] = 0;
//  q->tr[0][M2I] = -100000;
//  q->tr[0][M2D] = -100000;
//  q->tr[L][M2M] = 0;
//  q->tr[L][M2I] = -100000;
//  q->tr[L][M2D] = -100000;
//  q->Neff_M[0] = 99.999;  // Neff_av[0] is used for calculation of transition pseudocounts for the start state
//
//  // Set emission probabilities of zero'th (begin) state and L+1st (end) state to background probabilities
//  for (a = 0; a < 20; ++a)
//    q->f[0][a] = q->f[L + 1][a] = pb[a];
//
//  // Assign Neff_M[i] and calculate average over alignment, Neff_M[0]
//  if (use_global_weights == 1) {
//    for (i = 1; i <= L; ++i) {
//      float sum = 0.0f;
//      for (a = 0; a < 20; ++a)
//        if (q->f[i][a] > 1E-10)
//          sum -= q->f[i][a] * fast_log2(q->f[i][a]);
//      q->Neff_HMM += pow(2.0, sum);
//    }
//    q->Neff_HMM /= L;
//    float Nlim = fmax(10.0, q->Neff_HMM + 1.0);    // limiting Neff
//    float scale = log2((Nlim - q->Neff_HMM) / (Nlim - 1.0));  // for calculating Neff for those seqs with inserts at specific pos
//    for (i = 1; i <= L; ++i) {
//      float w_M = -1.0 / N_filtered;
//      for (k = 0; k < N_in; ++k)
//        if (in[k] && X[k][i] <= ANY)
//          w_M += wg[k];
//      if (w_M < 0)
//        q->Neff_M[i] = 1.0;
//      else
//        q->Neff_M[i] = Nlim - (Nlim - 1.0) * fpow2(scale * w_M);
////        fprintf(stderr,"M  i=%3i  ncol=---  Neff_M=%5.2f  Nlim=%5.2f  w_M=%5.3f  Neff_M=%5.2f\n",i,q->Neff_HMM,Nlim,w_M,q->Neff_M[i]);
//    }
//  } else {
//    for (i = 1; i <= L; ++i) {
//      q->Neff_HMM += Neff[i];
//      q->Neff_M[i] = Neff[i];
//      if (q->Neff_M[i] == 0) {
//        q->Neff_M[i] = 1;
//      }
//    }
//    q->Neff_HMM /= L;
//  }
//
//  return;
//}

/////////////////////////////////////////////////////////////////////////////////////
// Calculate transitions q->tr[i][a] (a=DM,DD) with pos-specific subalignments
/////////////////////////////////////////////////////////////////////////////////////
void Alignment::Transitions_from_I_state(HMM* q, char* in, const int maxres) {
  // Calculate position-dependent weights wi[k] for each i.
  // For calculation of weights in column i use sub-alignment
  // over sequences which have a INSERT in column i
  // and over columns where none of these sequences has an end gap.
  // This is done by calculating the arrays n[j][a] and rj[j] at each step i-1->i while letting i run from 1 to L.
  // n[j][a] = number of occurences of amino acid a at column j of the subalignment,
  //        => only columns with n[j][ENDGAP]=0 are contained in the subalignment!
  // If no sequences enter or leave the subalignment at the step i-1 -> i (i.e. change=0)
  // then the old values wi[k], Neff[i-1], and ncol are used for the new position i.
  // Index a can be an amino acid (0-19), ANY=20, GAP=21, or ENDGAP=22

  int k;                      // index of sequence
  int i, j;                    // position in alignment
  int a;                      // amino acid (0..19)
  int naa;                    // number of different amino acids
  int** n;  // n[j][a] = number of seq's with some residue at column i AND a at position j
  float wi[MAXSEQ];  // weight of sequence k in column i, calculated from subalignment i
  float Neff[maxres];         // diversity of subalignment i
  int nseqi;                  // number of sequences in subalignment i
  int ncol;                  // number of columns j that contribute to Neff[i]
  float fj[NAA + 3];            // to calculate entropy
  float sum;
  float Nlim = 0.0;                 // only for global weights
  float scale = 0.0;            // only for global weights

  // Global weights?
  //if (par.wg==1)
  if (1) {
    for (k = 0; k < N_in; ++k)
      wi[k] = wg[k];
    Nlim = fmax(10.0, q->Neff_HMM + 1.0);    // limiting Neff
    scale = flog2((Nlim - q->Neff_HMM) / (Nlim - 1.0));  // for calculating Neff for those seqs with inserts at specific pos
  }

  // Initialization
  n = new int*[L + 2];
  for (j = 1; j <= L; ++j)
    n[j] = new int[NAA + 3];

  //////////////////////////////////////////////////////////////////////////////////////////////
  // Main loop through alignment columns
  for (i = 1; i <= L; ++i)  // Calculate wi[k] at position i as well as Neff[i]
      {
    //if (par.wg==0) // local weights?
    if (0) {

      // Calculate n[j][a] and ri[j]
      nseqi = 0;
      for (k = 0; k < N_in; ++k) {
        if (in[k] && I[k][i] > 0) {
          if (nseqi == 0)  // Initialize only if inserts present! Otherwise O(L*L) even for single sequences!
              {
            // Initialization of n[j][a]
            for (j = 1; j <= L; ++j)
              for (a = 0; a < NAA + 3; ++a)
                n[j][a] = 0;
          }
          nseqi++;
          for (int j = 1; j <= L; ++j)
            n[j][(int) X[k][j]]++;
        }
      }  //end for (k)
      nseqs[i] = nseqi;

      // If there is no sequence in subalignment j ...
      if (nseqi == 0) {
        ncol = 0;
        Neff[i] = 0.0;  // effective number of sequence = 0!
        q->tr[i][I2M] = -100000;
        q->tr[i][I2I] = -100000;
        continue;
      }

      // update weights wi[k] and Neff[i]
//        if (1)
      {
        // Initialize weights and numbers of residues for subalignment i
        ncol = 0;
        for (k = 0; k < N_in; ++k)
          wi[k] = 1E-8;  // for pathological alignments all wi[k] can get 0;

        // sum wi[k] over all columns j and sequences k of subalignment
        for (j = 1; j <= L; ++j) {
          if (n[j][ENDGAP] > MAXENDGAPFRAC * nseqi)
            continue;
          naa = 0;
          for (a = 0; a < 20; ++a)
            if (n[j][a])
              naa++;
          if (naa == 0)
            continue;
          ncol++;
          for (k = 0; k < N_in; ++k) {
            if (in[k] && I[k][i] > 0 && X[k][j] < ANY) {
              if (!n[j][(int) X[k][j]]) {
                HH_LOG(ERROR) << "In " << __FILE__ << ":" << __LINE__ << ": " << __func__ << ":" << std::endl;
                HH_LOG(ERROR) << "\tIi=" << i << ": n[" << j << "][X[" << k << "][" << j << "]]=0! " <<
                    "(X[" << k << "][" << j << "]=" << X[k][j] <<")" << std::endl;
              }
              wi[k] += 1.0 / float(n[j][(int) X[k][j]] * naa);
            }
          }
        }

        // Check whether number of columns in subalignment is sufficient
        if (ncol >= NCOLMIN) {
          // Take global weights
          for (k = 0; k < N_in; ++k) {
            if (in[k] && I[k][i] > 0)
              wi[k] = wg[k];
            else
              wi[k] = 0.0;
          }
        }

        // Calculate Neff[i]
        Neff[i] = 0.0;
        for (j = 1; j <= L; ++j) {
          if (n[j][ENDGAP] > MAXENDGAPFRAC * nseqi)
            continue;
          for (a = 0; a < 20; ++a)
            fj[a] = 0;
          for (k = 0; k < N_in; ++k)
            if (in[k] && I[k][i] > 0 && X[k][j] < ANY)
              fj[(int) X[k][j]] += wi[k];
          NormalizeTo1(fj, NAA);
          for (a = 0; a < 20; ++a)
            if (fj[a] > 1E-10)
              Neff[i] -= fj[a] * fast_log2(fj[a]);
        }
        if (ncol > 0)
          Neff[i] = fpow2(Neff[i] / ncol);
        else
          Neff[i] = 1.0;

      }
      // Calculate transition probabilities from I state
      q->tr[i][I2M] = q->tr[i][I2I] = 0.0;
      for (k = 0; k < N_in; ++k)  //for all sequences
          {
        if (in[k] && I[k][i] > 0)  //current state is I
            {
          q->tr[i][I2M] += wi[k];
          q->tr[i][I2I] += wi[k] * (I[k][i] - 1);
        }
      }  // end for(k)
    }

    else  // fast global weights?
    {
      float w_I = -1.0 / N_filtered;
      ncol = 0;
      q->tr[i][I2M] = q->tr[i][I2I] = 0.0;
      // Calculate amino acid frequencies fj[a] from weights wg[k]
      for (k = 0; k < N_in; ++k)
        if (in[k] && I[k][i] > 0) {
          ncol++;
          w_I += wg[k];
          q->tr[i][I2M] += wi[k];
          q->tr[i][I2I] += wi[k] * (I[k][i] - 1);
        }
      if (ncol > 0) {
        if (w_I < 0)
          Neff[i] = 1.0;
        else
          Neff[i] = Nlim - (Nlim - 1.0) * fpow2(scale * w_I);
//            fprintf(stderr,"I  i=%3i  ncol=%3i  Neff_M=%5.2f  Nlim=%5.2f  w_I=%5.3f  Neff_I=%5.2f\n",i,ncol,q->Neff_HMM,Nlim,w_I,Neff[i]);
      } else {
        Neff[i] = 0.0;
        q->tr[i][I2M] = -100000;
        q->tr[i][I2I] = -100000;
        continue;
      }
    }

    // Normalize and take log
    sum = q->tr[i][I2M] + q->tr[i][I2I];
    q->tr[i][I2M] = flog2(q->tr[i][I2M] / sum);
    q->tr[i][I2I] = flog2(q->tr[i][I2I] / sum);

  }
  // end loop through alignment columns i
  //////////////////////////////////////////////////////////////////////////////////////////////

  // delete n[][]
  for (j = 1; j <= L; ++j)
    delete[] (n[j]);
  delete[] (n);

  q->tr[0][I2M] = 0;
  q->tr[0][I2I] = -100000;
  q->tr[L][I2M] = 0;
  q->tr[L][I2I] = -100000;
  q->Neff_I[0] = 99.999;

  // Assign Neff_I[i]
  for (i = 1; i <= L; ++i)  // Calculate wi[k] at position i as well as Neff[i] and Neff[i]
    q->Neff_I[i] = Neff[i];
  return;
}

/////////////////////////////////////////////////////////////////////////////////////
// Calculate transitions q->tr[i][a] (a=DM,DD) with pos-specific subalignments
/////////////////////////////////////////////////////////////////////////////////////
void Alignment::Transitions_from_D_state(HMM* q, char* in, const int maxres) {
  // Calculate position-dependent weights wi[k] for each i.
  // For calculation of weights in column i use sub-alignment
  // over sequences which have a DELETE in column i
  // and over columns where none of these sequences has an end gap.
  // This is done by updating the arrays n[j][a] and rj[j] at each step i-1->i while letting i run from 1 to L.
  // n[j][a] = number of occurences of index a at column j of the subalignment,
  //        => only columns with n[j][ENDGAP]=0 are contained in the subalignment!
  // If no sequences enter or leave the subalignment at the step i-1 -> i (i.e. change=0)
  // then the old values wi[k], Neff[i-1], and ncol are used for the new position i.
  // Index a can be an amino acid (0-19), ANY=20, GAP=21, or ENDGAP=22

  int k;                      // index of sequence
  int i, j;                    // position in alignment
  int a;                      // amino acid (0..19)
  int naa;                    // number of different amino acids
  int** n;  // n[j][a] = number of seq's with some residue at column i AND a at position j
  float wi[MAXSEQ];  // weight of sequence k in column i, calculated from subalignment i
  float Neff[maxres];         // diversity of subalignment i
  int nseqi = 0;      // number of sequences in subalignment i (for DEBUGGING)
  int ncol = 0;              // number of columns j that contribute to Neff[i]
  char change;  // has the set of sequences in subalignment changed? 0:no  1:yes
  float fj[NAA + 3];            // to calculate entropy
  float sum;
  float Nlim = 0.0;             // only for global weights
  float scale = 0.0;            // only for global weights

  // Global weights?
  //if (par.wg==1)
  if (1) {
    for (k = 0; k < N_in; ++k)
      wi[k] = wg[k];
    Nlim = fmax(10.0, q->Neff_HMM + 1.0);    // limiting Neff
    scale = flog2((Nlim - q->Neff_HMM) / (Nlim - 1.0));  // for calculating Neff for those seqs with dels at specific pos
  }

  // Initialization
  n = new int*[L + 2];
  for (j = 1; j <= L; ++j)
    n[j] = new int[NAA + 3];
  for (j = 1; j <= L; ++j)
    for (a = 0; a < NAA + 3; ++a)
      n[j][a] = 0;

  //////////////////////////////////////////////////////////////////////////////////////////////
  // Main loop through alignment columns
  for (i = 1; i <= L; ++i)  // Calculate wi[k] at position i as well as Neff[i]
      {
    //if (par.wg==0) // if local weights
    if (0) {
      change = 0;
      // Check all sequences k and update n[j][a] and ri[j] if necessary
      for (k = 0; k < N_in; ++k) {
        if (!in[k])
          continue;
        if (X[k][i - 1] != GAP && X[k][i] == GAP) {  // ... if sequence k was NOT included in i-1 and has to be included for column i
          change = 1;
          nseqi++;
          for (int j = 1; j <= L; ++j)
            n[j][(int) X[k][j]]++;
        } else if (X[k][i - 1] == GAP && X[k][i] != GAP) {  // ... if sequence k WAS included in i-1 and has to be thrown out for column i
          change = 1;
          nseqi--;
          for (int j = 1; j <= L; ++j)
            n[j][(int) X[k][j]]--;
        }
      }  //end for (k)
      nseqs[i] = nseqi;

      // If there is no sequence in subalignment j ...
      if (nseqi == 0) {
        ncol = 0;
        Neff[i] = 0.0;  // effective number of sequences = 0!
        q->tr[i][D2M] = -100000;
        q->tr[i][D2D] = -100000;
        continue;
      }

      // If subalignment changed: update weights wi[k] and Neff[i]
      if (change) {
        // Initialize weights and numbers of residues for subalignment i
        ncol = 0;
        for (k = 0; k < N_in; ++k)
          wi[k] = 1E-8;  // for pathological alignments all wi[k] can get 0;

        // sum wg[k][i] over all columns j and sequences k of subalignment
        for (j = 1; j <= L; ++j) {
          if (n[j][ENDGAP] > MAXENDGAPFRAC * nseqi)
            continue;
          naa = 0;
          for (a = 0; a < 20; ++a)
            if (n[j][a])
              naa++;
          if (naa == 0)
            continue;
          ncol++;
          for (k = 0; k < N_in; ++k) {
            if (in[k] && X[k][i] == GAP && X[k][j] < ANY) {
              if (!n[j][(int) X[k][j]]) {
                HH_LOG(ERROR) << "In " << __FILE__ << ":" << __LINE__ << ": " << __func__ << ":" << std::endl;
                HH_LOG(ERROR) << "\tDi="<<i<<": n["<<j<<"][X["<<k<<"]["<<j<<"]]=0! (X["<<k<<"]["<<j<<"]="<<X[k][j]<<")" << std::endl;
              }
              wi[k] += 1.0 / float(n[j][(int) X[k][j]] * naa);
            }
          }
        }

        // Check whether number of columns in subalignment is sufficient
        if (ncol < NCOLMIN) {
          // Take global weights
          for (k = 0; k < N_in; ++k) {
            if (in[k] && X[k][i] == GAP)
              wi[k] = wg[k];
            else
              wi[k] = 0.0;
          }
        }

        // Calculate Neff[i]
        Neff[i] = 0.0;
        for (j = 1; j <= L; ++j) {
          if (n[j][ENDGAP] > MAXENDGAPFRAC * nseqi)
            continue;
          for (a = 0; a < 20; ++a)
            fj[a] = 0;
          for (k = 0; k < N_in; ++k)
            if (in[k] && X[k][i] == GAP && X[k][j] < ANY)
              fj[(int) X[k][j]] += wi[k];
          NormalizeTo1(fj, NAA);
          for (a = 0; a < 20; ++a)
            if (fj[a] > 1E-10)
              Neff[i] -= fj[a] * fast_log2(fj[a]);
        }
        if (ncol > 0)
          Neff[i] = fpow2(Neff[i] / ncol);
        else
          Neff[i] = 1.0;

      }

      else  //no update was necessary; copy values for i-1
      {
        Neff[i] = Neff[i - 1];
      }

      // Calculate transition probabilities from D state
      q->tr[i][D2M] = q->tr[i][D2D] = 0.0;
      for (k = 0; k < N_in; ++k)  //for all sequences
          {
        if (in[k] && X[k][i] == GAP)            //current state is D
            {
          if (X[k][i + 1] == GAP)      //next state is D
            q->tr[i][D2D] += wi[k];
          else if (X[k][i + 1] <= ANY)  //next state is M
            q->tr[i][D2M] += wi[k];
        }
      }  // end for(k)
    }

    else  // fast global weights?
    {
      float w_D = -1.0 / N_filtered;
      ncol = 0;
      q->tr[i][D2M] = q->tr[i][D2D] = 0.0;
      // Calculate amino acid frequencies fj[a] from weights wg[k]
      for (k = 0; k < N_in; ++k)  //for all sequences
        if (in[k] && X[k][i] == GAP)            //current state is D
            {
          ncol++;
          w_D += wg[k];
          if (X[k][i + 1] == GAP)      //next state is D
            q->tr[i][D2D] += wi[k];
          else if (X[k][i + 1] <= ANY)  //next state is M
            q->tr[i][D2M] += wi[k];
        }
      if (ncol > 0) {
        if (w_D < 0)
          Neff[i] = 1.0;
        else
          Neff[i] = Nlim - (Nlim - 1.0) * fpow2(scale * w_D);
//            fprintf(stderr,"D  i=%3i  ncol=%3i  Neff_M=%5.2f  Nlim=%5.2f  w_D=%5.3f  Neff_D=%5.2f\n",i,ncol,q->Neff_HMM,Nlim,w_D,Neff[i]);
      } else {
        Neff[i] = 0.0;  // effective number of sequences = 0!
        q->tr[i][D2M] = -100000;
        q->tr[i][D2D] = -100000;
        continue;
      }
    }

    // Normalize and take log
    sum = q->tr[i][D2M] + q->tr[i][D2D];
    q->tr[i][D2M] = flog2(q->tr[i][D2M] / sum);
    q->tr[i][D2D] = flog2(q->tr[i][D2D] / sum);

  }
  // end loop through alignment columns i
  //////////////////////////////////////////////////////////////////////////////////////////////

  q->tr[0][D2M] = 0;
  q->tr[0][D2D] = -100000;
  q->Neff_D[0] = 99.999;

  // Assign Neff_D[i]
  for (i = 1; i <= L; ++i)
    q->Neff_D[i] = Neff[i];

  // delete n[][]
  for (j = 1; j <= L; ++j)
    delete[] (n[j]);
  delete[] (n);
  return;
}

/////////////////////////////////////////////////////////////////////////////////////
// Write alignment without insert states (lower case) to alignment file?
/////////////////////////////////////////////////////////////////////////////////////
void Alignment::WriteWithoutInsertsToFile(const char* alnfile,
                                          const char append) {
  HH_LOG(INFO) << "Writing alignment to " << alnfile << "\n";
  FILE* alnf;
  if (!append)
    alnf = fopen(alnfile, "w");
  else
    alnf = fopen(alnfile, "a");
  if (!alnf)
    OpenFileError(alnfile, __FILE__, __LINE__, __func__);
  // If alignment name is different from that of query: write name into commentary line
  if (strncmp(longname, sname[kfirst], DESCLEN - 1))
    fprintf(alnf, "#%s\n", longname);

  HH_LOG(INFO) << "Writing alignment to " << alnfile << "\n";

  // Write ss_ lines
  for (int k = 0; k < N_in; ++k)
    if (k == kss_pred || k == kss_conf || k == kss_dssp || k == ksa_dssp) {
      fprintf(alnf, ">%s\n", sname[k]);
      for (int i = 1; i <= L; ++i)
        fprintf(alnf, "%c", i2aa(X[k][i]));
      fprintf(alnf, "\n");
    }
  // Write other sequences
  for (int k = 0; k < N_in; ++k)
    if (!(k == kss_pred || k == kss_conf || k == kss_dssp || k == ksa_dssp)
        && (keep[k] || display[k] == 2))  // print if either in profile (keep[k]>0) or display is obligatory (display[k]==2)
        {
      fprintf(alnf, ">%s\n", sname[k]);
      for (int i = 1; i <= L; ++i)
        fprintf(alnf, "%c", i2aa(X[k][i]));
      fprintf(alnf, "\n");
    }
  fclose(alnf);
}

/////////////////////////////////////////////////////////////////////////////////////
// Write stored, filtered sequences WITH insert states (lower case) to alignment file?
/////////////////////////////////////////////////////////////////////////////////////
void Alignment::WriteToFile(const char* alnfile, const char append,
                            const char format[]) {
  std::stringstream out;
  WriteToFile(out, format);

  if (strcmp(alnfile, "stdout") == 0) {
    std::cout << out.str();
  } else {
    std::fstream outf;
    if (append)
      outf.open(alnfile, std::ios::out | std::ios::app);
    else
      outf.open(alnfile, std::ios::out);

    if (!outf.good())
      OpenFileError(alnfile, __FILE__, __LINE__, __func__);

    outf << out.str();

    outf.close();
  }
}

void Alignment::WriteToFile(std::stringstream& out, const char format[]) {

  if (!format || !strcmp(format, "a3m")) {
    // If alignment name is different from that of query: write name into commentary line
    if (strncmp(longname, sname[kfirst], DESCLEN - 1) || readCommentLine == '1')
      out << "#" << longname << std::endl;
    // Write ss_ lines
    for (int k = 0; k < N_in; ++k)
      if (k == kss_pred || k == kss_conf || k == kss_dssp || k == ksa_dssp)
        out << ">" << sname[k] << std::endl << seq[k] + 1 << std::endl;
    // Write other sequences
    for (int k = 0; k < N_in; ++k)
      if (!(k == kss_pred || k == kss_conf || k == kss_dssp || k == ksa_dssp)
          && (keep[k] || display[k] == 2))  // print if either in profile (keep[k]>0) or display is obligatory (display[k]==2)
        out << ">" << sname[k] << std::endl << seq[k] + 1 << std::endl;
  } else  // PSI-BLAST format
  {
    char line[LINELEN];
    char tmp_name[NAMELEN];
    for (int k = 0; k < N_in; ++k)  // skip sequences before kfirst!!
      if (!(k == kss_pred || k == kss_conf || k == kss_dssp || k == ksa_dssp)
          && (keep[k] || display[k] == 2))  // print if either in profile (keep[k]>0) or display is obligatory (display[k]==2)
          {
        strwrd(tmp_name, sname[k], NAMELEN);

        sprintf(line, "%-20.20s ", tmp_name);
        out << line;
        char* ptr = seq[k];
        for (; *ptr != '\0'; ptr++)
          if (*ptr == 45 || (*ptr >= 65 && *ptr <= 90))
            out << *ptr;
        out << std::endl;
      }
  }
}

/////////////////////////////////////////////////////////////////////////////////////
// Read a3m slave alignment of hit from file and merge into (query) master alignment
/////////////////////////////////////////////////////////////////////////////////////
void Alignment::MergeMasterSlave(Hit& hit, Alignment& Tali, char* ta3mfile,
                                 const int par_maxcol) {
  char* cur_seq = new char[par_maxcol];   // Sequence currently read in
  int maxcol = par_maxcol;
  int l, ll;        // position in unaligned template (T) sequence Tali.seq[l]
  int i;              // counts match states in query (Q) HMM
  int j;              // counts match states in T sequence Tali.seq[l]
  int h;              // position in aligned T sequence cur_seq[h]
  int k;              // sequence index
  char c;             //

  HH_LOG(DEBUG) << "Merging " << ta3mfile << " to query alignment\n";

  // Record imatch[j]
  int* imatch = new int[hit.j2 + 1];
  int step = hit.nsteps;
  for (j = hit.j1; j <= hit.j2; ++j) {
    // Advance to position of next T match state j
    while (hit.j[step] < j)
      step--;
    imatch[j] = hit.i[step];
    HH_LOG(DEBUG1) << "step=" << step << "  i=" << imatch[j] << " j="
                             << j << "    i1: " << hit.i1 << "  i2: " << hit.i2
                             << "  j1: " << hit.j1 << "  j2: " << hit.j2
                             << std::endl;
  }

  // Determine number of match states of Qali
  for (L = 0, l = 1; seq[kfirst][l] > '\0'; l++)
    if ((seq[kfirst][l] >= 'A' && seq[kfirst][l] <= 'Z')
        || seq[kfirst][l] == '-')
      L++;

  // For each sequence in T alignment: align to Qali
  for (k = 0; k < Tali.N_in; ++k) {
    if (!Tali.keep[k])
      continue;
    if (N_in >= MAXSEQ) {
      HH_LOG(WARNING) << "Warning in " << __FILE__ << ":" << __LINE__
                                << ": " << __func__ << ":" << std::endl;
      HH_LOG(WARNING)
          << "\tmaximum number of " << MAXSEQ
          << " sequences exceeded while reading " << ta3mfile
          << ". Skipping all following sequences of this MSA" << std::endl;
      break;
    }
    cur_seq[0] = ' ';     // 0'th position not used

    // Add the hit.i1-1 left end gaps to aligned sequence
    for (h = 1; h < hit.i1; h++)
      cur_seq[h] = '-';

    // Advance to match state hit.j1 of Tali.seq[k]
    for (j = 0, l = 1; (c = Tali.seq[k][l]) > '\0'; l++)
      if ((c >= 'A' && c <= 'Z') || c == '-')  // match state at position l?
        if ((++j) == hit.j1)
          break;       // yes: increment j. Reached hit,j1? yes: break

    if (j < hit.j1) {
      HH_LOG(ERROR) << "Error in " << __FILE__ << ":" << __LINE__
                              << ": " << __func__ << ":" << std::endl;
      HH_LOG(ERROR) << "\tdid not find " << hit.j1
                              << " match states in sequence " << k << " of "
                              << Tali.name << ". Sequence:\n" << Tali.seq[k]
                              << std::endl;
      exit(1);
    }

    // Write first match state to cur_seq
    int iprev = hit.i1;  // index of previous query match state
    int lprev = l;      // previous T match state in Tali.seq[k][l]
    cur_seq[h++] = Tali.seq[k][l];  // first column of alignment is Match-Match state

    // For each further match state j in alignment
    step = hit.nsteps;
    for (j = hit.j1 + 1; j <= hit.j2; ++j) {
      // Advance to position of next T match state j
      i = imatch[j];

      // Advance to position of next T match state j
      while ((c = Tali.seq[k][++l]) > '\0'
          && ((c >= 'a' && c <= 'z') || c == '.'))
        ;

      int di = i - iprev;  // number of Match states in Q between T match state j-1 and j
      int dl = l - lprev;  // 1 + number of inserted residues in T sequence between T match state j-1 and j
      if (di == 1) {
        // One Q match state for one T match state (treated as special case for speed reasons)
        // i:       i-1   i         di=1
        // Q:  XXXXXX.....XXXXXX
        // T:  YYYYYYyyyyyYYYYYY
        // j:       j-1   j
        // l:       lprev l         dl=6

        // Inserts in lower case
        for (ll = lprev + 1; ll < l; ll++)
          if (Tali.seq[k][ll] != '-' && Tali.seq[k][ll] != '.')
            cur_seq[h++] = lwrchr(Tali.seq[k][ll]);

        // Template Match state -> upper case
        cur_seq[h++] = Tali.seq[k][ll];
      } else if (di == 0) {
        // Gap in query: no Q match state for on T match state (special case for speed reasons)
        // i:       i-1   i-1       di=0
        // Q:  XXXXXX.....~~~XXX
        // T:  YYYYYYyyyyyYYYYYY
        // j:       j-1   j
        // l:       lprev l         dl=6

        // All T residues (including T match state) in lower case
        for (ll = lprev + 1; ll <= l; ll++)
          if (Tali.seq[k][ll] != '-' && Tali.seq[k][ll] != '.')
            cur_seq[h++] = lwrchr(Tali.seq[k][ll]);
      } else if (di >= dl) {
        // More Match states in Q than Inserts in the T sequence
        // => half T inserts y left, half right-aligned in uc, gaps to fill up
        // Number of T insert residues to be left-aligned: (int)(dl/2)
        // i:        iprev  i       di=7
        // Q:  XXXXXXXXXXXXXXXXXX
        // T:  YYYYYYYyyy-yyYYYYY
        // j:        j-1    j
        // l:        lprev  l       dl=6

        // Add left-bounded template residues
        for (ll = lprev + 1; ll <= lprev + (int) (dl / 2); ll++)
          cur_seq[h++] = uprchr(Tali.seq[k][ll]);

        // Add central gaps
        for (int gap = 1; gap <= di - dl; gap++)
          cur_seq[h++] = '-';

        // Add right-bounded residues
        for (; ll <= l; ll++)
          cur_seq[h++] = uprchr(Tali.seq[k][ll]);
      } else if (di < dl) {
        // Fewer Match states in Q than inserts in T sequence
        // => half of available space di for left- half for right-aligned T inserts, rest in lc
        // number of T inserts to be left-aligned in uc: (int)(di/2),
        // i:        iprev i       di=5
        // Q:  XXXXXXXXX.XXXXXXX
        // T:  YYYYYYYyyyyyYYYYY
        // j:        j-1   j
        // l:        lprev l       dl=6

        // Add left-bounded template residues
        for (ll = lprev + 1; ll <= lprev + (int) (di / 2); ll++)
          cur_seq[h++] = uprchr(Tali.seq[k][ll]);

        // Add central inserts
        for (int ins = 1; ins <= dl - di; ins++, ll++)
          if (Tali.seq[k][ll] != '-' && Tali.seq[k][ll] != '.')
            cur_seq[h++] = lwrchr(Tali.seq[k][ll]);

        // Add right-bounded residues
        for (; ll <= l; ll++)
          cur_seq[h++] = uprchr(Tali.seq[k][ll]);
      }
      HH_LOG(DEBUG3) << "i=" << i << " j=" << j << " l=" << l << " cur_seq="
                     << cur_seq << "\n";

      iprev = i;
      lprev = l;
      if (h >= maxcol - 1000)  // too few columns? Reserve double space
      {
        char* new_seq = new char[2 * maxcol];
        strncpy(new_seq, cur_seq, h);  //////// check: maxcol-1 ????
        delete[] (cur_seq);
        cur_seq = new_seq;
        maxcol *= 2;
      }
    }

    // Add the remaining gaps '-' to the end of the template sequence
    for (i = hit.i2 + 1; i <= L; ++i) {
      cur_seq[h++] = '-';

      // too few columns? Reserve double space
      if (h >= maxcol - 1000) {
        char* new_seq = new char[2 * maxcol];
        strncpy(new_seq, cur_seq, h);
        delete[] (cur_seq);
        cur_seq = new_seq;
        maxcol *= 2;
      }
    }
    cur_seq[h++] = '\0';

    keep[N_in] = display[N_in] = 1;
    seq[N_in] = new char[h];
    if (!seq[N_in])
      MemoryError("array for input sequences", __FILE__, __LINE__, __func__);
    strcpy(seq[N_in], cur_seq);
    X[N_in] = initX(h);
    if (!X[N_in])
      MemoryError("array for input sequences", __FILE__, __LINE__, __func__);
    I[N_in] = new short unsigned int[h];
    if (!I[N_in])
      MemoryError("array for input sequences", __FILE__, __LINE__, __func__);
    sname[N_in] = new char[strlen(Tali.sname[k]) + 1];
    if (!sname[N_in])
      MemoryError("array for input sequences", __FILE__, __LINE__, __func__);
    strcpy(sname[N_in], Tali.sname[k]);
    N_in++;

    HH_LOG(DEBUG1) << "k=" << k << " " << Tali.seq[k] << std::endl;
    HH_LOG(DEBUG1) << "Query " << seq[kfirst] << std::endl;
    HH_LOG(DEBUG1) << "k=" << k << " " << cur_seq << std::endl
                             << std::endl;

  }  // end for (k)

//   printf("N_in=%-5i  HMM=%s  with %i sequences\n",N_in,ta3mfile,Tali.N_filtered);

  delete[] cur_seq;
  delete[] imatch;
  delete[] ksort;
  ksort = NULL;  // if ksort already existed it will be to short for merged alignment
  delete[] first;
  first = NULL;  // if first already existed it will be to short for merged alignment
  delete[] last;
  last = NULL;  // if last  already existed it will be to short for merged alignment
  delete[] nres;
  nres = NULL;  // if nres  already existed it will be to short for merged alignment
}

/////////////////////////////////////////////////////////////////////////////////////
// Add a sequence to Qali
/////////////////////////////////////////////////////////////////////////////////////
void Alignment::AddSequence(char Xk[], int Ik[]) {
  int i;    // position in query and template
  if (L <= 0)
    InternalError("L is not set in AddSequence()", __FILE__, __LINE__,
                  __func__);
  X[N_in] = initX(L + 2);
  for (i = 0; i <= L + 1; ++i)
    X[N_in][i] = Xk[i];
  if (Ik == NULL)
    for (i = 0; i <= L + 1; ++i)
      I[N_in][i] = 0;
  else
    for (i = 0; i <= L + 1; ++i)
      I[N_in][i] = Ik[i];
  N_in++;

  delete[] ksort;
  ksort = NULL;  // if ksort already existed it will be to short for merged alignment
  delete[] first;
  first = NULL;  // if first already existed it will be to short for merged alignment
  delete[] last;
  last = NULL;  // if last  already existed it will be to short for merged alignment
}

/////////////////////////////////////////////////////////////////////////////////////
// Add secondary structure prediction to alignment (overwrite existing prediction)
/////////////////////////////////////////////////////////////////////////////////////
void Alignment::AddSSPrediction(char seq_pred[], char seq_conf[]) {
  unsigned int i;

  if ((int) strlen(seq_pred) != L + 1) {
    HH_LOG(WARNING) << "Could not add secondary struture prediction - unequal length!" << std::endl;
    return;
  }

  if (kss_pred < 0 || kss_conf < 0)   // At least one sequence is added
      {
    delete[] ksort;
    ksort = NULL;  // if ksort already existed it will be to short for merged alignment
    delete[] first;
    first = NULL;  // if first already existed it will be to short for merged alignment
    delete[] last;
    last = NULL;  // if last  already existed it will be to short for merged alignment
  }

  if (kss_pred < 0)  // No ss prediction exists
      {
    kss_pred = N_in;
    keep[N_in] = 0;
    display[N_in] = 1;
    seq[N_in] = new char[L + 2];
    strcpy(seq[N_in], seq_pred);
    X[N_in] = initX(L + 2);
    for (i = 0; i < strlen(seq_pred); ++i)
      X[N_in][i] = ss2i(seq_pred[i]);
    I[N_in] = new short unsigned int[L + 2];
    for (i = 0; i <= strlen(seq_pred); ++i)
      I[N_in][i] = 0;
    sname[N_in] = new char[50];
    strcpy(sname[N_in], "ss_pred PSIPRED predicted secondary structure");
    N_in++;
    N_ss++;
    n_display++;
  } else  // overwrite existing ss prediction
  {
    strcpy(seq[kss_pred], seq_pred);
    for (i = 0; i < strlen(seq_pred); ++i)
      X[kss_pred][i] = ss2i(seq_pred[i]);
  }

  if (kss_conf < 0)  // No ss prediction confidence exists
      {
    kss_conf = N_in;
    keep[N_in] = 0;
    display[N_in] = 1;
    seq[N_in] = new char[L + 2];
    strcpy(seq[N_in], seq_conf);
    X[N_in] = initX(L + 2);
    for (i = 0; i < strlen(seq_pred); ++i)
      X[N_in][i] = cf2i(seq_conf[i]);
    I[N_in] = new short unsigned int[L + 2];
    for (i = 0; i <= strlen(seq_pred); ++i)
      I[N_in][i] = 0;
    sname[N_in] = new char[35];
    strcpy(sname[N_in], "ss_conf PSIPRED confidence values");
    N_in++;
    N_ss++;
    n_display++;
  } else  // overwrite existing ss prediction confidence
  {
    strcpy(seq[kss_conf], seq_conf);
    for (i = 0; i < strlen(seq_pred); ++i)
      X[kss_conf][i] = cf2i(seq_conf[i]);
  }
}

/////////////////////////////////////////////////////////////////////////////////////
// Determine matrix of position-specific weights w[k][i] for multiple alignment
// Pos-specific weights are calculated like in "Amino_acid_frequencies_and_transitions_from_M_state()"
/////////////////////////////////////////////////////////////////////////////////////
void Alignment::GetPositionSpecificWeights(float* w[],
                                           char use_global_weights) {
  // Calculate position-dependent weights wi[k] for each i.
  // For calculation of weights in column i use sub-alignment
  // over sequences which have a *residue* in column i (no gap, no end gap)
  // and over columns where none of these sequences has an end gap.
  // This is done by updating the arrays n[j][a] at each step i-1->i while letting i run from 1 to L.
  // n[j][a] = number of occurences of amino acid a at column j of the subalignment,
  //        => only columns with n[j][ENDGAP]=0 are contained in the subalignment!
  // If no sequences enter or leave the subalignment at the step i-1 -> i (i.e. change=0)
  // then the old values w[k][i] and ncol are used for the new position i.
  // Index a can be an amino acid (0-19), ANY=20, GAP=21, or ENDGAP=22

  char* in = keep;  // to keep the code similar to Amino_acid_frequencies_and_transitions_from_M_state()
  int k;                      // index of sequence
  int i, j;                    // position in alignment
  int a;                      // amino acid (0..19)
  int naa;                    // number of different amino acids
  int** n;  // n[j][a] = number of seq's with some residue at column i AND a at position j
  int nseqi = 0;                // number of sequences in subalignment i
  int ncol = 0;              // number of columns j that contribute to Neff[i]
  char change;  // has the set of sequences in subalignment changed? 0:no  1:yes

  // Global weights?
  if (use_global_weights == 1) {
    for (k = 0; k < N_in; ++k)
      for (i = 1; i <= L; ++i)
        w[k][i] = wg[k];
  } else {

    // Initialization
    n = new int*[L + 2];
    for (j = 1; j <= L; ++j)
      n[j] = new int[NAA + 3];
    for (j = 1; j <= L; ++j)
      for (a = 0; a < NAA + 3; ++a)
        n[j][a] = 0;

    //////////////////////////////////////////////////////////////////////////////////////////////
    // Main loop through alignment columns
    for (i = 1; i <= L; ++i)  // Calculate w[k][i]
        {
      change = 0;
      // Check all sequences k and update n[j][a] and ri[j] if necessary
      for (k = 0; k < N_in; ++k) {
        if (!in[k])
          continue;
        if (X[k][i - 1] >= ANY && X[k][i] < ANY) {  // ... if sequence k was NOT included in i-1 and has to be included for column i
          change = 1;
          nseqi++;
          for (int j = 1; j <= L; ++j)
            n[j][(int) X[k][j]]++;
        } else if (X[k][i - 1] < ANY && X[k][i] >= ANY) {  // ... if sequence k WAS included in i-1 and has to be thrown out for column i
          change = 1;
          nseqi--;
          for (int j = 1; j <= L; ++j)
            n[j][(int) X[k][j]]--;
        }
      }  //end for (k)
      nseqs[i] = nseqi;

      // If subalignment changed: update weights w[k][i] and Neff[i]
      if (change) {
        // Initialize weights and numbers of residues for subalignment i
        ncol = 0;
        for (k = 0; k < N_in; ++k)
          w[k][i] = 0.0;

        // sum wi[k] over all columns j and sequences k of subalignment
        for (j = 1; j <= L; ++j) {
          //  do at least a fraction MAXENDGAPFRAC of sequences in subalignment contain an end gap in j?
          if (n[j][ENDGAP] > MAXENDGAPFRAC * nseqi)
            continue;
          naa = 0;
          for (a = 0; a < 20; ++a)
            if (n[j][a])
              naa++;
          if (naa == 0)
            continue;
          ncol++;
          for (k = 0; k < N_in; ++k) {
            if (in[k] && X[k][i] < ANY && X[k][j] < ANY) {
//                        if (!n[j][ (int)X[k][j]]) {fprintf(stderr,"Error in "<<par.argv[0]<<": Mi=%i: n[%i][X[%i]]=0! (X[%i]=%i)\n",i,j,k,k,X[k][j]);}
              w[k][i] += 1.0 / float(n[j][(int) X[k][j]] * naa);
            }
          }
        }

        // Check whether number of columns in subalignment is sufficient
        if (ncol < NCOLMIN)
          // Take global weights
          for (k = 0; k < N_in; ++k)
            if (in[k]) {
              if (X[k][i] < ANY)
                w[k][i] = wg[k];
              else
                w[k][i] = 0.0;
            }
      }
    }
    // end loop through alignment columns i
    //////////////////////////////////////////////////////////////////////////////////////////////

    // delete n[][]
    for (j = 1; j <= L; ++j)
      delete[] (n[j]);
    delete[] (n);

  }
  return;
}

// Set keep[] and display[] arrays to 0 to mark seqs as non-printable
void Alignment::MarkSeqsAsNonPrintable() {
  for (int k = 0; k < N_in; ++k)
    keep[k] = display[k] = 0;
}
