#include "hhutil.h"

/////////////////////////////////////////////////////////////////////////////////////
// Errors
/////////////////////////////////////////////////////////////////////////////////////
int FormatError(const char infile[], const char* file, const int line, const char* func, const char details[]) {
  HH_LOG(ERROR) << "In " << file << ":" << line << ": " << func << ":" << std::endl;
  HH_LOG(ERROR) << "\twrong format while reading file \'"<<infile<<". "<<details<<"\n";
  exit(1);
}

int OpenFileError(const char outfile[], const char* file, const int line, const char* func) {
  HH_LOG(ERROR) << "In " << file << ":" << line << ": " << func << ":" << std::endl;
  HH_LOG(ERROR) << "\tcould not open file \'" << outfile << "\'\n";
  exit(2);
}

int MemoryError(const char arrayname[], const char* file, const int line, const char* func) {
  HH_LOG(ERROR) << "In " << file << ":" << line << ": " << func << ":" << std::endl;
  HH_LOG(ERROR) << "\tCould not allocate memory for \'"<<arrayname<<"\'." << std::endl;
  HH_LOG(ERROR) << "\tDo you have >=4GB of RAM per core on your machine? " <<
      "Are your max memory size and stack sizes sufficient? " <<
      "(Check using '$ ulimit -a' under Linux and best set to 'unlimited')" << std::endl;
  exit(3);
}

int SyntaxError(const char* file, const int line, const char* func, const char details[]) {
  HH_LOG(ERROR) << "In " << file << ":" << line << ": " << func << ":" << std::endl;
  HH_LOG(ERROR) << "\ton command line: "<<details << std::endl;
  exit(4);
}

int InternalError(const char errstr[], const char* file, const int line, const char* func) {
  HH_LOG(ERROR) << "In " << file << ":" << line << ": " << func << ":" << std::endl;
  HH_LOG(ERROR) << "\t" << errstr << ". Please report this bug to the developers" << std::endl;
  exit(6);
}

/////////////////////////////////////////////////////////////////////////////////////
// Count number of lines in <file>
/////////////////////////////////////////////////////////////////////////////////////
int CountLinesInFile(const char* file) {
  std::unique_ptr<char[]> line_ptr(new char[LINELEN]);
  char* line = line_ptr.get();
  int numlines=0;
  char tmp_file[NAMELEN];
  strcpy(tmp_file, file);
  strcat(tmp_file, ".sizes");
  FILE* fin = fopen(tmp_file, "r");
  if (fin)
    {
      char* ptr=fgets(line, LINELEN, fin);
      numlines = strint(ptr);
      fclose(fin);
    } 
  else 
    {
      fin = fopen(file, "r");
      if (!fin) OpenFileError(file, __FILE__, __LINE__, __func__);
      while (fgets(line,LINELEN,fin)) numlines++; 
      fclose(fin);
    }
  return numlines;
}



void float_to_8_bit(float input, unsigned char& result) {
  const unsigned int normalization_111 = 939524096;
  const unsigned int complete_exp_mask = 2139095040;
  const unsigned int exp_mask = 125829120;
  const unsigned int mant_mask = 7864320;

  const unsigned int exp_shift = 19;
  const unsigned int mant_shift = 19;

  unsigned int in = *((unsigned int*)(&input));

  unsigned int e = in & complete_exp_mask;
  e -= normalization_111;
  e = e & exp_mask;
  e = e >> exp_shift;

  unsigned int m = in & mant_mask;
  m = m >> mant_shift;

  result = e | m;
}

void bit_8_to_float(unsigned char input, float& result) {
  const unsigned int normalization_111 = 939524096;
  const unsigned int exp_shift = 19;
  const unsigned int mant_shift = 19;

  unsigned int m = input & 15;
  m = m << mant_shift;

  unsigned int e = input & 240;
  e = e << exp_shift;
  e += normalization_111;

  unsigned int res = e | m;

  result = *((float*)(&res));
}

void float_to_16_bit(float input, unsigned short int& result) {
  const unsigned int normalization_63 = 536870912;
  const unsigned int complete_exp_mask = 2139095040;
  const unsigned int exp_mask = 528482304;
  const unsigned int mant_mask = 8380416;

  const unsigned int exp_shift = 13;
  const unsigned int mant_shift = 13;

  unsigned int in = *((unsigned int*)(&input));

  unsigned int e = in & complete_exp_mask;
  e -= normalization_63;
  e = e & exp_mask;
  e = e >> exp_shift;

  unsigned int m = in & mant_mask;
  m = m >> mant_shift;

  result = e | m;
}

void bit_16_to_float(unsigned short int input, float& result) {
  const unsigned int normalization_63 = 536870912;

  const unsigned int exp_shift = 13;
  const unsigned int mant_shift = 13;

  unsigned int m = input & 2046;
  m = m << mant_shift;

  unsigned int e = input & 129024;
  e = e << exp_shift;
  e += normalization_63;
  unsigned int res = e | m;

  result = *((float*)(&res));
}

void writeU16(std::ostream& file, unsigned short int val) {
  unsigned char bytes[2];

  // extract the individual bytes from our value
  bytes[1] = (val) & 0xFF;  // low byte
  bytes[0] = (val >> 8) & 0xFF;  // high byte

  // write those bytes to the file
  file.write( (char*)bytes, 2 );
}

void writeS16(std::ostream& file, signed short int val) {
  signed char bytes[2];

  // extract the individual bytes from our value
  bytes[1] = (val) & 0xFF;  // low byte
  bytes[0] = (val >> 8) & 0xFF;  // high byte

  // write those bytes to the file
  file.write( (char*)bytes, 2 );
}

