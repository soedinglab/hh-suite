#include "hhutil.h"

/////////////////////////////////////////////////////////////////////////////////////
// Errors
/////////////////////////////////////////////////////////////////////////////////////
int FormatError(const char infile[], const char details[])
{
  std::cerr<<"Error in "<<par.argv[0]<<": wrong format while reading file \'"<<infile<<". "<<details<<"\n";
  exit(1);
}

int OpenFileError(const char outfile[])
{
  std::cerr<<std::endl<<"Error in "<<par.argv[0]<<": could not open file \'"<<outfile<<"\'\n";
  exit(2);
}

int MemoryError(const char arrayname[])
{
  std::cerr<<"Error in "<<par.argv[0]<<": Could not allocate memory for \'"<<arrayname<<"\'.\nDo you have >=4GB of RAM per core on your machine? Are your max memory size and stack sizes sufficient? (Check using '$ ulimit -a' under Linux and best set to 'unlimited')"<<std::endl;
  exit(3);
}

int SyntaxError(const char details[])
{
  std::cerr<<"Error in "<<par.argv[0]<<" on command line: "<<details<<"\n";
  exit(4);
}

int InternalError(const char errstr[])
{
  std::cerr<<"Error in "<<par.argv[0]<<":  "<<errstr<<". Please report this bug to the developers\n";
  exit(6);
}


/////////////////////////////////////////////////////////////////////////////////////
//// Replace memalign by posix_memalign (Why? [JS])
/////////////////////////////////////////////////////////////////////////////////////
void *memalign(size_t boundary, size_t size, const char* what_for)
{
  void *pointer;
  if (posix_memalign(&pointer,boundary,size) != 0) {
    std::cerr<<"Error in "<<par.argv[0]<<": memalign could not allocate memory of " << size << " bytes, " << strerror(errno);
    if (what_for!=NULL)
      std::cerr<<"for "<<what_for;

    std::cerr<<".\nDo you have >=4GB of RAM per core on your machine? Are your max memory size and stack sizes sufficient? (Check using '$ ulimit -a' under Linux and best set to 'unlimited')"<<std::endl;
    exit(3);
  }

  return pointer;
}

/////////////////////////////////////////////////////////////////////////////////////
//// Execute system command
/////////////////////////////////////////////////////////////////////////////////////
void runSystem(std::string cmd, int v)
{
  if (v>2)
    std::cout << "Command: " << cmd << "!\n";
  int res = system(cmd.c_str());
  if (res!=0) 
    {
      std::cerr << std::endl << "Error in "<<par.argv[0]<<": Could not execute '" << cmd << "; "<<strerror(errno)<<std::endl;
      exit(1);
    }
    
}


/////////////////////////////////////////////////////////////////////////////////////
// Read up to n lines of outfile and write to screen (STDERR)
/////////////////////////////////////////////////////////////////////////////////////
void WriteToScreen(char* outfile, int n)
{
  char line[LINELEN]="";
  std::ifstream outf;
  outf.open(outfile, std::ios::in);
  if (!outf) {OpenFileError(outfile);}
  std::cout<<"\n";
  for(; n>0 && outf.getline(line,LINELEN); n--) std::cout<<line<<"\n";
  outf.close();
  std::cout<<"\n";
}

void WriteToScreen(char* outfile) {
	WriteToScreen(outfile,INT_MAX);
}

/////////////////////////////////////////////////////////////////////////////////////
// Read .hhdefaults file into array argv_conf (beginning at argv_conf[1])
/////////////////////////////////////////////////////////////////////////////////////
void ReadDefaultsFile(int& argc_conf, char** argv_conf, char* path)
{
  char line[LINELEN]="";
  char filename[NAMELEN];
  char* c_first;   //pointer to first character of argument string
  char* c;         //pointer to scan line read in for end of argument
  //  ifstream configf;
  FILE* configf=NULL;
  argc_conf=1;     //counts number of arguments read in

  // Open config file
  strcpy(filename,"./.hhdefaults");
  configf = fopen(filename,"r");
  if (!configf && path)
  {
    strcpy(filename,path);
    strcat(filename,".hhdefaults");
    configf = fopen(filename,"r");
  }
  if (!configf && getenv("HOME"))
  {
    strcpy(filename,getenv("HOME"));
    strcat(filename,"/.hhdefaults");
    configf = fopen(filename,"r");
  }
  if (!configf) return; // only webserver has no home directory => need no warning

  // Scan file until line 'program_nameANYTHING'
  while (fgets(line,LINELEN,configf))
    if (!strncmp(line,program_name,6)) break;
  // Found line 'program_nameANYTHING'?
  if (!strncmp(line,program_name,6))
    {
      // Read in options until end-of-file or empty line
      //while (fgets(line,LINELEN,configf) && strcmp(line,"\n"))
      while (fgets(line,LINELEN,configf))
        {
          // Analyze line
          c=line;
          do
            {
              // Find next word
              while (*c==' ' || *c=='\t') c++; //Advance until next non-white space
              if ((*c=='h' && *(c+1)=='h') || *c=='\0' || *c=='\n' || *c=='#' || *c==13) break;  //Is next word empty string? (char 13 needed for Windows!)
              c_first=c;
              while (*c!=' ' && *c!='\t'  && *c!='#' && *c!='\0' && *c!='\n' && *c!=13) c++; //Advance until next white space or '#' (char 13 needed for Windows!)
              if (*c=='\0' || *c=='\n' || *c=='#' || *c==13)         //Is end of line reached? (char 13 needed for Windows!)
                {
                  *c='\0';
                  argv_conf[argc_conf]=new char[strlen(c_first)+1];
                  strcpy(argv_conf[argc_conf++],c_first);
                  break;
                }
              *c='\0';
              argv_conf[argc_conf]=new char[strlen(c_first)+1];
              strcpy(argv_conf[argc_conf++],c_first);
              if (v>2) printf("Default argument: %s\n",c_first);
              c++;
            } while (1);
	  if (*c=='h' && *(c+1)=='h') break; // Next program found
        } //end read line
      if (v>=3)
        {
          std::cout<<"Arguments read in from .hhdefaults ("<<filename<<"):";
          for (int argc=1; argc<argc_conf; argc++) std::cout<<(argv_conf[argc][0]=='-'? " ":"")<<argv_conf[argc]<<" ";
          std::cout<<"\n";
        }
      else if (v>=3) std::cout<<"Read in "<<argc_conf<<" default arguments for "<<program_name<<" from "<<filename<<"\n";
    }
  else //found no line 'program_name   anything"
    {
      if (v>=3) std::cerr<<std::endl<<"WARNING: no default options for \'"<<program_name<<"\' found in "<<filename<<"\n";
      return; //no line 'program_name   anything' found
    }
  //   configf.close();
  fclose(configf);
}


/////////////////////////////////////////////////////////////////////////////////////
// Count the number of sequences "^>" in <file>
/////////////////////////////////////////////////////////////////////////////////////
int CountSeqsInFile(char* file, int& numseqs)
{
  char line[LINELEN]="";         // input line
  char tmp_file[NAMELEN];
  int LDB=0;
  numseqs=0;
  strcpy(tmp_file, file);
  strcat(tmp_file, ".sizes");
  FILE* fin = fopen(tmp_file, "r");
  if (fin)
    {
      char* ptr=fgets(line, LINELEN, fin);
      numseqs = strint(ptr);
      LDB = strint(ptr);
      fclose(fin);
     } 
  else 
    {
      fin = fopen(file, "r");
      if (!fin) OpenFileError(file);
      while (fgets(line,LINELEN,fin))
	{ 
	  if (line[0]=='>') numseqs++;
	  else LDB += strlen(line);
	}
      fclose(fin);
    }
  return LDB;
}


/////////////////////////////////////////////////////////////////////////////////////
// Count number of lines in <file>
/////////////////////////////////////////////////////////////////////////////////////
int CountLinesInFile(char* file)
{
  char line[LINELEN]="";         // input line
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
      if (!fin) OpenFileError(file);
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
