#include <iostream>   // cin, cout, cerr
#include <fstream>    // ofstream, ifstream 
#include <stdio.h>    // printf
#include <stdlib.h>   // exit
#include <string>     // strcmp, strstr
using std::cout;
using std::cerr;
using std::endl;
using std::ios;
using std::ifstream;
using std::ofstream;

// Removes the newline and other control characters at the end of a string (if present)
// and returns the new length of the string (-1 if str is NULL)
inline int chomp(char* str) 
{
  if (!str) return -1;
  int l;
  for (l=strlen(str)-1; l>=0 && str[l]<32; l--);
  str[++l]='\0';
  return l;
}

// Emulates the ifstream::getline method; similar to fgets(str,maxlen,FILE*), 
// but removes the newline at the end and returns NULL if at end of file or read error
inline char* fgetline(char* str, int maxlen, FILE* file) 
{
  if (fgets(str,maxlen,file)) 
    {
      chomp(str);
      return(str);
    }
  else return NULL;
}


int OpenFileError(char outfile[])
{
  cerr<<endl<<"Error: could not open file \'"<<outfile<<"\'\n"; 
  exit(2);
}

void Read(FILE* inf) 
{
  char line[1000];
  while (fgetline(line,1001,inf) && !(line[0]=='/' && line[1]=='/'))
    {
      if (!strcmp(line,"test")) printf("Found it: %s\n",line);
    }       
  fclose(inf);
}

int dummy_stack(int i, const int L) 
{
  char dummychar[L];
  dummychar[i]=i;
  int dummyint=dummychar[i];
  return dummyint;
}

int dummy_heap(int i, const int L) 
{
  char* dummychar=new(char[L]);
  dummychar[i]=i;
  int dummyint=dummychar[i];
  delete[] dummychar;
  return dummyint;
}


/////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{  
  const int L=100000000;
  int i;
  int dummysum=0;
    for (i=0; i<L-1; i++) dummysum+=dummy_stack(i,i+1);
    //  for (i=0; i<L-1; i++) dummysum+=dummy_heap(i,i+1);
  printf("%i\n",(int)dummysum);
  exit(0);
}
