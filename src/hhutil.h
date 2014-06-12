/*
 * hhutil.h
 *
 *  Created on: Mar 28, 2014
 *      Author: meiermark
 */

#include <errno.h>

#include <iostream>

#include "hhdecl.h"
#include "hhutil-inl.h"

#ifndef HHUTIL_H_
#define HHUTIL_H_


int FormatError(const char infile[], const char* file, const int line, const char* func, const char details[]="");

int OpenFileError(const char outfile[], const char* file, const int line, const char* func);

int MemoryError(const char arrayname[], const char* file, const int line, const char* func);

int SyntaxError(const char* file, const int line, const char* func, const char details[]="");

int InternalError(const char errstr[], const char* file, const int line, const char* func);

/////////////////////////////////////////////////////////////////////////////////////
//// Execute system command
/////////////////////////////////////////////////////////////////////////////////////
void runSystem(std::string cmd);


/////////////////////////////////////////////////////////////////////////////////////
// Read up to n lines of outfile and write to screen (STDERR)
/////////////////////////////////////////////////////////////////////////////////////
void WriteToScreen(char* outfile, int n);

void WriteToScreen(char* outfile);


/////////////////////////////////////////////////////////////////////////////////////
// Read .hhdefaults file into array argv_conf (beginning at argv_conf[1])
/////////////////////////////////////////////////////////////////////////////////////
void ReadDefaultsFile(int& argc_conf, char** argv_conf, char* path=NULL);


/////////////////////////////////////////////////////////////////////////////////////
// Count the number of sequences "^>" in <file>
/////////////////////////////////////////////////////////////////////////////////////
int CountSeqsInFile(char* file, int& numseqs);


/////////////////////////////////////////////////////////////////////////////////////
// Count number of lines in <file>
/////////////////////////////////////////////////////////////////////////////////////
int CountLinesInFile(char* file);


void float_to_8_bit(float input, unsigned char& result);

void bit_8_to_float(unsigned char input, float& result);

void float_to_16_bit(float input, unsigned short int& result);

void bit_16_to_float(unsigned short int input, float& result);

void writeU16(std::ostream& file, unsigned short int val);

#endif /* HHUTIL_H_ */
