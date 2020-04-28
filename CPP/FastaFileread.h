// --------------------------------------------------------------------------------------------------------------------------------------------------------------------
// fasta_fileread.h 
//
// Sam Humphrey  April 2020
//
//
// Functions for opening, reading and returning the reference sequences 
// --------------------------------------------------------------------------------------------------------------------------------------------------------------------

#ifndef fasta_filereadH
#define fasta_filereadH

#define FASTA_MAXLINE 128	/* Requires FASTA file lines to be <512 characters */

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string.h>
#include <fstream>
#include <ctype.h>
#include <vector>
#include <map>

typedef struct fastafile_s {
  FILE *fp;
  char  buffer[FASTA_MAXLINE];
} FASTAFILE;


extern FASTAFILE *OpenFASTA(char*);
extern int        ReadFASTA(FASTAFILE*, char**, char**, int*);
extern void       CloseFASTA(FASTAFILE*);
typedef std::map<std::string, int> Chr_Map;

void fileread(const std::string, std::vector<std::string>*, Chr_Map*, std::vector<int>*);
std::string remove_whitespace(const std::string);
std::string revcomp(const std::string);
int Nchecker(std::string);

#endif



















