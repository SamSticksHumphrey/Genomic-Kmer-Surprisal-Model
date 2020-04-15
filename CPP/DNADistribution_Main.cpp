// --------------------------------------------------------------------------------------------------------------------------------------------------------------------
// Distribution_DNA_main.cpp 
// Sam Humphrey, April 2020
// Produce the distrbution of all k-mers in the human reference genome
// 
// Open, Read and CloseFASTA codes weretaken from:::   
// http://www.cse.msu.edu/~yannisun/cse891/hmm-EM/fasta.c
//
// Requires number of processors to be equal to twice the number of chromosomes.
// --------------------------------------------------------------------------------------------------------------------------------------------------------------------

#include <mpi.h>
#include "sys/stat.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string.h>
#include <fstream>
#include <sstream>
#include <ctype.h>
#include <vector>
#include <algorithm>
#include <map>
#include <cmath>
#include <utility>
#include <iterator>
#include <limits>
#include <sys/time.h>
#include <sys/resource.h>
#include "mapreduce.h"
#include "keyvalue.h"

#define FASTA_MAXLINE 128   /* Requires FASTA file lines to be <512 characters */

class Userdata{
public:
	int me;
	int nprocs;
	std::string jobname;
	std::ifstream input_file;
	std::string input_filename;
	std::ofstream output_file;
	std::string output_filename;
	int k;
	std::string reference;
	std::string chr;
	int L; 
};

typedef struct fastafile_s {
  FILE *fp;
  char  buffer[FASTA_MAXLINE];
} FASTAFILE;


extern FASTAFILE *OpenFASTA(char*);
extern int        ReadFASTA(FASTAFILE*, char**, char**, int*);
extern void       CloseFASTA(FASTAFILE*);

void fileread(const std::string, std::string*, int*, std::string*, int, int);
std::string remove_whitespace (const std::string);
std::string revcomp (const std::string);
int Nchecker(std::string);
void split(int**, int**, const int, const int);

void MAP(int, MAPREDUCE_NS::KeyValue*, void*);
void REDUCE(char*, int, char*, int, int*, MAPREDUCE_NS::KeyValue*, void*);

// --------------------------------------------------------------------------------------------------------------------------------------------------------------------
// 		Main function  
// --------------------------------------------------------------------------------------------------------------------------------------------------------------------
int main(int narg, char **args){

	double tstart = MPI_Wtime();

	// Argument handling and MPI setup
	Userdata maindata;
	maindata.jobname = args[1];
	maindata.input_filename = std::string(args[2]);
	maindata.k = atoi(args[3]);
	
	MPI_Init(&narg,&args);
	MPI_Comm_rank(MPI_COMM_WORLD,&maindata.me);
	MPI_Comm_size(MPI_COMM_WORLD,&maindata.nprocs); 

    // // Output filename allocation
    std::string label = std::to_string(maindata.me);

	maindata.output_filename  =  std::string("Results/") + maindata.jobname + "_proc" + label + ".txt";
	if(maindata.me == 0){ std::cout << "Debug 1: check it runs. k = " << maindata.k << ", no of processors = " << maindata.nprocs << "\tme = " << label << std::endl; }

	if(maindata.me == 0) {std::cout << "Debug 2:" <<  maindata.input_filename << std::endl;}
	fileread(maindata.input_filename, &maindata.chr, &maindata.L, &maindata.reference, maindata.me, maindata.nprocs); // fileread is located in fasta_fileread.h


	std::cout << "Proc: " << maindata.me << "\tChr " << maindata.chr << "\tL = " << maindata.L << std::endl;  


	MAPREDUCE_NS::MapReduce *mr = new MAPREDUCE_NS::MapReduce(MPI_COMM_WORLD);	//create a mapreduce object
	mr->memsize = 1240;
	mr->verbosity = 1;
	mr->timer = 1;

	mr->map(maindata.nprocs, MAP, &maindata); // call MAP to add all k-mers into the map reduce object
	MPI_Barrier(MPI_COMM_WORLD);

	mr->collate(NULL);
	MPI_Barrier(MPI_COMM_WORLD);

	if(maindata.me == 0) {std::cout << "Debug 5:" <<  maindata.output_filename << std::endl;}
	maindata.output_file.open(maindata.output_filename.c_str());
	if(!(maindata.output_file)){std::cout << "Couldn't open  output file" << std::endl;  MPI_Finalize(); exit(0);}

	mr->reduce(REDUCE,&maindata);  // call REDUCE to output key \t nvalues 
	MPI_Barrier(MPI_COMM_WORLD);

	maindata.output_file.clear(); maindata.output_file.close();
	MPI_Barrier(MPI_COMM_WORLD);

	clock_t tend = MPI_Wtime();
	double time = (double) (tend-tstart);
	if(maindata.me == 0){std::cout << "Time for process = " << time << std::endl;}

	 MPI_Finalize(); 
	return(0);
}   // end main


// --------------------------------------------------------------------------------------------------------------------------------------------------------------------
// 		MapReduce Functions
// --------------------------------------------------------------------------------------------------------------------------------------------------------------------

void MAP(int itask, MAPREDUCE_NS::KeyValue* kv, void* ptr){		// MAP: requires 25 processors, one for each chromosome + reverse compliment 

	Userdata* _data = (Userdata*) ptr; 
	
	int index; 
	std::string chr;

	if(_data->me < _data->nprocs/2){ chr = _data->reference; }
	else{ chr = revcomp(_data->reference); }

	// for the length of the chromosome, take all k-mers and add to the Mapreduce object
	for(int i = 0; i < chr.length() - _data->k + 1; i++){	

		std::string tmp = chr.substr(i, _data->k);
		if(Nchecker(tmp) == 1){continue;}	// check for "N"s or "M"s etc... 

		kv->add( (char*) tmp.c_str(), sizeof(char)*tmp.length(), (char*) "NA", sizeof(char)*2 );
	}
}

// ----------------------------------------------- Output   ------------------------------------
void REDUCE(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, MAPREDUCE_NS::KeyValue *kv, void *ptr){ // REDUCE send "key \t nvalues" to the output file

    Userdata* _data = (Userdata*) ptr;
    std::string KEY(key);

    uint64_t nval = nvalues;
    _data->output_file << KEY.substr(0, (_data->k)) << "\t" << nval << std::endl;

}

// revcomp: Changes A <-> T  & C <-> G for the reverse  
std::string revcomp (const std::string kmer) {
        std::string revstring;
        for (int i = kmer.size() -1; i >= 0; i--) {
                char c = kmer[i];
                char revchar;
                switch (c) {
                        case 'g':
                        revchar = 'c'; break;

                        case 'G':
                        revchar = 'C'; break;

                        case 'a':
                        revchar = 't'; break;

                        case 'A':
                        revchar = 'T'; break;

                        case 't':
                        revchar = 'a'; break;

                        case 'T':
                        revchar = 'A'; break;

                        case 'c':
                        revchar = 'g'; break;

                        case 'C':
                        revchar = 'G'; break;

                        default:
                        revchar = 'N';
                }
                revstring += revchar;
        }  
        return (revstring.c_str());
}

// ---------------------------------------------------- Nchecker: Skips the N's and other incorrect Kmers  ----------------------------------
int Nchecker(std::string _seq){
        int flag = 0;
        std::string seq = _seq;
        for(std::string::iterator it = seq.begin(); it!= seq.end(); it++) {
                if( *it == 'N' || *it == 'n'){flag = 1;return flag;}
                else if( *it == 'M' || *it == 'm'){flag = 1;return flag;}
                else if( *it == 'R' || *it == 'r'){flag = 1;return flag;}  
                else if( *it == 'Z' || *it == 'z'){flag = 1;return flag;}
        }
  return flag;
}


/* Split the size of a file across processors --------------------------------------*/

void split(int** _start, int** _end, const int size, const int _nprocs){ 
        int nk = size / _nprocs;
        int remain = size % _nprocs;

        for(int i(0); i < _nprocs; i++) {
            int interval = nk;
            if(i<remain){interval++;}
            if(i==0){(*_start)[i] = 0;}
            else{(*_start)[i] = (*_end)[i-1] + 1;}
                
            (*_end)[i] = (*_start)[i] + interval;
        }
} // end split



// Remove whitespace ---------------------------------------------------------------------------------------
std::string remove_whitespace (std::string s) {

        std::string r = "";
        for (unsigned int i = 0; i < s.length(); i++) {
                char c = s[i];
                if (c != '\t' && c != '\n' && c != ' ' && c != '\0') {
                        r += c;
                }
        }
return(r);
}


// Open fasta file ---------------------------------------------------------------------------------------
FASTAFILE* OpenFASTA(char *seqfile){
        FASTAFILE *ffp;

        ffp = (FASTAFILE *) malloc(sizeof(FASTAFILE));
        ffp->fp = fopen(seqfile, "r");              /* Assume seqfile exists & readable!   */

        if (ffp->fp == NULL) { 
                free(ffp);
                std::cout << "Error 1: Couldn't open fasta file : " << seqfile << std::endl;
                exit(3);
        }

        if ((fgets(ffp->buffer, FASTA_MAXLINE, ffp->fp)) == NULL){
                free(ffp);
                std::cout << "Error 2: Couldn't open fasta file : " << seqfile << std::endl;
                exit(3);
        }

        //std::cout << "Fatsa open: ";
        return ffp;
} // end OpenFASTA




// Read fasta function ----------------------------------------------------
int ReadFASTA(FASTAFILE *ffp, char **ret_seq, char **ret_name, int *ret_L){
        char *s;
        char *name;
        char *seq;
        int   n;
        int   nalloc;

        /* Peek at the lookahead buffer; see if it appears to be a valid FASTA descline. */
        if (ffp->buffer[0] != '>') return 0;

        /* Parse out the name: the first non-whitespace token after the >*/
        s  = strtok(ffp->buffer+1, " \t\n");

        name = (char *) malloc(sizeof(char) * (strlen(s)+1));
        strcpy(name, s);

        seq = (char *) malloc(sizeof(char) * 128);     /* allocate seq in blocks of 128 residues */
        nalloc = 128;
        n = 0;

        while (fgets(ffp->buffer, FASTA_MAXLINE, ffp->fp)){
                if (ffp->buffer[0] == '>') break; /* a-ha, we've reached the next descline */
                for (s = ffp->buffer; *s != '\0'; s++){
                        if (! isalpha(*s)) continue;  /* accept any alphabetic character */
                        seq[n] = *s;                  /* store the character, bump length n */
                        n++;
                        if (nalloc == n){  
                        nalloc += 128;
                        seq = (char *) realloc(seq, sizeof(char) * nalloc);
                        }
                }
        }

        seq[n] = '\0';
        *ret_name = name;
        *ret_seq  = seq;
        *ret_L    = n;

        return 1;
}



// Close the fasta function   ---------------------------------------------------------------------------------------
void  CloseFASTA(FASTAFILE *ffp){
        fclose(ffp->fp);
        free(ffp);
}


// Read the reference file  ---------------------------------------------------------------------------------------
void fileread(const std::string fname, std::string* _chr, int* _L, std::string* _ref, int _me, int _nprocs){

        FASTAFILE* ffp;
        char* seq;
        char* chr_name;
        int L;

        int index;

        if(_me < _nprocs/2){index = _me;}
        else{index = _me - _nprocs/2; }

        int counter = 0;

        // for(int i = 0; i < _nprocs/2; i++){

        std::string fullname = fname;

        // open and read in the fasta file
        ffp = OpenFASTA((char*) fullname.c_str());


            while(ReadFASTA(ffp, &seq, &chr_name, &L)){

                // if(_me == 0){std::cout << "readfasta counter while: " << counter << std::endl;}

                if(counter == index){
                std::string name = std::string(chr_name);

                std::string tmp_seq(seq);
                tmp_seq = remove_whitespace(tmp_seq);

                *_ref = tmp_seq;
                *_L = L;
                *_chr = name;
                }
            counter++;
        
        }
        CloseFASTA(ffp);
        // }
} // end fileread function

