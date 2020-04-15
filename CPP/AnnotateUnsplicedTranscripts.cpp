// --------------------------------------------------------------------------------------------------------------------------------------------------------------------
// AnnotateUnsplicedTranscripts.cpp 
// Sam Humphrey, April 2020
// 
// We require the unspliced transcript sequences, which are too big to download directly using R. 
// So, here we: 
//  Read in transript files produced from extractTranscriptAnnotationFromGTF.R.
//  Take the transcript chromosome, strand, start and end points
//  Extract the full unspliced sequence from the reference
//  The mRNA file is then used to fill the annotation
//
// Note: reverse strand annotation is calculated in strand reverse strand specific coordinates, not forward strand.
// --------------------------------------------------------------------------------------------------------------------------------------------------------------------

#include "sys/stat.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string.h>
#include <fstream>
#include <sstream>
#include <ctype.h>
#include <vector>
#include <time.h>
#include <cmath>
#include <iterator>
#include <map>
#include <sys/time.h>
#include <sys/resource.h>
#include "FastaFileread.h"

// Class Transcript: contains all parameters for the transcript read in from the mRNA input file. 
class Transcript{
    private:    
        std::string chr;
        int start;
        int end;
        int forwardStart;
        int forwardEnd;
        int strand;
        std::string geneName;
        std::string geneID;
        std::string geneBiotype;
        std::string transID;
        std::string transName;
        std::string transBiotype;
        std::string proteinID;      
        std::string transSeq;

    public: 
        Transcript(){}; 
        Transcript(std::string&, Chr_Map* ,std::vector<int>&);
        ~Transcript(){}

        std::string get_transID(){return transID;}
        int get_chr(Chr_Map* _map){return _map->find(chr)->second;}
        std::string get_transSeq(){return transSeq;}
        int get_strand(){return strand;}
        int get_forwardStart(){return forwardStart;}
        int get_forwardEnd(){return forwardEnd;}

        std::ostream & operator<<(std::ostream&);

};

// ----------------------------------------------------------------------------------
// Transcript Class functions
// ----------------------------------------------------------------------------------
// Parameterised constructor
Transcript::Transcript(std::string& _str, Chr_Map* _map, std::vector<int>& _L){

        std::string RNAseq, CDSseq, AAseq;
        std::stringstream ss;
        ss << _str;    
        ss >> chr >> start >> end >> strand >> geneID >> geneName >> geneBiotype >> transID >> transName >> transBiotype >> proteinID >> RNAseq >> CDSseq >> AAseq;

        // Need to sort out the forward/reverse strand buisness, I've manually written all features to be from the respective start of the strand, so I 
        //  need to convert all stuff back to the forwrad strand for this analysis
        if(strand == 1){
            forwardStart = start;
            forwardEnd = end;
        } else if(strand == -1){
            int  chrEnd = _L[_map->find(chr)->second];
            forwardStart = chrEnd - end; 
            forwardEnd = chrEnd - start;
        }

}

// Overload the output operator 
std::ostream & Transcript::operator<<(std::ostream& os){

    os << chr << "\t" << start << "\t" << end << "\t" << strand << "\t" <<  geneName << "\t" << geneID << "\t" << geneBiotype << "\t" << transID << "\t" << transName << "\t" << transBiotype << "\t" << proteinID << "\t" << transSeq;
    return (os);
}

// --------------------------------------------------------------------------------------------------------------------------------------------------------------------
//      Main function  
// --------------------------------------------------------------------------------------------------------------------------------------------------------------------
int main(int narg, char **args){

    clock_t start = clock();

    // Argument handling
    std::string jobname = std::string(args[1]);
    std::string input_transcriptfilename = std::string(args[2]);
    std::string DNAfastq_filename = std::string(args[3]);
    std::cout << "Debug1: " << input_transcriptfilename << std::endl;
    
    // Open the input file & file check
    std::ifstream input_file(input_transcriptfilename.c_str());
    if(!(input_file)){std::cout << "Couldn't open input file" << std::endl; exit(0);}

    // Read in the genome fasta files into vectors 
    std::vector<std::string> reference; // Contains the reference sequence for each chromosome      // Contains the chromosome names (1-22, X, Y , MT)
    Chr_Map chr_map;                    // Map of chromosomes - indecies
    std::vector<int> L;                 // Contains chromosome size, to check that we're not going to go out of scope when looking for sequences0

    fileread(DNAfastq_filename, &reference, &chr_map, &L);
	std::cout << "Debug 2: Fileread Sucessful size = " << chr_map.size() << std::endl;

    // ----------------------------------------------------------------------------------------------------------------------------------------------------------------
    // Read in the exon region file to all threads
    // ----------------------------------------------------------------------------------------------------------------------------------------------------------------
    std::string output_fileName  = std::string(jobname) + ".txt";
    
    std::ofstream output_file(output_fileName.c_str());
    if(!(output_file)){std::cout << "Couldn't open output file" << std::endl; exit(0);}

    // ----------------------------------------------------------------------------------------------------------------------------------------------------------------
   	// While loop over the file to read in every line into a std::vector of Region Vtrans
    // ----------------------------------------------------------------------------------------------------------------------------------------------------------------
    int counter = 0;
    std::string line; 
    while(getline(input_file,line)){
        if(input_file.eof()){break;}                             // Break at the end of the file 
        
        if(input_file.fail() && !(input_file.eof())){           // If a line is incorrect
            input_file.clear();
            std::cout << "The datafile contains an error that was ignored\n" << std::endl;
            continue;
        }

        Transcript trans(line, &chr_map, L);                                 // Transcript parameterised constructor
        line.clear();                   

        std::string fullChr = reference[trans.get_chr(&chr_map)];                // get the full chromosome
        std::string seq = fullChr.substr(trans.get_forwardStart() - 1, (trans.get_forwardEnd() - trans.get_forwardStart() + 1)); // get the sequnence between the start and end points
        if(trans.get_strand() == -1){seq = revcomp(seq);}  // Deal with reverse strands 
        
        trans << output_file;
        output_file << seq << std::endl;

        counter++;
    } // End for loop over Vtrans to get all transcript ID's
    
    output_file.clear(); output_file.close();
    
    clock_t end = clock();

    std::cout << "Time to process " << counter << " transcripts = " << ((float) end - start)/CLOCKS_PER_SEC << "secs"<< std::endl; 

}


