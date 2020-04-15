// --------------------------------------------------------------------------------------------------------------------------------------------------------------------
// CodingSequenceDistribution.h
// Sam Humphrey, April 2020
// --------------------------------------------------------------------------------------------------------------------------------------------------------------------


#ifndef RNADIST_CDS_H
#define RNADIST_CDS_H

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
#include <iterator>
#include <sys/time.h>
#include <sys/resource.h>
#include "mapreduce.h"
#include "keyvalue.h"
#include "Blockmacros.h"

/* ------------------------------------------------------------
* TranscriptRegion class containing the information of the region around the Reg, 
* 	before tand after the region is mutated
* ------------------------------------------------------------*/
class Userdata{
	public:
		// Stuff for MPI
		int me;
		int nprocs;

		// Arguments
		int k;
		std::string jobname;
		int input_transcriptfilesize;
		int depth;

		std::map<std::string, int> chr_map;

		// file IO
		std::ifstream input_transcriptfile;
		std::string input_transcriptfilename;

		std::ifstream input_exonfile;
		std::string input_exonfilename;

		std::ofstream output_dist;
		std::string output_distname;

		std::ofstream output_exonExonJcts;
		std::string output_exonExonJctsname;
};


class Transcript{
	private:    
		std::string chr;
		int transStart;
		int transEnd;
		int strand;
		std::string geneID;
		std::string geneName;
		std::string geneBiotype;
		std::string transID;
		std::string transName;
		std::string transBiotype;
		std::string proteinID;     
		std::string RNAseq;
		std::string CDSseq;
		std::string AAseq;

	public:	
		Transcript(){};	
	    Transcript(std::string&);
		~Transcript(){}

		int get_chr(Userdata* _data){return _data->chr_map.find(chr)->second;}
		std::string get_transID(){return transID;}
		std::string get_CDSseq(){return CDSseq;}
		int get_strand(){return strand;}
		int get_transStart(){return transStart;}
		int get_transEnd(){return transEnd;}
	 	std::string print_value(int, int);
	 	std::ostream & operator<<(std::ostream&);	
};

class Exon {
	private:
  	        std::string transID;
			std::string chr;	        
			int exonStart;
	        int exonEnd;
	        int rank;	        
	        std::string exonID;

	public:	
		Exon(){};	
	    Exon(std::string&);
		~Exon(){}
	
		std::string get_chr(){return chr;}
		std::string get_transID(){return transID;}
		std::string get_exonID(){return exonID;}
		int get_exonWidth(){return (exonEnd - exonStart) + 1;}
		int get_exonStart(){return exonStart;}
		int get_exonEnd(){return exonEnd;}
		int get_rank(){return rank;}
};


class BoundInfo{
	public:
		std::string prevExonID;
		std::string nextExonID;
		int prevExonStart;
		int prevExonEnd;
		int nextExonStart;
		int nextExonEnd;
		int boundSiteTransCoord;
		std::string exonboundSeq;

		std::ostream & operator<<(std::ostream&);

		BoundInfo(){}
		~BoundInfo(){}
};

// ----------------------------------------------------------------------------------
// Transcript Class functions
// ----------------------------------------------------------------------------------
// Parameterised constructor
Transcript::Transcript(std::string& _str){

		std::string AAseq;
		std::stringstream ss;
		ss << _str;    
 		ss >> chr >> transStart >> transEnd >> strand >> geneID >> geneName >> geneBiotype >> transID >> transName >> transBiotype >> proteinID >> RNAseq >> CDSseq >> AAseq;
}

// Overload the output operator 
std::ostream & Transcript::operator<<(std::ostream& os){

 		os << chr << "\t" << transStart << "\t" << transEnd << "\t" << strand << "\t" << geneID << "\t" << geneName << "\t" << geneBiotype << "\t" << transID << "\t" << transName << "\t" << transBiotype << "\t" << proteinID;
 		return (os);
}

// ----------------------------------------------------------------------------------
// Exon Class functions
// ----------------------------------------------------------------------------------
// Parameterised constructor
Exon::Exon(std::string& _str){

		std::stringstream ss;
		ss << _str;    
		ss >> transID >> chr >> exonStart >> exonEnd >> rank >> exonID;
	}


// ----------------------------------------------------------------------------------
// General functions
// ----------------------------------------------------------------------------------
// Split the size of a file across processors 
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

int Nchecker(std::string _seq){
        int flag = 0;
        std::string seq = _seq;
        for(std::string::iterator it = seq.begin(); it!= seq.end(); it++) {
        	if(*it != 'A' && *it != 'C' && *it != 'G' && *it != 'T'){
        		flag = 1;
        	}
        }
  return flag;
}


void populateVpos(int _start, int _end, std::vector<int>* _Vpos){
	for(int i(_start); i <= _end; i++){
		_Vpos->push_back(i);
	}
}

// Overload the output operator 
std::ostream & BoundInfo::operator<<(std::ostream& os){ 

	os << "\t" <<  prevExonID << "\t" << nextExonID << "\t" << prevExonStart << "\t" << prevExonEnd << "\t" << nextExonStart << "\t" << nextExonEnd << "\t" << exonboundSeq << std::endl; 
	return(os);
} // end operator overload

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


// ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// MapReduce Functions
// ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// MAP: Read in the exon file and assign a map based on the transIDs 
void MAP_kmers(int itask, MAPREDUCE_NS::KeyValue* kv, void* ptr){

    Userdata* _data = (Userdata*) ptr; 

    if(_data->me == 0){ std::cout << "Debug MAP 1: check it runs" << std::endl; }

    int* START = new int[_data->nprocs];
    int* END   = new int[_data->nprocs];
    int size = _data->input_transcriptfilesize;
    split(&START, &END, size, _data->nprocs);

    std::cout << "Proc: " << _data->me << "\tstart = " << START[_data->me] << "\tend = " << END[_data->me] << std::endl; 

    std::string line;
    std::multimap<std::string, Exon> MAPexons;
    BoundInfo Binfo;

    int counterReg(0);
    int mapindex(0);
	int inputCounter(0);

    // While loop over the file to read in every line into a multimap of exons, with transID as the key
    while( getline(_data->input_exonfile,line) ){
        Exon exon(line);
        line.clear();
	    MAPexons.insert(std::pair<std::string, Exon>(exon.get_transID() , exon));

	    if(_data->chr_map.find(exon.get_chr()) == _data->chr_map.end()){
	    	_data->chr_map.insert(std::pair<std::string, int>(exon.get_chr(), mapindex));
	    	mapindex++;
	    }

    }

	if(_data->me == 0){ std::cout << "Debug MAP 2: Read in Exons size: " << MAPexons.size() << std::endl;}
	if(_data->me == 0){ 
		for(std::map<std::string, int>::iterator it = _data->chr_map.begin(); it != _data->chr_map.end(); ++it){
			std::cout << "chr = " << it->first << ", index = " << it->second << std::endl;
			}
		}

    // -------------------------------------------------------
    // While loop over each trancript in the transcript file
    // -------------------------------------------------------
    while( getline(_data->input_transcriptfile,line) ){

        // Break at the end of the file 
        if(_data->input_transcriptfile.eof()){break;}

        // counterReg is just a counter to make sure that each transcript is only read in once, the size of the transcript 
        // 	file is split by the number of processors into arrays of START and END points
        if(counterReg >= START[_data->me] && counterReg <= END[_data->me]){

            // If a line is incorrect
            if(_data->input_transcriptfile.fail() && !(_data->input_transcriptfile.eof())){
                _data->input_transcriptfile.clear();
                std::cout << "The input contains an error that was ignored: " << line << std::endl;
                continue;
            }

            // Create a varaiable trans of class transcript
            Transcript trans(line);          
	      	inputCounter++;


            // range is the associated iterator of all exons in the transcript
            std::pair< std::multimap<std::string, Exon>::iterator, std::multimap<std::string, Exon>::iterator> range;
            range = MAPexons.equal_range(trans.get_transID());

            int errorflag(0); 
            int rangeCounter(0);
            int seqCounter(0); // position of current sequence in transcript co-ordinates    
            int cdsFlag(0); // 0 = 5' UTR, 1 = cds, 2 = 3' UTR
            int firstExonRankflag(0);
            int firstExonRank(0);
            std::vector<int> Vpos; // genomic position of every nucleotide for the transcript
 			int rangeSize = MAPexons.count(trans.get_transID()); // number of exons
 			int cdsLength = trans.get_CDSseq().length(); // The annotation doesn't include the stop codon, however the sequence does.

 			if(rangeSize == 0){std::cout << "Error: no exons for: " << trans.get_transID() << " :: "<< line << std::endl; continue;}

	        // Populate the vector of all postions -> first chack if the transcript is a single exon,
 			// 	then iterate over the the exon extracting mucleotide locations and exon-exon boundaries
		        if(rangeSize == 1){
 					std::multimap<std::string, Exon>::iterator singleIT = MAPexons.find(trans.get_transID());
 					populateVpos(singleIT->second.get_exonStart(), singleIT->second.get_exonEnd(), &Vpos);

 					//std::cout << trans.get_transID() << " - single exon CDS region" << singleIT->second.get_exonStart() << " - " << singleIT->second.get_exonEnd() << std::endl; 
		         }
		        else{
			        for(int i(1); rangeCounter < rangeSize; i++){
			        	for(std::multimap<std::string, Exon>::iterator it = range.first; it != range.second; ++it){
			        		if(it->second.get_rank() == i){	
			        			populateVpos(it->second.get_exonStart(), it->second.get_exonEnd(), &Vpos);
			        			if(firstExonRankflag == 0){ firstExonRankflag++; firstExonRank = i;}
			        			rangeCounter++;
			        		}
			        	}
			        }

			
			// The annotation doesn't include the stop codon, however the sequence does.
			if(Vpos.size() == cdsLength || Vpos.size() == cdsLength - 3){}
				else{ std::cout << "Debug Vpos: " << trans.get_transID() << ", firstExonRank = " << firstExonRank << ", Vpos = " << Vpos.size() << " = " << trans.get_CDSseq().length() << std::endl;}
 			
 			// needs to be the first exon rank -1 since it gets ++'ed 
            int rankCounter(firstExonRank - 1);

	            // Iterate though the multimap, putting exons in order & finding boundaries
 				// Compare the exon start and end against the CDS region & then extract the +/- 100nt window surrounding the exon-exon boundary
	            do{
	            	// Errorflag checks that there is actually an exon associated with the next exon rank, if not it breaks out of the do-while.
	                if(errorflag == 1){ std::cout << "Missing Exon at proc = " << _data->me << " Trans = " << trans.get_transID() << " rank = " << rankCounter << " of " << MAPexons.count(trans.get_transID()) << std::endl; break; }

	                errorflag = 1;
	                rankCounter++;    

	                // Iterate over the range of exons associated with the transcript
	                for(std::multimap<std::string, Exon>::iterator it = range.first; it != range.second; ++it){

	                    //	Sort the exons in the right order according to the exon rank
	                    if(it->second.get_rank() == rankCounter){
							errorflag = 0;
	                    	
	                    	//  Add all exonWidths to the seqCounter to make sure indexes work later
	                    	if(cdsFlag == 0) {

	                    		cdsFlag = 1;
	                    		seqCounter += it->second.get_exonWidth();
	                    	
			                    Binfo.prevExonID = it->second.get_exonID();
			                    Binfo.prevExonStart = it->second.get_exonStart();
			                    Binfo.prevExonEnd = it->second.get_exonEnd();
			                    continue; // must continue here to avoid the next if statement
			                }
	                    	                	
	                    	// So we know the start of this next exon is coding so we need to put a boundary
	                    	if(cdsFlag == 1){

		                        Binfo.nextExonID = it->second.get_exonID(); 
		                        Binfo.nextExonStart = it->second.get_exonStart();
		                        Binfo.nextExonEnd = it->second.get_exonEnd();

		                        // Check if the boundary (inc depth) fits within the transcript. Have -3 from the total CDS length due to the presence of stop codons.
		                        if( seqCounter < _data->depth || seqCounter + _data->depth + _data->k - 1 > cdsLength){ Binfo.exonboundSeq = "NA"; } 
								else{ Binfo.exonboundSeq = remove_whitespace( trans.get_CDSseq().substr( seqCounter - _data->depth, 2*_data->depth + _data->k - 1) ); }

								// Check of there is a problem with the boundary size for some reason...
								if(Binfo.exonboundSeq.length() != 2*_data->depth + _data->k - 1 && Binfo.exonboundSeq.length() != 2){Binfo.exonboundSeq = "NA"; }

							    trans << _data->output_exonExonJcts;
		                        Binfo << _data->output_exonExonJcts;

			                    Binfo.prevExonID = it->second.get_exonID();
			                    Binfo.prevExonStart = it->second.get_exonStart();
			                    Binfo.prevExonEnd = it->second.get_exonEnd();	                    
			            
			                    // Boundary is at the end of the previous & start of the next exon
			                    seqCounter += it->second.get_exonWidth();
	                    	}
	                    	else{ break; }
	                	}
	                }

	            }while(rankCounter <= MAPexons.count(trans.get_transID()));
	        }

	            // Add boundaries to the mr object, including the start and end genomic postions to make sure they don't duplicate
				if( cdsLength < _data->k){ std::cout << trans.get_transID() << " CDS is less than " << _data->k << ", strand " << trans.get_strand() << std::endl;}
				else{
		                for(int i(0); i <= trans.get_CDSseq().length() - _data->k; i++){ 

		                    std::string key = trans.get_CDSseq().substr(i, _data->k);

		                    if(Nchecker(key) == 0){
		                    	int value [4]  = {trans.get_chr(_data), trans.get_strand(), Vpos[i], Vpos[ i + _data->k - 1]};
			                    kv->add( (char*) key.substr(0, _data->k).c_str(), sizeof(char)*_data->k, (char*) value, sizeof(int)*4);
			                }
		                } // end for adding the k-mers to the mr object
		            }		           

    } // end the while loop over the input

    counterReg++;	
	line.clear();
	
	}

	std::cout << "Proc = " << _data->me << ", total Transcripts = " << inputCounter << std::endl;

}	// end the MAP_kmers function



// --------------------------------------------------------------------------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------------------------------------------------------------------------
void REDUCE_output(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, MAPREDUCE_NS::KeyValue *kv, void *ptr){

    Userdata* _data = (Userdata*) ptr;


     // std::cout << "Proc = " << _data->me << "\tnvalues = "<< nvalues << "\tmultiValue" << multivalue << std::endl;
	char* value;

 	uint64_t nvalues_total;
	CHECK_FOR_BLOCKS(multivalue,valuebytes,nvalues,nvalues_total)
	BEGIN_BLOCK_LOOP(multivalue,valuebytes,nvalues)

	value = multivalue;
	int* Ivalue = (int*) value;

	if(nvalues == 1){
		_data->output_dist << std::string(key).substr(0,_data->k) << "\t" << nvalues << std::endl;	
	}else{

		std::vector<int*> Vvalues;
		Vvalues.push_back(Ivalue);

	    for(int i(0); i < nvalues - 1; i++){
	    
	    	value += valuebytes[i]; 
	     	Ivalue = (int*) value;    

			int flag = 0;
	       	for(std::vector<int*>::iterator it = Vvalues.begin(); it != Vvalues.end(); ++it){
	       		if( (*it)[0] == Ivalue[0] && (*it)[1] == Ivalue[1] && (*it)[2] == Ivalue[2] && (*it)[3] == Ivalue[3] ){flag = 1;}
	       	}
	       	if(flag == 0){Vvalues.push_back(Ivalue);}

	    }

	    _data->output_dist << std::string(key).substr(0,_data->k) << "\t" << Vvalues.size() << std::endl;		
	}

 	END_BLOCK_LOOP
                          
} // end the REDUCE_output function



#endif
