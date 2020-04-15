// --------------------------------------------------------------------------------------------------------------------------------------------------------------------
// UnsplicedDistribution_Main.cpp 
// Sam Humphrey, April 2020
// --------------------------------------------------------------------------------------------------------------------------------------------------------------------

#ifndef RNADist_Unspliced_H
#define RNADist_Unspliced_H

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

		std::ofstream output_exonIntronJcts;
		std::string output_exonIntronJctsname;
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
		std::string transSeq; 

	public:	
		Transcript(){};	
	    Transcript(std::string&);
		~Transcript(){}

		int get_chr(Userdata* _data){return _data->chr_map.find(chr)->second;}
		std::string get_transID(){return transID;}
		std::string get_transSeq(){return transSeq;}
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
		std::string exonID;
		std::string nextExonID;
		int exonStart;
		int exonEnd;
		std::string exon_5boundSeq;
		std::string exon_3boundSeq; 

		std::ostream & operator<<(std::ostream&);

		BoundInfo(){}
		~BoundInfo(){}
};


// ----------------------------------------------------------------------------------
// Transcript Class functions
// ----------------------------------------------------------------------------------
// Parameterised constructor
Transcript::Transcript(std::string& _str){

		std::stringstream ss;
		ss << _str;    
  		ss >> chr >> transStart >> transEnd >> strand >>  geneName >> geneID >> geneBiotype >> transID >> transName >> transBiotype >> proteinID >> transSeq;
    }

// Overload the output operator 
std::ostream & Transcript::operator<<(std::ostream& os){

 		os << chr << "\t" << strand << "\t" << geneID << "\t" << geneName << "\t" << geneBiotype << "\t" << transID << "\t" << transName << "\t" << transStart << "\t" << transEnd << "\t" << transBiotype << "\t" << proteinID;
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


void populateVpos(int _strand, int _start, int _end, std::vector<int>* _Vpos){

	if(_strand == 1){
		for(int i(_start); i <= _end; i++){
			_Vpos->push_back(i);
		}
	}else if(_strand == -1){
		for(int i(_end); i >= _start; i--){
			_Vpos->push_back(i);
		}
	}
}


// Overload the output operator 
std::ostream & BoundInfo::operator<<(std::ostream& os){ 

	os << "\t" <<  exonID << "\t" << exonStart << "\t" << exonEnd << "\t" << exon_5boundSeq << "\t" << exon_3boundSeq  << std::endl; 
	return(os);
} // end operator overload



// ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// MapReduce Functions
// ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// MAP: Read in the exon file and assign a map based on the transcriptIDs 
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
	// While loop over the file
	// -------------------------------------------------------
	while( getline(_data->input_transcriptfile,line) ){

		// Break at the end of the file 
		if(_data->input_transcriptfile.eof()){break;}
		if(counterReg >= START[_data->me] && counterReg <= END[_data->me]){

			// If a line is incorrect
			if(_data->input_transcriptfile.fail() && !(_data->input_transcriptfile.eof())){
				_data->input_transcriptfile.clear();
				std::cout << "The datafile contains an error that was ignored\n" << std::endl;
				continue;
			}

			Transcript trans(line);            

			std::pair< std::multimap<std::string, Exon>::iterator, std::multimap<std::string, Exon>::iterator > range;
			range = MAPexons.equal_range(trans.get_transID());
			int counter(0);
			int errorflag(0);
			int seqCounter(0); // position of current sequence in transcript co-ordinates    
          	int rangeSize = MAPexons.count(trans.get_transID()); // number of exons
 			
 			if(rangeSize > 1){

				// Iterate though the multimap, putting exons in order & finding boundaries
				do{

					if(errorflag == 1){
						std::cout << "Missing Exon at proc = " << _data->me << " Trans = " << trans.get_transID() << " rank = " << counter << " of " << MAPexons.count(trans.get_transID()) << std::endl;
						break;
					}

					errorflag = 1;
					counter++;
					
					int long transCoord_5bound = 0;
					int long transCoord_3bound = 0;
					int long transLength = trans.get_transSeq().length();

					for(std::multimap<std::string, Exon>::iterator it = range.first; it != range.second; ++it){

						// Sort the exons in the right order according to the exon rank
						// Then add the exon sequence to the overall trans sequence.
						if(it->second.get_rank() == counter){

						errorflag = 0;	
						Binfo.exonID = it->second.get_exonID();
						Binfo.exonStart = it->second.get_exonStart();
						Binfo.exonEnd = it->second.get_exonEnd();

						transCoord_5bound = Binfo.exonStart - trans.get_transStart();
						transCoord_3bound = Binfo.exonEnd - trans.get_transStart() + 1;

						if(transCoord_5bound < _data->depth || transCoord_5bound + _data->depth + _data->k - 1 > transLength){Binfo.exon_5boundSeq = "NA";}
						else{Binfo.exon_5boundSeq = trans.get_transSeq().substr(transCoord_5bound - _data->depth, 2*_data->depth + _data->k - 1);}

						if(transCoord_3bound < _data->depth || transCoord_3bound + _data->depth + _data->k - 1 > transLength){Binfo.exon_3boundSeq = "NA";}
						else{Binfo.exon_3boundSeq = trans.get_transSeq().substr(transCoord_3bound - _data->depth, 2*_data->depth + _data->k - 1);}
		
					
						trans << _data->output_exonIntronJcts;
						Binfo << _data->output_exonIntronJcts;

						}
					}
				}while(counter <= MAPexons.count(trans.get_transID()));
			}

#ifdef Dist
			int start, end;
			std::string seq = trans.get_transSeq();
			if(seq.length() < _data->k){continue;}
	            else{
	                for(int i(0); i < seq.length() - _data->k; i++){

	                    std::string key = seq.substr(i, _data->k);

	                    if(Nchecker(key) == 0){
	                        int value [4] = { trans.get_chr(_data), trans.get_strand(), trans.get_transStart() + i, trans.get_transStart() + i + _data->k - 1 };
		                    kv->add( (char*) key.substr(0, _data->k).c_str(), sizeof(char)*_data->k, (char*) value, sizeof(int)*4 );
		                }
	                } // end for adding the k-mers to the mr object
	            
	                line.clear();
				}	
#endif
			}
	counterReg++;
	} // end the while loop over the input

}// end the MAP_kmers function


// --------------------------------------------------------------------------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------------------------------------------------------------------------
void REDUCE_output(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, MAPREDUCE_NS::KeyValue *kv, void *ptr){

    Userdata* _data = (Userdata*) ptr;

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









