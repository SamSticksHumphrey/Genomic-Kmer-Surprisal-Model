#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string.h>
#include <fstream>
#include <sstream>
#include <ctype.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <map>
#include <iterator>


class Userdata{
	public:
		// Arguments
		int k;
		int depth;
		std::string jobname;

		// file IO
		std::ifstream input_boundFile;
		std::string input_boundFilename;

		std::ifstream input_distFile;
		std::string input_distFilename;

		std::ofstream output_file;
		std::string output_filename;

};

class Boundary{

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
	  	std::string prevExonID;
	 	std::string nextExonID;
	 	int prevExonStart;
	 	int prevExonEnd;
	 	int nextExonStart;
	 	int nextExonEnd;
 		std::string boundSeq;


	public:
		Boundary(){}
		Boundary(std::string);
		~Boundary(){}

		std::string get_boundSeq(){return boundSeq;}
		std::string print(){ return std::string(transID) + " " + prevExonID + " " + nextExonID + " " + boundSeq;}
		void set_boundSeq(std::string _seq){boundSeq = _seq;}

		std::ostream & operator<<(std::ostream&);
};

typedef std::map<std::string, int> MapDist;


std::string convertTripleToSingleAA(std::string);


// --------------------------------------------------------------------------------------------------------------------------------------------------------------------
//      Main function  
// --------------------------------------------------------------------------------------------------------------------------------------------------------------------
int main(int narg, char **args){

    // Argument handling 
    Userdata maindata;
    maindata.jobname = std::string(args[1]);
    maindata.k = atoi(args[2]);
    maindata.depth = atoi(args[3]);
    maindata.input_boundFilename = std::string(args[4]);
    maindata.input_distFilename = std::string(args[5]);

 	maindata.output_filename = std::string(maindata.jobname) + ".txt";
    
	// Open the input file for reading, the jct output & the error file, just in case... 
	std::cout << "Debug 1: k = " << maindata.k << "\tfilename = " <<  maindata.input_distFilename << std::endl;
    maindata.input_distFile.open(maindata.input_distFilename.c_str());
    if(!(maindata.input_distFile)){std::cout << "Couldn't open input distribution file" << std::endl; exit(0);}
  

    MapDist Mdist;

    std::stringstream ss;
    std::string line;

    int long totalOccurance(0);

    while( getline(maindata.input_distFile,line) ){
   		
   		int occurance;
    	std::string kmer;

    	ss << line;
    	ss >> kmer >> occurance;

    	Mdist.insert(std::pair<std::string, int>(kmer, occurance));

    	totalOccurance += occurance;

    	line.clear();
    	ss.clear(); ss.str("");
    }

	std::cout << std::endl << "Debug 2: totalOccurance = " <<  totalOccurance << std::endl;
	maindata.input_distFile.clear(); maindata.input_distFile.close();

	// Open the input file for reading, the jct output & the error file, just in case... 
	std::cout << std::endl << "Debug 3:" <<  maindata.input_boundFilename << std::endl;  
    maindata.input_boundFile.open(maindata.input_boundFilename.c_str());
    if(!(maindata.input_boundFile)){std::cout << "Couldn't open input boundary file: " << maindata.input_boundFilename << std::endl; exit(0);}

	std::cout << std::endl << "Debug 4:" <<  maindata.output_filename << std::endl;
    maindata.output_file.open(maindata.output_filename.c_str());
    if(!(maindata.output_file)){std::cout << "Couldn't open output file" << std::endl; exit(0);}

	while(getline(maindata.input_boundFile,line)){

			// If a line is incorrect
			if(maindata.input_boundFile.fail() && !(maindata.input_boundFile.eof())){
				maindata.input_boundFile.clear();
				std::cout << "The datafile contains an error that was ignored\n" << std::endl;
				continue;
			}

			Boundary bound(line);
			line.clear();

			std::vector<double> Vsurprisal;

			if(bound.get_boundSeq() == "NA"){ continue; }
			else{
				for(int i(0); i < bound.get_boundSeq().length() - maindata.k + 1; i++){

					std::string seq = bound.get_boundSeq().substr(i, maindata.k);

			        MapDist::iterator mapKmer = Mdist.find(convertTripleToSingleAA(seq));

			        if(mapKmer == Mdist.end()){ Vsurprisal.push_back(0); }
					else{ 
						
						int occur = mapKmer->second; 
						double surprisal = log2( totalOccurance ) - log2( occur );
						Vsurprisal.push_back(surprisal);
					}
				}

				bound.set_boundSeq(convertTripleToSingleAA(bound.get_boundSeq()));

				bound << maindata.output_file;

				for(std::vector<double>::iterator it = Vsurprisal.begin(); it != Vsurprisal.end(); ++it){
					maindata.output_file << "\t" << *it;
				}

				maindata.output_file << std::endl;

			}
	}

	maindata.output_file.clear(); maindata.output_file.close();
	std::cout << "All Done :)" << std::endl;
	return(0);
}


Boundary::Boundary(std::string _line){

	std::stringstream ss;
	ss << _line; 
	ss >> chr >> transStart >> transEnd >> strand >> geneID >> geneName >> geneBiotype >> transID >> transName >> transBiotype >> proteinID >>  prevExonID >> nextExonID >> prevExonStart >> prevExonEnd >> nextExonStart >> nextExonEnd >> boundSeq; 
}



std::ostream & Boundary::operator<<(std::ostream& os){

	os << chr << "\t" << transStart << "\t" << transEnd << "\t" << strand << "\t" << geneID << "\t" << geneName << "\t" << geneBiotype << "\t" << transID << "\t" << transName << "\t" << transBiotype << "\t" << proteinID << "\t" <<  prevExonID << "\t" << nextExonID << "\t" << prevExonStart << "\t" << prevExonEnd << "\t" << nextExonStart << "\t" << nextExonEnd << "\t" << boundSeq; 
	return(os);
}

std::string convertTripleToSingleAA(std::string _seq){
	std::string seqOut;
	int AAcounter = 0;
	for(std::string::iterator it = _seq.begin(); it != _seq.end(); it++){
		if(AAcounter % 3 == 0){
			char tmp = *it;
			seqOut.push_back(tmp);
		}
		AAcounter++;
	}
	return seqOut;
}






