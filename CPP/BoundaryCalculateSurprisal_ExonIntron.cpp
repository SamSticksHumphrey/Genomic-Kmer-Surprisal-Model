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
		std::ifstream input_file;
		std::string input_filename;

		std::ifstream input_distfile;
		std::string input_distfilename;

		std::ofstream output_file;
		std::string output_filename;

};

class Boundary{
	private:
		std::string chr;
		int strand;
		std::string geneID;
		std::string geneName;
		std::string geneBiotype;
		std::string transID;
		std::string transName;
		int transStart;
		int transEnd;
		std::string transBiotype;
		std::string proteinID;
	 	std::string exonID;
	 	std::string exonStart;
	 	std::string exonEnd;
	 	std::string exon5p_boundSeq;
	 	std::string exon3p_boundSeq;

	public:
		Boundary(){}
		Boundary(std::string);
		~Boundary(){}

		std::string get_exon3p_boundSeq(){return exon3p_boundSeq;}
		std::string get_exon5p_boundSeq(){return exon5p_boundSeq;}

		std::ostream & operator<<(std::ostream&);
};

typedef std::map<std::string, int> MapDist;

// --------------------------------------------------------------------------------------------------------------------------------------------------------------------
//      Main function  
// --------------------------------------------------------------------------------------------------------------------------------------------------------------------
int main(int narg, char **args){

    // Argument handling 
    Userdata maindata;
    maindata.jobname = std::string(args[1]);
    maindata.k = atoi(args[2]);
    maindata.depth = atoi(args[3]);
    maindata.input_filename = std::string(args[4]);
    maindata.input_distfilename = std::string(args[5]);

 	maindata.output_filename = std::string(maindata.jobname) + ".txt";
    
	// Open the input file for reading, the jct output & the error file, just in case... 
	std::cout << "Debug 1:" <<  maindata.input_distfilename << std::endl;
    maindata.input_distfile.open(maindata.input_distfilename.c_str());
    if(!(maindata.input_distfile)){std::cout << "Couldn't open input distribution file" << std::endl; exit(0);}
  

    MapDist Mdist;

    std::stringstream ss;
    std::string line;

    int long totalOccuance(0);

    while( getline(maindata.input_distfile,line) ){
   		
   		int occurance;
    	std::string kmer;

    	ss << line;
    	ss >> kmer >> occurance;

    	Mdist.insert(std::pair<std::string, int>(kmer, occurance));

    	if(occurance < 0 || totalOccuance < 0){
    		std::cout << "error, occurance , 0. occurance = " << occurance << "\ttotalOccuance = " << totalOccuance << "\t kmer = " << kmer << std::endl;
    	}

    	totalOccuance += (int long) occurance;

    	line.clear();
    	ss.clear(); ss.str("");
    }

	std::cout << std::endl << "Debug 2: totalOccuance = " <<  totalOccuance << std::endl;
	maindata.input_distfile.clear(); maindata.input_distfile.close();

	// Open the input file for reading, the jct output & the error file, just in case... 
	std::cout << std::endl << "Debug 3:" <<  maindata.input_filename << std::endl;  
    maindata.input_file.open(maindata.input_filename.c_str());
    if(!(maindata.input_file)){std::cout << "Couldn't open input file" << std::endl; exit(0);}

	std::cout << std::endl << "Debug 4:" <<  maindata.output_filename << std::endl;
    maindata.output_file.open(maindata.output_filename.c_str());
    if(!(maindata.output_file)){std::cout << "Couldn't open input file" << std::endl; exit(0);}

	while(getline(maindata.input_file,line)){

			// If a line is incorrect
			if(maindata.input_file.fail() && !(maindata.input_file.eof())){
				maindata.input_file.clear();
				std::cout << "The datafile contains an error that was ignored\n" << std::endl;
				continue;
			}

			Boundary bound(line);
			line.clear();

			std::vector<double> V3surprisal, V5surprisal;
			std::vector<std::string> VdiNuc_exon3p, VdiNuc_exon5p;

			bound << maindata.output_file;

			if( bound.get_exon5p_boundSeq() ==  "NA"){
				for(int i(0); i < 2*maindata.depth; i++){maindata.output_file << "\t" << -1;}
				for(int i(0); i < 2*maindata.depth; i++){maindata.output_file << "\t" << "NA";}
			} else {
				for(int i(0); i < bound.get_exon5p_boundSeq().length() - maindata.k + 1; i++){

					std::string seq = bound.get_exon5p_boundSeq().substr(i, maindata.k);
			        MapDist::iterator mapKmer = Mdist.find(seq);

			        if(mapKmer == Mdist.end()){ V5surprisal.push_back(0); }
					else{ 
						int occur = mapKmer->second; 
						double surprisal = log2( totalOccuance ) - log2( occur );
						V5surprisal.push_back(surprisal);
					}

					std::string diNuc_exon5p = bound.get_exon5p_boundSeq().substr(i, 2);
					VdiNuc_exon5p.push_back(diNuc_exon5p);
				}
	
				for(std::vector<double>::iterator it = V5surprisal.begin(); it != V5surprisal.end(); ++it){
					maindata.output_file << "\t" << *it;
				}
				
				for(std::vector<std::string>::iterator it = VdiNuc_exon5p.begin(); it != VdiNuc_exon5p.end(); ++it){
					maindata.output_file << "\t" << *it;
				}
				
			} 


			if( bound.get_exon3p_boundSeq() ==  "NA"){
				for(int i(0); i < 2*maindata.depth; i++){maindata.output_file << "\t" << -1;}
				for(int i(0); i < 2*maindata.depth; i++){maindata.output_file << "\t" << "NA";}
			} else {
				for(int i(0); i < bound.get_exon3p_boundSeq().length() - maindata.k + 1; i++){

					std::string seq = bound.get_exon3p_boundSeq().substr(i, maindata.k);
			        MapDist::iterator mapKmer = Mdist.find(seq);

			        if(mapKmer == Mdist.end()){ V3surprisal.push_back(0); }
					else{ 
						int occur = mapKmer->second; 
						double surprisal = log2( totalOccuance ) - log2( occur );
						V3surprisal.push_back(surprisal);
					}

					std::string diNuc_exon3p = bound.get_exon3p_boundSeq().substr(i, 2);
					VdiNuc_exon3p.push_back(diNuc_exon3p);
				}
	
				for(std::vector<double>::iterator it = V3surprisal.begin(); it != V3surprisal.end(); ++it){
					maindata.output_file << "\t" << *it;
				}
				
				for(std::vector<std::string>::iterator it = VdiNuc_exon3p.begin(); it != VdiNuc_exon3p.end(); ++it){
					maindata.output_file << "\t" << *it;
				}
				
			} 

			maindata.output_file << std::endl;
	}

	std::cout << "All Done :)" << std::endl;
	maindata.output_file.clear(); maindata.output_file.close();
	return(0);
}


Boundary::Boundary(std::string _line){
	std::stringstream ss;
	ss << _line;  
	ss >> chr >> strand >> geneID >> geneName >> geneBiotype >> transID >> transName >> transStart >> transEnd >> transBiotype >> proteinID >> exonID >> exonStart >> exonEnd >>  exon5p_boundSeq >>  exon3p_boundSeq;

	if(exon5p_boundSeq.length() == 0 || exon3p_boundSeq.length() == 0){std::cout << "Error with the read in of: " << _line << std::endl;}

}


std::ostream & Boundary::operator<<(std::ostream& os){

	os << chr << "\t" << strand << "\t" << geneID << "\t" << geneName << "\t" << geneBiotype << "\t" << transID << "\t" << transName << "\t" << transStart << "\t" << transEnd << "\t" << transBiotype << "\t" << proteinID << "\t" << exonID << "\t" << exonStart << "\t" << exonEnd << "\t" <<  exon5p_boundSeq << "\t" <<  exon3p_boundSeq;
	return(os);
}







