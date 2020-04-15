// --------------------------------------------------------------------------------------------------------------------------------------------------------------------
// AnnotateUnsplicedTranscripts.cpp 
// Sam Humphrey, April 2020
//
// Quick script to generate a k-mer, di-nucleotide lookup table
// --------------------------------------------------------------------------------------------------------------------------------------------------------------------

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <string>
#include "sys/stat.h"
#include <ctype.h>
#include <cmath>
#include <vector>
#include <iterator>
#include <sstream>
#include <fstream>
#include <map>
#include <utility>


class KmerDinucs{
	private:
		std::string kmer;
		int AA;
		int AC;
		int AG;
		int AT;
		int CA;
		int CC;
		int CG;
		int CT;
		int GA;
		int GC;
		int GG;
		int GT;
		int TA;
		int TC;
		int TG;
		int TT;

	public:
		KmerDinucs(){}
		KmerDinucs(std::string);
		~KmerDinucs(){}
		std::ostream & operator<<(std::ostream&);
};

char _int_to_base [4] = {'A', 'C', 'G', 'T'};
std::string decode_kmer_from_intval(int, int);

// -------- Main function  -------------------
int main(int narg, char **args){
	
	std::vector<KmerDinucs> vKmerDinucs;
	std::vector<std::string> vKmers;
	int k = atoi(args[1]);

    // Output filename allocation
    std::string kchar = std::to_string(k);
	std::string outfileName = "KmerDinucleotideLookupTable_k" + kchar + ".txt";

	std::ofstream outfile;
	outfile.open(outfileName);

 	outfile << "kmer" << "\t" << "AA" << "\t" << "AC" << "\t" << "AG" << "\t" << "AT" << "\t" << "CA" << "\t" << "CC" << "\t" << "CG" << "\t" << "CT" << "\t" << "GA" << "\t" << "GC" << "\t" << "GG" << "\t" << "GT" << "\t" << "TA" << "\t" << "TC" << "\t" << "TG" << "\t" << "TT" << std::endl;

	for(int i = 0; i < pow(4, k); i++){
		vKmers.push_back(decode_kmer_from_intval(i, k));
	}


	for(std::vector<std::string>::iterator it = vKmers.begin(); it != vKmers.end(); ++it){
		KmerDinucs out(*(it));

		out << outfile; 
	}

	outfile.clear();
	outfile.close();

}

// decode_kmer_from_intval: using intergers to define a kmer
std::string decode_kmer_from_intval(int intval, int kmer_length) {

	std::string kmer(kmer_length, ' ');
	for (int i = 1; i <= kmer_length; i++) {
        	int base_num = intval & 3ll;
        	kmer[kmer_length-i] = _int_to_base[base_num];
        	intval = intval >> 2;
        }
 return(kmer);
 }


KmerDinucs::KmerDinucs(std::string _kmer){

		kmer = _kmer;
		AA = 0;
		AC = 0;
		AG = 0;
		AT = 0;
		CA = 0;
		CC = 0;
		CG = 0;
		CT = 0;
		GA = 0;
		GC = 0;
		GG = 0;
		GT = 0;
		TA = 0;
		TC = 0;
		TG = 0;
		TT = 0;

	for(int i = 0; i < _kmer.length() - 1; i++){
		std::string dinuc = _kmer.substr(i, 2);

		if(dinuc == "AA"){ AA++;}
		else if(dinuc == "AC"){AC++;}
		else if(dinuc == "AG"){AG++;}
		else if(dinuc == "AT"){AT++;}
		else if(dinuc == "CA"){CA++;}
		else if(dinuc == "CC"){CC++;}
		else if(dinuc == "CG"){CG++;}
		else if(dinuc == "CT"){CT++;}
		else if(dinuc == "GA"){GA++;}
		else if(dinuc == "GC"){GC++;}
		else if(dinuc == "GG"){GG++;}
		else if(dinuc == "GT"){GT++;}
		else if(dinuc == "TA"){TA++;}
		else if(dinuc == "TC"){TC++;}
		else if(dinuc == "TG"){TG++;}
		else if(dinuc == "TT"){TT++;}
		else{std::cout << "Error; dinuc = " << dinuc << std::endl;}
	}

}

// Overload the output operator 
std::ostream & KmerDinucs::operator<<(std::ostream& os){
		os << kmer << "\t" << AA << "\t" << AC << "\t" << AG << "\t" << AT << "\t" << CA << "\t" << CC << "\t" << CG << "\t" << CT << "\t" << GA << "\t" << GC << "\t" << GG << "\t" << GT << "\t" << TA << "\t" << TC << "\t" << TG << "\t" << TT << std::endl;
 		return (os);
}
















