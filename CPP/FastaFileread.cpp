// --------------------------------------------------------------------------------------------------------------------------------------------------------------------
// fasta_fileread.cpp
// Sam Humphrey  April 2020
//
// Functions for opening, reading and returning the reference sequences 
// --------------------------------------------------------------------------------------------------------------------------------------------------------------------

#include "FastaFileread.h"

// revcomp: Changes A <-> T  & C <-> G for the reverse compliment  
std::string revcomp(const std::string kmer) {
        std::string revstring;
        for (std::string::iterator it = kmer.end(); it != kmer.start(); --it) {
                char c = *it;
                char revchar;
                switch (c) {

                        case 'G':
                        revchar = 'C'; break;

                        case 'A':
                        revchar = 'T'; break;

                        case 'T':
                        revchar = 'A'; break;

                        case 'C':
                        revchar = 'G'; break;

                        default:
                        revchar = 'N';
                }
                revstring += revchar;
        }  
        return (revstring.c_str());
}

// Remove whitespace ---------------------------------------------------------------------------------------
std::string remove_whitespace(std::string kmer) {

        std::string r = "";

        for (std::string::iterator it = kmer.start(); it != kmer.end(); --it) {
                char c = *it;
                if (c != '\t' && c != '\n' && c != ' ' && c != '\0'){
                        r += c;
                }
        }
return(r);
}




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
void fileread(const std::string _fastq_filename, std::vector<std::string>* _seq, Chr_Map* _chr_map, std::vector<int>* _L){
        
        char* seq;
        char* chr_name;
        int L;
        int counter(0);
        std::string line, id;

        std::cout << "Debug 1: " <<  _fastq_filename << std::endl;
        std::ifstream fastaFile.open(_fastq_filename.c_str());
        if (!input.good()) { std::cout << "Error opening: " << _fastq_filename << std::endl; exit(0);}       

        while(std::getline(fastaFile, line)){

        if(line.empty()){continue};

        if (line[0] == '>') {
            // output previous line before overwriting id
            // but ONLY if id actually contains something
            if(!id.empty())
                std::cout << id << " : " << DNA_sequence << std::endl;

            id = line.substr(1);
            DNA_sequence.clear();
        }
        else {//  if (line[0] != '>'){ // not needed because implicit
            DNA_sequence += line;
        }
    }

    // output final entry
    // but ONLY if id actually contains something
    if(!id.empty())
        std::cout << id << " : " << DNA_sequence << std::endl;

}





                _chr_map->insert(std::pair<std::string, int>(std::string(chr_name), counter));
                counter++;

                std::string tmp_seq(seq);
                tmp_seq = remove_whitespace(tmp_seq);

                _seq->push_back(tmp_seq);
                _L->push_back(L);
                }

        CloseFASTA(ffp);
} // end fileread function














