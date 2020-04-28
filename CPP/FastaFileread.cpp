// --------------------------------------------------------------------------------------------------------------------------------------------------------------------
// fasta_fileread.cpp
// Sam Humphrey  April 2020
//
// Open, Read and CloseFASTA codes were apapted from:::   
// http://www.cse.msu.edu/~yannisun/cse891/hmm-EM/fasta.c
//
// Functions for opening, reading and returning the reference sequences 
// --------------------------------------------------------------------------------------------------------------------------------------------------------------------

#include "FastaFileread.h"

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
void fileread(const std::string fname, std::vector<std::string>* _seq, Chr_Map* _chr_map, std::vector<int>* _L){

        FASTAFILE* ffp;
        char* seq;
        char* chr_name;
        int L;
        int counter(0);

        std::string fullname = fname;

        // open and read in the fasta file
        ffp = OpenFASTA((char*) fullname.c_str());

        while(ReadFASTA(ffp, &seq, &chr_name, &L)){

                _chr_map->insert(std::pair<std::string, int>(std::string(chr_name), counter));
                counter++;

                std::string tmp_seq(seq);
                tmp_seq = remove_whitespace(tmp_seq);

                _seq->push_back(tmp_seq);
                _L->push_back(L);
                }

        CloseFASTA(ffp);
} // end fileread function














