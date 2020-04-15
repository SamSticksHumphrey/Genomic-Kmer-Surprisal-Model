// --------------------------------------------------------------------------------------------------------------------------------------------------------------------
// CodingSequenceDistribution_Main.cpp 
// Sam Humphrey, April 2020
//
// This MapReduce MPI code outputs the distribution of all annotated coding sequence of size k.
// The code handles all arguments including the k-mer length in nucleotides, (k).
// 
// MapReduce:
// Create a lookup map of a transcriptID - exonID pair to the exon information, read in from the input file. 
// MAP the key = transcriptID against a CSV std::string of value = exonID, transcriptID, exonRank, exonWidth, exonStart, DNAsequence 
// Collate on the transcriptID
// REDUCE_ReMapAAmersGenomeCoords, clean up the multivalues into a vector, create a transcript & set it with the information from the 1st Exon.
//      If it's protein coding then stick the exons together and annotate the full coding sequence, meanwhile produce a vector of every nucleotide's position, 
//      within the CDS. Check for any errors, (CDS < k, or not % 3), and add to the mr object. key = Key object (geneID, k-mer start, k-mer end, k-mer seq), value = k-mer seq.
// Collate on the k-mer with cenomic annotation
// REDUCE_ReMapAAmers, remaps just the sequences, now independent of genomic location, to assure that the same sequence at the same location isn't repeated, key = k-mer seq, values isn't important.
// Collate on the k-mer sequences
// REDUCE_output, output the k-mer sequence against the number of values.
// --------------------------------------------------------------------------------------------------------------------------------------------------------------------

#include "CodingSequenceDistribution.h"

// --------------------------------------------------------------------------------------------------------------------------------------------------------------------
//      Main function  
// --------------------------------------------------------------------------------------------------------------------------------------------------------------------
int main(int narg, char **args){

        double tstart = MPI_Wtime();

        // Argument handling 
        Userdata maindata;
        maindata.jobname = std::string(args[1]);
        maindata.k = atoi(args[2]);
        maindata.depth = atoi(args[3]);
        maindata.input_transcriptfilename = std::string(args[4]);
        maindata.input_transcriptfilesize = atoi(args[5]);
        maindata.input_exonfilename = std::string(args[6]);

        // Stuff for MPI
        MPI_Init(&narg,&args);
        MPI_Comm_rank(MPI_COMM_WORLD,&maindata.me);
        MPI_Comm_size(MPI_COMM_WORLD,&maindata.nprocs); 

        // std::cout << maindata.me << "\t" << maindata.nprocs << std::endl;

        // // Output filename allocation
        std::string label = std::to_string(maindata.me);

        maindata.output_distname = std::string("Results/") + maindata.jobname + "_proc" + label + ".txt";
        maindata.output_exonExonJctsname = std::string("Results/") + maindata.jobname + "_Boundaries_proc" + label + ".txt";
        if(maindata.me == 0){ std::cout << "Debug 1: check it runs. k = " << maindata.k << "\tfilesize = " << maindata.input_transcriptfilesize << std::endl << std::endl; }

        // Open the input file for reading, the jct output & the error file, just in case... 
        maindata.input_transcriptfile.open(maindata.input_transcriptfilename.c_str());
        if(!(maindata.input_transcriptfile)){std::cout << "Couldn't open input file" << std::endl;  MPI_Finalize(); exit(0);}
        if(maindata.me == 0){ std::cout << "Debug 2:" <<  maindata.input_transcriptfilename << std::endl; } 

        // Open the input file for reading, the jct output & the error file, just in case... 
        maindata.input_exonfile.open(maindata.input_exonfilename.c_str());
        if(!(maindata.input_exonfile)){std::cout << "Couldn't open Annotation file" << std::endl;  MPI_Finalize(); exit(0);}
        if(maindata.me == 0) {std::cout << std::endl << "Debug 3:" <<  maindata.input_exonfilename << std::endl;}

        maindata.output_exonExonJcts.open(maindata.output_exonExonJctsname.c_str());
        if(!(maindata.output_exonExonJcts)){std::cout << "Couldn't open output EIJct file" << std::endl << std::endl;  MPI_Finalize(); exit(0);}
        if(maindata.me == 0) {std::cout << std::endl << "Debug 5:" <<  maindata.output_exonExonJctsname << std::endl ;}

        // Create a MapReduce object
        MAPREDUCE_NS::MapReduce *mr = new MAPREDUCE_NS::MapReduce(MPI_COMM_WORLD);
        mr->memsize = 5120;
        mr->verbosity = 1;
        mr->timer = 1;

        // MAP_kmers - read the transcript input file & add k-mer sequences to the mr object 
        mr->map(maindata.nprocs, MAP_kmers, &maindata);

        // Input & junction files is no needed anymore
        maindata.output_exonExonJcts.clear(); maindata.output_exonExonJcts.close();
        maindata.input_transcriptfile.clear(); maindata.input_transcriptfile.close(); 

        // Collate - k-mers into multivalues based on the group size
        mr->collate(NULL);
      
        // Open the output file
        if(maindata.me == 0) {std::cout << std::endl << "Debug 6:" <<  maindata.output_distname << std::endl;}
        maindata.output_dist.open(maindata.output_distname.c_str());
        if(!(maindata.output_dist)){std::cout << "Couldn't open output distribution file" << std::endl << std::endl;  MPI_Finalize(); exit(0);}

        // REDUCE_output - into the number of unique values
        mr->reduce(REDUCE_output,&maindata);  

        MPI_Barrier(MPI_COMM_WORLD);
        maindata.output_dist.clear(); maindata.output_dist.close();
        delete mr; 

        if(maindata.me == 0){std::cout << "All finished :) " << std::endl;}

        MPI_Finalize(); 
        return(0);
}   // end main




