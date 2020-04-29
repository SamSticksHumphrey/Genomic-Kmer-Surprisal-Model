# Genomic Information: k-mer Surprisal Model 
## Sam Humphrey, April 2020
## Email: sam.humphrey@postgrad.manchester.ac.uk

# C++ codes

## 1. Annotation

### AnnotateUnsplicedTranscripts.cpp

Arguments: 

* Output filename
* AllTranscript *.txt* file
* DNA *.fa* file

Output: The annotation of the transcript plus the full sequence from transcription start to transcription end.


## 2. Distributions and boundaries


All distribution codes are produced in C++ using MapReduce MPI (https://mapreduce.sandia.gov/) and open mpi (version 4.0.0). 

### DNADistribution_Main.cpp

Arguments: 

* Output filename
* DNA *.fa* file
* k 

Output: 

* all k-mers and the number of occurances in the DNA

Note: This code reads in each chromosome onto two different cores, one is then reverse complimented, hence it requires 2 x the number of chromosomes for the organism.


### CodingSequenceDistribution_Main.cpp

Arguments: 

* Output filename
* k
* boundary depth
* PCTranscript *.txt* file
* PCTranscript *.txt* file length
* PCannotations *.txt* file

Output files: 

* All k-mers and the number of occurances in the coding sequences
* All exon-exon boudary sequences (+/- depth) for coding sequences


### AminoAcidDistribution_Main.cpp

Arguments: 

* Output filename
* k
* boundary depth
* PCTranscript *.txt* file
* PCTranscript *.txt* file length
* PCannotations *.txt* file

Output files: 

* All k-mers and the number of occurances in the AA sequences
* All exon-exon boudary sequences (+/- depth) for AA sequences


### UnsplicedDistribution_Main.cpp

Arguments: 

* Output filename
* k
* boundary depth
* UnsplicedTranscripts *.txt* file
* UnsplicedTranscripts *.txt* file length
* AllExons *.txt* file

Outputs files: 

* All exon-intron boudary sequences (+/- depth) for all boundaries
* (optional as very large) All k-mers and the number of occurances in the pre-mRNA sequences


## 3. Calculate Surprisal

These C++ codes ouput the surprisal of the k-mer for the sequence string defined in the boudary file. They do not require MPI.

### BoundaryCalculateSurprisal_ExonIntron.cpp

Arguments: 

* Output filename
* k
* boundary depth
* Exon-Intron boundaries file (from *UnsplicedDistribution_Main.cpp*)
* Distribution file (can be DNA or Unspliced)

Outputs file: 

* Exon-intron boundary surprisal file


### BoundaryCalculateSurprisal_Nucleotides.cpp

Arguments: 

* Output filename
* k
* boundary depth
* Exon-exon boundaries file (from *CodingSequenceDistribution_Main.cpp*)
* Distribution file (can be DNA, Unspliced or coding sequences)

Outputs file: 

* Exon-exon boundary surprisal file for the distribution


### BoundaryCalculateSurprisal_AminoAcids.cpp

* Output filename
* k
* boundary depth
* Amino acid exon-exon boundaries file
* Amino acid distribution file

Outputs file: 

* Amino acid exon-exon boundary surprisal file


## 4. Figure Generation

### DinucleotideKmerLookup.cpp

Arguments: 

* k

Outputs file: 

* Lookup table for number of dinucleotides that occur within each k-mer














