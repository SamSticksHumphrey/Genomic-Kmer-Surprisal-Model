# K-mer surprisal model Codes
## Sam Humphrey, April 2020
## Email: sam.humphrey@postgrad.manchester.ac.uk

## Overview
This is a set of scripts used to generate the data used in the paper: *A model of k-mer surprisal to quantify local sequence information content surrounding splice regions.*

## Notes
* reverse strand exons/transcripts have been reversed.
* At the time this work was done, the current version of ensembl was version 97.
* All C++ codes can be compiled by updating the makefile with the location of the MapReduce MPI source files and typeing *make* in the CPP directory.
* DNA, CDS and AA stand for DNA, coding mRNA and amino acid sequences.

## Annotations
All data for use in these scripts can be downloaded from the Ensembl ftp server (ftp://ftp.ensembl.org/pub/). I downloaed DNA, cDNA, ncRNA, CDS, PEP and the gtf files for a species, example example for *H. sapiens*:

```
wget https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/cds/Homo_sapiens.GRCh38.cds.all.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz
wget https://ftp.ensembl.org/pub/current_gtf/homo_sapiens/Homo_sapiens.GRCh38.96.chr.gtf.gz
```

The primary assembly files were manually edited to remove all unlocalised sequences and put into a file named *Homo_sapiens.GRCh38.DNA.fa*.

```
grep -n ">" Homo_sapiens.GRCh38.dna.primary_assembly.fa
head -n lineNumber Homo_sapiens.GRCh38.dna.primary_assembly.fa > Homo_sapiens.GRCh38.DNA.fa
```

All references to the DNA sequnece uses this edited file, not the origially downloaded primary assembily.  

### extractTranscriptAnnotationsFromGTF.R

Arguments:

* A string of the species label
* *.gtf* file
* DNA *.fa* file
* cDNA *.fa* file
* CDS *.fa* file
* ncRNA *.fa* file
* AA *.fa* file

During the code - the function **reverseStrandCoords** swaps the coordinates for the reverse strand genes `rev_strand_coord = chromosome_length - for_strand_coord `. 

Outputs:

* *AllTranscript_species.txt* - Information about all transcripts, plus the full cDNA, coding mRNA and AA sequences for each transcript
* *PCTranscript_species.txt* - same, except protein coding sequences only.

* *AllExons_species.txt* - Just the transcript id, chromosome, start, end, exon rank and id for each exon.
* *PCannotations_species.txt* - Just the transcript id, chromosome, start, end, exon rank and exon id for each coding sequence (translation starts/ends are different between this and the exons file).


### AnnotateUnsplicedTranscripts.cpp

Arguments: 

* Output filename
* AllTranscript *.txt* file
* DNA *.fa* file

Output: The annotation of the transcript plus the full sequence from transcription start to transcription end.


## Distributions and Boundaries

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


## Calculate Surprisals

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


## Figure Generation

### DinucleotideKmerLookup.cpp

Arguments: 

* k

Outputs file: 

* Lookup table for number of dinucleotides that occur within each k-mer


### redundancyScript.R

Arguments: 

* Species label
* Distribution type
* Distribution file
* k_start
* k_end

Outputs file: 

* Redundancy for the distribution for all k between k_start and k_end

### redundancyScript_AA.R

Arguments: 

* Species label
* Distribution file
* k_start
* k_end

Outputs file: 

* Redundancy for AA sequences all k between k_start and k_end


### PaperFigures.Rmd

Arguments: 

* Redundancy files for DNA, CDS and AA
* Exon-intron surprisal files (from *BoundaryCalculateSurprisal_ExonIntron.cpp*)
* Exon-exon surprisal files (from *BoundaryCalculateSurprisal_Nucleotides.cpp* and *BoundaryCalculateSurprisal_AminoAcids.cpp*)
* DNA distributions for *H. sapiens* (k = 2 and 12)
* Dinucleotide lookup table (k = 12)


Outputs file: 

* Redundancy for AA sequences all k between k_start and k_end















