# K-mer surprisal model Codes
## Sam Humphrey, April 2020
## Email: sam.humphrey@postgrad.manchester.ac.uk

## R codes

## 1. Annotation

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


## 4. Figure Generation

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















