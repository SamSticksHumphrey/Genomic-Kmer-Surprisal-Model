# K-mer surprisal model Codes
## Sam Humphrey, April 2020
## Email: sam.humphrey@postgrad.manchester.ac.uk

This readme is in the process of being updated.


## Overview
This is a set of scripts used to generate the data used in the paper: *A model of k-mer surprisal to quantify local sequence information content surrounding splice regions.*

## Notes
* reverse strand exons/transcripts have been reversed.
* At the time this work was done, the current version of ensembl was version 97.

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

All references to the DNA uses this edited file, not the origially downloaded primary assembily.  

### extractTranscriptAnnotationsFromGTF.R
This R script takes in arguments of the species label plus the 6 files downloaded. All files are read into memory and the chromosome sizes were found. The function **reverseStrandCoords** swaps the coordinates for the reverse strand genes `rev_strand_coord = chromosome_length - for_strand_coord `. 

The script then outputs the set of transcripts, exons, protein coding transcripts and coding sequences, where 




## Distributions
## Calculate Surprisals
















