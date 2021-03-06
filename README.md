# Genomic Information: k-mer Surprisal Model
## Sam Humphrey, April 2020
## Email: sam.humphrey@postgrad.manchester.ac.uk

# Overview
This is a set of scripts used to generate the data used in the paper: *A model of k-mer surprisal to quantify local sequence information content surrounding splice regions.*

## General Notes
* Reverse strand exons/transcripts have been reversed.
* DNA, CDS and AA stand for DNA, coding mRNA and amino acid sequences.

## Raw data
All data for use in these scripts can be downloaded from the Ensembl ftp server (ftp://ftp.ensembl.org/pub/). I downloaded DNA, cDNA, ncRNA, CDS, PEP and the gtf files for a species, example example for *H. sapiens*:

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

## The project 

## 1. Annotation

Once the raw data was downloaded, the files were processed to extract the required information for downsteam analysis using *extractTranscriptAnnotationsFromGTF.R*. Then *AnnotateUnsplicedTranscripts.cpp* was used to produce the set of transcription start - end sequences for all transcripts. 


## 2. Distributions and boundaries

The DNA distribution can be found directly from the *DNA.fa* file. 

Since genes contain multiple transcripts k-mers will be repeated by the analysis if there are multiple protein coding transcripts for a gene. To correct for this, the *CodingSequenceDistribution_Main.cpp* and *AminoAcidDistribution_Main.cpp* map k-mers to genomic coordinates and only count each k-mer once. This method also provides the locations of exon-intron/exon-exon splice boundaries in the sequence, and hence these distribution codes aslo output the sequence window surrounding the splice site. Note: the correction is unnecessary for the *UnsplicedDistribution_Main.cpp*, but the boundaries are. 


## 3. Calculate Surprisal

The boundaries and distributions are then read into the *BoundaryCalculateSurprisal_ExonIntron.cpp*, *BoundaryCalculateSurprisal_Nucleotides.cpp* and *BoundaryCalculateSurprisal_AminoAcids.cpp* codes to calculate the surprisal of all k-mers in the boundary windows.


## 4. Figure Generation

For the figures, in addition to the files already produced - a lookup table to count the dinucleotides within each k-mers, was produced using *DinucleotideKmerLookup.cpp, and the redundancy scores and unique k-mer counts were generated by *redundancyScript.R*. These files were then used to generate all figures using *PaperFigures.Rmd*.














