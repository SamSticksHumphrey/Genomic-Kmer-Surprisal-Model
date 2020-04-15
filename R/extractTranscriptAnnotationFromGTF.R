#!/usr/bin/env Rscript

# ---------------------------------------------------------------------------
# getExonAnnotationFromGTF.R
# Sam Humphrey April 2020
#
# Read in the .gtf files and associated .fa sequences downloaded from ensembl (ftp://ftp.ensembl.org/pub/).
# Annotate the transcripts, exons and coding sequences for use in this script
# 
#
# ---------------------------------------------------------------------------

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} 

# ---------------------------------------------------------------------------
# Argument handling & librarys
# ---------------------------------------------------------------------------
species <- args[1]
gtfFile <- args[2]
DNAfile <- args[3]
cDNAfastaFile <- args[4]
CDSfastaFile <- args[5]
ncRNAfastaFile <- args[6]
AAfile <- args[7]


library(rtracklayer)
library(plyr)
library(seqinr)
library(stringr)

DNA <- read.fasta(DNAfile, as.string = TRUE, forceDNAtolower=FALSE)
DNAsize <- data.frame(sapply(DNA, nchar))
DNAsize$seqnames <- rownames(DNAsize)
colnames(DNAsize) <- c("size", "seqnames")

annotation.gtf <- import(gtfFile)
cDNA <- read.fasta(cDNAfastaFile, as.string = TRUE, forceDNAtolower=FALSE)
CDS <- read.fasta(CDSfastaFile, as.string = TRUE, forceDNAtolower=FALSE)
ncRNA <- read.fasta(ncRNAfastaFile, as.string = TRUE, forceDNAtolower=FALSE)
AA <- read.fasta(AAfile, as.string = TRUE, forceDNAtolower=FALSE)

reverseStrandCoords <- function(omic_type){
  df <- as.data.frame(subset(annotation.gtf, type == omic_type))
  print(paste0("number of ", omic_type, " = ", length(df[,1])))
  df <- merge(df, DNAsize, by = "seqnames", all = FALSE)
  df$revStart <- ifelse(df$strand == "+", df$start,  df$size - df$end)
  df$end <- ifelse(df$strand == "+", df$end, df$size - df$start)
  df$start <- df$revStart
  df$revStart <-NULL
  df$size <- NULL
  df$strand <- ifelse(df$strand == "+", 1, -1)
  df$gene_name <- word(df$gene_name, 1)		# remove any additional words assocated with the gene name many "(one of many)"s were  found in zebrafish
  
  return(df)
}

#genes <- reverseStrandCoords("gene")
transcripts <- reverseStrandCoords("transcript")
codingSequences <- reverseStrandCoords("CDS")
exons <- reverseStrandCoords("exon")

# The HomoSapiens, MusMusculus and DanioRerio transcript names have these dumb .1/.2 things that are unnecessary 
if(species  %in% c("HomoSapiens", "MusMusculus", "DanioRerio")){
	names(cDNA) <- as.character(read.table(text = names(cDNA), sep = ".")$V1)
	names(CDS) <- as.character(read.table(text = names(CDS), sep = ".")$V1)
	names(ncRNA) <- as.character(read.table(text = names(ncRNA), sep = ".")$V1)
	names(AA) <- as.character(read.table(text = names(AA), sep = ".")$V1)
}

if(species %in% c("Cerevisiae")){
  names(AA) <- as.character(read.table(text = names(AA), sep = ".")$V1)
}

cDNA.df <- data.frame(transcript_id = names(cDNA), RNAseq = as.character(cDNA))
ncRNA.df <- data.frame(transcript_id = names(ncRNA), RNAseq = as.character(ncRNA))
RNA.df <- rbind(cDNA.df, ncRNA.df)
CDS.df <- data.frame(transcript_id = names(CDS), CDSseq = as.character(CDS))
AA.df <- data.frame(protein_id = names(AA), AAseq = as.character(AA))

proteinID_transID <- unique(data.frame(protein_id = codingSequences$protein_id, transcript_id = codingSequences$transcript_id))

if(species %in% c("Cerevisiae", "Pombe")){
  transcripts <- transcripts[,c("seqnames", "start", "end", "strand", "gene_id", "gene_biotype", "transcript_id", "transcript_biotype")]
  } else { transcripts <- transcripts[,c("seqnames", "start", "end", "strand", "gene_id", "gene_name", "gene_biotype", "transcript_id", "transcript_name", "transcript_biotype", "tag")] }

if( species %in% c("Cerevisiae", "Pombe")){transcripts$gene_name <- transcripts$gene_id}
if( species %in% c("Cerevisiae", "Pombe")){transcripts$transcript_name <- transcripts$transcript_id}
    
transcripts <- merge(transcripts, proteinID_transID, by = "transcript_id", all.x = TRUE)
transcripts <- merge(transcripts, RNA.df, by = "transcript_id", all.x = TRUE)
transcripts <- merge(transcripts, CDS.df, by = "transcript_id", all.x = TRUE)
transcripts <- merge(transcripts, AA.df, by = "protein_id", all.x = TRUE)

# There is nothing we can do with retained introns so we should remove them here!
transcripts <- transcripts[!(transcripts$transcript_biotype == "retained_intron"),]
exons <- exons[,c("transcript_id", "seqnames", "start",  "end", "exon_number", "exon_id")]

# having tag = basic removed 16 of 19718 protein coding genes
if(species %in% c("Cerevisiae", "Pombe", "Drosophila", "DanioRerio")){
  PCtranscripts <- subset(transcripts, transcript_biotype == "protein_coding") 
  } else { PCtranscripts <- subset(transcripts, transcript_biotype == "protein_coding" & tag == "basic") }

cat(sprintf("number of PCtranscripts %d\n", length(PCtranscripts[,1])))

# Give codingSequences exonIDs
codingSequences <- codingSequences[c("transcript_id", "seqnames", "start",  "end", "exon_number")]
codingSequences <- merge(codingSequences, exons[,c("transcript_id", "seqnames", "exon_number", "exon_id")], by = c("transcript_id", "seqnames", "exon_number"), all.x = TRUE )

transcripts <- transcripts[,c("seqnames", "start", "end", "strand", "gene_id", "gene_name", "gene_biotype", "transcript_id", "transcript_name", "transcript_biotype", "protein_id", "RNAseq", "CDSseq", "AAseq")]
PCtranscripts <- PCtranscripts[,c("seqnames", "start", "end", "strand", "gene_id", "gene_name", "gene_biotype", "transcript_id", "transcript_name", "transcript_biotype", "protein_id", "RNAseq", "CDSseq", "AAseq")]
exons <- exons[,c("transcript_id", "seqnames", "start",  "end", "exon_number", "exon_id")]
codingSequences <- codingSequences[,c("transcript_id", "seqnames", "start", "end", "exon_number", "exon_id")]

write.table(as.data.frame(transcripts), paste0("/data/rnabiology/shumphrey/Project/1-Annotation/", species, "/allTranscripts_", species, ".txt"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
write.table(as.data.frame(exons), paste0("/data/rnabiology/shumphrey/Project/1-Annotation/", species, "/allExons_", species, ".txt"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
write.table(as.data.frame(PCtranscripts), paste0("/data/rnabiology/shumphrey/Project/1-Annotation/", species, "/PCtranscripts_", species, ".txt"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
write.table(as.data.frame(codingSequences), paste0("/data/rnabiology/shumphrey/Project/1-Annotation/", species, "/PCannotations_", species, ".txt"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

print(paste0("all done for ", species))







