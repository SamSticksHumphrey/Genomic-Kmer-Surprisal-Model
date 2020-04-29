# ------------------------------------------------------
# Top level makefile for the Distributions directory
# ------------------------------------------------------

# Declare variables
EXE_DIR =	./
MPICPP =	mpic++
MPICPP_FLAGS =	-g -O -I/locationOfMapReduceMPI/mrmpi/7Apr14/src/
USRLIB =	/locationOfMapReduceMPI/mrmpi/7Apr14/src/libmrmpi_mpicc.so

all: AnnotateUnsplicedTranscripts DNAdistribution UnsplicedDistribution CDSdistribution AAdistribution BoundaryCalculateSurprisal_ExonIntron BoundaryCalculateSurprisal_Nucleotides BoundaryCalculateSurprisal_AminoAcids DinucleotideKmerLookup

AnnotateUnsplicedTranscripts: FastaFileread.cpp AnnotateUnsplicedTranscripts.cpp
	g++ $(FASTAFILEREAD) FastaFileread.cpp AnnotateUnsplicedTranscripts.cpp -o $(EXE_DIR)annotate_unspliced_transcripts

DNAdistribution: DNADistribution_Main.cpp $(USRLIB)
	$(MPICPP) $(MPICPP_FLAGS) DNADistribution_Main.cpp $(USRLIB) -o $(EXE_DIR)dna_distribution

UnsplicedDistribution: UnsplicedDistribution_Main.cpp $(USRLIB)
	$(MPICPP) $(MPICPP_FLAGS) UnsplicedDistribution_Main.cpp $(USRLIB) -o $(EXE_DIR)unspliced_distribution

CDSdistribution: CodingSequenceDistribution_Main.cpp $(USRLIB)
	$(MPICPP) $(MPICPP_FLAGS) CodingSequenceDistribution_Main.cpp $(USRLIB) -o $(EXE_DIR)cds_distribution

AAdistribution: AminoAcidDistribution_Main.cpp $(USRLIB)
	$(MPICPP) $(MPICPP_FLAGS) AminoAcidDistribution_Main.cpp $(USRLIB) -o $(EXE_DIR)aa_distribution

BoundaryCalculateSurprisal_ExonIntron: BoundaryCalculateSurprisal_ExonIntron.cpp
	g++ BoundaryCalculateSurprisal_ExonIntron.cpp -o $(EXE_DIR)boundary_calculate_surprisal_exon_intron

BoundaryCalculateSurprisal_Nucleotides: BoundaryCalculateSurprisal_Nucleotides.cpp
	g++ BoundaryCalculateSurprisal_Nucleotides.cpp -o $(EXE_DIR)boundary_calculate_surprisal_nucleotides

BoundaryCalculateSurprisal_AminoAcids: BoundaryCalculateSurprisal_AminoAcids.cpp
	g++ BoundaryCalculateSurprisal_AminoAcids.cpp -o $(EXE_DIR)boundary_calculate_surprisal_amino_acids

DinucleotideKmerLookup: DinucleotideKmerLookup.cpp
	g++ DinucleotideKmerLookup.cpp -o $(EXE_DIR)dinucleotide_kmer_lookup


clean:
	rm $(EXE_DIR)annotate_unspliced_transcripts $(EXE_DIR)dna_distribution $(EXE_DIR)unspliced_distribution $(EXE_DIR)cds_distribution $(EXE_DIR)aa_distribution $(EXE_DIR)boundary_calculate_surprisal_exon_intron $(EXE_DIR)boundary_calculate_surprisal_nucleotides $(EXE_DIR)boundary_calculate_surprisal_amino_acids $(EXE_DIR)dinucleotide_kmer_lookup
