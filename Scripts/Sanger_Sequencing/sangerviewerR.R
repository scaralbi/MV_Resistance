# Load required packages
library(seqinr)
library(Biostrings)
library(ape)
library(ShortRead)
library(sangerseqR)

# Read the reference fasta sequence with annotations
ref_seq <- read.fasta("SP_prqR1_product.fasta")

# Read the Sanger sequences in ab1 format and convert to DNAStringSet
WT1 <- read.abif("/Volumes/Albi/Cambridge/PhD/Bioinformatics/prqRA/28042023_SP_prqR1_PCR_Products_Sequencing/WT1.ab1")
WT1_seq <- readDNAStringSet(WT1$seq)

WT2 <- read.abif("/Volumes/Albi/Cambridge/PhD/Bioinformatics/prqRA/28042023_SP_prqR1_PCR_Products_Sequencing/WT2.ab1")
WT2_seq <- readDNAStringSet(WT2$seq)

WT3 <- read.abif("/Volumes/Albi/Cambridge/PhD/Bioinformatics/prqRA/28042023_SP_prqR1_PCR_Products_Sequencing/WT3.ab1")
WT3_seq <- readDNAStringSet(WT3$seq)

# Align the Sanger sequences to the reference sequence
aligned_seq <- pairwiseAlignment(WT1_seq, WT2_seq, WT3_seq, ref_seq, type="global")

# Visualize the chromatogram traces
plotChromatograms(WT1_seq)

# Color and plot the quality scores
WT1_qual <- sangerQuality(WT1_seq)
plotQuality(WT1_qual, type = "color")

# Highlight the mismatches between reference and Sanger sequences
mismatches <- vmatchPattern(aligned_seq, ref_seq, max.mismatch = 1)
mismatch_sites <- matrix(as.integer(attributes(mismatches)$mismatch.loc), ncol = 2, byrow = TRUE)
mismatch_seqs <- cbind(as.character(WT1_seq), mismatch_sites)
