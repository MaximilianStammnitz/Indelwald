############ INDELWALD - HYBRID INDEL CALLING ##############
####################### Version 1.0 ########################

## STRINGENT FILTERING
## Last Update - 14/08/2015 ##
## mrs72 / Maximilian Stammnitz ##


## 1. Setting Up Work-Environment ## 
## ## ## ## ## ## ## ## ## ## ## ###

# Install and Load Packages
library(VennDiagram)
library(gridExtra)
library(stringr)
library(data.table)
library(GenomicRanges)

# Source Functions and Helper Files (~ 5 seconds)
source("1_Indelwald_functions.R")

# Load the Toy Dataset
# 1 Tumour-VCF of Platypus & Pindel, for Devil Chromosome 5 (Example)
# 1 Host-VCF of Platypus & Pindel, for Devil Chromosome 5
load(paste0(main.path,"/Toy/Toy1_Stringent.Rdata"))


## 2. Checking Caller's Raw Output  ## 
## ## ## ## ## ## ## ## ## ## ## ## ##

## Venn Overlap between Pindel and Platypus (~ 5 seconds)
raw.calls(pindel.vcf = tumour.pindel,
          platypus.vcf = tumour.platypus,
          chromosome = "5",
          sample = "Tumour_#1")

## Size Spectrum of Pindel and Platypus
raw.sizes(pindel.vcf = tumour.pindel,
          platypus.vcf = tumour.platypus,
          chromosome = "5",
          sample = "Tumour_#1")


## 3. Stringent Filtering  ## 
## ## ## ## ## ## ## ## ## ##

tumour.1 <- vector(mode="list", length=5)
names(tumour.1) <- c("overlap", "pindel-Q", 
                     "platypus-Q", "summary", "in-transcript")

## Take Platypus/Pindel-Intersection
tumour.1[[1]] <- overlap.vcfs(pindel.vcf = tumour.pindel,
                              platypus.vcf = tumour.platypus)

# Apply a Pindel-Quality Filter (default Threshold: 200)
tumour.1[[2]] <- pindel.filter(overlap = tumour.1$'overlap',
                               threshold = 200)

# Apply a Platypus-Quality Filter
tumour.1[[3]] <- platypus.filter(overlap = tumour.1$'pindel-Q',
                                 threshold = 200)

# Summarise: Match Contig-Positions to Chromosomal Positions
tumour.1[[4]] <- indel.summary(contigs.coordinates = contig.lengths,
                               pindel.vcf = tumour.1$'platypus-Q'[[1]],
                               platypus.vcf = tumour.1$'platypus-Q'[[2]])

# Have a look at the Summary-Table: Indel Positions, Types, BAFs
View(tumour.1$'summary')


## 4. Visualising Indels   ## 
## ## ## ## ## ## ## ## ## ##

# Indel-Density Plot on Chromosome
indel.density(hits = tumour.1$'summary',
              contigs.coordinates = contig.lengths,
              genes = devil.genes,
              chromosome = "5",
              sample = "Tumour_#1")

# Rainfall Plot
rainfall(hits = tumour.1$'summary',
         contigs.coordinates = contig.lengths,
         chromosome = "5",
         sample = "Tumour_#1",
         split = "Y")
                
# Rainfall Plot with Size-Type Specifications
rainfall.types(hits = tumour.1$'summary',
               chromosome = "5",
               sample = "Tumour_#1")

# Indel BAF Plot (default: with Platypus Calls only)
# Pindel BAFs are - due to the algorithm's read retrieval - constrained to lower BAF ratios.
# However, they can be visualised by flagging pindel="Y"
indel.baf(hits = tumour.1$'summary',
          contigs.coordinates = contig.lengths,
          chromosome = "5",
          sample = "Tumour_#1",
          pindel = "N")


## 5. Transcript Filtering  ##
## ## ## ## ## ## ## ## ## ###

# Which Indels are found in ENSEMBL-Transcriptlist?
tumour.1[[5]] <- transcript.hits(hits = tumour.1$'summary',
                                 chromosome = "5")

# Translate Transcripts/Genes to HGNC
tumour.1[[5]] <- hits.translate(hits = tumour.1$'in-transcript')

# Visualise these hits with one command
View(tumour.1$'in-transcript'[,c(12,14,1,3:7)])


## 6. Collect Stats  ##
## ## ## ## ## ## ## ##

final.results.plots(pindel.vcf = tumour.pindel,
                    platypus.vcf = tumour.platypus,
                    summary.list = tumour.1,
                    chromosome = "5",
                    sample = "Tumour_#1")