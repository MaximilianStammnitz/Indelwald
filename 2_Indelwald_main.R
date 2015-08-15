############ INDELWALD - HYBRID INDEL CALLING ##############
####################### Version 1.0 ########################

## CROSS-REFERENCING
## Last Update - 15/08/2015 ##
## mrs72 / Maximilian Stammnitz ##


## 1. Setting Up Work-Environment ## 
## ## ## ## ## ## ## ## ## ## ## ###

# Install and Load Packages
library(stringr)
library(GenomicRanges)
library(VennDiagram)
library(gridExtra)
library(cluster)

# Source Functions and Helper Files (~ 5 seconds)
source("2_Indelwald_functions.R")

# 6 filtered (Toy 1) VCFs, Devil Chromosome 5
# 6 unfiltered VCFs of Platypus and Pindel, Devil Chromosome 5
load(paste0(main.path,"/Toy/Toy2_Crossreference.Rdata"))


## 2. Build-Up Cross-Reference Tables  ## 
## ## ## ## ## ## ## ## ## ## ## ## ## ##

# Merge the Lists to one
prepared.sets <- cross.ref.prepare(filtered_high = filtered_200.both,
                                   platypus = raw_transcripts.platypus,
                                   pindel = raw_transcripts.pindel)

# Do the Cross-Referencing
cross.ref <- build.cross.ref(input = prepared.sets)


## 3. Shrink Set to Transcript Hits ##
## ## ## ## ## ## ## ## ## ## ## ## ##

cross.transcripts <- transcript.hits.cross(hits = cross.ref,
                                           devil.transcripts,
                                           devil.genes.translation,
                                           devil.transcripts.translation)

# Visualise one example
View(cross.transcripts$'Tumour1')


## 4. Visualise Overlap in Exons of Chromosome 5 ##
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ###

# Define a Lower Threshold for Calling (default: 4)
# Define "Type" of Overlaps - below: all options chosen once
# Output: Venn-Diagram & Hierarchical Clustering and Heatmap between Tumours/Germline

# Removing different Indels (e.g. by length) with shared Positions
cross.DFTD.venn(hits = cross.transcripts,
                threshold = 4,
                type = "cleared",
                titulo = "DFTD - Overlap (position-cleared Indels)")

# Unifying different Indels at the same Position
cross.DFTD.venn(hits = cross.transcripts,
                threshold = 4,
                type = "unified",
                titulo = "DFTD - Overlap (position-unified Indels)")

# All unique Indels
cross.DFTD.venn(hits = cross.transcripts,
                threshold = 4,
                type = "all",
                titulo = "DFTD - Overlap (all unique Indels)")