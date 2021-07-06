<b>####### INDELWALD - SCRIPTS FOR PROCESSING INDEL-CALLS ########</b>

<p>## First developed as part of a Computational Biology MPhil Internship in 2015 ##</p>
<p>## University of Cambridge, UK ##</p>
<p>## DAMTP / DVetMed / WTSI ##</p>

![Indelwald](./Indelwald_logo.png)


<b>############ 06-07-2021: UPDATE ##############</b>

Six years after its opening, I've finally added a new piece of R code to this repository. Using this, we generate indel variant spectra from input VCF files, in line with the PCAWG indel classifications specified in [Alexandrov et al., 2020](https://www.nature.com/articles/s41586-020-1943-3) ([COSMIC ID signature catalogues](https://cancer.sanger.ac.uk/signatures/id/)).

Note that, in order to produce your own indel spectra, you will need to provide the reference genome fasta file based on which your indel genotypes were called in the first place. My code then groups all insertion and deletion variants based on the (current, as of July 2021) 83 different different consensus types agreed upon by the PCAWG signature consortium, [based both on their lengths, sequence type and immediate sequence context](https://cancer.sanger.ac.uk/signatures/documents/4/PCAWG7_indel_classification_2017_12_08.xlsx). 

Disclaimer #1: This script and its associated results are yet unpublished. We won't pay your insurance on potentially false outputs. 
Disclaimer #2: I have stolen the color scheme from Alexandrov et al’s original [SigProfiler plotting](https://github.com/AlexandrovLab/SigProfilerPlotting/blob/master/sigProfilerPlotting/sigProfilerPlotting.py).

<b>############ 15-08-2015: UPDATE ##############</b>

The second Toy Example is now online, encompassing the crossreferencing steps.

Filtered and non-filtered sets of three anonymised hosts and tumours are provided as a binary RData-file: “Toy2_Crossreference.Rdata”. Users can now browse through all traits of Indelwald 1.0, in two step-by-step tutorials. These guide through plotting and data summary functions: 
- 1. Stringent Indel filtering (1_Indelwald_main.R)
- 2. Cross-referencing (2_Indelwald_main.R)

New plots include:
- Venn-Diagrams between the panel of hosts and each tumour
- Hierarchical clustering and heatmap for all samples, based on shared/unique indels present in exons.

<b>############ 14-08-2015 ##############</b>

The first toy example reaches daylight.

It involves all functions required to stringently filter and analyse
indel outputs from Pindel and Platypus. An RData-binary file was created, respectively encompassing 1 anonymised Pindel/Platypus host and tumour VCF-input - of Tasmanian devil chromosome 5. Just open “1_Indelwald_main.R”, and follow the step-by-step manual. Eventually, you should have the following output options:

- Platypus/Pindel Overlap Venn-Diagram
- Detected Indel Size Plot
- Indel-Density Plot
- Different Rainfall Plots
- BAF Plot
- Summary Statistics Plot

Note, that 5 R-packages need to be installed on your machine to get all functions running properly: VennDiagram, gridExtra, stringr, data.table, and GenomicRanges. The total zip-file comprises 20MB.

Enjoy Indelwald!
