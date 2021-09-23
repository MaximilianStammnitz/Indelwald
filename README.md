![Scripting](https://img.shields.io/badge/Language-R-yellow.svg) ![Copyright](https://img.shields.io/badge/Copyright-(c)_2021_Max\_Stammnitz\_@TCG\_Cambridge-green.svg)

Indelwald - Scripts for processing Indel calls
==============================================
![Indelwald](/Images/Indelwald_logo.png)

<b>############ 08-07-2021: UPDATE ##############</b>

Six years after its opening, I have added a [new piece of R code](/Scripts/3_Indelwald_spectrum.R) to the Indelwald repository. Using this, we generate indel variant spectra in line with the complex PCAWG classification scheme specified by [Alexandrov et al., 2020](https://www.nature.com/articles/s41586-020-1943-3) (see [COSMIC ID signature catalogues](https://cancer.sanger.ac.uk/signatures/id/)).

Note that, in order to produce your own indel spectra, you will need to specify:
* the reference genome fasta file based on which your alignments' indel calls were generated
* a VCF or dataframe object containing one indel per line, featuring the minimal set of columns "CHROM", "POS", "REF", "ALT"

My code then groups all of your short (< 80 bp) insertion and deletion variants based on the (current) 83 different indel types agreed upon by the PCAWG signature consortium – these reflect a consensus rule set [regarding variant lengths, sequence type and immediate sequence context](https://cancer.sanger.ac.uk/signatures/documents/4/PCAWG7_indel_classification_2017_12_08.xlsx). Plotted spectra look like this example:

![example](/Images/Example_spectrum.png)

In the above case, we see enrichments of single-T deletions or extensions at poly-T homopolymers (lengths ≥5 bp). These spikes are indicative of DNA polymerase slippage, which is particularly prominent in tissues with DNA mismatch repair (MMR) deficiency – commonly classified as COSMIC signatures ID1 and ID2.

Disclaimer #1: The script and its associated results are yet unpublished. It has nevertheless been extensively validated with indel calls derived from (i) human, (ii) dog and (iii) Tasmanian devil (cancer) genomes. In theory, this code should run smoothly for ANY species with a reference genome. If you do face a challenge in using the code or wish to provide general feedback, please get in touch directly via maxrupsta@gmail.com

Disclaimer #2: Color scheme is the same as the one used in the original PCAWG paper. The rest was written from scratch.

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

<b>############ 14-08-2015 ##############</b>

<p>## Developed as part of a Computational Biology MPhil Internship in 2015 ##</p>
<p>## University of Cambridge, UK ##</p>
<p>## DAMTP / DVetMed / WTSI ##</p>
