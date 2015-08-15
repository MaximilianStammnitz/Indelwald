<b>############ INDELWALD - HYBRID INDEL CALLING ##############</b>

<p>## Version 1.0 ##</p>
<p>## Developped as part of a Computational Biology MPhil Internship ##</p>
<p>## University of Cambridge, UK ##</p>
<p>## DAMTP / DVMS / WTSI ##</p>
<b>## Folder: Info-Files ##</b>


<b>############ 15-08-2015: UPDATE - INDELWALD 1.0 ONLINE (PART II) ##############</b>

The second Toy Example is now only, encompassing the crossreferencing steps.

Filtered and non-filtered sets of three anonymised hosts and tumours are provided as a binary RData-file: “Toy2_Crossreference.Rdata”. Users can now browse through all traits of Indelwald 1.0, in two step-by-step tutorials. These guide through plotting and data summary functions: 

- 1. Stringent Indel filtering (1_Indelwald_main.R)
- 2. Cross-referencing (2_Indelwald_main.R)

New plots include Venn-Diagrams between the panel of hosts and each tumour, as well as a hierarchical clustering and heatmap for all samples, based on shared/unique indels present in exons.

<b>############ 14-08-2015: UPDATE - INDELWALD 1.0 ONLINE (PART I) ##############</b>

The first toy example reaches daylight.

This involves all functions required to stringently filter and analyse
indel outputs from Pindel and Platypus. A RData-binary file was built, encompassing 1 anonymised
host and tumour VCF-input of Tasmanian Devil Chromosome 5. Just open “1_Indelwald_main.R”, and follow the step-by-step protocol. Eventually, you should have the following output options:
- Platypus/Pindel Overlap Venn-Diagram
- Detected Indel Size Plot
- Indel-Density Plot
- Different Rainfall Plots
- BAF Plot
- Summary Statistics Plot

Note, that 5 R-packages need to be installed on your machine to get all
functions running properly: VennDiagram, gridExtra, stringr, data.table,
and GenomicRanges.

Enjoy Indelwald!
