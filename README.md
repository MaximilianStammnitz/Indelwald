<b>############ INDELWALD - HYBRID INDEL CALLING ##############</b>

<p>## Version 1.0 ##</p>
<p>## Developped as part of a Computational Biology MPhil Internship ##</p>
<p>## University of Cambridge, UK ##</p>
<p>## DAMTP / DVMS / WTSI ##</p>


<b>############ 14-08-2015: UPDATE - INDELWALD 1.0 ONLINE (PART I) ##############</b>

The toy example reaches daylight.

This involves all functions required to stringently filter and analyse
indel outputs from Pindel and Platypus. A RData-binary file was built, encompassing 3 anonymised
host and tumour VCF-inputs of Tasmanian Devil Chromosome 5. Just open “1_Indelwald_main.R”, and follow the step-by-step protocol. Eventually, you should have the following output options:
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
