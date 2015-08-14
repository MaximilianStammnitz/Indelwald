<b>############ INDELWALD - HYBRID INDEL CALLING ##############</b>

<p>## Version 1.0 ##</p>
<p>## Developped as part of a Computational Biology MPhil Internship ##</p>
<p>## by Maximilian R. Stammnitz ##</p>
<p>## University of Cambridge, UK ##</p>


<b>############ 14-08-2015: UPDATE - INDELWALD 1.0 PART I ONLINE ##############</b>

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


<b>############ 14-08-2015: UPDATE - FILE UPLOAD IN VERY SHORT. ##############</b>

Portions of the raw input of our toy example, and associated file/function names of Indelwald, comprise strictly confidential informations on a recently submitted article on the origin of these samples (Pye et al.). Our group was hoping until the last that this paper would be released before the deadline of the Indelwald-associated thesis ('Indels in Tasmanian Devils'), but unfortunately it is ongoingly reviewed. In order to share the code nevertheless, a careful renaming (and some revision) has been undertaken in the last two days. Non-conflicting code will be released later today, on <b> Friday, 14th August</b>.
