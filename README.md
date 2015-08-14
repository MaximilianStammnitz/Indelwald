<b>############ 14-08-2015: UPDATE - INDELWALD 1.0 PART I ONLINE ##############</b>

Finally, the toy example reaches daylight.

This involves all functions required to stringently filter and analyse
indel outputs from Pindel and Platypus. A RData-binary file was built, encompassing 3 anonymised
host and tumour VCF-inputs of Tasmanian Devil Chromosome 5. Just open “1_Indelwald_main.R”, and follow the step-by-step protocol. Eventually, you should have the following output options:
- Platypus/Pindel Overlap Venn-Diagram
- Detected Indel Size Plot
- Indel-Density Plot
- Rainfall Plot
- BAF Plot
- Summary Statistics Plot

Note, that 5 R-packages need to be installed on your machine to get all
functions running properly: VennDiagram, gridExtra, stringr, data.table,
and GenomicRanges.

Enjoy Indelwald!


<b>############ 14-08-2015: UPDATE - FILE UPLOAD IN VERY SHORT. ##############</b>

Portions of the raw input of our toy example, and associated file/function names of Indelwald, comprise strictly confidential informations on a recently submitted article on the origin of these samples (Pye et al.). Our group was hoping until the last that this paper would be released before the deadline of the Indelwald-associated thesis ('Indels in Tasmanian Devils'), but unfortunately it is ongoingly reviewed. In order to share the code nevertheless, a careful renaming (and some revision) has been undertaken in the last two days. Non-conflicting code will be released later today, on <b> Friday, 14th August</b>.


<b>############ INDELWALD - HYBRID INDEL CALLING ##############</b>

<p>## Version 1.0 ##</p>
<p>## Developped as part of a Computational Biology MPhil Internship ##</p>
<p>## by Maximilian R. Stammnitz ##</p>
<p>## University of Cambridge, UK ##</p>

<p>This small framework consists of a few components for analysing DFTD indels:</p>
1. An R-file containing all functions for post-processing and displaying Pindel and Platypus calls
2. An R-file example, guiding through the main options of output lists and plots
3. A small Platypus output example VCF
4. A small Pindel output example VCF
5. A table with > 36,000 sorted supercontigs of the Devil Reference Genome 7.1, their respective sizes and chromosomal coordinates
6. A table with all biomaRt-translated ENSEMBL gene and transcript IDs
7. A table comprising all annotated devil exons and introns, chromosomal coordinate, strand orientation, ENSEMBL gene and transcript ID
8. A CSV-table with all annotated COSMIC cancer gene consensus genes and their characteristics
