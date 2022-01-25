![Scripting](https://img.shields.io/badge/Language-R-yellow.svg) ![Copyright](https://img.shields.io/badge/Copyright-(c)_2021_Max\_Stammnitz\_@TCG\_Cambridge-green.svg)

Indelwald – from indel calls to indel spectra
=============================================
![Indelwald](/Images/Indelwald_logo.png)


Using the functions in [Indelwald](/Indelwald.R), we generate short insertion and deletion – also termed indels (ID) – spectra in line with the complex PCAWG classification scheme first specified by [Alexandrov et al., 2020](https://www.nature.com/articles/s41586-020-1943-3) (see [COSMIC ID signature catalogues](https://cancer.sanger.ac.uk/signatures/id/)).

Note that, in order to produce your own indel spectra, you will need to specify:
* the reference genome fasta file based on which your alignments' indel calls were generated
* a VCF or dataframe object containing one indel per line, featuring the minimal set of columns "CHROM", "POS", "REF", "ALT"

My code then groups all of your short (< 80 bp) insertion and deletion variants based on the (current) 83 different indel types agreed upon by the PCAWG signature consortium – these reflect a consensus rule set [regarding variant lengths, sequence type and immediate sequence context](https://cancer.sanger.ac.uk/signatures/documents/4/PCAWG7_indel_classification_2017_12_08.xlsx). Plotted spectra look like this example:

![example](/Images/Example_spectrum.png)

In the above case, we see enrichments of single-T deletions or extensions at poly-T homopolymers (lengths ≥5 bp). These spikes are indicative of DNA polymerase slippage, which is particularly prominent in tissues with DNA mismatch repair (MMR) deficiency – commonly classified as COSMIC signatures ID1 and ID2.

The script and its associated results are yet unpublished/peer-reviewed. It has nevertheless already been extensively benchmarked with indel calls against the following reference genomes (special thanks to Adrian Baez-Ortega, Wellcome Sanger Institute): Human, Black-and-white colobus, Cat, Cow, Dog, Ferret, Giraffe, Harbour porpoise, Horse, Lion, Mouse, Naked mole-rat, Rabbit, Rat, Ringle-tailed lemur, Tiger and – last but not least – Tasmanian devil. In theory, this code should run smoothly for ANY species with a reference genome. 

Find it helpful or you simply enjoy making Indelwald spectra against your brand-new, awesome reference genome? Then why not surf the wave with your substitution calls? – have a look at its sister library [SubstitutionSafari](https://github.com/MaximilianStammnitz/SubstitutionSafari)! If you do face a challenge in using the code or wish to provide general feedback, please get in touch directly via maxrupsta@gmail.com
