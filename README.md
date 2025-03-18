![Scripting](https://img.shields.io/badge/Language-R-yellow.svg) ![Copyright](https://img.shields.io/badge/Copyright-(c)_2021_Max\_Stammnitz\_@TCG\_Cambridge-green.svg)

Indelwald – from indel calls to indel spectra
=============================================

Using the functions in [Indelwald](/Indelwald.R), we generate short insertion and deletion – also termed indels (ID) – spectra in line with the complex PCAWG classification scheme first specified by [Alexandrov et al., 2020](https://www.nature.com/articles/s41586-020-1943-3) (see [COSMIC ID signature catalogues](https://cancer.sanger.ac.uk/signatures/id/)).

Note that, in order to produce your own indel spectra, you will need to specify:
* the reference genome fasta file based on which your alignments' indel calls were generated
* a VCF or dataframe object containing one indel per line, featuring the minimal set of columns "CHROM", "POS", "REF", "ALT"

My code then groups all of your short (< 80 bp) insertion and deletion variants based on the (current) 83 different indel types agreed upon by the PCAWG signature consortium – these reflect a consensus rule set [regarding variant lengths, sequence type and immediate sequence context](https://cancer.sanger.ac.uk/signatures/documents/4/PCAWG7_indel_classification_2021_08_31.xlsx). Plotted spectra look like this example:

![example](/Images/Example_spectrum.png)

In the above case, we see enrichments of single-T deletions or extensions at poly-T homopolymers (lengths ≥5 bp). These spikes are indicative of DNA polymerase slippage, which is particularly prominent in tissues with DNA mismatch repair (MMR) deficiency – commonly classified as COSMIC signatures ID1 and ID2.

We have extensively benchmarked this script with indel calls from the Wellcome Sanger Institute's cross-species mutation rate project ([Cagan, Baez-Ortega _et al._ 2022, Nature](https://www.nature.com/articles/s41586-022-04618-z)), featuring the following species: Human, Black-and-white colobus, Cat, Cow, Dog, Ferret, Giraffe, Harbour porpoise, Horse, Lion, Mouse, Naked mole-rat, Rabbit, Rat, Ringle-tailed lemur and Tiger. In theory, this code should run smoothly for ANY species with a reference genome.

Find it helpful or you simply enjoy making Indelwald spectra against your brand-new, awesome reference genome? Then why not surf the wave with your substitution calls? – have a look at its sister library [SubstitutionSafari](https://github.com/MaximilianStammnitz/SubstitutionSafari)! If you do face a challenge in using the code or wish to provide general feedback, please get in touch directly via maxrupsta@gmail.com

---

## Citation

If you can make good use of these functions in your work, I would be grateful for your citation of our associated preprint: **The evolution of two transmissible cancers in Tasmanian devils ([Stammnitz _et al._ 2023, Science 380:6642](https://doi.org/10.1126/science.abq6453))**
