# AC_phage_analysis
Scripts used in the following manuscript: Stop codon recoding is widespread in diverse phage lineages and has the potential to regulate translation of late stage and lytic genes

doi: https://doi.org/10.1101/2021.08.26.457843


The purpose of the script get_CD.py is to calculate coding density (CD) in standard code, code 4, and code 15. It works by summing the phage gene lengths and dividing by the phage genome length. This script requires that you have already predicted genes using prodigal in single mode with the following genetic codes: code 4 (TGA recoded), code 11 (standard code), and code 15 (TAG recoded). This script will provide a genetic code prediction, which is informed by the coding density in standard and alternative code. See associated manuscript for more details.

Be sure to check genetic code predictions, as many factors other than genetic code choice can influence the coding density of a contig!

Phage genetic code prediction becomes more accurate when the phage genome is complete, or near complete. Best practices are to check genetic code assignments by manual assesment: If the alternative genetic code results in more contiguous operon structure, reduced strand switching, correct-length genes (as checked by blastp against NCBI database), and does not result in gene fusions (as checked by blastp against NCBI database) the phage can be confirmed as a alt-coded. 





