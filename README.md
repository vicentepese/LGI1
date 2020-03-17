# LGI1
Genome-wide association study of paraneoplastic syndrom patients associated with LGI1. This repository contains the basic pipeline for GWAS analysis of data. The following steps were taking, using non-binary GWAS files as inputs to the pipeline. 

- Binarization of files 
- Parsing of common SNPs within controls and cases 
- Compute intersection between common SNPs across controls and cases
- Filter individual files based on common SNPs 
- Merge files into a master file 
- Inheritance by descendence computationn
    - Create list of patients with high PI score for future exclusion
- Quality control computation
    - Filter out patients and samples with genome or sample missingness lower than 10%
    - Minimum Allele Frequency threshold of 5%
    - Filter out MHC area to avoid over-representation
- Add phenotype based on clinical data 
- Compute PCA and plot PC1 vs PC2 for batch effect visual analysis
- Re-filter data:
    - Filter out patients and samples with genome or sample missingness lower than 10%
    - Minimum Allele Frequency threshold of 5%
- Patient matching based on euclidian distance
- Post-matching PCA and PC1 vs PC2 plot 
- Manual removal of non-matching subjects
- Logistic regression computation for association analysis with matched subjects and exclusin of patients based on IBD computation
- Manhattan plot for visual examination