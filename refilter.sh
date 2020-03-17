#!/bin/bash 

# Go to directory 
cd Data/GWAS 

# Rewrite files with 10% missingness of genotype and sample, and remove MHC region in CHR6 
plink --bfile LGI1 --remove ../../Resources/excludeID.txt \
    --no-fid --no-sex --no-parents --not-chr 25,26 \
    --maf 0.05 --geno 0.1 --mind 0.1 \
    --make-bed --out filtLGI1
mv filtLGI1* ../
cp ../modLGI1.fam ../filtLGI1.fam