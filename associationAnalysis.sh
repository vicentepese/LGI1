#!/bin/bash

cd Data
# plink --bfile filtLGI1 --keep ../Resources/filtPostMatchPatList.txt --linear --covar ../Resources/LGI1.eigenvec --allow-no-sex --out tst
cp modLGI1.fam filtLGI1.fam
plink --bfile filtLGI1 --keep ../Resources/filtPostMatchPatList.txt \
    --covar ../Resources/covariates.txt --covar-name PC1, PC2, PC3, PC4, HLAA_allele1, HLAA_allele2, HLAB_allele1, HLAB_allele2, HLAC_allele1, HLAC_allele2 \
    --logistic --allow-no-sex --out filtLGI1
mv filtLGI1.assoc.logistic ../Resources/filtLGI1.assoc.logistic
