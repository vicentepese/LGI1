#!/bin/bash

cd Data
# plink --bfile filtLGI1 --keep ../Resources/filtPostMatchPatList.txt --linear --covar ../Resources/LGI1.eigenvec --allow-no-sex --out tst
cp modLGI1.fam filtLGI1.fam
plink --bfile filtLGI1 --keep ../Resources/filtPostMatchPatList.txt \
    --covar ../Resources/PCAtst.eigenvec --covar-name PC1,PC2,PC3,PC4 \
    --logistic --allow-no-sex --out filtLGI1

