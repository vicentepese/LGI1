#!/bin/bash

# Go to directory 
cd Data

# Count number of subjects 
nsub=$(awk 'END{print NR}' filtLGI1.fam)

# Compute PCA
mv modLGI1.fam filtLGI1.fam
plink --bfile filtLGI1 --pca 20 --out LGI1

# Move PCs to Resources
mv LGI1.eigen* ../Resources/
rm LGI1.*
