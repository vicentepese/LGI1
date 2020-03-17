#!/bin/bash 

cd Data
plink --bfile filtLGI1 --keep ../Resources/patList.txt --make-bed --out temp
plink --bfile temp --pca 20 --out postLGI1
mv postLGI1.eigen* ../Resources/
rm -r postLGI1*
rm -r temp*