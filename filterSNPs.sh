#!/bin/bash

########## CONTROLS ############

# Copy resources
cp Resources/commonSNPs.txt Data/GWAS/Controls/Binaries/commonSNPs.txt

# Go to binaries directory
cd Data/GWAS/Controls/Binaries/

files=($(ls))
for file in ${files[@]} ; do 
    if [[ $file == *".bed"* ]] ; then 
        echo "Filtering SNPs in " $file
        IFS='.' read -a strarr <<< "$file"
        plink --bfile ${strarr[0]} --no-sex --no-pheno --no-fid --no-parents --extract commonSNPs.txt --make-bed --out ../FiltBinaries/${strarr[0]}
        mv ${strarr[0]}.log ../FiltBinaries/Log/${strarr[0]}.log
    fi 
done
rm commonSNPs.txt

cd ../FiltBinaries/
files=($(ls))
for file in ${files[@]} ; do 
    if [[ $file == *".bim"* ]] ; then 
        echo $file
        awk 'END{print NR}' $file
    fi 
done 

########### CASES ##############

cd ../../../../
cp Resources/commonSNPs.txt Data/GWAS/Cases/Binaries/commonSNPs.txt

# Go to binaries directory
cd Data/GWAS/Cases/Binaries/

files=($(ls))
for file in ${files[@]} ; do 
    if [[ $file == *".bed"* ]] ; then 
        echo "Filtering SNPs in " $file
        IFS='.' read -a strarr <<< "$file"
        plink --bfile ${strarr[0]} --no-sex --no-pheno --no-fid --no-parents --extract commonSNPs.txt --make-bed --out ../FiltBinaries/${strarr[0]}
        mv ${strarr[0]}.log ../FiltBinaries/Log/${strarr[0]}.log
    fi 
done
rm commonSNPs.txt

cd ../FiltBinaries/
files=($(ls))
for file in ${files[@]} ; do 
    if [[ $file == *".bim"* ]] ; then 
        echo $file
        awk 'END{print NR}' $file
    fi 
done 