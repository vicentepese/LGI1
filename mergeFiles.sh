#!/bin/bash 

# Go to GWAS directory 
cd Data/GWAS

# Create temporary folder
mkdir tmp

# Copy cases
cp Controls/FiltBinaries/* tmp

# Copy controls 
cp Cases/FiltBinaries/* Cases/
cd Cases
files=($(ls *Plates*))
for file in ${files[@]} ; do
    IFS='.' read -a str <<< "$file"
    mv $file LGI1_cases.${str[1]}    
done
cp -r LGI1_cases* ../tmp 
cd ..

# Copy list of files 
cp ../../Resources/mergeList.txt tmp

# Merge 
cd tmp
plink --bfile LGI1_cases --merge-list mergeList.txt --out LGI1
cd ..

# Move file to GWAS
mv tmp/LGI1.* ./
cp LGI1.* ../filtLGI1*

# Delete temporary folder
rm -R tmp
