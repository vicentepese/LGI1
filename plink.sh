#!/bin/bash

########### CONTROLS #############
cd Data/GWAS/Controls
controls="LGI1_controls" 

# Convert to binary 
files=$(ls)
if [[ $files == *"LGI1_controls.ped"* ]] 
then 
    if [[ $files != *"BLGI1_controls.bed"* ]]
    then 
        echo "Creating duplicate vars"
        plink --file $controls --no-sex --no-parents --no-fid --list-duplicate-vars --noweb --out BLGI1_controls
        echo "Converting LGI1_controls.ped to binary"
        plink --file $controls --no-fid --no-sex --no-parents --noweb --make-bed --out BLGI1_controls
        mv BLGI1_controls.log Log/BLGI1_controls.log
        echo "Conversion to binary was successful"
        plink --bfile 
    else
        echo "LGI1_controls already convervet to binary"
    fi 
else
    echo "LGI1_controls.ped not found. Please run getControls.py to create the .ped file"
fi 

# Remove duplicates
echo "Removing duplicates of LGI1_controls"
plink --file LGI1_controls --no-sex --no-parents --no-fid --list-duplicate-vars suppress-first --noweb --out BLGI1_controls
mv BLGI1_controls.dupvar Log/BLGI1_controls.dupvar
echo "Duplicates removed"



# Perform quality controls 
files=$(ls)
log=$(ls)
if [[ $files == *"BLGI1_controls.bed"* ]] && [[ $log =~ *"BLGI1_controls.fmendel"* ]]
then 
    echo "Performing quality control"
    echo "Missing genotypes"
    plink --bfile BLGI1_controls --no-fid --no-sex --no-parents --noweb --missing --out BLGI1_controls
    echo "Missingness of case/control status"
    plink --bfile BLGI1_controls --no-fid --no-sex --no-parents --noweb --test-missing --out BLGI1_controls
    echo "Haplotype-based test for non-random missing genotype data"
    plink --bfile BLGI1_controls --no-fid --no-sex --no-parents --noweb --test-mishap  --out BLGI1_controls
    echo "Hardy-Weinberg Equilibrium"
    plink --bfile BLGI1_controls --no-fid --no-sex --no-parents --noweb --hardy  --out BLGI1_controls
    echo "Allele frequency"
    plink --bfile BLGI1_controls --no-fid --no-sex --no-parents --noweb --freq  --out BLGI1_controls
    echo "Mendel error"
    plink --bfile BLGI1_controls --no-fid --no-sex --no-parents --noweb --mendel  --out BLGI1_controls

else
    echo "File not converted to binary. Please convert to binary before performing quality control."
fi 

########### CASES #############
cd ~/Documents/LGI1
cd Data/GWAS/Cases
cases="LGI1_cases" 

# Convert to binary 
files=$(ls)
if [[ $files == *"LGI1_cases.ped"* ]] 
then 
    if [[ $files != *"BLGI1_cases.bed"* ]]
    then 
        echo "Converting LGI1_cases.ped to binary"
        plink --file $cases --no-fid --no-sex --no-parents --noweb --make-bed --out BLGI1_cases
        echo ""
        mv BLGI1_cases.log Log/BLGI1_cases.log
        echo "Conversion to binary was successful"
    else
        echo "LGI1_cases already convervet to binary"
    fi 
else
    echo "LGI1_cases.ped not found. Please run getCases.py to create the .ped file"
fi 

# Perform quality controls 
files=$(ls)
log=$(ls)
if [[ $files == *"BLGI1_cases.bed"* ]]
then 
    echo "Performing quality control"
    echo "Missing genotypes"
    plink --bfile BLGI1_cases --no-fid --no-sex --no-parents --noweb --missing --out Log/BLGI1_cases
    echo "Missingness of case/control status"
    plink --bfile BLGI1_cases --no-fid --no-sex --no-parents --noweb --test-missing --out Log/BLGI1_cases
    echo "Haplotype-based test for non-random missing genotype data"
    plink --bfile BLGI1_cases --no-fid --no-sex --no-parents --noweb --test-mishap  --out Log/BLGI1_cases
    echo "Hardy-Weinberg Equilibrium"
    plink --bfile BLGI1_cases --no-fid --no-sex --no-parents --noweb --hardy  --out Log/BLGI1_cases
    echo "Allele frequency"
    plink --bfile BLGI1_cases --no-fid --no-sex --no-parents --noweb --freq  --out Log/BLGI1_cases
    echo "Mendel error"
    plink --bfile BLGI1_cases --no-fid --no-sex --no-parents --noweb --mendel  --out Log/BLGI1_cases

else
    echo "File not converted to binary. Please convert to binary before performing quality control."
fi 

# plink --file $controls --no-fid --no-sex --no-parents --noweb --make-bed --out BLGI1_controls

# Missing genotypes
#plink --bfile BLGI1_controls --no-fid --no-sex --no-parents --noweb --missing

# Allele frequency statistics
# plink --bfile BForward --freq --out freq_stat --noweb

#Sort MAF
# sort --key=5 -nr freq_stat.frq > sortFreq_stat.frq
# awk '$5 < 0.05 && $5 > 0 {print $0}' sortFreq_stat.frq > filtFreq_stat.frq 

# Missing statistics
# plink --bfile BForward --missing --out miss_stat --noweb

#  Association analysis 
# plink --bfile BForward --assoc --out ../../Output/ForwardAssoc --noweb

plink --file BLGI1_cases --no-sex --no-fid --no-parents --merge LGI1_controls.ped LGI1_controls.map --make-bed tst 
plink --file LGI1_cases --no-sex --no-fid --no-parents --make-bed --out BLGI1_cases
plink --file LGI1_controls --no-sex --no-fid --no-parents --make-bed --out BLGI1_controls
plink --bfile BLGI1_cases --bmerge BLGI1_controls.bed BLGI1_controls.bim BLGI1_controls.fam --make-bed --out merge