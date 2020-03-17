#!/bin/bash 

########## CONTROLS ###########

# - No missalignment in controls. No "copys". This section takes each .ped and .map file and 
#       converst it to a binary that is then moved to the Binaries folder
# - Logs are moved to Logs directory

cd Data/GWAS/Controls/OG
files=($(ls))
for file in ${files[@]} ; do
    if [[ $file == *'.ped'* ]] ; then
        echo "Converting $file to binary"
        IFS='.' read -a strarr <<< "$file"
        plink --file ${strarr[0]} --no-sex --no-pheno --no-fid --no-parents --noweb --make-bed --out ../Binaries/${strarr[0]} >> ../Binaries/${strarr[0]}.log
        mv ../Binaries/${strarr[0]}.log ../Binaries/Log/${strarr[0]}.log
    fi  
done 


########## CASES ###########

# There is are copys of cases in the .ped file. 
# - .ped files are filtered and copies are rejected. 
# - Files are converted to binaries, sent to the ../Binaries directory and .log to the ../Log directory

cd ../../Cases/OG

# Filter data
files=$ls
if [[$files != *'_filt.ped'*]] ; then 
    files=($(ls))
    for file in ${files[@]} ; do 
        if [[ $file == *'.ped'* ]] ; then 
            echo "Removing duplicate patients of " $file 
            IFS='.' read -a strarr <<< "$file"  
            awk '$3 !~ /Copy/ { print $0 ; }' $file > ${strarr[0]}"_filt.ped"
            cp ${strarr[0]}.map ${strarr[0]}_filt.map
            echo "Duplicates removed"
        fi 
    done
else 
    echo "Cases already filtered"
fi 

files=($(ls))
for file in ${files[@]} ; do
    echo $file 
    if [[ $file == *'.ped'* ]] && [[ $file == *"filt"* ]] ; then
        echo "Converting $file to binary"
        IFS='.' read -a strarr <<< "$file"
        plink --file ${strarr[0]} --no-sex --no-pheno --no-fid --no-parents --noweb --make-bed --out ../Binaries/${strarr[0]} >> ../Binaries/${strarr[0]}.log
        mv ../Binaries/${strarr[0]}.log ../Binaries/Log/${strarr[0]}.log
    fi 
done 

