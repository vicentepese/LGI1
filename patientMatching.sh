#!/bin/bash 

# Compute distance between cases 
cd Data 
plink --bfile filtLGI1 --distance square --out distance
mv *distance* ../Resources/