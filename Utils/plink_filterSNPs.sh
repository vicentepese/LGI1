#!/bin/bash

mv ../Resources/notCmnrsIDs.txt ../Data/GWAS/Controls/
cd ../Data/GWAS/Controls
plink --file LGI1_controls --no-fid --no-sex --no-parents --exclude notComnrsIDs.txt --make-bed --noweb --out tst 
plink --bfile tst --recode --tab --out pedtst
rm notComnrsIDs.txt

plink --file LGI1_cases --merge LGI1_controls.ped LGI1_controls.map --recode --out merge