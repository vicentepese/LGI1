#!/bin/bash
cd Data/GWAS/

king -b LGI1.bed --ibd
mv king.seg ../../Resources/LGI1.genome
rm -r *king*
