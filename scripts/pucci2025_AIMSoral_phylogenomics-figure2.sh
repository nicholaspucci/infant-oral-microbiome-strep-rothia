#! /bin/bash

for name in `awk '{print $1}' genus_phylo_analysis.txt`

do

#Run anvi'o (v8.0) module to create a phylogenomics tree based on concatenated amino acid sequences of single-copy core genes from 33 and 63 Streptococcus and Rothia 
#genomes/MAGs, respectively

anvi-gen-phylogenomic-tree -f concatenated_AAseqs-$name.fa \
                           -o phylogenomic_tree-$name.txt	

done

#Generated tree files where further customized in iTOL: Interactive Tree Of Life (Letunic and Bork (2024) Nucleic Acids Res doi: 10.1093/nar/gkae268)) and Adobe Illustrator (v28.7.5, ©1984-2024 Adobe)
