#!/bin/bash -ue
Rscript /Users/putriramadani/Documents/GitHub/measles-rubella-vaccinedesign/bin/combine_epitopes.R \
    --bcell=accession2_accession2_NP_740664.1_bcell_epitopes.csv \
    --tcell-i=accession2_accession2_NP_740664.1_tcell_i_epitopes.csv \
    --tcell-ii=accession2_accession2_NP_740664.1_tcell_ii_epitopes.csv \
    --protein-type=accession2 \
    --output=accession2_combined_epitopes.csv
