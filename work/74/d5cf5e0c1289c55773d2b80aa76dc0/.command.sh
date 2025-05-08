#!/bin/bash -ue
Rscript /Users/putriramadani/Documents/GitHub/measles-rubella-vaccinedesign/bin/combine_epitopes.R \
    --bcell=accession1_accession1_NP_054712.1_bcell_epitopes.csv \
    --tcell-i=accession1_accession1_NP_054712.1_tcell_i_epitopes.csv \
    --tcell-ii=accession1_accession1_NP_054712.1_tcell_ii_epitopes.csv \
    --protein-type=accession1 \
    --output=accession1_combined_epitopes.csv
