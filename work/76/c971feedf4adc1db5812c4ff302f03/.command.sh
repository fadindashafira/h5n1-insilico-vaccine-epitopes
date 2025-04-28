#!/bin/bash -ue
Rscript /Users/putriramadani/Documents/GitHub/h5n1-insilico-vaccine-epitopes/bin/combine_epitopes.R \
    --bcell=neuraminidase_neuraminidase_BAL61230.1_bcell_epitopes.csv \
    --tcell-i=neuraminidase_neuraminidase_BAL61230.1_tcell_i_epitopes.csv \
    --tcell-ii=neuraminidase_neuraminidase_BAL61230.1_tcell_ii_epitopes.csv \
    --protein-type=neuraminidase \
    --output=neuraminidase_combined_epitopes.csv
