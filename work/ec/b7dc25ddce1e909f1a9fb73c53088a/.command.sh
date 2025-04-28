#!/bin/bash -ue
Rscript /Users/putriramadani/Documents/GitHub/h5n1-insilico-vaccine-epitopes/bin/combine_epitopes.R \
    --bcell=hemagglutinin_hemagglutinin_BAL61222.1_bcell_epitopes.csv \
    --tcell-i=hemagglutinin_hemagglutinin_BAL61222.1_tcell_i_epitopes.csv \
    --tcell-ii=hemagglutinin_hemagglutinin_BAL61222.1_tcell_ii_epitopes.csv \
    --protein-type=hemagglutinin \
    --output=hemagglutinin_combined_epitopes.csv
