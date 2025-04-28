#!/bin/bash -ue
python /Users/putriramadani/Documents/GitHub/h5n1-insilico-vaccine-epitopes/bin/predict_tcell_ii.py \
    --fasta=neuraminidase_BAL61230.1.fasta \
    --protein-type=neuraminidase \
    --method=NetMHCIIpan \
    --threshold=500 \
    --length=15 \
    --alleles='HLA-DRB1*03:01,HLA-DRB1*04:01,HLA-DRB1*01:01,HLA-DRB1*07:01,HLA-DRB1*15:01,HLA-DRB1*13:01,HLA-DRB1*11:01' \
    --output=neuraminidase_neuraminidase_BAL61230.1_tcell_ii_epitopes.csv
