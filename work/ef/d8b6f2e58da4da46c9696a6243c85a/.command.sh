#!/bin/bash -ue
python /Users/putriramadani/Documents/GitHub/automated-insilico-vaccinedesign/bin/predict_tcell_ii.py \
    --fasta=accession2_NP_740664.1.fasta \
    --protein-type=accession2 \
    --method=netmhciipan \
    --threshold=500 \
    --length=15 \
    --alleles='HLA-DRB1*03:01,HLA-DRB1*04:01,HLA-DRB1*01:01,HLA-DRB1*07:01,HLA-DRB1*15:01,HLA-DRB1*13:01,HLA-DRB1*11:01' \
    --output=accession2_accession2_NP_740664.1_tcell_ii_epitopes.csv
