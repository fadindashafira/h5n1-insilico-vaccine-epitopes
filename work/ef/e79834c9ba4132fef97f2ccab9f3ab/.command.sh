#!/bin/bash -ue
python /Users/putriramadani/Documents/GitHub/h5n1-insilico-vaccine-epitopes/bin/predict_tcell_ii.py \
    --fasta=hemagglutinin_BAL61222.1.fasta \
    --protein-type=hemagglutinin \
    --method=netmhciipan \
    --threshold=500 \
    --length=15 \
    --alleles='HLA-DRB1*03:01,HLA-DRB1*04:01,HLA-DRB1*01:01,HLA-DRB1*07:01,HLA-DRB1*15:01,HLA-DRB1*13:01,HLA-DRB1*11:01' \
    --output=hemagglutinin_hemagglutinin_BAL61222.1_tcell_ii_epitopes.csv
