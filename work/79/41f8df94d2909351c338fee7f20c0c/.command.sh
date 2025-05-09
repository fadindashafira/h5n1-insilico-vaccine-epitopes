#!/bin/bash -ue
python /Users/putriramadani/Documents/GitHub/h5n1-insilico-vaccine-epitopes/bin/predict_tcell_i.py \
    --fasta=hemagglutinin_BAL61222.1.fasta \
    --protein-type=hemagglutinin \
    --method=NetMHCpan \
    --threshold=500 \
    --length=9 \
    --alleles='HLA-A*02:01,HLA-B*07:02,HLA-B*35:01,HLA-A*11:01,HLA-A*24:02,HLA-A*01:01,HLA-C*07:01' \
    --output=hemagglutinin_hemagglutinin_BAL61222.1_tcell_i_epitopes.csv
