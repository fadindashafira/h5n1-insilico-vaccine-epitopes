#!/bin/bash -ue
python /Users/putriramadani/Documents/GitHub/measles-rubella-vaccinedesign/bin/predict_tcell_i.py \
    --fasta=accession2_NP_740664.1.fasta \
    --protein-type=accession2 \
    --method=netmhcpan \
    --threshold=500 \
    --length=9 \
    --alleles='HLA-A*02:01,HLA-B*07:02,HLA-B*35:01,HLA-A*11:01,HLA-A*24:02,HLA-A*01:01,HLA-C*07:01' \
    --output=accession2_accession2_NP_740664.1_tcell_i_epitopes.csv
