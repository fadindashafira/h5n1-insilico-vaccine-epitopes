#!/bin/bash -ue
python /Users/putriramadani/Documents/GitHub/measles-rubella-vaccinedesign/bin/predict_tcell_ii.py \
    --fasta=accession1_ABK40530.1.fasta \
    --protein-type=accession1 \
    --method=netmhciipan \
    --threshold=500 \
    --length=15 \
    --alleles='HLA-DRB1*03:01,HLA-DRB1*04:01,HLA-DRB1*01:01,HLA-DRB1*07:01,HLA-DRB1*15:01,HLA-DRB1*13:01,HLA-DRB1*11:01' \
    --output=accession1_accession1_ABK40530.1_tcell_ii_epitopes.csv
