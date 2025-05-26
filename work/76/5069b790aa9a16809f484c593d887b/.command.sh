#!/bin/bash -ue
# Add error handling
set -e
python /Users/putriramadani/Documents/GitHub/automated-insilico-vaccinedesign/bin/predict_bcell.py \
    --fasta=accession2_NP_740664.1.fasta \
    --protein-type=accession2 \
    --method=Bepipred \
    --threshold=0.5 \
    --window-size=9 \
    --output=accession2_accession2_NP_740664.1_bcell_epitopes.csv || touch accession2_accession2_NP_740664.1_bcell_epitopes.csv
