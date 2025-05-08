#!/bin/bash -ue
# Add error handling
set -e
python /Users/putriramadani/Documents/GitHub/h5n1-insilico-vaccine-epitopes/bin/predict_bcell.py \
    --fasta=neuraminidase_BAL61230.1.fasta \
    --protein-type=neuraminidase \
    --method=Bepipred \
    --threshold=0.5 \
    --window-size=9 \
    --output=neuraminidase_neuraminidase_BAL61230.1_bcell_epitopes.csv || touch neuraminidase_neuraminidase_BAL61230.1_bcell_epitopes.csv
