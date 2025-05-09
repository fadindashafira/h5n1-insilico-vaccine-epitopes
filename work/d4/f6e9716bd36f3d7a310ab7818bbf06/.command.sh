#!/bin/bash -ue
python /Users/putriramadani/Documents/GitHub/h5n1-insilico-vaccine-epitopes/bin/predict_bcell.py \
    --fasta=neuraminidase_BAL61230.1.fasta \
    --protein-type=neuraminidase \
    --method=Bepipred-2.0 \
    --threshold=0.5 \
    --window-size=9 \
    --output=neuraminidase_neuraminidase_BAL61230.1_bcell_epitopes.csv
