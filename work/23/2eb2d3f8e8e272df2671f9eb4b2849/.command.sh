#!/bin/bash -ue
python /Users/putriramadani/Documents/GitHub/h5n1-insilico-vaccine-epitopes/bin/predict_bcell.py \
    --fasta=hemagglutinin_BAL61222.1.fasta \
    --protein-type=hemagglutinin \
    --method=Bepipred-2.0 \
    --threshold=0.5 \
    --window-size=9 \
    --output=hemagglutinin_hemagglutinin_BAL61222.1_bcell_epitopes.csv
