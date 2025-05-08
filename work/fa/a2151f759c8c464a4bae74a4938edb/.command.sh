#!/bin/bash -ue
# Add error handling
set -e
python /Users/putriramadani/Documents/GitHub/h5n1-insilico-vaccine-epitopes/bin/predict_bcell.py \
    --fasta=hemagglutinin_BAL61222.1.fasta \
    --protein-type=hemagglutinin \
    --method=Bepipred \
    --threshold=0.5 \
    --window-size=9 \
    --output=hemagglutinin_hemagglutinin_BAL61222.1_bcell_epitopes.csv || touch hemagglutinin_hemagglutinin_BAL61222.1_bcell_epitopes.csv
