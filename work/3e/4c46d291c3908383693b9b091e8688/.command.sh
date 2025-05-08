#!/bin/bash -ue
# Add error handling
set -e
python /Users/putriramadani/Documents/GitHub/measles-rubella-vaccinedesign/bin/predict_bcell.py \
    --fasta=accession1_ABK40530.1.fasta \
    --protein-type=accession1 \
    --method=Bepipred \
    --threshold=0.5 \
    --window-size=9 \
    --output=accession1_accession1_ABK40530.1_bcell_epitopes.csv || touch accession1_accession1_ABK40530.1_bcell_epitopes.csv
