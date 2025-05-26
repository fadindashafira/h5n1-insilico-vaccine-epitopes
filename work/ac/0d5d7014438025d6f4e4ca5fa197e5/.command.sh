#!/bin/bash -ue
python /Users/putriramadani/Documents/GitHub/automated-insilico-vaccinedesign/bin/evaluate_vaccine.py \
    --vaccine=accession1_vaccine_construct.fasta \
    --protein-type=accession1 \
    --linker=GPGPG \
    --iedb-api-url=http://tools-api.iedb.org/tools_api/ \
    --output-evaluation=accession1_vaccine_evaluation.txt \
    --output-properties=accession1_vaccine_properties.csv \
    --output-colabfold=accession1_vaccine_for_colabfold.fasta
