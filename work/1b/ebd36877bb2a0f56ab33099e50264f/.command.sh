#!/bin/bash -ue
python /Users/putriramadani/Documents/GitHub/automated-insilico-vaccinedesign/bin/evaluate_vaccine.py \
    --vaccine=accession2_vaccine_construct.fasta \
    --protein-type=accession2 \
    --linker=GPGPG \
    --iedb-api-url=http://tools-api.iedb.org/tools_api/ \
    --output-evaluation=accession2_vaccine_evaluation.txt \
    --output-properties=accession2_vaccine_properties.csv \
    --output-colabfold=accession2_vaccine_for_colabfold.fasta
