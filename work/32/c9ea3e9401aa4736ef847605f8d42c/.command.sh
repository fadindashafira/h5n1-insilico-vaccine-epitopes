#!/bin/bash -ue
python /Users/putriramadani/Documents/GitHub/automated-insilico-vaccinedesign/bin/evaluate_vaccine.py \
    --vaccine=combined_vaccine_construct.fasta \
    --protein-type=combined \
    --linker=GPGPG \
    --iedb-api-url=http://tools-api.iedb.org/tools_api/ \
    --output-evaluation=combined_vaccine_evaluation.txt \
    --output-properties=combined_vaccine_properties.csv \
    --output-colabfold=combined_vaccine_for_colabfold.fasta
