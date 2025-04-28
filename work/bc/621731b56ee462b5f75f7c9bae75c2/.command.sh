#!/bin/bash -ue
python /Users/putriramadani/Documents/GitHub/h5n1-insilico-vaccine-epitopes/bin/evaluate_vaccine.py \
    --vaccine=neuraminidase_vaccine_construct.fasta \
    --protein-type=neuraminidase \
    --linker=GPGPG \
    --iedb-api-url=http://tools-api.iedb.org/tools_api/ \
    --output-evaluation=neuraminidase_vaccine_evaluation.txt \
    --output-properties=neuraminidase_vaccine_properties.csv \
    --output-colabfold=neuraminidase_vaccine_for_colabfold.fasta
