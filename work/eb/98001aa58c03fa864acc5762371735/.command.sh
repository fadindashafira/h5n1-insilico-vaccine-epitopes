#!/bin/bash -ue
python /Users/putriramadani/Documents/GitHub/h5n1-insilico-vaccine-epitopes/bin/evaluate_vaccine.py \
    --vaccine=hemagglutinin_vaccine_construct.fasta \
    --protein-type=hemagglutinin \
    --linker=GPGPG \
    --iedb-api-url=http://tools-api.iedb.org/tools_api/ \
    --output-evaluation=hemagglutinin_vaccine_evaluation.txt \
    --output-properties=hemagglutinin_vaccine_properties.csv \
    --output-colabfold=hemagglutinin_vaccine_for_colabfold.fasta
