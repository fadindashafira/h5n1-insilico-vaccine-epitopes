#!/bin/bash -ue
# First, concatenate the epitope files
head -n 1 neuraminidase_combined_epitopes.csv > combined_epitopes_all.csv
for f in neuraminidase_combined_epitopes.csv hemagglutinin_combined_epitopes.csv; do
    tail -n +2 $f >> combined_epitopes_all.csv
done

# Now run with the combined file
python /Users/putriramadani/Documents/GitHub/h5n1-insilico-vaccine-epitopes/bin/design_vaccine.py             --combined-epitopes=combined_epitopes_all.csv             --protein-type=combined             --linker=GPGPG             --max-epitopes=10             --min-epitopes=3             --leading-seq=M             --trailing-seq=             --similarity-threshold=0.7             --output-fasta=combined_vaccine_construct.fasta             --output-report=combined_vaccine_report.html
