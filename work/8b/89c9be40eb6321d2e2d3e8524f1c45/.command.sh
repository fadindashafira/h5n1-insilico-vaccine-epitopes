#!/bin/bash -ue
python /Users/putriramadani/Documents/GitHub/automated-insilico-vaccinedesign/bin/screen_epitopes.py \
  --bcell accession2_accession2_NP_740664.1_bcell_epitopes.csv \
  --tcelli accession2_accession2_NP_740664.1_tcell_i_epitopes.csv \
  --tcellii accession2_accession2_NP_740664.1_tcell_ii_epitopes.csv \
  --output accession2_epitope_screening.csv
