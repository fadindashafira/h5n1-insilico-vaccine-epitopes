#!/bin/bash -ue
python /Users/putriramadani/Documents/GitHub/automated-insilico-vaccinedesign/bin/screen_epitopes.py \
  --bcell accession1_accession1_ABK40530.1_bcell_epitopes.csv \
  --tcelli accession1_accession1_ABK40530.1_tcell_i_epitopes.csv \
  --tcellii accession1_accession1_ABK40530.1_tcell_ii_epitopes.csv \
  --output accession1_epitope_screening.csv
