#!/bin/bash

# Run RECCS
cd ../build
cmake .. && make
./reccs ../data/cit_hepph.tsv -v -c ../data/cit_hepph_leiden_0.01.tsv

# Run stats on the output
cd ..
python3 extlib/stats.py -i build/output.tsv -e data/cit_hepph_leiden_0.01.tsv -o end.csv

# Run stats on the reference
python3 extlib/stats.py -i data/cit_hepph.tsv -e data/cit_hepph_leiden_0.01.tsv -o ref.csv

# Run the evaluation script
python3 eval/check_outputs.py -s end.csv -d end_degree_sequences.json -rs ref.csv -rd ref_degree_sequences.json
