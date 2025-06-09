#!/bin/bash

EVAL_DIR="eval/eval_dir"
STATS_OUTPUT=$EVAL_DIR/end.csv
REF_OUTPUT=$EVAL_DIR/ref.csv
REF_DEGREE_SEQ=$EVAL_DIR/ref_degree_sequences.json
END_DEGREE_SEQ=$EVAL_DIR/end_degree_sequences.json

# Run RECCS
cd ../build
cmake .. && make
./reccs ../data/cit_hepph.tsv -v -c ../data/cit_hepph_leiden_0.01.tsv

cd ..

# Create the evaluation directory if it doesn't exist
mkdir -p $EVAL_DIR

# Run stats on the output
python3 extlib/stats.py -i build/output.tsv -e data/cit_hepph_leiden_0.01.tsv -o $STATS_OUTPUT

# Run stats on the reference
python3 extlib/stats.py -i data/cit_hepph.tsv -e data/cit_hepph_leiden_0.01.tsv -o $REF_OUTPUT

# Run the evaluation script
python3 eval/check_outputs.py -s $STATS_OUTPUT \
                               -d $END_DEGREE_SEQ \
                               -rs $REF_OUTPUT \
                               -rd $REF_DEGREE_SEQ
