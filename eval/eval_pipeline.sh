#!/bin/bash

NETWORK=$(realpath ../data/cen.tsv)
CLUSTERING=$(realpath ../data/cen_0.01_cm.tsv)

EVAL_DIR="eval/eval_dir"
STATS_OUTPUT=$EVAL_DIR/end.csv
REF_OUTPUT=$EVAL_DIR/ref.csv
REF_DEGREE_SEQ=$EVAL_DIR/ref_degree_sequences.json
END_DEGREE_SEQ=$EVAL_DIR/end_degree_sequences.json
SBM_OUTPUT=$EVAL_DIR/sbm.csv
SBM_DEGREE_SEQ=$EVAL_DIR/sbm_degree_sequences.json
PLOT_OUTPUT=$EVAL_DIR/RECCS_v_SBM.png

# Run RECCS
cd ../build
cmake .. && make
./reccs $NETWORK -c $CLUSTERING -o output.tsv -v > test.txt 2> err.txt

# Get the most recent temp{timestamp}/ directory
TEMP_DIR=build/$(ls -td temp*/ | head -n 1)
CLUSTERED_SBM_OUTPUT=$TEMP_DIR/clustered_sbm/syn_sbm.tsv

cd ..

# Create the evaluation directory if it doesn't exist
mkdir -p $EVAL_DIR

# Run stats on the output
python3 extlib/stats.py -i build/output.tsv -e $CLUSTERING -o $STATS_OUTPUT

# Run stats on the SBM output
python3 extlib/stats.py -i $CLUSTERED_SBM_OUTPUT -e $CLUSTERING -o $SBM_OUTPUT

# Run stats on the reference
python3 extlib/stats.py -i $NETWORK -e $CLUSTERING -o $REF_OUTPUT

# Run the evaluation script
python3 eval/check_outputs.py -s $STATS_OUTPUT \
                              -d $END_DEGREE_SEQ \
                              -rs $REF_OUTPUT \
                              -rd $REF_DEGREE_SEQ \
                              -ss $SBM_OUTPUT \
                              -sd $SBM_DEGREE_SEQ \
                              -p $PLOT_OUTPUT
