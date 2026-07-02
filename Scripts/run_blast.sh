#!/bin/bash

QUERY=$1
DATABASE=$2
NAME=$3

echo "Running BLAST..."

blastn \
-query "$QUERY" \
-db "$DATABASE" \
-evalue 1e-20 \
-max_target_seqs 10 \
-max_hsps 20 \
-outfmt 6 \
-out Results/blast/${NAME}.tsv

blastn \
-query "$QUERY" \
-db "$DATABASE" \
-evalue 1e-20 \
-max_target_seqs 10 \
-max_hsps 20 \
-out Results/alignments/${NAME}_alignment.txt

blastn \
-query "$QUERY" \
-db "$DATABASE" \
-evalue 1e-20 \
-max_target_seqs 10 \
-max_hsps 20 \
-outfmt 5 \
-out Results/xml/${NAME}.xml

echo ""
echo "BLAST completed!"