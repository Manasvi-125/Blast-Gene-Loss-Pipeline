#!/bin/bash

QUERY=$1
DB=$2
PREFIX=$3

echo "======================================"
echo "Running BLAST..."
echo "======================================"

bash Scripts/run_blast.sh "$QUERY" "$DB" "$PREFIX"

echo
echo "======================================"
echo "Running Gene Loss Analysis..."
echo "======================================"

python3 Scripts/blast_gene_loss_analysis.py Results/xml/${PREFIX}.xml > Results/xml/${PREFIX}_report.txt

echo
echo "Done!"
echo
echo "XML:      Results/xml/${PREFIX}.xml"
echo "CSV:      Results/xml/${PREFIX}.summary.csv"
echo "REPORT:   Results/xml/${PREFIX}_report.txt"