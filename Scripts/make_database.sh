#!/bin/bash
if [ $# -ne 2 ]; then
    echo "Usage:"
    echo "bash Scripts/make_database.sh <genome.fna> <output_database_prefix>"
    exit 1
fi

GENOME=$1
DBNAME=$2

makeblastdb \
-in "$GENOME" \
-dbtype nucl \
-out "$DBNAME"

echo ""
echo "Database created successfully!"