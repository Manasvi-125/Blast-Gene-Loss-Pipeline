#!/bin/bash

GENOME=$1
DBNAME=$2

makeblastdb \
-in "$GENOME" \
-dbtype nucl \
-out "$DBNAME"

echo ""
echo "Database created successfully!"