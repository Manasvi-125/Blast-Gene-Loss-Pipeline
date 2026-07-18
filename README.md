# BLAST Gene Loss Pipeline

A command-line pipeline developed to assist in the validation of gene loss by analysing BLAST alignments for coding sequence disruptions.

## Features

- Performs nucleotide BLAST searches using BLAST+
- Parses BLAST XML output automatically
- Detects:
  - Query and subject indels
  - Frameshift-causing insertions and deletions
  - Premature stop codons
- Generates a summary CSV file
- Generates a detailed text report for each BLAST analysis
- Supports automated execution using shell scripts

## Repository Structure

```
Blast-Gene-Loss/
│
├── Scripts/
│   ├── blast_gene_loss_analysis.py
│   ├── blast_indel_analysis.py
│   ├── make_database.sh
│   ├── run_blast.sh
│   └── run_pipeline.sh
│
├── example_data/
│   └── example_query.fasta
│
└── .gitignore
```

## Requirements

- Python 3
- BLAST+ (NCBI)
- Linux / WSL

## Running the Pipeline

1. Build a BLAST database.

```bash
bash Scripts/make_database.sh
```

2. Run BLAST.

```bash
bash Scripts/run_blast.sh <query.fasta> <database> <output_prefix>
```

3. Run the complete pipeline.

```bash
bash Scripts/run_pipeline.sh
```

## Outputs

The pipeline generates:

- BLAST XML alignment
- Detailed text report
- Summary CSV containing:
  - Identity
  - Alignment length
  - Query and subject coordinates
  - Indels
  - Frameshift detection
  - Premature stop codons
  - Stop codon positions
  - Stop codon mutations

## Example

An example FASTA query sequence is included in the `example_data` folder for testing.

## Status

**This project is under active development, and additional features and improvements are being added.**
