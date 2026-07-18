"""
Gene Loss Detection Pipeline (BLASTN XML)
-------------------------------------------
Automates manual BLAST-alignment inspection for comparative genomics /
gene loss studies (e.g. WDR93, RGS4). Parses BLASTN XML with Biopython,
reports indels, builds codons in the QUERY's reading frame, screens for
premature stop codons, and produces a conservative gene-loss summary
per HSP.

Design note: the functions below that already worked (query_anchor,
find_gaps, get_first_codon_offset, build_codons, check_stop_codons) are
left untouched. Everything new is additive — small functions wrapping
each report section, plus an orchestration layer (process_hsp /
process_alignment / process_record / main) that replaces the old flat
script body without changing what gets printed or in what order.
"""

import csv
import argparse
import contextlib
from pathlib import Path
from Bio.Blast import NCBIXML

# When False: print only Gap Analysis, Stop Codon Analysis, Statistics,
# and the Gene Loss Summary — fast screening across many HSPs/species.
# When True: also print the full alignment view and codon-by-codon table
# — for manual inspection / validation against NCBI BLAST.
VERBOSE = False


# --------------------------------------------------
# Function: Convert subject to query-anchored dots
# --------------------------------------------------
def query_anchor(query, subject):
    anchored = ""

    for q, s in zip(query, subject):
        if q == s:
            anchored += "."
        elif q == "-" or s == "-":
            anchored += "-"
        else:
            anchored += s

    return anchored


# --------------------------------------------------
# Function: Find continuous gaps
# Returns (start_position, length)
# --------------------------------------------------
def find_gaps(sequence):

    gaps = []
    in_gap = False
    start = 0

    for i, base in enumerate(sequence):

        if base == "-" and not in_gap:
            start = i
            in_gap = True

        elif base != "-" and in_gap:
            gaps.append((start, i - start))
            in_gap = False

    if in_gap:
        gaps.append((start, len(sequence) - start))

    return gaps


# --------------------------------------------------
# Codon Analysis Functions
# --------------------------------------------------
def get_first_codon_offset(query_start):
    """
    Returns the number of nucleotides to skip
    before reaching the first complete codon.

    Biological reasoning: hsp.query_start is a coordinate in the
    ORIGINAL query CDS, not necessarily the start of a codon. If the
    CDS reading frame began upstream of this HSP, we need to know how
    many bases into this HSP the next codon boundary falls, so query
    position 14 (mid-codon) correctly resolves to first-codon-at-16,
    not a spurious frame that would make every downstream stop-codon
    call wrong.
    """

    return (3 - ((query_start - 1) % 3)) % 3


def build_codons(query_seq, subject_seq, query_start):
    """
    Builds codons using the QUERY reading frame and keeps track of
    query nucleotide positions for every codon.
    """

    offset = get_first_codon_offset(query_start)

    query_seq = query_seq[offset:]
    subject_seq = subject_seq[offset:]

    codons = []
    insertions = []

    q_codon = ""
    s_codon = ""
    pending_insertion = ""

    current_query_pos = query_start + offset
    codon_start = None

    i = 0

    while i < len(query_seq):

        q = query_seq[i]
        s = subject_seq[i]

        # Query gap -> insertion in subject
        if q == "-":
            pending_insertion += s
            i += 1
            continue

        if pending_insertion:
            insertions.append({
                "codon_number": len(codons) + 1,
                "insert_seq": pending_insertion,
            })
            pending_insertion = ""

        if codon_start is None:
            codon_start = current_query_pos

        q_codon += q
        s_codon += s

        current_query_pos += 1

        if len(q_codon) == 3:

            codons.append({
                "query_codon": q_codon,
                "subject_codon": s_codon,
                "query_start": codon_start,
                "query_end": current_query_pos - 1,
            })

            q_codon = ""
            s_codon = ""
            codon_start = None

        i += 1

    if pending_insertion:
        insertions.append({
            "codon_number": len(codons) + 1,
            "insert_seq": pending_insertion,
        })

    return codons, insertions

# --------------------------------------------------
# Stop Codon Functions
# --------------------------------------------------

STOP_CODONS = {"TAA", "TAG", "TGA"}


def check_stop_codons(codons):
    """
    Checks subject codons for premature stop codons.
    """

    print("\nPremature Stop Codon Analysis")
    print("-" * 40)

    found_stop = False
    stop_details = []

    for number, codon in enumerate(codons, start=1):

        query_codon = codon["query_codon"]
        subject_codon = codon["subject_codon"]

        if "-" in subject_codon:
            continue

        if subject_codon.upper() in STOP_CODONS:

            print("\n>>> PREMATURE STOP CODON DETECTED <<<")

            print(f"Codon Number    : {number}")
            print(f"Query Position  : {codon['query_start']}-{codon['query_end']}")
            print(f"Query Codon     : {query_codon}")
            print(f"Subject Codon   : {subject_codon}")
            print(f"Mutation        : {query_codon} -> {subject_codon}")
            print("Interpretation  : Possible premature termination of protein")
            stop_details.append({
                "position": f"{codon['query_start']}..{codon['query_end']}",
                "mutation": f"{query_codon}->{subject_codon}"
        })
       
       
            found_stop = True

    if not found_stop:
        print("No premature stop codons detected.")

    return found_stop, stop_details


# --------------------------------------------------
# NEW: Gap reporting helpers
# --------------------------------------------------
def analyze_gaps(hsp):

    query_gap_positions = find_gaps(hsp.query)
    subject_gap_positions = find_gaps(hsp.sbjct)

    query_gap_list = []
    subject_gap_list = []

    frameshift_found = False

    for start, length in query_gap_positions:

        sequence = hsp.sbjct[start:start+length].replace("-", "")

        query_gap_list.append({
            "start": start + 1,
            "end": start + length,
            "length": length,
            "sequence": sequence
        })

        if length % 3 != 0:
            frameshift_found = True

    for start, length in subject_gap_positions:

        sequence = hsp.query[start:start+length].replace("-", "")

        subject_gap_list.append({
            "start": start + 1,
            "end": start + length,
            "length": length,
            "sequence": sequence
        })

        if length % 3 != 0:
            frameshift_found = True

    return query_gap_list, subject_gap_list, frameshift_found


def print_indel_section(label, gap_list):

    print(f"{label} Indels")
    print("-" * 40)

    if len(gap_list) == 0:
        print("None")
        print()
        return

    for gap in gap_list:

        print(f"Position             : {gap['start']}-{gap['end']}")
        print(f"Length               : {gap['length']} bp")
        print(f"Sequence             : {gap['sequence']}")

        if gap["length"] % 3 == 0:
            print("Possible Frameshift  : NO")
        else:
            print("Possible Frameshift  : YES")

        print()


# --------------------------------------------------
# NEW: Alignment view + codon table printing
# --------------------------------------------------
def print_alignment_view(hsp):
    """Prints the query-anchored dot-format alignment in 60-column
    blocks, exactly as before."""
    anchored = query_anchor(hsp.query, hsp.sbjct)

    print("Legend:")
    print("  . = same as query")
    print("  A/C/G/T = mismatch")
    print("  - = insertion/deletion")
    print()

    for i in range(0, len(hsp.query), 60):
        print("Query          :", hsp.query[i:i + 60])
        print("               ", hsp.match[i:i + 60])
        print("Subject (dots) :", anchored[i:i + 60])
        print("Subject (DNA)  :", hsp.sbjct[i:i + 60])
        print()


def print_codon_table(codons, insertions, hsp):
    """Prints the first-codon-position note and full codon table,
    including any subject-only insertions right before the codon they
    precede (see build_codons docstring for why these are separate
    from the codon pairs rather than folded into one of them)."""
    print(
        f"\nFirst complete codon starts at query position "
        f"{hsp.query_start + get_first_codon_offset(hsp.query_start)}"
    )
    print("\nCODONS")
    print("-" * 40)

    insertions_by_codon = {}
    for ins in insertions:
        insertions_by_codon.setdefault(ins["codon_number"], []).append(ins["insert_seq"])

    for number, codon in enumerate(codons, start=1):

        q = codon["query_codon"]
        s = codon["subject_codon"]
        for insert_seq in insertions_by_codon.get(number, []):
            print(f"Insertion detected in subject before codon {number}: {insert_seq}")
            print()

        print(f"Codon {number}")
        print(f"Query   : {q}")
        print(f"Subject : {s}")

        if "-" in s:
            print(f"Deletion detected in subject codon ({s})")

        print()

    # Trailing insertion after the last real codon (alignment ends on
    # a query gap)
    trailing_codon_number = len(codons) + 1
    for insert_seq in insertions_by_codon.get(trailing_codon_number, []):
        print(f"Insertion detected in subject after codon {len(codons)}: {insert_seq}")
        print()


# --------------------------------------------------
# NEW: Statistics + gene loss summary
# --------------------------------------------------
def print_statistics(hsp):
    """Prints the alignment statistics block and returns percent identity
    so the gene loss summary doesn't need to recompute it."""
    identity = (hsp.identities / hsp.align_length) * 100

    print("Statistics")
    print("-" * 70)
    print(f"Alignment Length : {hsp.align_length}")
    print(f"Identity         : {identity:.2f}%")
    print(f"E-value          : {hsp.expect}")
    print(f"Bit Score        : {hsp.bits}")
    print(f"Gap Characters   : {hsp.gaps}")
    print(f"Mismatches       : {hsp.align_length - hsp.identities}")
    print(f"\nQuery Coordinates   : {hsp.query_start} - {hsp.query_end}")
    print(f"Subject Coordinates : {hsp.sbjct_start} - {hsp.sbjct_end}")

    return identity


def print_gene_loss_summary(identity, query_gap_list, subject_gap_list,
                             insertion_count, frameshift_found, stop_found, hsp):
    """
    Prints the conservative gene-loss summary block.

    Biological reasoning: this deliberately never claims gene loss is
    confirmed. Frameshifts and premature stops are each necessary but
    not sufficient evidence on their own — assembly errors, sequencing
    gaps, or alignment artifacts can produce the same signature. The
    interpretation ladder here only strengthens when BOTH lines of
    evidence agree, and always recommends manual verification short of
    that.

    Subject-only insertions are reported here too: they're bases the
    subject has that the query doesn't, which is just as relevant to
    a gene-loss call as a deletion is (e.g. a large in-frame insertion
    could still disrupt a functional domain even without shifting the
    reading frame).

    Note on naming: this heuristic only flags a "candidate" — a gap
    whose length isn't divisible by 3. It does not prove a biologically
    significant frameshift, since multiple indels in the same HSP can
    combine or cancel out their effect on the reading frame later in
    the alignment. That's why the CSV field is "frameshift_candidate"
    rather than "potential_frameshift" or a bare boolean with no
    context.

    Returns a flat dict of everything printed, for CSV aggregation
    across HSPs/subjects.
    """
    print("\n")
    print("=" * 60)
    print("GENE LOSS SUMMARY")
    print("=" * 60)

    print(f"Identity                : {identity:.2f}%")
    print(f"Subject Indels          : {len(subject_gap_list)}")
    print(f"Query Indels            : {len(query_gap_list)}")
    print(f"Subject-only Insertions : {insertion_count}")

    print(f"Frameshift Candidate    : {'YES' if frameshift_found else 'NO'}")
    print(f"Premature Stop Codon    : {'YES' if stop_found else 'NO'}")

    print()

    if frameshift_found and stop_found:
        interpretation = "Strong evidence for coding sequence disruption."
    elif frameshift_found:
        interpretation = (
            "Frameshift candidate detected (gap length not divisible by 3). "
            "Verify manually."
        )
    elif stop_found:
        interpretation = "Premature stop codon detected. Verify manually."
    else:
        interpretation = "No obvious evidence of coding sequence disruption."

    print(f"Interpretation : {interpretation}")
    query_gap_text = "; ".join(
        f"{g['start']}-{g['end']} ({g['length']} bp): {g['sequence']}"
        for g in query_gap_list
    )

    subject_gap_text = "; ".join(
        f"{g['start']}-{g['end']} ({g['length']} bp): {g['sequence']}"
        for g in subject_gap_list
    )
    return {
        "identity_pct": round(identity, 2),
        "alignment_length": hsp.align_length,
        "query_start": hsp.query_start,
        "query_end": hsp.query_end,
        "subject_start": hsp.sbjct_start,
        "subject_end": hsp.sbjct_end,
        "subject_indels": len(subject_gap_list),
        "query_indels": len(query_gap_list),
        "subject_only_insertions": insertion_count,
        "frameshift_candidate": frameshift_found,
        "premature_stop_codon": stop_found,
        "interpretation": interpretation,
        "query_gap_details": query_gap_text,
        "subject_gap_details": subject_gap_text,
    
    }


# --------------------------------------------------
# NEW: Orchestration (replaces the old flat script body)
# --------------------------------------------------
def process_hsp(hsp, hsp_number):
    """Runs the full report for a single HSP, in the same order as the
    original script: gaps -> alignment view -> codons -> stop codons ->
    statistics -> gene loss summary. Returns the summary row dict for
    CSV aggregation."""
    print(f"\nHSP {hsp_number}")
    print("\nAlignment")
    print("-" * 70)

    query_gap_list, subject_gap_list, frameshift_found = analyze_gaps(hsp)

    print("Gap Analysis")
    print("-" * 70)
    print(f"Number of Query Indels   : {len(query_gap_list)}")
    print(f"Number of Subject Indels : {len(subject_gap_list)}")

    print_indel_section("Subject", subject_gap_list)

    if VERBOSE:
        print_alignment_view(hsp)

    # Codon building always runs (stop-codon screening needs it even in
    # compact mode) — only the printed table is gated by VERBOSE.
    codons, insertions = build_codons(hsp.query, hsp.sbjct, hsp.query_start)

    if VERBOSE:
        print_codon_table(codons, insertions, hsp)

    stop_found, stop_details = check_stop_codons(codons)

    # Biological reasoning: stop-codon detection assumes each codon is
    # correctly in-frame. If an earlier indel already shifted the frame,
    # every codon from that point on is read against the wrong triplet
    # boundaries — so a "no premature stop found" result (or one that
    # IS found) may not reflect the true translated sequence. Flag this
    # explicitly rather than let a frameshift silently undermine a
    # stop-codon call downstream.
    if frameshift_found:
        print(
            "\nCaution: a frameshift was detected earlier in this HSP — "
            "downstream codons (including any stop-codon calls above) may "
            "not reflect the true reading frame. Interpret with caution."
        )

    print()
    print_indel_section("Query", query_gap_list)

    identity = print_statistics(hsp)
    
    summary_row = print_gene_loss_summary(
        identity, query_gap_list, subject_gap_list,
        len(insertions), frameshift_found, stop_found, hsp
    )
    summary_row["hsp_number"] = hsp_number
    summary_row["stop_positions"] = "; ".join(
        s["position"] for s in stop_details
    )

    summary_row["stop_mutations"] = "; ".join(
        s["mutation"] for s in stop_details
)

    return summary_row


def process_alignment(alignment):
    """Prints the SUBJECT header and processes every HSP under one
    alignment (one query-vs-one-subject-sequence hit). Returns the list
    of summary rows for all HSPs in this alignment, tagged with the
    subject title for CSV aggregation."""
    print("\n" + "=" * 70)
    print("SUBJECT")
    print("=" * 70)
    print(f"Subject : {alignment.title}")

    rows = []
    for hsp_number, hsp in enumerate(alignment.hsps, start=1):
        row = process_hsp(hsp, hsp_number)
        row["subject"] = alignment.title
        rows.append(row)

    return rows


def process_record(blast_record):
    """Prints the BLAST summary and processes every alignment (hit) for
    one query record. Returns the combined list of summary rows across
    every HSP in every alignment."""
    print("=" * 60)
    print("BLAST SUMMARY")
    print("=" * 60)

    print("Query:")
    print(blast_record.query)

    print("\nQuery Length:")
    print(blast_record.query_length)

    print("\nTotal Hits:")
    print(len(blast_record.alignments))

    all_rows = []
    for alignment in blast_record.alignments:
        all_rows.extend(process_alignment(alignment))

    return all_rows


def write_summary_csv(xml_file, rows):
    """
    Writes a cleaner per-HSP summary table.
    """

    import re
    from pathlib import Path
    import csv

    if not rows:
        print("No HSPs processed — no summary CSV written.")
        return

    csv_path = Path(xml_file).with_suffix(".summary.csv")

    headers = [
        "HSP",
        "Chromosome",
        "Identity (%)",
        "Alignment (bp)",
        "Query Range",
        "Subject Range",
        "Query Deletions",
        "Subject Deletions",
        "Insertions",
        "Frameshift",
        "Stop Codon",
        "Stop Position",
        "Stop Mutation",
        "Interpretation"
    ]

    with open(csv_path, "w", newline="") as f:

        writer = csv.writer(f)
        writer.writerow(headers)

        for row in rows:

            subject = row["subject"]

            match = re.search(r'chromosome\s+([A-Za-z0-9]+)', subject, re.IGNORECASE)

            if match:
                chromosome = "Chr " + match.group(1)
            else:
                chromosome = subject.split()[0]

            writer.writerow([
                row["hsp_number"],
                chromosome,
                f'{row["identity_pct"]:.2f}',
                row["alignment_length"],
                f'{row["query_start"]}-{row["query_end"]}',
                f'{row["subject_start"]}-{row["subject_end"]}',
                row["query_gap_details"],
                row["subject_gap_details"],
                row["subject_only_insertions"],
                "YES" if row["frameshift_candidate"] else "NO",
                "YES" if row["premature_stop_codon"] else "NO",
                row["stop_positions"],
                row["stop_mutations"],
                row["interpretation"]
            ])

    print(f"Summary CSV written to: {csv_path}")


def main(xml_file, output_file=None):
    """
    Entry point. Uses NCBIXML.read() as before (single-query XML) — kept
    as-is since that matches your current one-query-per-file workflow.
    Swapping to NCBIXML.parse() is a small, isolated change for later
    when "support multiple XML files in one run" gets built.

    If output_file is given, the full text report is redirected there
    instead of the terminal — useful once VERBOSE=True output gets long
    enough that scrolling the terminal isn't practical. The CSV summary
    and the final status messages always go to the real terminal, even
    when the report itself is redirected to a file, so you always see
    confirmation of what was written and where.
    """
    with open(xml_file) as result_handle:
        blast_record = NCBIXML.read(result_handle)

    if output_file:
        with open(output_file, "w") as f, contextlib.redirect_stdout(f):
            summary_rows = process_record(blast_record)
        print(f"Report written to: {output_file}")
    else:
        summary_rows = process_record(blast_record)

    write_summary_csv(xml_file, summary_rows)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Gene loss detection pipeline (BLASTN XML)."
    )
    parser.add_argument(
        "xml_file", nargs="?", default="Results/xml/Ostrich_vs_Anser.xml",
        help="Path to the BLASTN XML result file "
             "(default: Results/xml/Ostrich_vs_Anser.xml)",
    )
    parser.add_argument(
        "output_file", nargs="?", default=None,
        help="Optional path to write the text report to a file instead "
             "of the terminal (e.g. output.txt)",
    )
    args = parser.parse_args()

    main(args.xml_file, args.output_file)
