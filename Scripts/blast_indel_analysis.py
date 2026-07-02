from Bio.Blast import NCBIXML


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
# Read BLAST XML
# --------------------------------------------------
xml_file = "Results/xml/Ostrich_vs_Anser.xml"

with open(xml_file) as result_handle:
    blast_record = NCBIXML.read(result_handle)


# --------------------------------------------------
# BLAST Summary
# --------------------------------------------------
print("=" * 60)
print("BLAST SUMMARY")
print("=" * 60)

print("Query:")
print(blast_record.query)

print("\nQuery Length:")
print(blast_record.query_length)

print("\nTotal Hits:")
print(len(blast_record.alignments))


# --------------------------------------------------
# Loop through hits
# --------------------------------------------------
for alignment in blast_record.alignments:

    print("\n" + "=" * 70)
    print("SUBJECT")
    print("=" * 70)

    print(alignment.title)

    # ----------------------------------------------
    # Loop through each HSP
    # ----------------------------------------------
    for hsp in alignment.hsps:

        print("\nAlignment")
        print("-" * 70)
    
        # ------------------------------------------
        # Gap analysis
        # ------------------------------------------
        query_gap_list = find_gaps(hsp.query)
        subject_gap_list = find_gaps(hsp.sbjct)

        print("Gap Analysis")
        print("-" * 70)

        print(f"Number of Query Indels   : {len(query_gap_list)}")
        print(f"Number of Subject Indels : {len(subject_gap_list)}")

    
        # ------------------------------------------
        # Subject gaps
        # ------------------------------------------
        print("Subject Indels")
        print("-" * 30)

        if len(subject_gap_list) == 0:
            print("None")

        for start, length in subject_gap_list:

            print(f"Start Position : {start + 1}")
            print(f"Length         : {length} bp")

            if length % 3 == 0:
                print("Possible Frameshift : NO")
            else:
                print("Possible Frameshift : YES (based on gap length)")

            print()

        # ------------------------------------------
        # Print Alignment
        # ------------------------------------------
        anchored = query_anchor(hsp.query, hsp.sbjct)
        print("Legend:")
        print("  . = same as query")
        print("  A/C/G/T = mismatch")
        print("  - = insertion/deletion")
        print()

        for i in range(0, len(hsp.query), 60):

            print("Query          :", hsp.query[i:i+60])
            print("               ", hsp.match[i:i+60])
            print("Subject (dots) :", anchored[i:i+60])
            print("Subject (DNA)  :", hsp.sbjct[i:i+60])
            print()
        # ------------------------------------------
        # Query gaps
        # ------------------------------------------
        print("\nQuery Indels")
        print("-" * 30)

        if len(query_gap_list) == 0:
            print("None")

        for start, length in query_gap_list:

            print(f"Start Position : {start + 1}")
            print(f"Length         : {length} bp")

            if length % 3 == 0:
                print("Possible Frameshift : NO")
            else:
                print("Possible Frameshift : YES")

            print()
        # ------------------------------------------
        # Alignment statistics
        # ------------------------------------------
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