# Title: Isolating rows of interest from raw RNA-seq counts matrix
# Author: Dr. Alice M. Godden

import csv

# File paths
raw_counts_file = 'rawcounts.csv' # rnaseq raw counts file
target_ids_file = 'mir200' # list of gene ids targeted by mirna
output_file = 'mir_matched_counts_200.csv' # output of raw counts for matched gene

# Read target miRNA IDs from file into a list
target_mirna_ids = []
with open(target_ids_file, 'r') as f:
    for line in f:
        target_mirna_ids.append(line.strip())  # Assuming each line contains one miRNA ID

# Set to store matched rows
matched_rows = []

# Read raw counts file and collect rows matching the target miRNA IDs
with open(raw_counts_file, 'r') as f:
    reader = csv.reader(f)
    header = next(reader)  # Read header
    matched_rows.append(header)  # Store header in matched_rows
    for row in reader:
        miRNA = row[0].strip()  # Assuming miRNA ID is in the first column
        if miRNA in target_mirna_ids:
            matched_rows.append(row)

# Write matched rows to output CSV file
with open(output_file, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerows(matched_rows)

print(f"Matching rows saved to {output_file}.")
