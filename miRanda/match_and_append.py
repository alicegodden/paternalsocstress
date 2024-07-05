# Title: match and append genes and miRNAs
# Author: Dr. Alice M. Godden

import pandas as pd

# Load matchy.txt into a DataFrame
matchy_df = pd.read_csv("matchy.txt", delimiter='\t')

# Load matched_genes_mirnas.txt into a DataFrame
matched_genes_mirnas_df = pd.read_csv("matched_genes_mirnas.txt", delimiter='\t')

# Create a dictionary for easy lookup from matchy_df
matchy_dict = matchy_df.set_index('initial_alias').T.to_dict()

# Initialize a list to store the new rows with appended information
new_rows = []

# Iterate through each row in matched_genes_mirnas_df
for index, row in matched_genes_mirnas_df.iterrows():
    initial_alias = row['Seq2']
    if initial_alias in matchy_dict:
        # If a match is found, append the additional information
        match_info = matchy_dict[initial_alias]
        new_row = list(row) + list(match_info.values())
    else:
        # If no match is found, append empty strings for the additional columns
        new_row = list(row) + [""] * (len(matchy_df.columns) - 1)
    new_rows.append(new_row)

# Create a new DataFrame with the combined information
columns = list(matched_genes_mirnas_df.columns) + list(matchy_df.columns[1:])
combined_df = pd.DataFrame(new_rows, columns=columns)

# Write the combined DataFrame to a new file
combined_df.to_csv("matched_annotated_miranda_output.txt", sep='\t', index=False)
