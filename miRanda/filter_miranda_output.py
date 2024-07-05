# Title: filter_mirands_output.py
# Author: Dr. Alice M. Godden

def read_sig_transcripts(file_path):
    with open(file_path, 'r') as file:
        return set(line.strip() for line in file)

def filter_miranda_output(sig_transcripts, miranda_file, output_file):
    with open(miranda_file, 'r') as infile, open(output_file, 'w') as outfile:
        header = infile.readline()  # Assuming the first line is the header
        outfile.write(header)  # Write header to output file

        for line in infile:
            parts = line.strip().split('\t')
            if len(parts) > 1 and parts[1] in sig_transcripts:
                outfile.write(line)

if __name__ == "__main__":
    sig_transcripts_file = "sig_transcripts"  # Replace with your actual file name
    miranda_output_file = "MirandaOutput.tab"  # Replace with your actual file name
    matched_genes_mirnas_file = "matched_genes_mirnas.txt"

    sig_transcripts = read_sig_transcripts(sig_transcripts_file)
    filter_miranda_output(sig_transcripts, miranda_output_file, matched_genes_mirnas_file)
    print(f"Matched rows saved to {matched_genes_mirnas_file}")
