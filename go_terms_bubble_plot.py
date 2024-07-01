# Title: Plotting GO terms bubble plots
# Subtitle: Use of ShinyGO .csv outputs to plot custom bubble plots
# Author: Dr. Alice M. Godden

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Read data from CSV file
df = pd.read_csv('enrichment.csv')

# Extract columns
enrichment_fdr = df['Enrichment FDR'].tolist()
n_genes = df['nGenes'].tolist()
go_termsFE = df['Fold Enrichment'].tolist()
descriptions = df['Pathway'].tolist()

# Convert the descriptions list to a NumPy array
descriptions = np.array(descriptions)
# Replace space with newline character
descriptions = [desc.replace(' ', '\n', 1) for desc in descriptions]

# Calculate the bubble sizes
bubble_sizes = np.array(n_genes) * 50  # Scaling up for better visualization adjust to less or more than 100 percent to change size down or up

# Create a bubble plot
plt.figure(figsize=(12, 8)) # 10 x 6 width v height 12 x 8

# Reverse the y-axis to have each term in a row
plt.gca().invert_yaxis()

# Bubble plot with x and y data
sc = plt.scatter(go_termsFE, np.arange(len(descriptions)), s=bubble_sizes, c=enrichment_fdr, cmap='magma', alpha=0.7)

# Add labels and title
plt.xlabel('Fold Enrichment', fontweight='bold', fontsize=14)
plt.ylabel('GO Biological Processes', fontweight='bold', fontsize=14)
plt.title('RNA-seq:High v Low social stress embryonic transcriptome', fontweight='bold', fontsize=16, loc='center')
# RNA-seq Testes temperature v control GO Biological processes
# GO Biological processes'
# KEGG pathways

# Set x-ticks and labels to go_terms
plt.xticks(np.arange(2, 8, 2))
plt.yticks(np.arange(len(descriptions)), descriptions, fontweight='bold')

# Set x-ticks and labels to go_terms
#min_go_term = min(go_termsFE)
#max_go_term = max(go_termsFE)
#plt.xticks(np.arange(min_go_term, max_go_term + 1, 1))  # Adjust the step size as needed


# Create a colorbar
colorbar = plt.colorbar(sc, label='FDR', orientation='horizontal')
colorbar.ax.xaxis.set_label_position('bottom')

# Create a legend with bubble sizes
sizes = [10, 20, 30] # Adjust the sizes as needed #adjust s=size to make bubbles fit in box


legend_labels = [
    plt.scatter([], [], s=size * 50, c='gray', alpha=0.5, label=f'{size} genes') for size in sizes #adjust s=size to make bubbles fit in box
]

# Display the legend with the plot
legend = plt.legend(handles=legend_labels, loc='best', bbox_to_anchor=(1, 1), title='NGenes')
legend.get_title().set_fontweight('bold')

# Show the plot
plt.grid(True, alpha=0.25)
plt.tight_layout()
plt.savefig('bubble_plot_GO_BP_socstressRNASEQ.png', dpi=600)
plt.show()

