# Title: Barchart plotting
# Author: Dr. Alice M. Godden

import matplotlib.pyplot as plt
import pandas as pd

# Data
data = {
    "miRNA": [
        "dre-miR-129-5p",
        "dre-miR-184",
        "dre-miR-181a-5-3p",
        "dre-miR-183-5p",
        "dre-miR-193b-5p",
        "dre-miR-10b-5p",
        "dre-miR-200a-5p"
    ],
    "num_matching_sigDE_genes": [140, 30, 22, 72, 27, 35, 47]
}

# Create DataFrame
df = pd.DataFrame(data)

# Plot
plt.figure(figsize=(10, 6))

# Use magma colormap for colors
colors = plt.cm.magma(df['num_matching_sigDE_genes'] / max(df['num_matching_sigDE_genes']))

# Bar plot with colored bars
bars = plt.bar(df['miRNA'], df['num_matching_sigDE_genes'], color=colors, zorder=2)

# Adding gridlines
plt.grid(alpha=0.25, zorder=0)

plt.xlabel('miRNA', fontweight='bold')
plt.ylabel('Num of matching sig differentially expressed genes', fontweight='bold')
plt.title('Number of matching significantly differentially expressed genes', fontweight='bold')
plt.xticks(rotation=45, ha='right', fontweight='bold')
plt.yticks(fontweight='bold')
plt.tight_layout()

# Save the plot with 600 ppi resolution
plt.savefig("mirna_bar_chart_600ppi.png", dpi=600)

plt.show()
