import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Load the data
file_path = "Data/variant_analysis/sll1180.xlsx"
reloaded_data = pd.read_excel(file_path, header=None)
reloaded_data.columns = reloaded_data.iloc[0]
reloaded_data.drop(0, inplace=True)
reloaded_data.reset_index(drop=True, inplace=True)

# Warmer colors for the nucleotides
warmer_colors = {
    'A': '#FFC125',
    'T': '#4575B4',
    'G': '#D73027',
    'C': '#91BFDB',
}

locus_dataframes = []
coverage_values_all_loci = []

# Iterate over all loci columns in steps of 5
for i in range(0, reloaded_data.shape[1]-1, 5):
    # Extract data for the current locus
    locus_data = reloaded_data.iloc[:, i:i+4]
    
    # Filter out non-numeric rows
    locus_df = locus_data[~locus_data[locus_data.columns[0]].str.contains(r'[a-zA-Z]')].copy()
    
    # Extract percentage values and convert to float
    for col in locus_df.columns:
        locus_df[col] = locus_df[col].str.extract(r'(\d+.\d+)%')[0].astype(float)
    
    # Rename columns to match the nucleotide
    locus_df.columns = reloaded_data.iloc[1, i:i+4].values
    locus_df.reset_index(drop=True, inplace=True)
    
    # Extract and sum the coverage values for the current locus
    coverage_data = reloaded_data.iloc[:, i:i+4][~reloaded_data['locus_id'].str.contains(r'[a-zA-Z]')]
    total_coverage = coverage_data.apply(lambda x: x.str.extractall(r'\((\d+)\)').astype(float).sum().sum(), axis=1)
    
    # Store the processed dataframe and coverage values
    locus_dataframes.append(locus_df)
    coverage_values_all_loci.append(total_coverage.tolist())

# Plotting grouped stacked bar chart for all loci
fig, ax = plt.subplots(figsize=(15, 6))

# Width of a bar 
width = 0.2

# Position of bars on x axis
r = np.arange(len(locus_dataframes[0]))

for i, (locus_df, coverage_values) in enumerate(zip(locus_dataframes, coverage_values_all_loci)):
    # Bottom values for stacking
    bottoms = np.array([0] * len(locus_df))
    
    for column in locus_df.columns[1:]:
        ax.bar(r + i*width, locus_df[column], width=width, edgecolor='black', bottom=bottoms, label=column, color=warmer_colors[column])
        bottoms += locus_df[column].values

    # Adding the total coverage values as text on top of every bar
    for j, rect in enumerate(ax.containers[i]):
        height = rect.get_height()
        width = rect.get_x() + rect.get_width() / 2
        ax.text(width, height + 5, f'Coverage=\n{int(coverage_values[j])}', ha='center', va='bottom', fontsize=8, color='black')

# Setting labels, title, and legend
ax.set_xlabel('Strains', fontweight='bold')
ax.set_ylabel('Percentage', fontweight='bold')
ax.set_xticks(r + 1.5*width)  # Adjust the x-ticks to be centered
ax.set_xticklabels(locus_dataframes[0]['locus_id'])
ax.legend(title="Bases", loc="upper left", bbox_to_anchor=(1,1))
plt.tight_layout()
plt.show()
