import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import seaborn as sns
import matplotlib.colors as mcolors
import matplotlib.colorbar as mcolorbar
import scipy.stats as stats

mpl.rcParams['font.family'] = 'Helvetica'
# Load the CSV file into a DataFrame
df = pd.read_csv('../Data/variant_analysis/NZCP073017_allmutations.csv')

# Process the 'Variant Frequency' column to ensure it is correctly formatted and numeric
df["Variant Frequency"] = df["Variant Frequency"].str.split(',').str[0].str.replace('%', '').astype(float) / 100

# Drop rows with NaN values in 'Variant Frequency' column
df = df.dropna(subset=['Variant Frequency'])

# Ensure 'Variant Frequency' column is numeric (already ensured by conversion above)
df['Variant Frequency'] = pd.to_numeric(df['Variant Frequency'], errors='coerce')

# Drop rows where 'Variant Frequency' couldn't be converted to a number
df = df.dropna(subset=['Variant Frequency'])

# Remove duplicate rows
df = df.drop_duplicates()
# Define bins for the histogram
bins = [0, 0.33333333333, 0.66666666667, 1.0]

# Define colors for different strains
colors = sns.color_palette("tab10", len(df['Strain'].unique()))
colormap = sns.color_palette("coolwarm", len(bins)-1)
# Define the colors for different strains
wt_colors = sns.dark_palette("Blue", 2)
mvR_colors = sns.light_palette("Red", 8)


# Define the order of strains
strain_order = ['wt_Nixon', 'mvR1_Nixon', 'mvR2_Nixon', 'mvR3_Nixon', 'mvR6_Nixon',
                'wt_Howe', 'mvR9_Howe', 'mvR10_Howe', 'mvR11_Howe', 'mvR12_Howe']

# Ensure the strains are ordered correctly in the DataFrame
df['Strain'] = pd.Categorical(df['Strain'], categories=strain_order, ordered=True)
df = df.sort_values('Strain')


# Calculate the fractions of low to high frequency variants for each strain
fraction_data = []
p_values = []

for strain in strain_order:
    strain_data = df[df['Strain'] == strain]['Variant Frequency']
    low_freq_count = (strain_data < 0.33333333333).sum()
    high_freq_count = (strain_data > 0.66666666667).sum()    
    total_freq_count = len(strain_data)

    if high_freq_count == 0:  # To avoid division by zero
        ratio = np.nan
    else:
        ratio = low_freq_count * 100 / total_freq_count
        
    fraction_data.append(ratio)

    # Perform Mann-Whitney U test for low vs high frequency
    low_freq = strain_data[strain_data <= 0.33333333333]
    high_freq = strain_data[strain_data >= 0.66666666667]
    if len(low_freq) > 0 and len(high_freq) > 0:
        u_stat, p_value = stats.mannwhitneyu(low_freq, high_freq)
    else:
        p_value = np.nan  # If either category is empty, set p-value to NaN
    p_values.append(p_value)




# Create a figure and axis object with adjusted height ratios
fig, axes = plt.subplots(3, 5, figsize=(10, 9), gridspec_kw={'height_ratios': [1, 1, 1]})
axes = axes.flatten()

# Plot histograms for each strain
for i, strain in enumerate(strain_order):
    ax = axes[i]
    strain_data = df[df['Strain'] == strain]
    
    # Plot histogram and get the patches
    n, bins, patches = ax.hist(strain_data["Variant Frequency"], bins=bins, alpha=0.7, edgecolor='black', linewidth=1)
    
    # Color the bars according to the colormap
    for j in range(len(patches)):
        patches[j].set_facecolor(colormap[j])
    
    ax.set_title(strain)
    ax.set_xlabel('Variant frequency')
    ax.set_ylabel('Number of variants')
    ax.set_ylim(0, 140)

# Remove the extra axes
for i in range(len(strain_order), len(axes)):
    fig.delaxes(axes[i])

# Create the strain bar plot in the last row spanning all columns
strainbar_ax = fig.add_subplot(3, 1, 3)

bar_colors = []
for strain in strain_order:
    if strain.startswith('wt'):
        bar_colors.append(wt_colors[0]) if strain == 'wt_Nixon' else wt_colors[1]
    else:
        bar_colors.append(mvR_colors[0]) if strain.startswith('mvR1') else mvR_colors[1]

strainbar_ax.bar(strain_order, fraction_data, color=bar_colors, edgecolor='black', linewidth=1)
strainbar_ax.set_xlabel('Strain')
strainbar_ax.set_ylabel('Low frequency variants (%)')
strainbar_ax.set_xticks(np.arange(len(strain_order)))  # Set the x-axis ticks
strainbar_ax.set_xticklabels(strain_order, rotation=45)  # Set the x-axis tick labels with rotation
strainbar_ax.set_ylim(30, 70)  # Adjust the y-axis limits

# Remove unwanted ticks and labels on both axes
strainbar_ax.tick_params(axis='x', which='both', bottom=True, top=False)
strainbar_ax.tick_params(axis='y', which='both', left=True, right=False)
strainbar_ax.spines['top'].set_visible(False)
strainbar_ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig("../Figures/WGS/variant_frequency_fraction_comparison.png", dpi=600, transparent=False, bbox_inches='tight')
plt.show()