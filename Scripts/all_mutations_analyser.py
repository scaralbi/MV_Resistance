import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import seaborn as sns
import matplotlib.colors as mcolors
import matplotlib.colorbar as mcolorbar
import scipy.stats as stats
import matplotlib.gridspec as gridspec
import string


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

#bins = [0, 0.25, 0.5, 0.75, 1.0]
# bins = [0, 0.2, 0.4, 0.6, 0.8, 1.0]

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
lowfraction_data = []
highfraction_data = []
low_threshold = 0.33
high_threshold = 0.66

p_values = []

for strain in strain_order:
    strain_data = df[df['Strain'] == strain]['Variant Frequency']
    # low_freq_count = (strain_data < 0.33333333333).sum()
    # high_freq_count = (strain_data > 0.66666666667).sum()    
    low_freq_count = (strain_data <= low_threshold).sum()
    high_freq_count = (strain_data > high_threshold).sum()  
    total_freq_count = len(strain_data)

    if high_freq_count == 0:  # To avoid division by zero
        ratio = np.nan
    else:
        low_ratio = low_freq_count * 100 / total_freq_count
        high_ratio = high_freq_count * 100 / total_freq_count

        
    lowfraction_data.append(low_ratio)
    highfraction_data.append(high_ratio)

    # # Perform Mann-Whitney U test for low vs high frequency
    # low_freq = strain_data[strain_data <= 0.33333333333]
    # high_freq = strain_data[strain_data >= 0.66666666667]
    # if len(low_freq) > 0 and len(high_freq) > 0:
    #     u_stat, p_value = stats.mannwhitneyu(low_freq, high_freq)
    # else:
    #     p_value = np.nan  # If either category is empty, set p-value to NaN
    # p_values.append(p_value)


# Create a GridSpec object
gs = gridspec.GridSpec(3, 10)

# Create a figure
fig = plt.figure(figsize=(10, 8))

# Create axes for the histograms
axes = [plt.subplot(gs[i // 5, 2*(i % 5):2*(i % 5 + 1)]) for i in range(10)]  # Adjust this line

# Plot histograms for each strain
for i, strain in enumerate(strain_order):
    ax = axes[i]
    strain_data = df[df['Strain'] == strain]
    
    # Plot histogram and get the patches
    n, bins, patches = ax.hist(strain_data["Variant Frequency"], bins=bins, alpha=0.7, edgecolor='black', linewidth=1)
    
    # Calculate total number of mutations
    total_mutations = int(n.sum())
    
    # Add total number of mutations to the plot
    ax.text(0.95, 0.95, f'Total = {total_mutations}', transform=ax.transAxes, fontsize=10, ha='right', va='top')
    
    # Color the bars according to the colormap
    for j in range(len(patches)):
        patches[j].set_facecolor(colormap[j])
    
    ax.set_title(strain)
    ax.set_xlabel('Variant frequency')
    ax.set_ylabel('Number of variants')
    ax.set_xticks(bins)
    ax.tick_params(axis='both', labelsize=8)
    ax.set_ylim(0, 130)
 # Add letter to the subplot
    ax.text(-0.2, 1.12, string.ascii_lowercase[i], transform=ax.transAxes, fontsize=12, fontweight='bold', va='top', ha='left')

# Create the strain bar plots in the last row
strainbar_ax1 = plt.subplot(gs[2, :5])
strainbar_ax2 = plt.subplot(gs[2, 5:])
# Add letters to the strain bar plots
strainbar_ax1.text(-0.05, 1.12, string.ascii_lowercase[len(strain_order)], transform=strainbar_ax1.transAxes, fontsize=12, fontweight='bold', va='top', ha='left')
strainbar_ax2.text(-0.05, 1.12, string.ascii_lowercase[len(strain_order) + 1], transform=strainbar_ax2.transAxes, fontsize=12, fontweight='bold', va='top', ha='left')

# Create a list of edge widths
edge_widths = [2 if strain.startswith('wt') else 1 for strain in strain_order]
# Plot the bars and horizontal lines for the low frequency variants
for i, (strain, data) in enumerate(zip(strain_order, lowfraction_data)):
    edge_width = 2 if strain.startswith('wt') else 1
    strainbar_ax1.bar(i, data, color=colormap[0], edgecolor='black', linewidth=edge_width)
    if strain.startswith('wt'):
        strainbar_ax1.plot([i, i + 4.4], [data, data], color='black', linestyle='dashed')

strainbar_ax1.set_xlabel('Strain')
strainbar_ax1.set_ylabel(f'Fraction of total variants \n with frequency < {low_threshold} (%)')
strainbar_ax1.set_xticks(np.arange(len(strain_order)))  # Set the x-axis ticks
strainbar_ax1.set_xticklabels(strain_order, rotation=45, ha='right')  # Set the x-axis tick labels with rotation
strainbar_ax1.set_title('Low frequency variants')

# strainbar_ax.set_ylim(30, 70)  # Adjust the y-axis limits

# Remove unwanted ticks and labels on both axes
strainbar_ax1.tick_params(axis='x', which='both', bottom=True, top=False)
strainbar_ax1.tick_params(axis='y', which='both', left=True, right=False)
strainbar_ax1.spines['top'].set_visible(False)
strainbar_ax1.spines['right'].set_visible(False)


# Plot the bars and horizontal lines for the high frequency variants
for i, (strain, data) in enumerate(zip(strain_order, highfraction_data)):
    edge_width = 2 if strain.startswith('wt') else 1
    strainbar_ax2.bar(i, data, color=colormap[-1], edgecolor='black', linewidth=edge_width)
    if strain.startswith('wt'):
        strainbar_ax2.plot([i, i + 4.4], [data, data], color='black', linestyle='dashed')
        
strainbar_ax2.set_xlabel('Strain')
strainbar_ax2.set_ylabel(f'Fraction of total variants \n with frequency > {high_threshold} (%)')
strainbar_ax2.set_xticks(np.arange(len(strain_order)))  # Set the x-axis ticks
strainbar_ax2.set_xticklabels(strain_order, rotation=45, ha='right')  # Set the x-axis tick labels with rotation
strainbar_ax2.set_title('High frequency variants')

# strainbar_ax.set_ylim(30, 70)  # Adjust the y-axis limits

# Remove unwanted ticks and labels on both axes
strainbar_ax2.tick_params(axis='x', which='both', bottom=True, top=False)
strainbar_ax2.tick_params(axis='y', which='both', left=True, right=False)
strainbar_ax2.spines['top'].set_visible(False)
strainbar_ax2.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig("../Figures/WGS/variant_frequency_fraction_comparison_03-06.png", dpi=600, transparent=False, bbox_inches='tight')
plt.show()