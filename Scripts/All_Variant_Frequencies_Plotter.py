import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import seaborn as sns
import matplotlib.colors as mcolors
import matplotlib.colorbar as mcolorbar
import scipy.stats as stats
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
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
    low_freq_count = (strain_data < 0.25).sum()
    high_freq_count = (strain_data >= 0.25).sum()
    
    if high_freq_count == 0:  # To avoid division by zero
        ratio = np.nan
    else:
        ratio = low_freq_count / high_freq_count
        
    fraction_data.append(ratio)

    # Perform Mann-Whitney U test for low vs high frequency
    low_freq = strain_data[strain_data < 0.25]
    high_freq = strain_data[strain_data >= 0.25]
    if len(low_freq) > 0 and len(high_freq) > 0:
        u_stat, p_value = stats.mannwhitneyu(low_freq, high_freq)
    else:
        p_value = np.nan  # If either category is empty, set p-value to NaN
    p_values.append(p_value)

# Plotting
plt.figure(figsize=(12, 6))

bar_colors = sns.color_palette("tab10", len(fraction_data))
bars = plt.bar(strain_order, fraction_data, color=bar_colors)

plt.xlabel('Strain')
plt.ylabel('Fraction of low to high frequency variants')
plt.title('Fraction of low to high frequency variants across strains')
plt.xticks(rotation=45)
plt.tight_layout()

plt.savefig("../Figures/WGS/variant_frequency_fraction_comparison.png", dpi=600, transparent=False, bbox_inches='tight')
plt.show()
