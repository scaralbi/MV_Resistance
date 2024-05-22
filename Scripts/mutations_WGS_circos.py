from Bio import SeqIO
from pycirclize import Circos
from pycirclize.parser import Gff
from pycirclize.utils import ColorCycler
from matplotlib.patches import Patch
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


# Load the .xlsx file
xls = pd.ExcelFile("WGS_SNPs.xlsx")
# Load mutation data
mutations = pd.read_csv("NC000911_filtered.csv")
# Get the names of the sheets
sheet_names = xls.sheet_names

# Identify the sheets that end with '_all'
all_sheets = [name for name in sheet_names if name.endswith('_filtered')]

# Load the data from these sheets into a dictionary
all_data = {sheet: pd.read_excel(xls, sheet) for sheet in all_sheets}

# Filter the data and create a new column 'Mutation'
mutation_number = 1
mutation_numbers = {}
for sheet, data in all_data.items():
    data = data[(data['CDS'] != 'non coding') & (data['Protein Effect'].notna())]
    data['Mutation'] = data['CDS'] + " (" + data['Codon Change'].astype(str) + " -> " + data['Amino Acid Change'].astype(str) + ")"
    data['Mutation Number'] = 0
    for i, row in data.iterrows():
        key = (row['Minimum'], row['CDS'], row['Amino Acid Change'], row['Protein Effect'])
        if key not in mutation_numbers:
            mutation_numbers[key] = mutation_number
            mutation_number += 1
        data.loc[i, 'Mutation Number'] = mutation_numbers[key]
    all_data[sheet] = data

# Load the .gff file
gbk = Gff("/Data/genome_annotations/NC_000911.gff")



# Get the data for the current sheet
data_NC00011_all = all_data['NC000911_filtered']

# Initialize the Circos plot
circos = Circos(sectors={gbk.name: gbk.range_size})
sector = circos.get_sector(gbk.name)

# Plot outer track with xticks
major_ticks_interval = 1000000
minor_ticks_interval = 100000
outer_track = sector.add_track((75, 77))
outer_track.axis(fc="lightgrey")
outer_track.xticks_by_interval(
    major_ticks_interval, label_formatter=lambda v: f"{v/ 10 ** 6:.1f} Mb"
)
outer_track.xticks_by_interval(minor_ticks_interval, tick_length=1, show_label=False)


# Extract unique strains from the first dataset
strains_NC00011_all = data_NC00011_all['Strain'].unique()


# Define colors for each strain
colors = ['tomato', 'skyblue', 'magenta', 'green', 'purple', 'brown', 'orange']
strain_color = {strain: color for strain, color in zip(strains_NC00011_all, colors)}


# Plot a track for each strain
for i, strain in enumerate(strains_NC00011_all):
    strain_data = data_NC00011_all[data_NC00011_all['Strain'] == strain]
    radial_range = (69 - i*6, 73 - i*6)
    strain_track = sector.add_track(radial_range, r_pad_ratio=0.1)
    strain_track.axis(fc=strain_color[strain],alpha=0.4)
    for index, mutation in strain_data.iterrows():
        strain_track.rect(mutation['Minimum']-25000, mutation['Minimum']+25000, fc="white", ec=strain_color[strain], alpha=0.8, lw=1)
        #strain_track.scatter([mutation['Minimum']], [50], color="black", s=10, marker=".", linewidths =0.5, ec="black")
        strain_track.text(str(mutation['Mutation Number']), mutation['Minimum'], orientation="vertical", color="black", size=8)

 
# Add the outermost track with labels for mutations
label_track = sector.add_track((82, 85))  # Adjust the radial range as needed


# Group mutations by 'CDS' and select the first row from each group
grouped_mutations = mutations.groupby('product').first().reset_index()

# Now create the labels and pos_list from this grouped data
labels = [f"{row['product']}" for _, row in grouped_mutations.iterrows()]
pos_list = grouped_mutations["Minimum"].values.tolist()


# Plot outer xticks (labels)
label_track.xticks(
    pos_list,
    labels,
    outer=True,   # Make sure labels are plotted outside
    tick_length=1,  # Set tick_length to 0 so ticks won't show up
    label_margin=1,
    label_orientation="vertical",  # Set label orientation to vertical
    label_size= 8,
)

# Create a new figure for the table
fig, ax = plt.subplots(figsize=(10, 8))

# Prepare the table data
table_data = []
for i, (index, row) in enumerate(data_NC00011_all.iterrows()):
    table_data.append([row['Mutation Number'], row['CDS'], row['Codon Change'], row['Amino Acid Change'], row['Variant Frequency']])

# Create the table
table = ax.table(cellText=table_data, colLabels=["Mutation Number", "CDS", "Codon Change", "Amino Acid Change", "Variant Frequency"], loc='center')

# Hide axes
ax.axis('off')

# Save the Circos plot
fig_circos = circos.plotfig()
_ = circos.ax.legend(
    handles=[
        Patch(color="gray", label="NC_000911", alpha=0.6, edgecolor="black"),
        Patch(color="tomato", label="mvR1", alpha=0.6, edgecolor="black"),
        Patch(color="skyblue", label="mvR2", alpha=0.6, edgecolor="black"),
        Patch(color="magenta", label="mvR3", alpha=0.6, edgecolor="black"),
        Patch(color="green", label="mvR6", alpha=0.6, edgecolor="black"),
        Patch(color="purple", label="mvR9", alpha=0.6, edgecolor="black"),

    ],
    bbox_to_anchor=(0.5, 0.5),
    loc="center",
    ncols=1,
    fontsize=14,
)


fig_circos.savefig("NC00011_all_circos_plot.png", dpi=600)

# Save the table figure
fig.savefig("NC00011_all_mutation_table.png", dpi=600)

# Display the table figure
plt.show()
