import pandas as pd
from Bio import SeqIO
import numpy as np
from pycirclize import Circos
from pycirclize.parser import Gff
from pycirclize.utils import ColorCycler
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import os
import re

    
def wrap_text(text, max_length):
    """
    Wrap the text every max_length characters.

    Parameters:
        text (str): The text to wrap.
        max_length (int): The maximum line length.

    Returns:
        wrapped_text (str): The wrapped text.
    """
    return '\n'.join(text[i:i+max_length] for i in range(0, len(text), max_length))



def min_percentage(s):
    """Return the smallest percentage value in a string."""
    if pd.isna(s):
        return 0
    if not isinstance(s, str):
        s = str(s)
    # Extract all values from the string
    values = re.findall(r"(\d+(?:\.\d+)?%?)", s)
    
    # Convert to float, and if the value was a percentage, divide by 100
    floats = [float(x.rstrip('%')) if '%' in x else float(x) * 100 for x in values]
    
    return min(floats) if floats else 0


def adjust_positions(positions, min_distance):
    """
    Adjust positions to ensure a minimum distance between each position so that mutated gene labels do not overlap

    Parameters:
        positions (list): List of positions, sorted in ascending order
        min_distance (int): Minimum distance that should be between each position

    Returns:
        new_positions (list): List of adjusted positions
    """
    new_positions = positions.copy()
    for i in range(len(new_positions) - 1):
        # Check if distance to next position is less than min_distance
        if new_positions[i + 1] - new_positions[i] < min_distance:
            # Adjust next position
            new_positions[i + 1] = new_positions[i] + min_distance

    # need to iterate until no more adjustments are needed,
    # because an adjustment can cause a position to encroach on its next neighbor.
    # If any adjustment is made, then the loop is run again.
    if new_positions != positions:
        new_positions = adjust_positions(new_positions, min_distance)

    return new_positions

def plot_mutation_circos(ref):

    # Load mutation data
    all_data = pd.read_csv(f"../Data/variant_analysis/{ref}_filtered_mutations.csv")

    # Load the .gff file
    gbk = Gff(f"../Data/genome_annotations/{ref}.gff")

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
    all_strains = all_data['Strain'].unique()
    print(all_strains)

    # Define colors for each strain
    colors = ['tomato', 'skyblue', 'magenta', 'green', 'purple', 'brown', 'orange', 'lightgreen', 'cadetblue']
    strain_color = {strain: color for strain, color in zip(all_strains, colors)}

    # Group mutations by 'CDS' and select the first row from each group
    # grouped_mutations = all_data.groupby('product').first().reset_index()

    grouped_mutations = all_data[all_data['Protein Effect'].notna() & all_data['Protein Effect'].ne('')].groupby(['locus_tag', 'Minimum', 'Change']).first().reset_index()

    # Create a mapping of strains and frequencies per mutation
    # mutation_strain_map = all_data.groupby(['locus_tag', 'Minimum', 'Change']).agg({'Strain': lambda x: ', '.join(x), 'Variant Frequency': lambda x: ', '.join(x)}).to_dict(orient='index')

    mutation_strain_map = all_data.groupby(['locus_tag', 'Minimum', 'Change']).agg({
        'Strain': lambda x: ', '.join(x), 
        'Variant Frequency': lambda x: ', '.join(map(str, x))
    }).to_dict(orient='index')



    # Group by 'locus_tag' and get the first entry for each unique gene
    unique_grouped_mutations = grouped_mutations.groupby('locus_tag').first()

    labels = [f"{row['product']}" for _, row in unique_grouped_mutations.iterrows()]
    pos_list = unique_grouped_mutations["Minimum"].values.tolist()

    # Initialize new_labels and pos_list_CDSonly as empty lists

    new_labels = []
    pos_list_CDSonly = []

    # Step 1: Define mutation_details before the loop
    mutation_details = []

    for i, (index, row) in enumerate(grouped_mutations.iterrows(), start=1):
        freq = min_percentage(row['Variant Frequency'])
        if freq >= 75:
            label = f"{i}) {row['locus_tag']}" if row['Protein Effect'] not in [None, ''] else ''
            label_table = f"{i}) {row['locus_tag']}"
            new_labels.append(label)
            pos_list_CDSonly.append(row['Minimum'])  # Add this position to pos_list_CDSonly
            # Append mutation details to the dataframe
            # Step 2: Append mutation details to mutation_details
            amino_acid_change = row['Amino Acid Change']
            if not amino_acid_change:
                amino_acid_change = "fs"

            mutation_strain_freq = mutation_strain_map.get((row['locus_tag'], row['Minimum'], row['Change']), {})
            mutation_details.append({
                "Mutation": label_table,
                "Gene": row['product'],
                "Position": row['Minimum'],
                "Change": amino_acid_change,
                "Strain": wrap_text(mutation_strain_freq.get('Strain', ''), 15),
                "Variant Frequency": wrap_text(mutation_strain_freq.get('Variant Frequency', ''), 15)
            })
            
        else:
            new_labels.append('')
            pos_list_CDSonly.append(row['Minimum'])  # Add this position to pos_list_CDSonly


    # Plot a track for each strain
    for i, strain in enumerate(all_strains):
        strain_data = all_data[all_data['Strain'] == strain]
        radial_range = (69 - i*5, 73 - i*5)
        strain_track = sector.add_track(radial_range, r_pad_ratio=0.1)
        strain_track.axis(fc=strain_color[strain],alpha=0.4)
        for index, mutation in strain_data.iterrows():
            strain_track.rect(mutation['Minimum']-0, mutation['Minimum']+1000, fc="white", ec=strain_color[strain], alpha=0.8, lw=1)
            #strain_track.scatter([mutation['Minimum']
            # ], [50], color="black", s=10, marker=".", linewidths =0.5, ec="black")

    
    # Add the outermost track with labels for mutations
    label_track = sector.add_track((82, 85))  # Adjust the radial range as needed

    # Use new_labels instead of labels
    pos_list = grouped_mutations["Minimum"].values.tolist()

    # #Add mutation numbers at each location onto reads track
    # for i in range(len(grouped_mutations)):
    #     label_track.text('.', grouped_mutations['Minimum'].iloc[i], orientation="vertical", color="black", size=22)

    # Adjust positions to ensure minimum distance of usually works well with 25000
#    pos_list_CDSonly = adjust_positions(sorted(pos_list), 500)

    # Plot outer xticks (labels)
    label_track.xticks(
        pos_list,
        new_labels,
        outer=True,   # Make sure labels are plotted outside
        tick_length=0,  # Set tick_length to 0 so ticks won't show up
        label_margin=1,
        label_orientation="vertical",  # Set label orientation to vertical
        label_size= 10,
    )


    # Save the Circos plot
    fig_circos = circos.plotfig()
    _ = circos.ax.legend(
        handles=[
            Patch(color="gray", label=f"{ref}", alpha=0.6, edgecolor="black"),
            Patch(color="tomato", label="mvR1", alpha=0.6, edgecolor="black"),
            Patch(color="skyblue", label="mvR2", alpha=0.6, edgecolor="black"),
            Patch(color="magenta", label="mvR3", alpha=0.6, edgecolor="black"),
            Patch(color="green", label="mvR6", alpha=0.6, edgecolor="black"),
            Patch(color="purple", label="mvR9", alpha=0.6, edgecolor="black"),
            Patch(color="brown", label="mvR10", alpha=0.6, edgecolor="black"),
            Patch(color="orange", label="mvR11", alpha=0.6, edgecolor="black"),
            Patch(color="lightgreen", label="mvR12", alpha=0.6, edgecolor="black"),
        ],
        bbox_to_anchor=(0.5, 0.5),
        loc="center",
        ncols=1,
        fontsize=10,
    )

    fig_circos.savefig(f"../Figures/WGS/{ref}_mutants_vs_WT_genome_views_New3.svg", dpi=600)


    # Step 3: Use mutation_details to generate the table
    headers = list(mutation_details[1].keys())  # Convert dict_keys to list
    cell_text = [list(d.values()) for d in mutation_details]  # Get cell values
    
    fig = plt.figure(figsize=(15,10), dpi=300)
    ax = plt.subplot()

    # Determine maximum text length in each column and use it as width
    column_widths = [
        max(len(str(x)) for x in column) * 8 if header in ["Amino Acid Change", "product"] else max(len(str(x)) for x in column) * 1
        for column, header in zip(zip(*([headers] + cell_text)), headers)  # Adding header to cell_text for each column
    ]

    # Calculate positions based on column widths
    positions = np.cumsum([0] + column_widths)  # Here we include the rightmost border

    ncols = len(headers)
    nrows = len(cell_text)

    print(f"Mutations found:{nrows}")

    # Add 1 for headers and adjust limits according to column widths
    ax.set_xlim(0, sum(column_widths))
    ax.set_ylim(0, nrows + 1)

    # Add table's main text
    for i in range(nrows):
        for j, text in enumerate(cell_text[i]):
            ax.annotate(
                xy=(positions[j] + 0.5 * column_widths[j], nrows - i - 0.5),  # Adjust Y position for text to align it in the center of the bounding box
                text=text,
                ha='center',
                va='center',       
                fontsize=10  # Make header text larger
            )

    # Add column names

    for index, header in enumerate(headers):
        ax.annotate(
            xy=(positions[index] + 0.5 * column_widths[index], nrows + 0.5),  # Adjust Y position for header text to align it in the center of the bounding box - for 20 mutation change 1 to 0.5
            text=header,
            ha='center',
            va='center',
            weight='heavy',  # Make header text bold
            fontsize=12  # Make header text larger
        )

    # Add dividing lines
    ax.plot([ax.get_xlim()[0], ax.get_xlim()[1]], [nrows, nrows], lw=1.5, color='black', marker='', zorder=4)
    ax.plot([ax.get_xlim()[0], ax.get_xlim()[1]], [0, 0], lw=1.5, color='black', marker='', zorder=4)
    for x in range(1, nrows):
        ax.plot([ax.get_xlim()[0], ax.get_xlim()[1]], [x, x], lw=1.15, color='gray', ls='--', zorder=3 , marker='')

    # Add vertical grid lines but exclude the outermost ones
    for j in range(1, ncols):  # Excluding outermost vertical lines
        ax.plot([positions[j], positions[j]], [0, nrows], color='gray', linewidth=1.15, ls='--', zorder=3 , marker='')

    ax.set_axis_off()
    plt.savefig(
        f"..Data/Figures/WGS/Table_Mutants2_{ref}.png",
        dpi=900,
        transparent=False,
        bbox_inches='tight'
    )
    print(f"Table Saved succesfully as: ../Figures/WGS/Table2_{ref}.png")

    return

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Filter mutations and plot a circos diagram.')
    parser.add_argument('-ref', help='The reference genome.')
    args = parser.parse_args()

    plot_mutation_circos(args.ref)


