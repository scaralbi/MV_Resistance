from pycirclize import Circos
from pycirclize.parser import Gff
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import argparse
import matplotlib.gridspec as gridspec
from matplotlib.font_manager import FontProperties
from matplotlib.table import Table
import matplotlib as mpl


# Set the default font to 'Helvetica'
plt.rcParams['font.family'] = 'Helvetica'


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

    # We need to iterate until no more adjustments are needed,
    # because an adjustment can cause a position to encroach on its next neighbor.
    # If we make any adjustments, we'll run the loop again.
    if new_positions != positions:
        new_positions = adjust_positions(new_positions, min_distance)

    return new_positions


def main():
    parser = argparse.ArgumentParser(description='Draw a circular diagram of a genome with mutations.')
    parser.add_argument('-strain', type=str, required=True, help='Strain to visualize.')
    parser.add_argument('-ref', type=str, required=True, help='Reference genome.')
    args = parser.parse_args()

    #Parameters
    Strain = args.strain #the strain from which the DNA sequencing reads are obtained (has to be the same name as found in the csv files)
    Reference = args.ref #the accession number of the reference geome sequencing the reads have been aligned to (no underscores)
    coverage_threshold = 200#only highlight the bad coverage regions (below the desired coverage threshol) otherwise it takes too long to draw
    identity_threshold = 0.99#only highlight the unconserved genomic regions(below the desired identity threshold) otherwise it takes too long to draw

    # Load your Coverage data from the CSV file
    data = pd.read_csv(f"Data/reads_alignments/{Strain}_{Reference}.csv")
    coverage_data = data[['Position', 'Coverage']]
    identity_data = data[['Position', 'Identity']]
    #FIlter coverage and identiy to highlight only problematic regions
    filtered_coverage = coverage_data[coverage_data['Coverage'] < coverage_threshold]
    filtered_identity = identity_data[identity_data['Identity'] < identity_threshold]

    print("Coverage data loaded")

    # Load mutation data
    mutations = pd.read_csv(f"Data/variant_analysis/{Reference}_allmutations.csv")
    strain_specific_mutations = mutations.loc[mutations['Strain'] == Strain]
    print("Mutations succesfully loaded")

    # Correct the position data
    data["Position"] = data["Position"] - 1
    strain_specific_mutations["Minimum"] = strain_specific_mutations["Minimum"] - 1

    # Load Genbank file
    gbk = Gff(f"Data/genome_annotations/{Reference}.gff")
    print("Genome annotation file loaded")
    circos = Circos(sectors={gbk.name: gbk.range_size})

    sector = circos.get_sector(gbk.name)

    # Check if there are any positions in the data greater than the genome size
    excess_positions = data[data["Position"] >= gbk.range_size]

    if not excess_positions.empty:
        print("Positions outside range:\n", excess_positions)
    else:
        print("All positions are within the genome range.")

    # Plot outer track with xticks
    major_ticks_interval = 1000000
    minor_ticks_interval = 100000
    outer_track = sector.add_track((62, 64))
    outer_track.axis(fc="lightgrey")
    outer_track.xticks_by_interval(
        major_ticks_interval, label_formatter=lambda v: f"{v/ 10 ** 6:.1f} Mb"
    )
    outer_track.xticks_by_interval(minor_ticks_interval, tick_length=0.5, show_label=False)

    # Plot Forward CDS, Reverse CDS, rRNA, tRNA
    print("Drawing genomic annotations...")

    ref_track = sector.add_track((62, 64))
    ref_track.rect(1, gbk.range_size, fc="lightgrey", ec="black", alpha=0.7)

    f_cds_track = sector.add_track((54, 61), r_pad_ratio=0.1)
    f_cds_track.genomic_features(gbk.extract_features("CDS", target_strand=1), fc="tomato")

    r_cds_track = sector.add_track((47, 54), r_pad_ratio=0.1)
    r_cds_track.genomic_features(gbk.extract_features("CDS", target_strand=-1), fc="skyblue")

    rrna_track = sector.add_track((40, 47), r_pad_ratio=0.1)
    rrna_track.genomic_features(gbk.extract_features("rRNA"), fc="magenta")

    trna_track = sector.add_track((40, 47), r_pad_ratio=0.1)
    trna_track.genomic_features(gbk.extract_features("tRNA"), color="green", lw=0.1)


    # # Add the coverage track
    coverage_track = sector.add_track((71, 75), r_pad_ratio=0.1)
    coverage_track.rect(1, gbk.range_size, fc="white", ec="black")
    coverage_track.fill_between(filtered_coverage["Position"].values, np.log(filtered_coverage["Coverage"].values + 1), fc="plum", ec="white", lw=0.1) 

    # # Add the identity track
    identity_track = sector.add_track((75, 79), r_pad_ratio=0.1)
    identity_track.rect(1, gbk.range_size, fc="seagreen", ec="white", alpha=0.5)
    identity_track.fill_between(filtered_identity["Position"].values, filtered_identity["Identity"].values, fc="seagreen", ec="white", lw=0.1) 


    # Add the outermost track with labels for mutations
    label_track = sector.add_track((80, 82))  # Adjust the radial range as needed
    label_track.rect(1, gbk.range_size, fc="black", ec="mediumseagreen", alpha=1)

    print("Locating mutations...")

    grouped_mutations = strain_specific_mutations[strain_specific_mutations['Protein Effect'].notna() & strain_specific_mutations['Protein Effect'].ne('')].groupby(['locus_tag', 'Minimum', 'Change']).first().reset_index()



    # Adding mutation numbers to labels and creating a dictionary to store mutation details
    mutation_details = {}
    new_labels = []

    for i, (index, row) in enumerate(grouped_mutations.iterrows(), start=1):
        #label = f"{i}) {row['locus_tag']}, {row['Protein Effect']}" if row['Protein Effect'] not in [None, ''] else ''
        label = f"{i}) {row['locus_tag']}" if row['Protein Effect'] not in [None, ''] else ''
        label_table = f"{i}) {row['locus_tag']}"
        new_labels.append(label)
        # Store mutation details for the table
        mutation_details[i] = {
            'Gene': label_table,
            'Product': row['product'],
            'Protein\nEffect': row['Protein Effect'],
            'Amino Acid\nChange': row['Amino Acid Change'],
            'Variant\nFrequency': row['Variant Frequency']
        }

    # Use new_labels instead of labels
    pos_list = grouped_mutations["Minimum"].values.tolist()

    #Add mutation numbers at each location onto reads track
    for i in range(len(grouped_mutations)):
        label_track.text('-', grouped_mutations['Minimum'].iloc[i], orientation="vertical", color="white", size=22)

    # Adjust positions to ensure minimum distance of usually works well with 25000
    pos_list_CDSonly = adjust_positions(sorted(pos_list), 25000)

    # Plot outer xticks (labels)
    label_track.xticks(
        pos_list_CDSonly,
        new_labels,
        outer=True,   # Make sure labels are plotted outside
        tick_length=0,  # Set tick_length to 0 so ticks won't show up
        label_margin=0.5,
        label_orientation="vertical",  # Set label orientation to vertical
        label_size= 8,
    )

    fig = circos.plotfig()
    # Add legend
    handles = [
        Line2D([0], [0], color='lightgrey', linewidth=8, label=f"Reference:{Reference}", mec='black', mew=1.5),
        Patch(color="tomato", label="Forward CDS"),
        Patch(color="skyblue", label="Reverse CDS"),
        Patch(color="magenta", label="rRNA"),
        Patch(color="green", label="tRNA"),
        Line2D([0], [0], color='black', linewidth=8, label=f"WGS reads from\nstrain: {Strain}", mec='mediumseagreen', mew=1.5),
        Patch(color="plum", label="Coverage", alpha=1),
        Patch(color="seagreen", label="Identity", alpha=1),
    ]
    legend = fig.legend(handles=handles, loc="center", fontsize=8)

    # Find indices of labels to be bolded
    to_bold = [i for i, text_obj in enumerate(legend.texts) if text_obj.get_text() in [f"{Reference}", f"{Strain}"]]

    # Apply bold font style
    for i in to_bold:
        legend.texts[i].set_weight('bold')

    fig.savefig(f"Figures/WGS/{Strain}_vs_{Reference}_genomeview.png", dpi=900)
    print("Genome Figure saved succesfully")

    ## Draw Table
    headers = list(mutation_details[1].keys())  # Convert dict_keys to list
    cell_text = [list(d.values()) for d in mutation_details.values()]  # Get cell values

    print("Printing table...")

    fig = plt.figure(figsize=(12,10), dpi=300)
    ax = plt.subplot()

    # Determine maximum text length in each column and use it as width
    column_widths = [
        max(len(str(x)) for x in column) *1
        for column in zip(*([headers] + cell_text))  # Adding header to cell_text for each column
    ]

    # Normalize widths (to make sure they sum up to the number of columns)
    column_widths = [w / sum(column_widths) * len(headers) for w in column_widths]

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
            xy=(positions[index] + 0.5 * column_widths[index], nrows + 1),  # Adjust Y position for header text to align it in the center of the bounding box - for 20 mutation change 1 to 0.5
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
        f"Figures/WGS/Table_{Strain}_{Reference}.png",
        dpi=900,
        transparent=False,
        bbox_inches='tight'
    )

    print(f"Table Saved succesfully as: Figures/WGS/Table_{Strain}_{Reference}.png")


if __name__ == "__main__":
    main()


