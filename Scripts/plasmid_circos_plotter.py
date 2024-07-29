import argparse
import pandas as pd
from pycirclize import Circos
from pycirclize.parser import Gff
import re
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import matplotlib


# Set the font globally to Helvetica
matplotlib.rcParams['font.family'] = 'Helvetica'


def wrap_text(text, max_length):
    """Wrap the text every max_length characters."""
    return '\n'.join(text[i:i+max_length] for i in range(0, len(text), max_length))

def adjust_positions(positions, min_distance):
    """Adjust positions to ensure a minimum distance between each position."""
    new_positions = positions.copy()
    for i in range(len(new_positions) - 1):
        if new_positions[i + 1] - new_positions[i] < min_distance:
            new_positions[i + 1] = new_positions[i] + min_distance
    if new_positions != positions:
        new_positions = adjust_positions(new_positions, min_distance)
    return new_positions


def min_percentage(s):
    """Return the smallest percentage value in a string."""
    if pd.isna(s):
        return 0
    if not isinstance(s, str):
        s = str(s)
    values = re.findall(r"(\d+(?:\.\d+)?%?)", s)
    floats = [float(x.rstrip('%')) if '%' in x else float(x) * 100 for x in values]
    return min(floats) if floats else 0

def load_data(mutation_file, annotation_file):
    """Load mutation data and GFF file."""
    all_data = pd.read_csv(mutation_file)
    gff = Gff(annotation_file)
    print(f"Loaded mutation data:\n{all_data.head()}")
    return all_data, gff

def process_mutation_data(all_data):
    """Process mutation data to extract relevant information."""
    all_strains = all_data['Strain'].unique()
    print(f"Strains found: {all_strains}")
    # colors = ['purple', 'brown', 'orange', 'lightgreen']
    colors = ['tomato', 'skyblue', 'magenta', 'green', 'purple', 'brown', 'orange', 'lightgreen', 'cadetblue']
    strain_color = {strain: color for strain, color in zip(all_strains, colors)}
    
    # Adjusted filtering criteria
    grouped_mutations = all_data.groupby(['locus_tag', 'Minimum', 'Change']).first().reset_index()
    mutation_strain_map = all_data.groupby(['locus_tag', 'Minimum', 'Change']).agg({
        'Strain': lambda x: ', '.join(x), 
        'Variant Frequency': lambda x: ', '.join(map(str, x))
    }).to_dict(orient='index')
    print(f"Grouped mutations:\n{grouped_mutations.head()}")
    return all_strains, strain_color, grouped_mutations, mutation_strain_map

def create_mutation_details(grouped_mutations, mutation_strain_map):
    """Create a list of mutation details."""
    new_labels, pos_list_CDSonly = [], []
    mutation_details = []
    for _, row in grouped_mutations.iterrows():
        freq = min_percentage(row['Variant Frequency'])
        if int(freq) >= 75:
            amino_acid_change = row['Amino Acid Change'] if row['Amino Acid Change'] else "frameshift"
            # label = f"{row['product']}, ({row['Strain']}), {amino_acid_change}, {row['CDS Codon Number']}, {row['Minimum']}" if row['Protein Effect'] not in [None, ''] else ''
            label = f"{row['locus_tag']}, ({amino_acid_change}, {row['CDS Codon Number']})" if row['Protein Effect'] not in [None, ''] else ''
            # new_labels.append(row['locus_tag'])
            new_labels.append(label)
            pos_list_CDSonly.append(row['Minimum'])
            mutation_strain_freq = mutation_strain_map.get((row['locus_tag'], row['Minimum'], row['Change']), {})
            mutation_details.append({
                "Label": label,
                "Description": row['product'],
                "AA\nNumber": pos_list_CDSonly,
                "Change": amino_acid_change,
                "Strain(s)": wrap_text(mutation_strain_freq.get('Strain', ''), 15),
                "Variant Frequency": wrap_text(mutation_strain_freq.get('Variant Frequency', ''), 15)
            })
        else:
            new_labels.append('')
            pos_list_CDSonly.append(row['Minimum'])
    print(f"Generated {len(new_labels)} labels and {len(pos_list_CDSonly)} positions")
    return new_labels, pos_list_CDSonly, mutation_details

def save_plotted_mutations(plotted_mutations, output_csv):
    """Save plotted mutation details to a CSV file."""
    mutation_df = pd.DataFrame(plotted_mutations)
    mutation_df.to_csv(output_csv, index=False)
    print(f"Plotted mutations saved to {output_csv}")


def plot_circos(gff, all_data, all_strains, strain_color, grouped_mutations, new_labels, pos_list_CDSonly, output_file):
    """Create and save a Circos plot."""
    circos = Circos(sectors={gff.name: gff.range_size})
    sector = circos.get_sector(gff.name)
    
    
    # Plot outer track with xticks
    major_ticks_interval = 15000
    minor_ticks_interval = 1500
    outer_track = sector.add_track((90, 91))
    outer_track.axis(fc="lightgrey")
    outer_track.xticks_by_interval(
        major_ticks_interval, label_formatter=lambda v: f"{v/ 10 ** 3:.0f} kbp"
    )
    outer_track.xticks_by_interval(minor_ticks_interval, tick_length=1, show_label=False)



    cds_track = sector.add_track((80,90))
    cds_track.axis(fc="#E5E4E2", ec="black")

    f_cds_feats = gff.extract_features("CDS", target_strand=1)
    cds_track.genomic_features(f_cds_feats, plotstyle="arrow", r_lim=(85, 90), fc="white", ec="black", lw=0.5)
    r_cds_feats = gff.extract_features("CDS", target_strand=-1)
    cds_track.genomic_features(r_cds_feats, plotstyle="arrow", r_lim=(80, 85), fc="white", ec="black", lw=0.5)

    # pos_list, labels = [], []
    # for feat in gff.extract_features("CDS"):
        
    #     start, end = int(str(feat.location.end)), int(str(feat.location.start))
    #     pos = (start + end) / 2
    #     label = feat.qualifiers.get("product", [""])[0]
    #     if label == "thiswasempty" or label.startswith("thiswashypothetical"):
    #         continue
    #     if len(label) > 50:
    #         label = label[:50] + "..."
    #     pos_list.append(pos)
    #     labels.append(label)

    plotted_mutations = []

    for i, strain in enumerate(all_strains):
        strain_data = all_data[all_data['Strain'] == strain]
        radial_range = (75 - i*5, 80 - i*5)
        strain_track = sector.add_track(radial_range, r_pad_ratio=0.1)
        strain_track.axis(fc=strain_color[strain], alpha=0.5, ec="black", lw=1)
        
        for _, mutation in strain_data.iterrows(): 
            mutation_detail = {
            "locus_tag": mutation['locus_tag'],
            "Strain": strain,
            "Minimum": mutation['Minimum'],
            "Variant Frequency": mutation['Variant Frequency'],
            "Amino Acid Change": mutation['Amino Acid Change'],
            "Product": mutation['product'],
            "protein_id": mutation['protein_id'],
            "CDS": mutation['CDS'],
            "Change": mutation['Change'],
            "Coverage": mutation['Coverage'],
            }       
            # if min_percentage(mutation['Variant Frequency']) >= 75 and pd.notna(mutation['Protein Effect']) and mutation['Protein Effect'] != 'None':
            #     strain_track.rect(mutation['Minimum'], mutation['Minimum']+1, fc="black", ec="black", lw=2)
            #     plotted_mutations.append(mutation_detail)
            #     # if mutation['Minimum'] not in pos_list_CDSonly:
            #     #     pos_list_CDSonly.append(mutation['Minimum'])
            #     #     new_labels.append(mutation['locus_tag'])
            # # else:
            # #     strain_track.rect(mutation['Minimum'], mutation['Minimum']+1, fc="white", ec="white", lw=1)
  
            if min_percentage(mutation['Variant Frequency']) >= 75 and pd.notna(mutation['Protein Effect']) and mutation['Protein Effect'] != 'None':
                strain_track.rect(mutation['Minimum'], mutation['Minimum']+1, fc="black", ec="black", lw=2)
                plotted_mutations.append(mutation_detail)
                if mutation['Minimum'] not in pos_list_CDSonly:
                    pos_list_CDSonly.append(mutation['Minimum'])
                    new_labels.append(mutation['locus_tag'])
            # elif 50 <= min_percentage(mutation['Variant Frequency']) < 75 and pd.notna(mutation['Protein Effect']) and mutation['Protein Effect'] != 'None':
            #     strain_track.rect(mutation['Minimum'], mutation['Minimum']+1, fc="lightgrey", ec="lightgrey", lw=1, alpha=0.2)
            #     # plotted_mutations.append(mutation_detail)
            # elif 1 <= min_percentage(mutation['Variant Frequency']) < 50 and pd.notna(mutation['Protein Effect']) and mutation['Protein Effect'] != 'None':
            #     strain_track.rect(mutation['Minimum'], mutation['Minimum']+1, fc="white", ec="white", lw=1, alpha=0.2)
            # #     # plotted_mutations.append(mutation_detail)




    label_track = sector.add_track((95, 100))

    # # for i in range(len(grouped_mutations)):
    # #     label_track.text('.', grouped_mutations['Minimum'].iloc[i], orientation="vertical", color="black", size=22)

    print(f"Final positions list length: {len(pos_list_CDSonly)}, labels list length: {len(new_labels)}")
    # pos_list_CDSonly = adjust_positions(sorted(pos_list_CDSonly), 100)
    label_track.xticks(pos_list_CDSonly, new_labels, outer=True, tick_length=0, label_margin=1, label_orientation="vertical", label_size=9)

    fig_circos = circos.plotfig()
    _ = circos.ax.legend(
        handles=[
            Patch(color="gray", label=f"{gff.name}", alpha=0.6, edgecolor="black", linewidth=0.5),
            # Patch(color="purple", label="mvR09", alpha=0.8, edgecolor="black", linewidth=2),
            # Patch(color="brown", label="mvR10", alpha=0.8, edgecolor="black", linewidth=2),
            # Patch(color="orange", label="mvR11", alpha=0.8, edgecolor="black", linewidth=2),
            # Patch(color="lightgreen", label="mvR12", alpha=0.8, edgecolor="black", linewidth=2),
            # Patch(color="gray", label="Synechocystis \n sp. PCC 6803", alpha=0.6, edgecolor="black"),
            Patch(color="tomato", label="mvR01_Nixon", alpha=0.6, edgecolor="black", linewidth=0.5),
            Patch(color="skyblue", label="mvR02_Nixon", alpha=0.6, edgecolor="black", linewidth=0.5),
            Patch(color="magenta", label="mvR03_Nixon", alpha=0.6, edgecolor="black", linewidth=0.5),
            Patch(color="green", label="mvR06_Nixon", alpha=0.6, edgecolor="black", linewidth=0.5),
            Patch(color="purple", label="mvR09_Howe", alpha=0.6, edgecolor="black", linewidth=0.5),
            Patch(color="brown", label="mvR10_Howe", alpha=0.6, edgecolor="black", linewidth=0.5),
            Patch(color="orange", label="mvR11_Howe", alpha=0.6, edgecolor="black", linewidth=0.5),
            Patch(color="lightgreen", label="mvR12_Howe", alpha=0.6, edgecolor="black", linewidth=0.5),      
        ],
        bbox_to_anchor=(0.5, 0.5),
        loc="center",
        ncols=1,
        fontsize=9,
    )

    fig_circos.savefig(output_file, dpi=600)
    print(f"Figure saved successfully as: {output_file}")
    return plotted_mutations
def main(ref):
    mutation_file = f"../Data/variant_analysis/{ref}_filtered_mutations.csv"
    annotation_file = f"../Data/genome_annotations/{ref}.gff"
    output_file = f"../Figures/WGS/{ref}_plasmid_circos2.svg"
    output_csv = f"../Data/variant_analysis/{ref}_mutation_details2.csv"

    all_data, gff = load_data(mutation_file, annotation_file)
    all_strains, strain_color, grouped_mutations, mutation_strain_map = process_mutation_data(all_data)
    new_labels, pos_list_CDSonly, mutation_details = create_mutation_details(grouped_mutations, mutation_strain_map)
    plotted_mutations = plot_circos(gff, all_data, all_strains, strain_color, grouped_mutations, new_labels, pos_list_CDSonly, output_file)
    save_plotted_mutations(plotted_mutations, output_csv)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate a Circos plot for a given reference accession.")
    parser.add_argument("-ref", required=True, help="Reference accession")
    args = parser.parse_args()
    main(args.ref)