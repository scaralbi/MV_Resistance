import argparse
import pandas as pd
from pycirclize import Circos
from pycirclize.parser import Gff
import re
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

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
    colors = ['purple', 'brown', 'orange', 'lightgreen']
    strain_color = {strain: color for strain, color in zip(all_strains, colors)}
    
    # Adjusted filtering criteria
    grouped_mutations = all_data.groupby(['old_locus_tag', 'Minimum', 'Change']).first().reset_index()
    mutation_strain_map = all_data.groupby(['old_locus_tag', 'Minimum', 'Change']).agg({
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
            label = f"{row['old_locus_tag']}" if row['Protein Effect'] not in [None, ''] else ''
            new_labels.append(label)
            pos_list_CDSonly.append(row['Minimum'])
            amino_acid_change = row['Amino Acid Change'] if row['Amino Acid Change'] else "fs"
            mutation_strain_freq = mutation_strain_map.get((row['old_locus_tag'], row['Minimum'], row['Change']), {})
            mutation_details.append({
                "Gene": label,
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

def plot_circos(gff, all_data, all_strains, strain_color, grouped_mutations, new_labels, pos_list_CDSonly, output_file):
    """Create and save a Circos plot."""
    circos = Circos(sectors={gff.name: gff.range_size})
    sector = circos.get_sector(gff.name)
    cds_track = sector.add_track((90, 100))
    cds_track.axis(fc="#E5E4E2", ec="black")

    f_cds_feats = gff.extract_features("CDS", target_strand=1)
    cds_track.genomic_features(f_cds_feats, plotstyle="arrow", r_lim=(95, 100), fc="salmon", ec="black", lw=0.5)
    r_cds_feats = gff.extract_features("CDS", target_strand=-1)
    cds_track.genomic_features(r_cds_feats, plotstyle="arrow", r_lim=(90, 95), fc="skyblue", ec="black", lw=0.5)

    pos_list, labels = [], []
    for feat in gff.extract_features("CDS"):
        
        start, end = int(str(feat.location.end)), int(str(feat.location.start))
        pos = (start + end) / 2
        label = feat.qualifiers.get("product", [""])[0]
        if label == "" or label.startswith("hypothetical"):
            continue
        if len(label) > 30:
            label = label[:30] + "..."
        pos_list.append(pos)
        labels.append(label)

    for i, strain in enumerate(all_strains):
        strain_data = all_data[all_data['Strain'] == strain]
        radial_range = (75 - i*5, 80 - i*5)
        strain_track = sector.add_track(radial_range, r_pad_ratio=0.1)
        strain_track.axis(fc=strain_color[strain], alpha=0.4, ec="black", lw=1)
        
        for _, mutation in strain_data.iterrows():
            if min_percentage(mutation['Variant Frequency']) >= 75:
                strain_track.rect(mutation['Minimum']-0, mutation['Minimum']+1, fc="black", ec="black", lw=2)
            else:
                strain_track.rect(mutation['Minimum']-0, mutation['Minimum']+1, fc="white", ec="white", alpha=0, lw=1)

    label_track = sector.add_track((105, 110))

    for i in range(len(grouped_mutations)):
        label_track.text('.', grouped_mutations['Minimum'].iloc[i], orientation="vertical", color="black", size=22)

    print(f"Final positions list length: {len(pos_list_CDSonly)}, labels list length: {len(labels)}")
    pos_list_CDSonly = adjust_positions(sorted(pos_list), 20)
    label_track.xticks(pos_list_CDSonly, labels, outer=True, tick_length=1, label_margin=1, label_orientation="vertical", label_size=9)

    fig_circos = circos.plotfig()
    _ = circos.ax.legend(
        handles=[
            Patch(color="gray", label=f"{gff.name}", alpha=0.8, edgecolor="black", linewidth=2),
            Patch(color="purple", label="mvR09", alpha=0.8, edgecolor="black", linewidth=2),
            Patch(color="brown", label="mvR10", alpha=0.8, edgecolor="black", linewidth=2),
            Patch(color="orange", label="mvR11", alpha=0.8, edgecolor="black", linewidth=2),
            Patch(color="lightgreen", label="mvR12", alpha=0.8, edgecolor="black", linewidth=2),
        ],
        bbox_to_anchor=(0.5, 0.5),
        loc="center",
        ncols=1,
        fontsize=10,
    )

    fig_circos.savefig(output_file, dpi=600)
    print(f"Figure saved successfully as: {output_file}")

def main(ref):
    mutation_file = f"../Data/variant_analysis/{ref}_filtered_mutations.csv"
    annotation_file = f"../Data/genome_annotations/{ref}.gff"
    output_file = f"../Figures/WGS/{ref}_plasmid_circos.png"

    all_data, gff = load_data(mutation_file, annotation_file)
    all_strains, strain_color, grouped_mutations, mutation_strain_map = process_mutation_data(all_data)
    new_labels, pos_list_CDSonly, mutation_details = create_mutation_details(grouped_mutations, mutation_strain_map)
    plot_circos(gff, all_data, all_strains, strain_color, grouped_mutations, new_labels, pos_list_CDSonly, output_file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate a Circos plot for a given reference accession.")
    parser.add_argument("-ref", required=True, help="Reference accession")
    args = parser.parse_args()
    main(args.ref)

