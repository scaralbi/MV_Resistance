import os
import glob
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Align.Applications import MuscleCommandline
from Bio.SeqRecord import SeqRecord
from Bio.Seq import MutableSeq
import matplotlib.patches as mpatches
from matplotlib.path import Path
import dnaplotlib as dpl
from dnaplotlib.colors import colorset

# Read AB1 FILES
def read_abi_files(file_list):
    return [SeqIO.read(file, "abi") for file in file_list]

# Read GenBANK FILES
def read_gb_file(file_path):
    record = SeqIO.read(file_path, "genbank")
    keywords = record.annotations.get("keywords")
    if keywords:
        record.id = keywords[0]

    start_position = int(record.features[0].location.start)
    return record, start_position

    def genbank_features_to_dpl_glyphs(reference_record):
        glyphs = []
        for feature in reference_record.features:
            if feature.type in ["promoter"]:
                start = int(feature.location.start)
                end = int(feature.location.end)
                direction = 1 if feature.strand >= 0 else -1
                glyph = {'type': 'Promoter', 'name': feature.qualifiers.get("label", [""])[0], 'start': start,     'end': end, 'fwd': direction == 1}
                glyphs.append(glyph)
        return glyphs


def plot_alignment(alignment, reference_record, filename="alignment_plot_3.png", start_position=0):
    num_sequences = len(alignment)
    seq_length = len(alignment[0].seq)

    # Define the colors for the nucleotides
    colors = dict(zip("ACGTacgt", sns.color_palette("Set2", 4) * 2))
    colors["-"] = "white"

    # Sort the sequences based on their ids
    wt_seq_ids = sorted([seq.id for seq in alignment if seq.id.startswith("WT")])
    mv_seq_ids = sorted([seq.id for seq in alignment if seq.id.startswith("MV")])
    reference_seq_id = reference_record.id
    alignment.sort(key=lambda seq_record: ([reference_seq_id] + wt_seq_ids + mv_seq_ids).index(seq_record.id))

    # Create the plot
    plt.figure(figsize=(20, num_sequences))
    matplotlib.rcParams['font.family'] = "Arial"


    for i, seq in enumerate(alignment[::-1]):  # Change this line
        if seq.id == reference_record.id:
            seq.id = reference_record.id  # Change this line

        for j, aa in enumerate(seq.seq):
            y_position = num_sequences - i - 1  # Add this line
            column = [s.seq[j] for s in alignment]
            consensus_all = len(set(column)) == 1
            consensus_no_ref = len(set([s.seq[j] for s in alignment if s.id != reference_record.id])) == 1  # Change this line

            if seq.id == reference_record.id:  
                if consensus_all:
                    base_color = "gray"
                else:
                    base_color = colors[aa]
            else:
                if consensus_all:
                    base_color = "gray"
                elif consensus_no_ref:
                    if '-' in column:
                        base_color = "lightgray"
                    else:
                        base_color = colors[aa]
                else:
                    # Check if there is a mismatch between any WT sequence and any MV sequence
                    wt_bases = [s.seq[j] for s in alignment if s.id.startswith("WT")]
                    mv_bases = [s.seq[j] for s in alignment if s.id.startswith("MV")]
                    mismatch = not np.any([wt_base == mv_base for wt_base in wt_bases for mv_base in mv_bases])

                    if mismatch or ('-' in wt_bases and '-' in mv_bases):
                        base_color = colors[aa]
                    elif '-' in column:
                        base_color = "lightgray"
                    else:
                        base_color = colors[aa]

            rect = plt.Rectangle((j, y_position), 1, 0.4, facecolor=base_color, lw=1)  # Change this line
            plt.gca().add_patch(rect)

            seq_name = seq.id.replace('.ab1', '')
            name_color = 'blue' if seq_name.startswith('WT') else 'red'
            if seq.id == reference_seq_id:
                name_color = 'black'
            plt.gca().annotate(seq_name, (-0.5, y_position + 0.2), fontsize=12, ha='right', va='center', color=name_color)  # Change this line

    plt.yticks(range(num_sequences)[::-1], ['' for _ in alignment], fontsize=12)
    
    
    
    # Modify the x-axis numbering
    xticks = np.arange(0, alignment.get_alignment_length(), 100)
    new_xtick_labels = [str(x + start_position) for x in xticks]
    plt.xticks(xticks, new_xtick_labels, fontsize=10, rotation=90)
    plt.grid(axis='y', linestyle='--')
    plt.gca().xaxis.tick_top()
    plt.xlim(-1, len(reference_record.seq))
    plt.ylim(-0.6, num_sequences - 0.4)
    

    # Add legend for base colors
    legend_elements = [plt.Rectangle((0, 0), 1, 1, facecolor=colors[base], edgecolor='black', label=base,  lw=1)
                    for base in 'ACGT'] + \
                    [plt.Rectangle((0, 0), 1, 1, facecolor='gray', edgecolor='black', label='Consensus\n(with ref)',  lw=1),  # Modify this line
                    plt.Rectangle((0, 0), 1, 1, facecolor='lightgray', edgecolor='black', label='Consensus\n(Sanger only)',  lw=1),  # Modify this line
                    plt.Rectangle((0, 0), 1, 1, label='Gap beween\nref and Sanger', edgecolor='black', facecolor='white', linestyle='-', lw=1)]  # Modify this line
    plt.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(1.15, 1))  # Modify this line

    # Add feature annotations to the reference sequence
    ax = plt.gca()
    y_offset = -0.6 
    # Save the plot
    plt.savefig(filename, dpi=300)

    for feature in reference_record.features:
        if feature.type in ["CDS", "promoter", "RBS", "misc_feature"]:
            start = feature.location.start
            end = feature.location.end
            label = feature.qualifiers.get("gene", [""])[0] or feature.qualifiers.get("product", [""])[0] or feature.qualifiers.get("label", [""])[0]
            line = plt.Line2D((start, end), (y_offset, y_offset), lw=1, color='k')  # Change this line
            ax.add_line(line)
            ax.annotate(label, (start + (end - start) / 2, y_offset), color="k", fontsize=10, ha="center", va="top")  # Change this line

    plt.ylim(num_sequences - 0.4, -0.8)

    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    plt.show()




### Define Read AB1 
# Update the read_abi_files function
def read_abi_files(file_list):
    seq_records = []
    for file in file_list:
        seq_record = SeqIO.read(file, "abi")
        seq_record.id = os.path.splitext(os.path.basename(file))[0]
        seq_records.append(seq_record)
    return seq_records



### Script
# Read all .ab1 files in the 'ab1_files' directory
file_list = glob.glob("ab1_files/*.ab1")

# Read all gb files in the 'ab1_files' directory
gb_file = "/Volumes/Albi/Cambridge/PhD/Bioinformatics/prqRA/28042023_SP_prqR1_PCR_Products_Sequencing/SP_prqR1_product.gb"
reference_record, start_position = read_gb_file(gb_file)

seq_records = read_abi_files(file_list)

# Add the reference record to seq_records
seq_records.append(reference_record)

# Sort seq_records based on whether the name starts with "WT" or "MV"
seq_records.sort(key=lambda x: x.id.startswith("WT"), reverse=True)

# Write all SeqRecord objects to a FASTA file
SeqIO.write(seq_records, "input_sequences.fasta", "fasta")

# Run MUSCLE alignment
muscle_cline = MuscleCommandline(cmd="/Users/albi/bin/muscle3.8.31_i86darwin64", input="input_sequences.fasta", out="aligned_sequences.aln", clwstrict=True)
stdout, stderr = muscle_cline()

# Read the aligned sequences
aligned_seqs = SeqIO.parse("aligned_sequences.aln", "clustal")
alignment = MultipleSeqAlignment(aligned_seqs)

# Plot the alignment
plot_alignment(alignment, reference_record, start_position=start_position)

