import os
import glob
import numpy as np
import matplotlib
import csv
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.gridspec as gridspec
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio.Align import MultipleSeqAlignment
from Bio.Align.Applications import MuscleCommandline
from Bio.SeqRecord import SeqRecord
from Bio.Seq import MutableSeq
import matplotlib.patches as mpatches
import dnaplotlib as dpl
from dnaplotlib import DNARenderer
import pySBOL2


# Function definitions

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



def plot_dna_design(ax, design):
    dr = dpl.DNARenderer(scale=0.7, linewidth=1.5)
    start, end = dr.renderDNA(ax, design, dr.trace_part_renderers())
    ax.set_xlim([start_position, end])
    ax.set_ylim([-5, 35])
    ax.set_aspect('auto')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.axis('off')

def plot_alignment(ax, alignment, reference_record, start_position=0, ax_dna=None):
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

    # Disable spines and ticks for the subplot
    ax.set_frame_on(False)
    ax.xaxis.set_ticks([])
    ax.yaxis.set_ticks([])

    for i, seq in enumerate(alignment[::-1]):
        if seq.id == reference_record.id:
            seq.id = reference_record.id

        for j, aa in enumerate(seq.seq):
            y_position = num_sequences - i - 1
            column = [s.seq[j] for s in alignment]
            consensus_all = len(set(column)) == 1
            consensus_no_ref = len(set([s.seq[j] for s in alignment if s.id != reference_record.id])) == 1

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
                    wt_bases = [s.seq[j] for s in alignment if s.id.startswith("WT")]
                    mv_bases = [s.seq[j] for s in alignment if s.id.startswith("MV")]
                    mismatch = not np.any([wt_base == mv_base for wt_base in wt_bases for mv_base in mv_bases])

                    if mismatch or ('-' in wt_bases and '-' in mv_bases):
                        base_color = colors[aa]
                    elif '-' in column:
                        base_color = "lightgray"
                    else:
                        base_color = colors[aa]

            rect = plt.Rectangle((j, y_position), 1, 0.4, facecolor=base_color, lw=1)
            ax.add_patch(rect)

            seq_name = seq.id.replace('.ab1', '')
            name_color = 'blue' if seq_name.startswith('WT') else 'red'
            if seq.id == reference_seq_id:
                name_color = 'black'
            ax.annotate(seq_name, (-0.5, y_position + 0.2), fontsize=14, ha='right', va='center', color=name_color)
    ax.set_yticks(range(num_sequences)[::-1])
    ax.set_yticklabels(['' for _ in alignment], fontsize=12)

    # Modify the x-axis numbering
    xticks = np.arange(0, alignment.get_alignment_length(), 100)
    new_xtick_labels = [str(x + start_position) for x in xticks]
    ax.set_xticks(xticks)
    ax.set_xticklabels(new_xtick_labels, fontsize=12, rotation=0)
    ax.set_xlim(-1, len(reference_record.seq))
    ax.set_ylim(num_sequences - 0.4, -0.8)
    ax.grid(axis='y', linestyle='--')
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position('top')
    ax.xaxis.set_label_coords(0.5, 1.15)
    ax.xaxis.set_tick_params(width=0.5, direction='out')
    ax.xaxis.set_tick_params(pad=15)

 # Add legend for base colors
    legend_elements = [plt.Rectangle((0, 0), 1, 1, facecolor=colors[base], edgecolor='black', label=base,  lw=1)
                    for base in 'ACGT'] + \
                    [plt.Rectangle((0, 0), 1, 1, facecolor='gray', edgecolor='black', label='Consensus\n(with ref)',  lw=1),
                    plt.Rectangle((0, 0), 1, 1, facecolor='lightgray', edgecolor='black', label='Consensus\n(Sanger only)',  lw=1),
                    plt.Rectangle((0, 0), 1, 1, label='Gap', edgecolor='black', facecolor='white', linestyle='-', lw=1)]
    ax.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(1.2, 0.9))


####Update the read_abi_files function
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

# Color maps for formatting
col_map = {}
col_map['red']     = (0.95, 0.30, 0.25)
col_map['green']   = (0.38, 0.82, 0.32)
col_map['blue']    = (0.38, 0.65, 0.87)
col_map['orange']  = (1.00, 0.75, 0.17)
col_map['purple']  = (0.55, 0.35, 0.64)
col_map['yellow']  = (0.98, 0.97, 0.35)
col_map['grey']    = (0.70, 0.70, 0.70)
col_map['dark_grey'] = (0.60, 0.60, 0.60)
col_map['light_grey'] = (0.9, 0.9, 0.9)

# CDS formatting options
opt_CDS1 = {'label':'slr0895 (prqR)', 'label_style':'italic', 'label_y_offset':-7, 'color':col_map['orange'], 'label_size':14}
opt_CDS2 = {'label':'slr0896 (prqA)', 'label_style':'italic', 'label_y_offset':-7, 'color':col_map['purple'], 'label_size':14}


design = [
    {'type': 'EmptySpace', 'name': 'sp1', 'start': 1, 'end': 115, 'fwd': True, 'opts': {'color': col_map['light_grey']}},
    {'type': 'Promoter', 'name': 'P1', 'start': 115, 'end': 160, 'fwd': True, 'opts': {'label':'PprqR', 'color': col_map['green'], 'label_y_offset':-7, 'label_size':14}},
    {'type': 'RBS', 'name': 'RBS1', 'start': 171, 'end': 175, 'fwd': True, 'opts': {'label':'RBS', 'color': col_map['blue'], 'label_y_offset':-7, 'label_size':14}},
    {'type': 'CDS', 'name': 'CDS1', 'start': 185, 'end': 748, 'fwd': True, 'opts': opt_CDS1},
    {'type': 'Terminator', 'name': 'term', 'start': 755, 'end': 785, 'fwd': True, 'opts': {'label':'intergenic region', 'color': col_map['purple'], 'label_y_offset':-7, 'label_size':14}},
    {'type': 'UserDefined', 'name': 'ins', 'start': 809, 'end': 829, 'fwd': True, 'opts': {'color': col_map['red']}},
    {'type': 'CDS', 'name': 'CDS2', 'start': 882, 'end': 927, 'fwd': True, 'opts': opt_CDS2},
]


# Create the plot
num_sequences = len(alignment)
alignment_length = alignment.get_alignment_length()
matplotlib.rcParams['font.family'] = "Arial"
matplotlib.rcParams['font.size'] = 14  # Increase the font size as needed

fig = plt.figure(figsize=(16, 9))

# Set up the gridspec layout
gs = gridspec.GridSpec(2, 1, height_ratios=[2, 10], hspace=0)
ax_dna = plt.subplot(gs[0])
ax_alignment = plt.subplot(gs[1])

# Plot the custom design
plot_dna_design(ax_dna, design)

# Plot the alignment
plot_alignment(ax_alignment, alignment, reference_record, start_position=0)
plt.tight_layout()
plt.subplots_adjust(top=0.99)

# Save the plot
plt.savefig("alignment_plot_dnaplotlib.png", dpi=600)



