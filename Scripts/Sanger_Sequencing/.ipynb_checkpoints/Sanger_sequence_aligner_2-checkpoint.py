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
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio import Phylo
import graphviz
from collections import defaultdict

#Parameter
start = 70
end = 200

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
    dr = dpl.DNARenderer(scale=0.1, linewidth=1.5)
    start, end = dr.renderDNA(ax, design, dr.trace_part_renderers())
    ax.set_xlim([start_position, end])
    ax.set_ylim([-5, 35])
    ax.set_aspect('auto')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.axis('off')

def map_trace_to_alignment(aligned_seq, trace):
    gap_positions = [i for i, char in enumerate(aligned_seq) if char == "-"]
    # Make sure trace and aligned_seq have the same length
    trace = trace[:len(aligned_seq)]
    mapped_trace = np.insert(trace, gap_positions, np.nan)
    return mapped_trace


# Function to plot the chromatogram traces
def plot_chromatogram(ax, trace_data, aligned_seq, start, end, y_position):
    mapped_trace_data9 = map_trace_to_alignment(aligned_seq, np.array(trace_data["DATA9"]))
    mapped_trace_data10 = map_trace_to_alignment(aligned_seq, np.array(trace_data["DATA10"]))
    mapped_trace_data11 = map_trace_to_alignment(aligned_seq, np.array(trace_data["DATA11"]))
    mapped_trace_data12 = map_trace_to_alignment(aligned_seq, np.array(trace_data["DATA12"]))
    ax.plot(range(start, end), mapped_trace_data9[start:end] + y_position*3, color="blue")
    ax.plot(range(start, end), mapped_trace_data10[start:end] + y_position*3, color="red")
    ax.plot(range(start, end), mapped_trace_data11[start:end] + y_position*3, color="green")
    ax.plot(range(start, end), mapped_trace_data12[start:end] + y_position*3, color="yellow")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.grid(False)
    ax.spines["left"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.set_facecolor("none")


def plot_alignment(fig, gs, alignment, reference_record, region=None):
    num_sequences = len(alignment)
    seq_length = len(alignment[0].seq)
    if region is None:
        region = (0, seq_length)
    start, end = region

    colors = dict(zip("ACGTacgt", sns.color_palette("bright", 4) * 2))
    colors["-"] = "white"

    wt_seq_ids = sorted([seq.id for seq in alignment if seq.id.startswith("WT")])
    mv_seq_ids = sorted([seq.id for seq in alignment if seq.id.startswith("MV")])
    reference_seq_id = reference_record.id
    alignment.sort(key=lambda seq_record: ([reference_seq_id] + wt_seq_ids + mv_seq_ids).index(seq_record.id))
    
    
    # Create a list to store all subplot axes
    axes = []
    for i, seq in enumerate(alignment):  # enumerate through alignment
        ax = fig.add_subplot(gs[i+1, :])  # Add a subplot for each row in the GridSpec
        if seq.id != reference_record.id:  # seq is now defined as an element in alignment
            if "abif_raw" in seq.annotations:
                trace = defaultdict(list)
                print("Parsing Chromatogram Traces in ABI files...")
                for c in ["DATA9", "DATA10", "DATA11", "DATA12"]:
                    trace[c] = seq.annotations["abif_raw"][c][::5]
                print("Plotting Chromatogram Traces...")
                plot_chromatogram(ax, trace, str(seq.seq), start, end, 0)
                #ax.set_xticks([]) #Here
                #ax.set_yticks([])
                #ax.grid(False)
                #ax.spines["left"].set_visible(False)
                #ax.spines["right"].set_visible(False)
                #ax.spines["top"].set_visible(False)
                #ax.spines["bottom"].set_visible(False)
               # ax.set_facecolor("none") #New
            for j, aa in enumerate(str(seq.seq)[start:end]):
                y_position = num_sequences - i - 1
                column = [str(s.seq)[j + start] for s in alignment]
                consensus_all = len(set(column)) == 1
                consensus_no_ref = len(set([str(s.seq)[j + start] for s in alignment if s.id != reference_record.id])) == 1
                base_color = colors[aa]
                rectangle_color = "grey" if consensus_all else base_color
                rectangle_alpha = 0.6  # adjust as needed
                print("Drawing Conservation Boxes...")
                ax.add_patch(plt.Rectangle((j-0.5, y_position-3), 1, 2, facecolor=rectangle_color, alpha=rectangle_alpha, edgecolor="black"))
                print("Drawing DNA Letters...A..C....G......T....")
                ax.text(j, y_position, aa, fontsize=10, color=base_color, ha='center', va='center')
                seq_name = seq.id.replace('.ab1', '')
                name_color = 'blue' if seq_name.startswith('WT') else 'red'
                if seq.id == reference_seq_id:
                    name_color = 'black'
                ax.annotate(seq_name, (-1, y_position), fontsize=14, ha='right', va='center', color=name_color)
                ax.set_xticks([]) #Here
                ax.set_yticks([])
                ax.grid(False)
                ax.spines["left"].set_visible(False)
                ax.spines["right"].set_visible(False)
                ax.spines["top"].set_visible(False)
                ax.spines["bottom"].set_visible(False)
                ax.set_facecolor("none") #New
                ax.axis('off')
        if seq.id == reference_record.id:  # seq is now defined as an element in alignment
            for j, aa in enumerate(str(seq.seq)[start:end]):
                y_position = num_sequences - i - 1
                column = [str(s.seq)[j + start] for s in alignment]
                consensus_all = len(set(column)) == 1
                consensus_no_ref = len(set([str(s.seq)[j + start] for s in alignment if s.id != reference_record.id])) == 1
                base_color = colors[aa]
                rectangle_color = "grey" if consensus_all else base_color
                rectangle_alpha = 1  # adjust as needed
                print("Drawing Conservation Boxes in Reference Sequence...")
                ax.add_patch(plt.Rectangle((j-0.5, y_position-3), 1, 2, facecolor=rectangle_color, alpha=rectangle_alpha, edgecolor="black"))
                print("Drawing Reference Sequence...A..C....G......T....")
                ax.text(j, y_position, aa, fontsize=10, color=base_color, ha='center', va='center')
                seq_name = seq.id.replace('.ab1', '')
                name_color = 'blue' if seq_name.startswith('WT') else 'red'
                if seq.id == reference_seq_id:
                    name_color = 'black'
                ax.annotate(seq_name, (-1, y_position), fontsize=14, ha='right', va='center', color=name_color)
                ax.set_xticks([]) #Here
                ax.set_yticks([])
                ax.grid(False)
                ax.spines["left"].set_visible(False)
                ax.spines["right"].set_visible(False)
                ax.spines["top"].set_visible(False)
                ax.spines["bottom"].set_visible(False)
                ax.set_facecolor("none") #New
                ax.axis('off')        
           # ax.set_yticks(range(num_sequences)[::-1])
           # ax.set_yticklabels(['' for _ in alignment], fontsize=12)


        # Modify the x-axis numbering
        xticks = np.arange(0, end - start, 10)
        new_xtick_labels = [str(x + start) for x in xticks]
        ax.set_xticks(xticks)
        ax.set_xticklabels(new_xtick_labels, fontsize=12, rotation=0)
        ax.set_xlim(-1, end - start)
        ax.set_ylim(num_sequences - 0.4, -0.8)
        ax.xaxis.tick_top()
        ax.xaxis.set_label_position('top')
        ax.xaxis.set_label_coords(0.5, 1.15)
        ax.xaxis.set_tick_params(width=0.5, direction='out')
        ax.xaxis.set_tick_params(pad=15)
        ax.set_yticks([])
        ax.grid(False)
        ax.spines["left"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        ax.set_facecolor("none")


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
print("Reading gb files...")
gb_file = "/Volumes/Albi/Cambridge/PhD/Bioinformatics/prqRA/28042023_SP_prqR1_PCR_Products_Sequencing/SP_prqR1_product_NC_000911.gb"
reference_record, start_position = read_gb_file(gb_file)

print("Reading ABI files...")
seq_records = read_abi_files(file_list)

# Add the reference record to seq_records
seq_records.append(reference_record)

# Sort seq_records based on whether the name starts with "WT" or "MV"
seq_records.sort(key=lambda x: x.id.startswith("WT"), reverse=True)

# Write all SeqRecord objects to a FASTA file
SeqIO.write(seq_records, "input_sequences.fasta", "fasta")

# Run MUSCLE alignment
print("Running MUSCLE alignment...")
muscle_cline = MuscleCommandline(cmd="/Users/albi/bin/muscle3.8.31_i86darwin64", input="input_sequences.fasta", out="aligned_sequences.aln", clwstrict=True)
stdout, stderr = muscle_cline()

# Read the aligned sequences
aligned_seqs = SeqIO.parse("aligned_sequences.aln", "clustal")
alignment = MultipleSeqAlignment(aligned_seqs)



## DNAPLOTLIB SPECS
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
opt_CDS1 = {'label':'slr0895', 'label_style':'italic', 'label_y_offset':-7, 'color':col_map['orange'], 'label_size':14}
opt_CDS2 = {'label':'slr0896 (prqA)', 'label_style':'italic', 'label_y_offset':-7, 'color':col_map['purple'], 'label_size':14}


design = [
    {'type': 'UserDefined', 'name': 'PBS', 'start': 70, 'end': 100, 'fwd': True, 'opts': {'color': 'white'}},
    {'type': 'Promoter', 'name': 'P1', 'start': 115, 'end': 160, 'fwd': True, 'opts': {'color': col_map['green']}},
    {'type': 'RBS', 'name': 'RBS1', 'start': 168, 'end': 172, 'fwd': True, 'opts': {'label':'RBS', 'color': col_map['blue'], 'label_y_offset':-7, 'label_size':14}},
    {'type': 'CDS', 'name': 'CDS1', 'start': 180, 'end': 200, 'fwd': True, 'opts': opt_CDS1},
]


# Create the plot
num_sequences = len(alignment)
alignment_length = alignment.get_alignment_length()
matplotlib.rcParams['font.family'] = "Arial"
matplotlib.rcParams['font.size'] = 16  # Increase the font size as needed

fig = plt.figure(figsize=(26, 14))

# Set up the gridspec layout
# Calculate the number of rows for the grid
num_rows = num_sequences + 1   # plus one for the reference plot

# Set the height ratios for the grid cells. The first cell (for the reference plot)
# will have a height of 1, and all other cells will have a height of 3.
height_ratios = [1] + [6] * num_sequences 

# Create the GridSpec object with the calculated number of rows and height ratios
gs = gridspec.GridSpec(num_rows, 1, height_ratios=height_ratios, wspace=0, hspace=0)

# Assign the first and second subplots to ax_dna and ax_alignment
ax_dna = plt.subplot(gs[0, :])
ax_alignment = plt.subplot(gs[1, :])

# Plot the alignment
print("Drawing the alignment...")
plot_alignment(fig, gs, alignment, reference_record, region=(start, end))

#plt.subplots_adjust(top=1)

# Save the plot
print("..Saving...")
plt.savefig("alignment_plot_dnaplotlib_phylotree_zoom.png", dpi=600)
print("..Completed succesfully...")



