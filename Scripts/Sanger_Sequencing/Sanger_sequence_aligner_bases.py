import os
import matplotlib.pyplot as plt
from Bio import SeqIO, AlignIO
from Bio.Align.Applications import MuscleCommandline
from collections import defaultdict

ab1_folder = "/Volumes/Albi/Cambridge/PhD/Bioinformatics/prqRA/28042023_SP_prqR1_PCR_Products_Sequencing/ab1_files"
start = 50
end = 200

# Read all AB1 files and store them in a list
seq_records = [SeqIO.read(os.path.join(ab1_folder, f), "abi") for f in os.listdir(ab1_folder) if f.endswith(".ab1")]

# Write all SeqRecord objects to a FASTA file
SeqIO.write(seq_records, "input_sequences.fasta", "fasta")

# Run MUSCLE alignment
muscle_cline = MuscleCommandline(cmd="/Users/albi/bin/muscle3.8.31_i86darwin64", input="input_sequences.fasta", out="aligned_sequences.aln", clwstrict=True)
stdout, stderr = muscle_cline()

# Read the aligned sequences
alignment = AlignIO.read("aligned_sequences.aln", "clustal")

# Function to plot the chromatogram traces
def plot_chromatogram(ax, trace_data, start, end):
    ax.plot(trace_data["DATA9"][start:end], color="blue")
    ax.plot(trace_data["DATA10"][start:end], color="red")
    ax.plot(trace_data["DATA11"][start:end], color="green")
    ax.plot(trace_data["DATA12"][start:end], color="yellow")
    ax.set_yticks([])

# Plot the aligned sequences and chromatogram traces separately
fig, axes = plt.subplots(len(seq_records) * 2, 1, figsize=(5, len(seq_records) * 2), sharex=True)

for i, record in enumerate(seq_records):
    # Get the traces from the record
    channels = ["DATA9", "DATA10", "DATA11", "DATA12"]
    trace = defaultdict(list)
    for c in channels:
        trace[c] = record.annotations["abif_raw"][c][::5]  # Grab every 5th value

    # Get the aligned sequence
    aligned_seq = alignment[i]

    # Plot the traces
    plot_chromatogram(axes[i * 2], trace, start, end)

    # Plot the aligned sequence
    for pos, base in enumerate(aligned_seq[start:end]):
        axes[i * 2 + 1].text(pos, 0.5, base, fontsize=4, ha='center', va='bottom', color='black')
    axes[i * 2 + 1].set_yticks([])


plt.xlabel("Position")
plt.tight_layout()
plt.show()
