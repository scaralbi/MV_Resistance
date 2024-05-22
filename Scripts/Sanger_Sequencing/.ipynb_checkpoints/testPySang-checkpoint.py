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
import PySanger
from collections import defaultdict
from pysanger import * 

# Read GenBANK FILES
def read_gb_file(file_path):
    record = SeqIO.read(file_path, "genbank")
    keywords = record.annotations.get("keywords")
    if keywords:
        record.id = keywords[0]

    start_position = int(record.features[0].location.start)
    return record, start_position




gb_file="/Volumes/Albi/Cambridge/PhD/Bioinformatics/prqRA/28042023_SP_prqR1_PCR_Products_Sequencing/SP_prqR1_product_NC_000911.gb"
reference_record, start_position = read_gb_file(gb_file)

abidata=abi_to_dict(filename='/Volumes/Albi/Cambridge/PhD/Bioinformatics/prqRA/28042023_SP_prqR1_PCR_Products_Sequencing/ab1_files/WT1.ab1')  
fseq, rseq = generate_consensusseq(abidata)  
fig        = visualize(abidata, template=gb_file, region="aligned") 
fig.savefig("testPySanger.png", bbox_inches="tight") 


