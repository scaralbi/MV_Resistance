#!/usr/bin/env python3

import sangerseq_viewer.sangerseq_viewer as sang
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patches as mpatches
import dnaplotlib as dpl
from dnaplotlib import DNARenderer
import pySBOL2
import graphviz
from collections import defaultdict
import matplotlib.gridspec as gridspec
import matplotlib.image as mpimg


wdir = '/Volumes/Albi/Cambridge/PhD/Bioinformatics/prqRA/28042023_SP_prqR1_PCR_Products_Sequencing'
ref_dir = wdir+'/gb_files'
ref_data = ref_dir+'/prqR_PCR_product.gb'
seq_dir = wdir+'/ab1_files'
out_dir = wdir+'/figs'

X0 = 700
X1 = 850



fig = sang.view_sanger(ref_data, seq_dir, output=out_dir+'/zoomedincds.png', start = X0, end=X1)

# Adjust the layout
plt.tight_layout()

# Show or save the figure
plt.show()






