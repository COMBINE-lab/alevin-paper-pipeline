from __future__ import print_function
from vpolo.alevin import parser
import gzip
import pandas as pd
import numpy as np
from collections import defaultdict

mappings = defaultdict(set)
for (_, row) in pd.read_csv("./mappings.txt", sep="\n", header=None).iterrows():
    toks = row[0].strip().split('\t')
    ref_cb = toks[0].split("_")[0]
    query_cb = toks[1].split("_")[0]
    mappings[ref_cb].add(query_cb)
    
print(len(mappings))

rna_counts = parser.read_quants_bin("/mnt/scratch5/laraib/alevin2/paperSimAviFinal/em/")

# generate prior based on mapping
rna_peak_counts = np.zeros( (len(mappings), rna_counts.shape[1]) )
rna_peak_counts_row_names = []
ctr = 0

for (rna_cb, atac_cbs) in mappings.items():
    num_atac_cbs = len(atac_cbs)
    frac = 1.0 / num_atac_cbs 

    rna_peak_counts_row_names.append(rna_cb)
    for atac_cb in atac_cbs:
        rna_peak_counts[ctr] += rna_counts.loc[atac_cb] * frac
        
    ctr += 1
    print("\r Done {:.2}".format( ctr / len(mappings)), end="")
    
rna_peak_counts = pd.DataFrame(rna_peak_counts)
rna_peak_counts.index = rna_peak_counts_row_names
rna_peak_counts.columns = rna_counts.columns

with open("prior/quants_mat_cols.txt", 'w') as f:
    for cb in rna_peak_counts.columns:
        f.write(cb+"\n")

# tr_cbs = pd.read_csv("./minnow_human/alevin/true_cell_names.txt", header=None)
# qt_cbs = pd.read_csv("./minnow_human/alevin/quants_mat_rows.txt", header=None)
# ct = pd.concat([qt_cbs, tr_cbs], axis=1)
# ct.columns = ['tr_cb', 'qt_cb']
# cb_map = ct.set_index('qt_cb').to_dict()['tr_cb']

with open("prior/quants_mat_rows.txt", 'w') as f:
    for cb in rna_peak_counts.index:
        f.write( cb+"\n")

np.savetxt("prior/quants_mat.csv", rna_peak_counts, delimiter=",", fmt='%f')
