"""
Pipeline to optimise admixture graph of the known topology
"""
from ete3 import Tree

import os
from os import path
from multiprocess import Pool
from functools import partial

from tree_cov import get_populations
from utils import *
from pandas import read_csv
import numpy as np

from migadmi_opt import decomposition_of_variance, migadmi

# --------------------------

# Admixture type
kabuli_tur = False
if kabuli_tur:
    adm_type = 'kabuli_from_tur'
else:
    adm_type = 'kabuli_from_uzb'


# Files and paths
path_data = 'data/'
path_data_in = path_data + 'in/'
path_data_out = path_data + f'out_{adm_type}/'
path_decomposition = path_data + f'decomposition_{adm_type}/'

if not os.path.exists(path_data_out):
    os.mkdir(path_data_out)
if not os.path.exists(path_decomposition):
    os.mkdir(path_decomposition)

# Initial Tree
file_tree = 'tree_init.nwk'
# Admixture chain
file_admixture = f'{adm_type}.txt'

# --------------------------
# Read tree
tree = Tree(path_data + file_tree)

# Read or create admixture
admixtures = read_admixture(path_data + file_admixture)


# Steps to optimize admixture in a sequence
admixture_steps = [[0, 1], [2, 3], [4, 5]]
# admixture_steps = [[0, 1], [2, 3, 4, 5]]

# --------------------------
# Optimisation in windows
#
n_wnd = 167  # number of windows
n_thr = 4

def opt(iwnd):
    file_output = f'{path_data_out}wnd_{iwnd}.txt'
    # if path.exists(file_output):
    #     return
    d_mx = read_csv(f'{path_data_in}wnd_{iwnd + 1}.txt', sep='\t', index_col=0)

    variables, variance_decomposition, weight_sets = migadmi(tree=tree,
                                                             admixtures=admixtures,
                                                             admixture_steps=admixture_steps,
                                                             dist_matrix=d_mx, alpha=1)


    # with open(file_output, 'w') as f:
    #     for k, v in variables.items():
    #         f.write(f'{k}\t{v}\n')

    return variance_decomposition

#

with Pool(n_thr) as workers:
    pmap = workers.map
    decomp_all = pmap(opt, range(n_wnd))


# -------------------------
# Save (Write) variance

decomp_names = [[v1, v2] for v1, v2, _ in decomp_all[0]]
decomp_var = np.array([[v for _, _, v in decomp_all[i]] for i in range(len(decomp_all))])
decomp_var = decomp_var.transpose()

with open(f'{path_decomposition}decomposition.txt', 'w') as f:
    for i, v in enumerate(decomp_var):
        # print(corr + estim)
        f.write(f'{decomp_names[i][0]}\t{decomp_names[i][1]}\t')
        f.write('\t'.join([str(x) for x in v]))
        f.write('\n')

