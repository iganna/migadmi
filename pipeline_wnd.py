"""
Pipeline to optimise admixture graph of the known topology
"""
import numpy as np
from ete3 import Tree

import os
from multiprocess import Pool
from functools import partial

from tree_cov import get_populations
from optimization import opt_wnd
from utils import *


# --------------------------
# Files and paths
path_data = 'data/'
path_data_in = path_data + 'in/'
path_data_out = path_data + 'out/'

if not os.path.exists(path_data_in):
    os.mkdir(path_data_in)
if not os.path.exists(path_data_out):
    os.mkdir(path_data_out)

# Initial Tree
file_tree = 'tree_init.nwk'

# Admixture
kabuli_tur = True
if kabuli_tur:
    file_admixture = 'kabuli_from_tur.txt'
else:
    file_admixture = 'kabuli_from_uzb.txt'

# --------------------------
# Read tree
tree = Tree(path_data + file_tree)
pop_names = [node.name for node in tree.traverse("levelorder") if node.is_leaf()]
print(tree)
popset_init, variables_init = get_populations(tree, pop_names)


# Read or create admixture
pop_admix, pop_sources_names = read_admixture(path_data + file_admixture)

# # If you want to create admixtures without a file,
# # please create (1) list of names on mixed populations
# # (2) the corresponding list of sources for each mixed populations.
# # See code in the following example:
# from data.example import example_admixtures
# pop_admix, pop_sources_names = example_admixtures()

# Indexes of source populations
pop_sources = [[i for i, p in enumerate(pop_names + pop_admix) if p in pop]
                    for pop in pop_sources_names]

# Steps to optimize admixture in a sequence
admixture_steps = [[0, 1], [2, 3], [4, 5]]
# admixture_steps = [[0, 1], [2, 3, 4, 5]]

# Check possibility of admixtures
check_admixtures(pop_names, pop_sources, admixture_steps)

# --------------------------
# Optimisation in windows
#
n_wnd = 167  # number of windows
n_thr = 4
with Pool(n_thr) as workers:
    pmap = workers.map
    opt = partial(opt_wnd, pop_names=pop_names,
                  pop_admix=pop_admix,
                  pop_sources=pop_sources,
                  admixture_steps=admixture_steps,
                  popset_init=popset_init,
                  variables_init=variables_init,
                  path_data_in=path_data_in,
                  path_data_out=path_data_out)

    pmap(opt, range(n_wnd))





