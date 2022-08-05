"""
Pipeline to optimise admixture graph of the known topology
"""
from ete3 import Tree

import os
from os import path
from multiprocess import Pool
from functools import partial

from tree_cov import get_populations
from optimisation import opt_wnd
from utils import *
from pandas import read_csv
import numpy as np

from migadmi_opt import decomposition_of_variance, migadmi_old

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
pop_names = [node.name for node in tree.traverse("levelorder") if node.is_leaf()]
print(tree)
popset_init, variables_init = get_populations(tree, pop_names)


# Read or create admixture
pop_admix, pop_sources_names = read_admixture_old(path_data + file_admixture)

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

def opt(iwnd):
    file_output = f'{path_data_out}wnd_{iwnd}.txt'
    if path.exists(file_output):
        return
    d_mx = read_csv(f'{path_data_in}wnd_{iwnd + 1}.txt', sep='\t', index_col=0)
    res = migadmi_old(pop_names=pop_names,
                      popset_init=popset_init,
                      variables_init=variables_init,
                      pop_admix=pop_admix,
                      admixture_steps=admixture_steps,
                      dist_matrix=d_mx * 2,
                      pop_sources=pop_sources)

    with open(file_output, 'w') as f:
        for k, v in res.items():
            # print(corr + estim)
            f.write(f'{k}\t{v}\n')


# opt(0)

with Pool(n_thr) as workers:
    pmap = workers.map
    pmap(opt, range(n_wnd))


# -------------------------
# Calculate variance
# -------------------------


v_diag, weight_sets = decomposition_of_variance(pop_names=pop_names,
            popset_init=popset_init,
            variables_init=variables_init,
            pop_admix=pop_admix,
            admixture_steps=admixture_steps,
            pop_sources=pop_sources)


decomp_all = []
for iwnd in range(n_wnd):
    print(iwnd)
    variables = dict()
    with open(f'{path_data_out}wnd_{iwnd}.txt', 'r') as f:
        for line in f:
            # print(line)
            k, v = line.split()
            variables[k] = float(v)

    pop_all = pop_names + pop_admix
    decomp_wnd = []
    for iwset, wset in enumerate(weight_sets):

        if len(wset) == 0:
            continue

        v_tmp_0 = dict(variables)
        for w_tmp in wset:
            v_tmp_0[str(w_tmp)] = 0
        v_values_0 = [v.subs(v_tmp_0) for v in v_diag]

        variances = []
        for w in wset:
            v_tmp = dict(v_tmp_0)
            v_tmp[str(w)] = variables[str(w)]

            # v_values_w += [[v.subs(v_tmp) - v_values_0[i] for i,v in enumerate(v_diag)]]

            variances += [v_diag[iwset + len(pop_names)].subs(v_tmp) - v_values_0[iwset + len(pop_names)]]
        variances += [v_values_0[iwset + len(pop_names)]]

        percent = [v / sum(variances) for v in variances]
        for i in range(len(wset)):
            decomp_wnd += [[pop_admix[iwset], pop_all[pop_sources[iwset][i]], percent[i]]]
        decomp_wnd += [[pop_admix[iwset], pop_admix[iwset], percent[len(wset)]]]

    decomp_all += [decomp_wnd]

decomp_names = [[v1, v2] for v1, v2, _ in decomp_all[0]]
decomp_var = np.array([[v for _, _, v in decomp_all[i]] for i in range(len(decomp_all))])
decomp_var = decomp_var.transpose()



with open(f'{path_decomposition}decomposition.txt', 'w') as f:
    for i, v in enumerate(decomp_var):
        # print(corr + estim)
        f.write(f'{decomp_names[i][0]}\t{decomp_names[i][1]}\t')
        f.write('\t'.join([str(x) for x in v]))
        f.write('\n')