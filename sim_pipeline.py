"""
Pipeline to optimise admixture graph of the known topology
"""
from ete3 import Tree

import os
from multiprocess import Pool
from functools import partial

from tree_cov import get_populations
from migadmi_opt import *
from utils import *
from pandas import read_csv




# --------------------------
# Files and paths
path_data = 'simulation/sim_data/'
path_data_in = path_data + 'in/'
path_data_out = path_data + 'out/'

path_data = 'simulation/sim_data/'
path_data_in = path_data + 'in2/'
path_data_out = path_data + 'out2/'


if not os.path.exists(path_data_out):
    os.mkdir(path_data_out)
n_sim = 100

for i_sim in range(1, n_sim+1):
    # Initial Tree
    file_tree = path_data + 'sim_tree.nwk'
    # file_dist = path_data_in + 'd01.txt'
    # file_dist = path_data_in + 'd02.txt'
    file_dist = '{}s{}_dist.txt'.format(path_data_in, i_sim)
    file_admixture = '{}s{}_adm.txt'.format(path_data_in, i_sim)
    file_estim = '{}s{}_estim.txt'.format(path_data_out, i_sim)


    file_dist = path_data_in + 'd01.txt'
    file_admixture = path_data_in + 'admix01.txt'
    file_estim = '{}s{}_estim.txt'.format(path_data_out, i_sim)


    # --------------------------
    # Optimise just a tree
    # Read tree
    tree = Tree(file_tree)
    pop_names = [node.name for node in tree.traverse("levelorder") if node.is_leaf()]
    print(tree)
    popset_init, variables_init = get_populations(tree, pop_names)

    # Read distance matrix
    dist_real = read_csv(file_dist, sep='\t', index_col=0)

    pop_admix, pop_sources_names = read_admixture_old(file_admixture)
    # Indexes of source populations
    pop_sources = [[i for i, p in enumerate(pop_names + pop_admix) if p in pop]
                        for pop in pop_sources_names]

    # Steps to optimize admixture in a sequence
    admixture_steps = [[0,1]]
    # admixture_steps = [[0], [1], [2]]
    # admixture_steps = [[0, 1, 2]]
    # Check possibility of admixtures
    check_admixtures(pop_names, pop_sources, admixture_steps)

    # --------------------------
    # Optimisation
    variables_estimated = migadmi_old(pop_names,
                                      popset_init,
                                      variables_init,
                                      pop_admix,
                                      admixture_steps,
                                      dist_real,
                                      pop_sources=pop_sources)

    # Save results
    with open(file_estim, 'w') as file:
        for k, v in variables_estimated.items():
            file.write('{}\t{}\n'.format(k,v))



