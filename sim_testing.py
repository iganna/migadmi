"""
Pipeline to optimise admixture graph of the known topology
"""
from ete3 import Tree

import os
from multiprocess import Pool
from functools import partial
from itertools import product

from tree_cov import get_populations
from migadmi_opt import *
from utils import *
from pandas import read_csv, DataFrame
import sys


# --------------------------
# Files and paths
path_data = 'simulation/sim_data/'
# path_data_in = path_data + 'in2/'

path_data_out = path_data + 'out_testing/'
if not os.path.exists(path_data_out):
    os.mkdir(path_data_out)



def simulation(i_sim):

    np.random.seed(seed=i_sim + 239)

    n_adm = 2
    pop_admix = ['admix' + str(i + 1) for i in range(n_adm)]


    # Initial Tree
    file_tree = path_data + 'sim_tree.nwk'
    # file_dist = '{}s{}_dist.txt'.format(path_data_in, i_sim)
    # file_admixture = '{}s{}_adm.txt'.format(path_data_in, i_sim)


    # --------------------------
    # Structure of the dataset
    # Read tree
    tree = Tree(file_tree)
    pop_names = [node.name for node in tree.traverse("levelorder") if node.is_leaf()]
    pop_names_all = pop_names + pop_admix
    print(tree)
    popset_init, variables_init = get_populations(tree, pop_names)
    n_source = len(pop_names)

    # Events
    # Indexes of source populations
    pop_sources = [list(np.sort(np.random.choice(n_source + i, np.random.choice(2)+2, replace=False)))
     for i in range(n_adm)]

    with open('{}test{}_pop_sources.txt'.format(path_data_out, i_sim), 'w') as f:
        f.write('{}\n'.format(pop_sources))

    pop_sources_names = [[pop_names_all[i] for i in idxs]
        for idxs in pop_sources]

    # Admixture steps
    admixture_steps = [list(range(n_adm))]
    # Check possibility of admixtures
    check_admixtures(pop_names, pop_sources, admixture_steps)

    # --------------------------
    # Get distance matrix
    idx_pop = list(range(n_adm))
    popset = popset_init
    variables = variables_init
    weight_sets = []
    weights_used = []
    var_mix_all_init = []

    for i in idx_pop:
        popset, variables, var_mix = create_admixture_alpha(popset,
                                                            variables, pop_sources[i])
        weight_sets, weights_used = define_used_weights(weight_sets, weights_used, variables)
        var_mix_all_init += [var_mix]

    # matrix of distances from the root (covariance matrix)
    v_mx = get_cov(popset, variables)
    # matrix of variances between leaves
    d_mx = get_dist(v_mx)

    # Create random values for parameters
    variables_created = dict()
    for v in variables:
        if is_alpha(v) or is_branch(v):
            variables_created[v] = np.random.uniform() * 0.9 + 0.1

    for weights in weight_sets:
        n_w = len(weights)
        tmp = np.random.uniform(0, 1, n_w)
        tmp /= sum(tmp)
        for i, w in enumerate(weights):
            variables_created[w] = tmp[i]

    # Calculate distance matrix
    variance_real = d_mx.subs(variables_created)
    variance_real = DataFrame(np.array(variance_real), columns=pop_names_all,
                   index=pop_names_all)
    dist_real = variance_real ** (1/2)  # why square root?

    np.savetxt(X=dist_real, fname='{}test{}_dist.txt'.format(path_data_out, i_sim),
               delimiter='\t')

    f = open('{}test{}_variables_created.txt'.format(path_data_out, i_sim), 'w')
    for k, v in variables_created.items():
        # print(corr + estim)
        f.write(f'{k}:{v}\n')
    f.close()


    # --------------------------
    # Simulate frequencies
    n_sim = 100

    # increments are only for branches! weights and admixture parameters are the same
    x_incr = dict()
    for v in variables:
        if not is_branch(v):
            continue
        # the scale parameter is not variance it is a square root of variance
        x_incr[v] = np.random.normal(loc=0,
                                     scale=variables_created[v] ** (1/2),
                                     size=n_sim-1)

    x = []
    x_in = []
    for perc_init in list(range(1, n_sim)):
        f_init = perc_init / n_sim
        x_init = np.log((1-f_init) / f_init)

        # increment of X along the branches: could be positive or negative
        x_increment = dict(variables_created)
        for v in variables:
            if not is_branch(v):
                continue
            x_increment[v] = x_incr[v][perc_init-1]

        x_mx = np.array(v_mx.subs(x_increment)) + x_init
        x_end = np.array(list(x_mx.diagonal()))

        x += [x_end]
        x_in += [f_init]

    x_mx = np.array(x).transpose()

    variance_sim = np.array([[np.array(x_mx[i] - x_mx[j], dtype=np.float64).std() ** 2
                  for j in range(len(x_mx))]
                 for i in range(len(x_mx))])

    dist_sim = DataFrame(np.array(variance_sim ** (1/2)), columns=pop_names_all,
                   index=pop_names_all)


    # --------------------------------------------------
    # Estimation
    # admixture_steps = [[0], [1]]
    popdisp_part = partial(migadmi_old, pop_names=pop_names,
                           popset_init=popset_init,
                           variables_init=variables_init,
                           pop_admix=pop_admix,
                           admixture_steps=admixture_steps,
                           pop_sources=pop_sources)
    # Estimation
    variables_estim_real = popdisp_part(dist_matrix=dist_real * 2)
    variables_estim_sim = popdisp_part(dist_matrix=dist_sim)

    delta_w_real = [np.mean([abs(variables_estim_real[w] - variables_created[w])
                             for w in weights]) for weights in weight_sets]
    delta_w_sim = [np.mean([abs(variables_estim_sim[w] - variables_created[w])
                            for w in weights]) for weights in weight_sets]



    with open('{}test{}.txt'.format(path_data_out, i_sim), 'w') as f:
        f.write('{}\t{}\t{}'.format(i_sim, delta_w_real, delta_w_sim))


n_sim = 100
n_thr = 30
with Pool(n_thr) as workers:
    pmap = workers.map
    pmap(simulation, range(n_sim))




