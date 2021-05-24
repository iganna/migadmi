"""
This module contain functions to optimise parameters
"""

import numpy as np
from scipy import optimize
from scipy.stats import norm
from pandas import read_csv
from utils import showl

from itertools import combinations, chain
from functools import partial

from tree_cov import eval_mx, eval_vect, is_branch, is_weight, is_alpha, \
    create_admixture_alpha, define_used_weights, get_cov, get_dist


def likelihood(params, wnd, d_mx, variables_d, alpha):
    """
    Computes the likelihood value
    :param params: values of parameters
    :param wnd: observed matrix of differences
    :param d_mx: symbolic matrix of differences between populations
    :param variables_d: symbols for parameters
    :param alpha: parameter for the Dirichlet prior
    :return:
    """

    n_pop = wnd.shape[0]
    d_values = eval_mx(d_mx, variables_d, params).tolist()

    var_tol = 10 ** -10

    ll_value = 0
    for i, j in combinations(range(n_pop), 2):
        variance = d_values[i][j]
        if np.abs(variance) < var_tol:
            variance = var_tol
        ll_value += norm.logpdf(x=wnd.iloc[i, j],
                                loc=0, scale=float(variance ** (1/2)))

    # # Ridge/Lasso
    # for i, v in enumerate(variables_d):
    #     if is_weight(v):
    #         ll_value -= (params[i] ** 2) * 0.1


    # # Paired
    # idx_weights = [i for i, v in enumerate(variables_d) if is_weight(v)]
    # for i, j in combinations(idx_weights, 2):
    #     ll_value -= (params[i] * params[j]) * 0.5


    # Dirichlet
    idx_weights = [i for i, v in enumerate(variables_d) if is_weight(v)]
    for i in idx_weights:
        if params[i] > 10 ** -5:
            ll_value += (alpha - 1) * np.log(params[i])
        else:
            ll_value += (alpha - 1) * np.log(10 ** -5)

    return -ll_value


def param_estim(d_mx, variables_d, weight_sets, var_mix_all, wnd, alpha, opt_method='SLSQP'):
    """

    :param d_mx:
    :param variables_d:
    :param weight_sets:
    :param var_mix_all:
    :param wnd:
    :return:
    """

    bnds = [(0, float('inf')) if is_branch(v) else (0, 1)
            for v in variables_d]
    bnds = [(0.1, 1) if is_alpha(v) else bnds[i] for i, v in enumerate(variables_d)]

    cons_w = []
    for i in range(len(weight_sets)):
        if len(weight_sets[i]) == 0:
            continue
        cons_w += [[1 if v in weight_sets[i] else 0 for v in variables_d]]
    A = cons_w
    lbnd = [1] * len(cons_w)
    upbnd = [1] * len(cons_w)

    linear_constraint = optimize.LinearConstraint(A, lbnd, upbnd)

    obj_func = partial(likelihood, wnd=wnd, d_mx=d_mx, variables_d=variables_d, alpha=alpha)

    max_fun = 1000
    for i in [2, 4, 6, 8]:
        x0 = [i + 1 if is_branch(v) else 0.5
              for v in variables_d]
        if len(weight_sets) == 0:
            res = optimize.minimize(obj_func, x0, method=opt_method, bounds=bnds)
        else:
            res = optimize.minimize(obj_func, x0, method=opt_method, bounds=bnds,
                                    constraints=linear_constraint)
        if res.fun < max_fun:
            max_fun = res.fun
            max_res = res

    res = max_res

    variances_all = [eval_vect(var_mix, variables_d, res.x) for var_mix in var_mix_all]
    variances_val = [[float(v / sum(variances)) for v in variances] for variances in variances_all]

    return res.x, [item for sublist in variances_val for item in sublist]


def opt_wnd(pop_names, popset_init,
            pop_admix,
            pop_sources,
            admixture_steps,
            variables_init,
            path_data_in,
            path_data_out,
            pref_wnd_init='wnd_',
            pref_out_res='res_',
            pref_out_var='var_',
            opt_method='SLSQP',
            alpha=0.9,
            iwnd=None):
    """
    Optimise all admixtures in one window
    :param iwnd: id of the window
    :param pop_names: names of initial populations
    :param pop_admix: names of mixed populations
    :param pop_sources: sources of admixtures
    :param admixture_steps: sequence of admixtures
    :param popset_init: initial
    :param variables_init:
    :param path_data_in: path to the folder with data
    :param path_data_out: path to the folder for the output
    :param pref_wnd_init: prefix for the input files with the observed distance matrix
    :param pref_out_res: prefix for the output file with parameter estimates
    :param pref_out_var: prefix for the output file with variances
    :param opt_method: optimisation method
    :return: -
    """

    if (iwnd is None):
        swnd = ''
    else:
        swnd = str(iwnd + 1)


    variables_estimated = dict()
    for istep in range(len(admixture_steps)):

        idx_pop = [item for i in range(istep + 1) for item in admixture_steps[i]]
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

        v_mx = get_cov(popset, variables)
        d_mx = get_dist(v_mx)

        # remove already estimated variables
        d_mx = d_mx.subs(variables_estimated)
        var_mix_all = [[v.subs(variables_estimated) for v in var_mix]
                       for var_mix in var_mix_all_init]

        variables_d = [v for v in variables if v not in variables_estimated.keys()]
        weight_sets = [[w for w in elem if w not in variables_estimated.keys()] for elem in weight_sets]

        # --------------------------
        # Optimisation


        wnd_init = read_csv(path_data_in + pref_wnd_init + swnd + '.txt', sep='\t', index_col=0)

        pop_current = pop_names + [p for i, p in enumerate(pop_admix) if i in idx_pop]
        wnd = wnd_init.loc[pop_current, pop_current]

        # print(wnd)
        # print(d_mx)


        param_val, variances = param_estim(d_mx, variables_d, weight_sets, var_mix_all, wnd,
                                           opt_method=opt_method,
                                           alpha=alpha)

        param_val[param_val < 10 ** -6] = 0


        showl(list(zip(variables_d, param_val)))
        showl(variances)

        variables_estimated.update(dict(zip(variables_d, param_val)))

        # If you estimate the first tree without admixtures - do not estimate
        # parameters around the root
        if istep == 0 and len(admixture_steps[0]) == 0:
            variables_estimated.pop(variables[0])
            variables_estimated.pop(variables[1])

    # np.savetxt(X=list(variables_estimated.values()), fname=path_data_out + pref_out_res + str(iwnd + 1) + '.txt', delimiter='\t')

    pop_all = pop_names + pop_admix
    alpha_corresp = list(chain(*[[[pop_admix[i], pop_all[j]] for j in pop_sources[i]]
                      for i in range(len(pop_sources))]))
    alpha_extims = [[x, y] for x, y in variables_estimated.items()
                    if is_weight(x)]

    f = open(path_data_out + pref_out_res + swnd + '.txt', 'w')
    for corr, estim in list(zip(alpha_corresp, alpha_extims)):
        # print(corr + estim)
        f.write('{} {}: {}: {}\n'.format(corr[0], corr[1], estim[0], estim[1]))
    f.close()


    np.savetxt(X=variances, fname=path_data_out + pref_out_var + swnd + '.txt', delimiter='\t')
