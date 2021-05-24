"""
This module containf functions to work with variance-covariance matrices,
based on phylogenetic tree
"""
import numpy as np
from ete3 import Tree
import sympy as sym
from itertools import product, combinations_with_replacement


def define_used_weights(weight_sets, weights_used, variables):
    """
    Check which weights have already been defined
    :param weight_sets:
    :param weights_used:
    :param variables:
    :return:
    """
    weight_sets += [[v for v in variables
                     if is_weight(v) and v not in weights_used]]
    weights_used += weight_sets[-1]

    return weight_sets, weights_used


def is_branch(v):
    """
    Check wether a variable is branch length or weight parameter
    :param v: variable
    :return: True/False
    """
    return str(v)[0] == 't'


def is_alpha(v):
    """
    Check wether a variable is branch length or weight parameter
    :param v: variable
    :return: True/False
    """
    return str(v)[0] == 'a'


def is_weight(v):
    """
    Check wether a variable is branch length or weight parameter
    :param v: variable
    :return: True/False
    """
    return str(v)[0] == 'w'


def creare_variable(variables, s):
    """
    Create new variable with prefix s
    :param variables: list of previous variables
    :param s: prefix
    :return:
    """

    if s not in ['w', 't', 'a']:
        print('Warning: variable with an unknown prefix')

    if len(variables) == 0:
        t = sym.Symbol(s + '0')
    else:
        free_id = int(str(variables[-1])[1:]) + 1
        t = sym.Symbol(s + str(free_id))

    return t, variables + [t]


def is_root_identifiable(mx, variables):
    """
    It is hardcoded, that two first variables in variable list are branches
    from the root towards two clades.
    However, in distance matrix is nome cases they can be found
    only in a combination 't0 + t1'.
    If so, than one of these parameters is not identifiable.
    :param mx:
    :param variables:
    :return:
    """

    for v in variables[2:]:
        if is_branch(v):
            mx = mx.subs(v, 0)

    w, _ = creare_variable(variables, 'w')
    t, _ = creare_variable(variables, 't')
    mx = mx.subs(variables[0], w * t)
    mx = mx.subs(variables[1], (1-w) * t)
    mx.simplify()
    # print(mx)

    return w in mx.free_symbols

    # for i, j in combinations_with_replacement(range(mx.shape[0]), 2):
    #     if (mx[i, j] != 0) and (mx[i, j] != (variables[0] + variables[1])):
    #         return True
    # return False


# def remove_root_param


def n_undefined(variables_d, popset):
    """
    Checks the identification of the problem
    :param variables_d:
    :param popset:
    :return:
    """
    n_pop = len(popset)
    return len(variables_d) - (n_pop*(n_pop-1)/2)


def get_dist(v):
    """
    Get distance matrix from the matrix of variances
    :param v: covariance matrix
    :return:
    """
    n_pop = v.shape[0]
    d = sym.Matrix(np.zeros((n_pop, n_pop)))
    for i, j in product(range(n_pop), repeat=2):
        d[i, j] = v[i, i] + v[j, j] - 2 * v[i, j]

    return d


def get_populations(tree, pop_names):
    """
    Get variances for populations
    :param tree: binary tree
    :param pop_names:  names of populations at leaves
    :return:
    """
    variables = []

    n_pop = len(tree)
    # lam = sym.Symbol('lam')
    # variables += [lam]
    popset = [0] * n_pop

    for node in tree.traverse("levelorder"):
        if node.is_leaf():
            continue
        children0 = [node_tmp.name
                     for node_tmp in node.children[0].traverse("levelorder") if node_tmp.is_leaf()]
        children1 = [node_tmp.name
                     for node_tmp in node.children[1].traverse("levelorder") if node_tmp.is_leaf()]
        idx0 = [i
                for i, x in enumerate(pop_names)
                for j, y in enumerate(children0) if x == y]

        idx1 = [i
                for i, x in enumerate(pop_names)
                for j, y in enumerate(children1) if x == y]


        t0, variables = creare_variable(variables, 't')
        t1, variables = creare_variable(variables, 't')

        for i in idx0:
            popset[i] += t0
        for i in idx1:
            popset[i] += t1

    return popset, variables


def create_admixture(popset, variables, idx_sources=None):
    """
    Introduce parameters for mixtures
    """
    n_pop = len(popset)
    if idx_sources is None:
        idx_sources = range(n_pop)

    # New population
    y, variables = creare_variable(variables, 't')
    variances = []

    for idx in idx_sources:
        w, variables = creare_variable(variables, 'w')
        y += w * popset[idx]
        variances += [w**2 * popset[idx]]

    return popset + [y], variables, variances


def create_admixture_alpha(popset, variables, idx_sources=None):
    """
    Introduce parameters for mixtures via alpha
    """
    n_pop = len(popset)
    if idx_sources is None:
        idx_sources = range(n_pop)

    # New population
    y, variables = creare_variable(variables, 't')
    y_remain = y
    a, variables = creare_variable(variables, 'a')
    variances = []
    var_own = get_own_variance(popset, variables)

    variables_new = []

    for idx in idx_sources:
        w, variables = creare_variable(variables, 'w')
        variables_new += [w]
        # y += w * popset[idx] + w * var_own[idx] * (a - 1)
        y += w * (popset[idx] - var_own[idx]) + a * w * var_own[idx]
        # y.simplify()
        variances += [w**2 * (popset[idx] - var_own[idx]) + a**2 * w**2 * var_own[idx]]

    if len(variables_new) == 1:
        w = variables[-1]
        variables = variables[:-1]
        y = y.subs(w, 1)
        variances = [v.subs(w, 1) for v in variances]


    return popset + [y], variables, variances + [y_remain]


def create_admixture2(popset_init, variables, idx_sources=None):
    """
    Admixture with own variance for each population
    :param popset:
    :param variables:
    :param idx_sources:
    :return:
    """
    popset = list(popset_init)
    n_pop = len(popset)
    if idx_sources is None:
        idx_sources = range(n_pop)

    # New population

    y, variables = creare_variable(variables, 't')
    variables_new = []

    for idx in idx_sources:
        w, variables = creare_variable(variables, 'w')
        variables_new += [w]
        y += w * popset[idx]

    # Set sum to 1 assumption
    variables = variables[:-1]
    y = y.subs(variables_new[-1], 1 - sum(variables_new[:-1]))
    y.simplify()

    # Add own variance for each population
    for i in range(n_pop):
        if i not in idx_sources:
            continue
        t, variables = creare_variable(variables, 't')
        popset[i] += t


    return popset + [y], variables


def create_admixture2_equality(popset_init, variables, idx_sources=None):
    """
    Admixture with own variance for each population
    :param popset:
    :param variables:
    :param idx_sources:
    :return:
    """
    popset = list(popset_init)
    n_pop = len(popset)
    if idx_sources is None:
        idx_sources = range(n_pop)

    # New population

    y, variables = creare_variable(variables, 't')
    variables_new = []

    for idx in idx_sources:
        w, variables = creare_variable(variables, 'w')
        variables_new += [w]
        y += w * popset[idx]



    # Add own variance for each population
    for i in range(n_pop):
        if i not in idx_sources:
            continue
        t, variables = creare_variable(variables, 't')
        popset[i] += t


    return popset + [y], variables


def create_admixture2_equality_var(popset_init, variables, idx_sources=None):
    """
    Admixture with own variance for each population
    :param popset:
    :param variables:
    :param idx_sources:
    :return:
    """
    popset = list(popset_init)
    n_pop = len(popset)
    if idx_sources is None:
        idx_sources = range(n_pop)

    # New population

    y, variables = creare_variable(variables, 't')
    variables_new = []
    variances = []

    for idx in idx_sources:
        w, variables = creare_variable(variables, 'w')
        variables_new += [w]
        y += w * popset[idx]

        variances +=  [ popset[idx]* w ** 2]



    # Add own variance for each population
    for i in range(n_pop):
        if i not in idx_sources:
            continue
        t, variables = creare_variable(variables, 't')
        popset[i] += t


    return popset + [y], variables, variances


def pop_cov(pop1, pop2, variables):
    """
    Covariance between two populations
    :param pop1:
    :param pop2:
    :param variables:
    :return:
    """
    cov_val = 0
    pop_product = (pop1) * (pop2)

    # pop_product = pop_product.subs('lam', 0)
    for v1 in variables:
        pop_product_tmp = pop_product.subs(v1, 1)
        for v0 in variables:
            if v1 == v0:
                continue
            # if str(v0)[0] == 'w':

            if not is_branch(v0):
                continue
            pop_product_tmp = pop_product_tmp.subs(v0, 0)
        cov_val += pop_product_tmp * v1
    return cov_val


def get_cov(popset, variables):
    """
    Get covariance matrix
    :param popset:
    :param variables:
    :return:
    """
    n_pop = len(popset)

    v_mx = sym.Matrix(np.zeros((n_pop, n_pop)))

    for i, j in combinations_with_replacement(range(n_pop), 2):
        tmp = pop_cov(popset[i], popset[j], variables)
        v_mx[i, j] = tmp
        v_mx[j, i] = tmp

    return v_mx


def eval_mx(mx, variables, params):
    """
    Evaluate matrix by parameters
    :param mx: matrix
    :param variables: symbol names of parameters
    :param params: values of parameters
    :return:
    """

    d = dict(zip(variables, params))
    mx = mx.evalf(subs=d)

    # for v, p in zip(variables, params):
    #     mx = mx.subs(v, p)
    return mx


def eval_vect(vect, variables, params):
    """
    Evaluate vector by parameters
    :param vect: vector
    :param variables: symbol names of parameters
    :param params: values of parameters
    :return:
    """
    vect_val = []
    for i in range(len(vect)):
        tmp = vect[i]
        for v, p in zip(variables, params):
            # print(v, tmp)
            tmp = tmp.subs(v, p)
        vect_val += [tmp]
    return vect_val


def get_own_variance(popset, variables):
    """
    Get own variance of a population
    :param popset:
    :return:
    """
    n_pop = len(popset)

    popset_no_weights = list(popset)
    for i in range(n_pop):
        for v in variables:
            if is_branch(v):
                continue
            popset_no_weights[i] = popset_no_weights[i].subs(v, 0)



    var_own = [];
    for i in range(n_pop):
        sym_own = popset_no_weights[i].free_symbols
        sym_other = set()
        for j in range(n_pop):
            if j == i:
                continue
            sym_other |= popset_no_weights[j].free_symbols

        var_own += list(sym_own - sym_other)

    return var_own




