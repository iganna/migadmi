"""
Example to run the migadmi
"""

from ete3 import Tree
from pandas import read_csv

from utils import read_admixture
from migadmi_opt import migadmi

# Setup
path_data = 'example/'


# ------- Tree -------
file_tree = 'tree_init.nwk'
tree = Tree(path_data + file_tree)
# Get names of the base source populations
pop_names = [node.name for node in tree.traverse("levelorder") if node.is_leaf()]

# ------- Distance matrix -------

d_mx = read_csv(f'{path_data}wnd_40.txt', sep='\t', index_col=0)

# # User can also provide the allele frequencies:


# ------- Admixture events -------

file_admixture = 'kabuli_from_tur.txt'
admixtures = read_admixture(path_data + file_admixture)

# # If you want to create admixtures without a file,
# # please create (1) list of names on mixed populations
# # (2) the corresponding list of sources for each mixed populations.
# # See code in the following example:
# admixtures = {'ETHI_d': ['LEB_d', 'TUR_d', 'IND_d'],
#               'MOR_d': ['LEB_d', 'TUR_d'],
#               'TUR_k': ['TUR_d'],
#               'UZB_k': ['UZB_d', 'TUR_k'],
#               'MOR_k': ['MOR_d', 'TUR_k'],
#               'LEB_k': ['LEB_d', 'TUR_k']}


# Steps to optimize admixture in a sequence
admixture_steps = [[0, 1], [2, 3], [4, 5]]

# ------- Optimisation -------

res = migadmi(pop_names=pop_names,
              tree=tree,
              admixtures=admixtures,
              admixture_steps=admixture_steps,
              dist_matrix=d_mx)

# ------- Get results -------