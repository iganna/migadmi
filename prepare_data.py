"""
Prepare windows of distance marices
"""

import numpy as np
import glob
import os

from pandas import DataFrame
from itertools import combinations


# ------ Setup ------

path_blocks = 'data/in/'
path_estim = 'data/estim/'
file_snps = 'data/snps_names.txt'

if not os.path.exists(path_blocks):
    # Create a new directory because it does not exist
    os.makedirs(path_blocks)

# Setup for blocks of SNPs

s_step = 10**6
s_window = 3*s_step

# Setup for frequencies

pop_region = ['ETHI', 'LEB', 'TUR',  'UZB', 'IND', 'MOR']
pop_type = ['desi', 'kabuli']
movement_type = 'routes'
tol = 0.01


# ------ Get SNPs coordinates ------

snps_names = []
with open(file_snps, 'r') as file:
    for line in file:  # reading each line
        for word in line.split():  # reading each word
            snps_names += [word]

snps_pos = []
snps_chr = []
for sn in snps_names:
    if (sn[0] == 'C') and (sn[1] == 'a'):

        tmp = sn.split('.')
        chr = int(tmp[0][2])
        pos = int(tmp[1].split('_')[0])
        snps_pos += [pos]
        snps_chr += [chr]
    else:
        snps_pos += [0]
        snps_chr += [0]

snps_chr = np.array(snps_chr)
snps_pos = np.array(snps_pos)


# ------ Get blocks ------

blocks = []
block_coord = []  #coordinates
for i_chr in range(1,(max(snps_chr)+1)):
    if i_chr == 0:
        break
    idx_chr = (snps_chr == i_chr)
    pos_max = max(snps_pos[idx_chr])
    pos = snps_pos[idx_chr][0]

    while (pos <= pos_max):
        tmp = idx_chr & (snps_pos >= pos) & (snps_pos < pos+s_window)
        pos_block = np.where(tmp)[0]
        pos = pos + s_step

        if len(pos_block) < 10:
            continue
        blocks += [pos_block]
        block_coord += [[i_chr, round(pos - s_step/2)]]

with open(f'{path_blocks}block_coord.txt', 'w') as f:
    for i_chr, pos in block_coord:
    f.write(f'{i_chr}\t{pos}\n')


# ------ Read Frequencies ------

ilrs = []
pop_names = []
for p_region in pop_region:
    for p_type in pop_type:
        pattern = f'{path_estim}{p_region}_{movement_type}_{p_type}*'
        files = glob.glob(pattern)
        if len(files) == 0:
            continue
        pop_names += [f'{p_region}_{p_type[0]}']
        freqs = np.array([[float(x) for x in open(f, 'r')] for f in files])

        freqs[freqs >= 1 - tol] = 1 - tol
        freqs[freqs < tol] = tol

        ilr = np.log((1 - freqs) / freqs)  # ilr
        ilr = ilr.mean(axis=0)  # mean over mcmc iterations

        ilrs += [ilr]
n_pop = len(pop_names)
ilrs = np.array(ilrs)


# ------ Normalised ilr balances ------

ilrs_mean = np.mean(ilrs, axis=0)
m_fa = 1 / (1 + np.exp(ilrs_mean))  # geometric mean version
m_fa = np.sqrt( m_fa * (1 - m_fa))
ilrs_norm = ilrs * m_fa  # multiplication, yes, it's true by formulas


# ------ Calculate distance matrices in windows ------

for ib, b in enumerate(blocks):
    freqs_blocks = ilrs_norm[:, b]
    mx = np.zeros((n_pop, n_pop))
    for i, j in list(combinations(range(n_pop),2)):
        mx[i, j] = mx[j, i] = np.linalg.norm(freqs_blocks[i] - freqs_blocks[j]) / np.sqrt(len(b))
    mx = DataFrame(mx)
    mx.index = pop_names
    mx.columns = pop_names
    mx.to_csv(f'{path_blocks}wnd_{ib+1}.txt', sep = '\t')


