"""
Utilities
"""
import re
from itertools import chain

def showl(a):
    """
    Print list by lines
    :param a:
    :return:
    """
    for tmp in a:
        print(tmp)


def read_admixture_old(file):
    """
    Read admixture file by a pattern
    :param file: each row in the file represent an admixture,
    first word in the line represents the mixed populations,
    the rest words - sources
    :return: list of mixed population and list of their sources
    """

    pattern = r'\b\w+\b'

    with open(file, 'r') as f:
        lines = f.readlines()

    pop_admix, pop_sources_names = zip(
        *[[admix[0], admix[1:len(admix)]] for admix in [re.findall(pattern, line) for line in lines]])

    return list(pop_admix), list(pop_sources_names)


def read_admixture(file):
    """
    Read admixture file by a pattern
    :param file: each row in the file represent an admixture,
    first word in the line represents the mixed populations,
    the rest words - sources
    :return: list of mixed population and list of their sources
    """

    pattern = r'\b\w+\b'

    with open(file, 'r') as f:
        lines = f.readlines()

    pop_admix, pop_sources_names = zip(
        *[[admix[0], admix[1:len(admix)]] for admix in [re.findall(pattern, line) for line in lines]])

    admixtures = {i:x for i, x in zip(pop_admix, pop_sources_names)}

    return admixtures


def check_admixtures(pop_names, pop_sources, admixture_steps, echo=False):
    """
    :param pop_names:
    :param pop_sources:
    :param admixture_steps:
    :return:
    """

    sources_tmp = set(range(len(pop_names)))
    for istep in sum(admixture_steps, []):
        sources_step = set(pop_sources[istep])
        if len(sources_step - sources_tmp) > 0:
            raise TypeError('Sources and Mixtures are not consistent')

        sources_tmp |= set([istep + len(pop_names)])

        if not echo:
            continue
        print(sources_tmp)
        print(sources_step)
        print(istep + len(pop_names))

