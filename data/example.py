def example_admixtures(kabuli_tur=True):
    # Initial hypotheses
    pop_admix = ['ETHI_d', 'MOR_d']
    pop_sources_names = [['LEB_d', 'TUR_d', 'IND_d'],
                         ['LEB_d', 'TUR_d']]

    # Kabuli source: two alternatives

    if kabuli_tur:
        pop_admix += ['TUR_k', 'UZB_k']
        pop_sources_names += [['TUR_d'],
                              ['UZB_d', 'TUR_k']]
    else:
        pop_admix += ['UZB_k', 'TUR_k']
        pop_sources_names += [['UZB_d'],
                              ['TUR_d', 'UZB_k']]

    # additional Kabulis
    pop_admix += ['MOR_k', 'LEB_k']
    pop_sources_names += [['MOR_d', 'TUR_k'],
                          ['LEB_d', 'TUR_k']]

    return (pop_admix, pop_sources_names)