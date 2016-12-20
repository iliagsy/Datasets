# coding: utf-8

BASE_DIR = '/home/gaoshiyu/Datasets'

SPECIES = ['human', 'fly', 'worm']

_fn_fmt = '/cut_%s_gene.csv'
HM_FILE, FLY_FILE, WORM_FILE = map(lambda s: BASE_DIR + _fn_fmt % s,
                                   SPECIES)

NPY_FILE = dict(zip(
    SPECIES,
    map(lambda s: BASE_DIR + _fn_fmt.replace('.csv', '.npy') % s,
        SPECIES)
))

GRAPH_FILE = dict(zip(
    SPECIES,
    map(lambda s: BASE_DIR + _fn_fmt.replace('cut_', 'graph_') % s,
        SPECIES)
))
