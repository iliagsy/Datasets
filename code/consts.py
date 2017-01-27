# coding: utf-8
import numpy as np

BASE_DIR = '/Users/gaoshiyu/desktop/LB/Datasets'
BASE_DIR1 = '/Volumes/GSY_SS/LB/'

SPECIES = ['fly', 'worm', 'human']
SPECIES_REL = ['worm', 'fly', 'human']

ORI_GENE_FILE = dict(zip(
    SPECIES,
    map(lambda s: BASE_DIR + '/LB-{}_gene.csv'.format(s), SPECIES)
))

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

GRAPH_FILE_NPY = dict(zip(
    SPECIES,
    map(lambda s: GRAPH_FILE[s].replace('.csv', '.npy'),
        SPECIES)
))

GRAPH_FILE_GZ = dict(zip(
    SPECIES,
    map(lambda sp: GRAPH_FILE[sp] + '.gz',
        SPECIES)
))

FINAL_FILE = dict(zip(
    SPECIES,
    map(lambda s: (BASE_DIR
                   + _fn_fmt.replace('cut_', 'final_').replace('.csv', '') % s),
        SPECIES)
))

def FINAL_GRAPH(species):
    return np.load(FINAL_FILE[species]+'.npz')['arr']

_GRAPH_PRF = BASE_DIR + '/profile_graph.json'
GRAPH_PRF = dict(zip(
    SPECIES,
    map(lambda s: (BASE_DIR
                   + '/profile_graph_%s.json' % s),
        SPECIES)
))

ZIDX_FILE = dict(zip(
    SPECIES,
    map(lambda s: BASE_DIR + '/zidx_{}'.format(s),
        SPECIES)
))

ZIDX = {}
for sp in SPECIES:
    ZIDX[sp] = np.load(ZIDX_FILE[sp]+'.npy').tolist()

REL_FILE = BASE_DIR + '/LB-Modencode.merged.orth20120611_wfh_comm_all.csv'

def REL_MAT_FILE(sp1, sp2):
    assert sp1 in SPECIES and sp2 in SPECIES
    return BASE_DIR + '/rel_mat_{}_{}'.format(sp1, sp2)
