# coding: utf-8
import numpy as np
from numpy import *
from datetime import *

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
    map(lambda s: BASE_DIR
                  + _fn_fmt.replace('cut_', 'graph_').replace('.csv', '') % s,
        SPECIES)
))

GRAPH_FILE_NPY = dict(zip(
    SPECIES,
    map(lambda s: GRAPH_FILE[s] + '.npy',
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

def _GRAPH_PRF():
    return BASE_DIR + '/profile_graph_{}.json'.format(datetime.now().date())

def GRAPH_PRF(species):
    assert species in SPECIES
    return BASE_DIR + '/profile_graph_{}_{}.json'.format(species, datetime.now().date())

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
    return BASE_DIR + '/rel_mat_{}_{}'.format(sp1, sp2)

def REL_MAT(t):
    sp1, sp2 = map(SPECIES.__getitem__, t)
    try:
        return load(REL_MAT_FILE(sp1, sp2)+'.npz')['arr']
    except:
        sp1,sp2 = sp2,sp1
        return load(REL_MAT_FILE(sp1, sp2)+'.npz')['arr'].T
