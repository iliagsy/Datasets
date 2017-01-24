from numpy import *
from consts import GRAPH_FILE, SPECIES, NPY_FILE, BASE_DIR, ZIDX_FILE

from Alg.pearson import Pearson as _P


def get_graph_matrix(data_arr, savefn=None):
# a wrapper
    arr = _P.PearsonMat0(data_arr)
    for r in range(arr.shape[0]):
        arr[r, r] = 0.
    if savefn:
        save(savefn, arr)
    return arr


def gen_graph(species):
    assert species.lower() in SPECIES
    data_arr = remove_zero(species)
    res_arr = get_graph_matrix(data_arr,
                               savefn=GRAPH_FILE[species])


def remove_zero(species):
    data_arr = load(NPY_FILE[species.lower()])
    zidx_ = add.reduce(data_arr != 0., axis=1) == 0.
    zidx = arange(data_arr.shape[0])[zidx_]
    save(ZIDX_FILE[species], zidx)
    return data_arr[logical_not(zidx_), :]


def test_zero(species):
    assert species.lower() in SPECIES
    data_arr = load(NPY_FILE[species.lower()])
    print sum(add.reduce(data_arr != 0., axis=1) == 0.) / float(data_arr.shape[0])
# 0.013 worm
# 0.003 human


if __name__ == '__main__':
    for sp in SPECIES:
        gen_graph(sp)
