import json
from numpy import *
from consts import GRAPH_FILE, SPECIES, NPY_FILE, BASE_DIR

from Alg.pearson import Pearson as _P


def get_graph_matrix(data_arr, savefn=None):
# a wrapper
    arr = _P.PearsonMat(data_arr)
    for r in range(arr.shape[0]):
        arr[r, r] = 0.
    if savefn:
        savetxt(savefn, arr)
    return arr


def gen_graph(species):
    assert species.lower() in SPECIES
    data_arr = load(NPY_FILE[species.lower()])
    res_arr = get_pearson_matrix(data_arr,
                                 savefn=GRAPH_FILE[species] + test_str)


if __name__ == '__main__':
    for sp in SPECIES:
        gen_graph(sp)
