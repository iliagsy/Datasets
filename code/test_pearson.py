# coding: utf-8
# Done! Do not repeat this operation!

from numpy import *
from Alg.pearson import Pearson as _P
from consts import BASE_DIR, NPY_FILE, GRAPH_FILE_NPY, ZIDX, GRAPH_FILE, SPECIES
import json

def test_pearson():
    A = load(NPY_FILE['human'])
    A = A[:5, :]

    print A
    savetxt('p1.txt', _P.PearsonMat1(A))
    savetxt('p0.txt', _P.PearsonMat0(A))
    savetxt('p.txt', _P.PearsonMat(A))


def remove_zero(sp):
    arr = load(GRAPH_FILE_NPY[sp])
    idx = array(sorted(list(
        set(range(arr.shape[0])) - set(ZIDX[sp])
    )))
    arr_ = arr[idx, :][:, idx]
    del arr
    assert all(isfinite(arr_))
    save(GRAPH_FILE_NPY[sp].replace('.npy', '')+'.new', arr_)
    return arr_


if __name__ == '__main__':
    skip = []
    for sp in SPECIES:
        if sp in skip: continue
        remove_zero(sp)
