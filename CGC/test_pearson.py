# coding: utf-8
# Done! Do not repeat this operation!

from numpy import *
import numpy as np
from Alg.pearson import Pearson as _P
from consts import (BASE_DIR, NPY_FILE, GRAPH_FILE_NPY, ZIDX,
                    GRAPH_FILE, SPECIES, FINAL_GRAPH)
import json

def test_pearson():
    A = random.rand(100, 20) * 100 + 15
    p1 = _P.PearsonMat1(A)
    p0 = _P.PearsonMat0(A)
    p = _P.PearsonMat(A)
    print sum(p0>1), sum(p>1), sum(p1>1)  # 20+, 0, 0
    print np.max(abs(p1-p)), np.max(abs(p0-p))  # <10e-14


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


def test_Var():
    from scipy.stats import describe
    from scipy import stats
    A = random.rand(20)
    print _P.Var(A)
    print np.var(A)
    print describe(A).variance


def plotPearsonMat():
    import matplotlib.pyplot as plt
    for sp in SPECIES:
        arr = load(GRAPH_FILE_NPY[sp])
        plt.hist(arr.flatten(), bins=int(2/0.001), histtype='bar')
        plt.savefig(BASE_DIR + '/histPearsonMatrix_{}.pdf'.format(sp))


if __name__ == '__main__':
    test_pearson()
