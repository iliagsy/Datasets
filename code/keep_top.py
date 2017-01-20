# coding: utf-8
import numpy as np
import json
from numpy import *

from consts import GRAPH_FILE_GZ, SPECIES, FINAL_FILE, GRAPH_PRF

def keep_top(species, PM=None, dl=[5], fnl=None, ret_arr=True):
    # PM: 2darray
    # if ret_arr = False, only return metadata
    if fnl is None:
        fnl = []
        for i in range(len(dl)):
            fnl.append(None)
    else:
        assert len(dl) == len(fnl)

    prfs = []
    if PM is None:
        PM = loadtxt(GRAPH_FILE_GZ[species])  # Pearson Matrix
    nV = PM.shape[0]
    for d,savefn in zip(dl, fnl):
        # PM_f: PM-filter -- bool type
        PM_f = ndarray(shape=PM.shape, dtype=bool)
        for i in range(PM.shape[0]):
            sidx = abs(PM[i]).argsort()[::-1]  # sorted index array
            PM_f[i, sidx[:d]] = True
            PM_f[i, sidx[d:]] = False
            PM_f[i, abs(PM[i]) < 10**(-14)] = False
        nE = sum(logical_or(PM_f.T, PM_f))
        prf = {
            'd': d,
            'species': species,
            'nE': nE,
            'nV': nV,
            'nE_per_V': float(nE) / nV
        }
        if not ret_arr and not savefn:
            prfs.append(prf); del PM_f; continue

        # calc new PM
        PM = where(logical_or(PM_f.T, PM_f),
                   PM,
                   zeros_like(PM))   # array-wise `c ? a : b`
        del PM_f
        # wrap up
        if savefn:
            savetxt(savefn, PM); save(savefn, PM)
        prf.update({'PM': PM})
        prfs.append(prf)

    return prfs


def gen_graph_profile():
    skip =  []
    for sp in SPECIES:
        if sp in skip: continue
        prf = keep_top(sp, dl=[5,8,10], ret_arr=False)
        with open(GRAPH_PRF[sp], 'w') as fh:
            fh.write(json.dumps(prf, indent=4))


if __name__ == '__main__':
    gen_graph_profile()
