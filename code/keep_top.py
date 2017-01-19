# coding: utf-8
import numpy as np
from numpy import *

from consts import GRAPH_FILE_GZ, SPECIES, FINAL_FILE, GRAPH_PRF

def keep_top(species, PM=None, d=5, savefn=None, ret_arr=True):
    # PM: 2darray
    # if ret_arr = False, only return metadata
    if PM is None:
        PM = loadtxt(GRAPH_FILE_GZ[species])  # Pearson Matrix
    # keep top d
    # PM_f: PM-filter -- bool type
    PM_f = ndarray(shape=PM.shape, dtype=bool)
    for i in range(PM.shape[0]):
        sidx = PM[i].argsort()[::-1]  # sorted index array
        PM_f[i, sidx[:d]] = True
        PM_f[i, sidx[d:]] = False
    # calc new PM, with OR-complemented PM_f as condition
    PM = where(logical_or(PM_f.transpose(), PM_f),
               PM,
               zeros_like(PM))   # array-wise `c ? a : b`
    # wrap up
    if savefn:
        savetxt(savefn, PM)
        save(savefn, PM)
    nE = sum(abs(PM) > 10**(-14)); nV = PM.shape[0]
    res = {
        'd': d,
        'species': species,
        'nE': nE,
        'nV': nV,
        'nE_per_V': float(nE) / nV
    }
    if not ret_arr:
        return res
    res.update({'PM': PM})
    return res


if __name__ == '__main__':
    with open(GRAPH_PRF, 'w') as fh:
        fh.write('%5s%10s%10s%10s%10s\n' % ('d','species','nE','nV','nE_per_V'))
    for sp in SPECIES:
        PM = loadtxt(GRAPH_FILE_GZ[sp])  # only load once for different d
        for d in (5, 8, 10):
            data = keep_top(sp, d=d, PM=PM, ret_arr=False)
            prf = map(data.get, ['d','species','nE','nV','nE_per_V'])
            s = '%5d%10s%10d%10d%10.4f\n' % tuple(prf)  # 不可以是list
            with open(GRAPH_PRF, 'a') as fh:
                fh.write(s)
        del PM
