# coding: utf-8
import numpy as np
from numpy import *

from consts import GRAPH_FILE_GZ, SPECIES, FINAL_FILE, GRAPH_PRF

def keep_top(species, PM=None, d=5, savefn=None, ret_arr=True):
    # PM: 2darray
    # if ret_arr==False: only return metadata
    if PM is None:
        PM = loadtxt(GRAPH_FILE_GZ[species])  # Pearson Matrix
    # keep top d
    # PM_f: PM filtered -- bool type
    PM_f = ndarray(shape=PM.shape, dtype=bool)
    for i in range(PM.shape[0]):
        sidx = PM[i].argsort()[::-1]  # sorted index array
        PM_f[i, sidx[:d]] = True
        PM_f[i, sidx[d:]] = False
    # PM_b: PM binary, i.e. complemented PM_f
    PM_b = logical_or(PM_f.transpose(), PM_f)
    # calc new PM, with PM_b as condition
    # numpy.where(condition[, x, y]) -- Return elements, either from x or y, depending on condition.
    PM = where(PM_b, PM, zeros_like(PM))  # array-wise `c ? a : b`
    nE = sum(abs(PM) > 10**(-14)); nV = PM.shape[0]
    res = {
        'd': d,
        'species': species,
        'nE': nE,
        'nV': nV,
        'nE_per_V': float(nE) / nV
    }
    if not ret_arr and savefn is None:
        return res
    # wrap up
    if savefn:
        savetxt(savefn, PM)
        save(savefn, PM)
    res.update({'PM': PM})
    return res

if __name__ == '__main__':
    res_l = []
    for sp in SPECIES:
        PM = loadtxt(GRAPH_FILE_GZ[sp])  # only load once for different d
        for d in (5, 8, 10):
            data = keep_top(sp, d=d, PM=PM, ret_arr=False)
            data = map(data.get, ['d','species','nE','nV','nE_per_V'])
            s = '%4d%10s%8d%8d%8.4f' % tuple(data)  # 不可以是list
            res_l.append(s)
    with open(GRAPH_PRF, 'w') as fh:
        fh.write('%4s%10s%8s%8s%8s\n' % ('d','species','nE','nV','nE_per_V'))
        fh.write('\n'.join(res_l))
