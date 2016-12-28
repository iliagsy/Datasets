# coding: utf-8
import numpy as np
from numpy import *

from consts import GRAPH_FILE_GZ, SPECIES, FINAL_FILE, GRAPH_PRF

def keep_top(species, PM=None, d=5, savefn=None, ret_arr=True):
    # PM: 2darray
    if PM is None:
        PM = loadtxt(GRAPH_FILE_GZ[species])  # Pearson Matrix
    # keep top d
    # PM_f: PM filtered -- bool type
    PM_f = ndarray(shape=PM.shape, dtype=bool)
    for i in range(PM.shape[0]):
        kv_arr_ = array(sorted(
            zip(range(len(PM[i])), PM[i].tolist()),
            key=lambda t: t[1],
            reverse=True
        ))
        kv_arr_[d:, 1] = 0
        kv_arr = array(sorted(kv_arr_.tolist(), key=lambda t: t[0]))
        val_arr = kv_arr[:, 1]
        PM_f[i] = (val_arr != 0)
    # PM_b: PM binary, i.e. complemented PM_f
    PM_b = logical_or(PM_f.transpose(), PM_f)
    nE = sum(PM_b) / 2; nV = PM.shape[0]
    res = {
        'd': d,
        'species': species,
        'nE': nE,
        'nV': nV,
        'nE_per_V': float(nE) / nV
    }
    if not ret_arr and savefn is None:
        return res
    # calc new PM, with PM_b as condition
    # numpy.where(condition[, x, y]) -- Return elements, either from x or y, depending on condition.
    PM = where(PM_b, PM, zeros_like(PM))
    # wrap up
    if savefn:
        savetxt(savefn, PM)
        save(savefn, PM)
    res.update({'PM': PM})
    return res

for d in (5, 8, 10):
    for sp in SPECIES:
        res_l = []
        data = keep_top(sp, d=d, ret_arr=False)
        data = map(data.get, ['d','species','nE','nV','nE_per_V'])
        s = '%4d%10s%8d%8d%8.4f' % tuple(data)
        res_l.append(s)
    with open(GRAPH_PRF, 'w') as fh:
        fh.write('\n'.join(res_l))
