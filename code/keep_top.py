# coding: utf-8
import numpy as np
from numpy import *

from consts import GRAPH_FILE_GZ, SPECIES, FINAL_FILE

def keep_top(species, d=5, savefn=None):
    PM = loadtxt(GRAPH_FILE_GZ[species])  # Pearson Matrix
    for i in range(PM.shape[0]):
        kv_arr_ = array(sorted(
            zip(range(len(PM[i])), PM[i].tolist()),
            key=lambda t: t[1],
            reverse=True
        ))
        kv_arr_[d:] = 0
        kv_arr = array(sorted(kv_arr_.tolist(), key=lambda t: t[0]))
        val_arr = kv_arr[:, 1]
        PM[i] = val_arr
    if savefn:
        savetxt(savefn, PM)
        save(savefn, PM)
    return PM

# for sp in SPECIES:
sp = ''
keep_top(sp, savefn=FINAL_FILE[sp])
