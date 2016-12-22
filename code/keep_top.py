# coding: utf-8
from numpy import *

from consts import GRAPH_FILE_GZ, SPECIES, FINAL_FILE

def keep_top(species, d=5, savefn=None):
    PM = loadtxt(GRAPH_FILE_GZ[species])  # Pearson Matrix
    for i in range(PM.shape[0]):
        kv_lst_ = sorted(
            zip(range(len(PM[i])), PM[i].tolist()),
            key=lambda t: t[1],
            reverse=True
        )
        kv_lst_[d:] = 0
        kv_lst = sorted(kv_lst_, key=lambda t: t[0])
        val_lst = [t[1] for t in kv_lst]
        PM[i] = array(val_lst)
    if savefn:
        savetxt(savefn, PM)
        save(savefn, PM)
    return PM

for sp in SPECIES:
    keep_top(sp, savefn=FINAL_FILE[sp])
