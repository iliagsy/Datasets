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
        PM_f = keep_top_d(PM, d)
        nE = sum(PM_f!=0)
        prf = {
            'd': d,
            'species': species,
            'nE': nE,
            'nV': nV,
            'nE_per_V': float(nE) / nV
        }
        if not ret_arr and not savefn:
            prfs.append(prf); continue
        if savefn:
            savetxt(savefn, PM_f); save(savefn, PM_f)
        prf.update({'PM': PM_f})
        prfs.append(prf)
    return prfs


def keep_top_d(G, d):
    for r in range(G.shape[0]):
        idxs = abs(G[r]).argsort()[::-1]
        G[r, idxs[d:]] = 0.
        G[r, abs(G[r]) < 10**(-14)] = 0.
    for r in range(G.shape[0]):
        for c in range(r+1, G.shape[0]):
            u,l = G[r,c], G[c,r]
            if u != 0 and l == 0:
                G[c,r] = u
            elif l != 0 and u == 0:
                G[r,c] = l
    return G


def gen_graph_profile():
    skip = []
    for sp in SPECIES:
        if sp in skip: continue
        prf = keep_top(sp, dl=[5,8,10], ret_arr=False)
        with open(GRAPH_PRF[sp], 'w') as fh:
            fh.write(json.dumps(prf, indent=4))


if __name__ == '__main__':
    gen_graph_profile()
