# coding: utf-8
from numpy import *
import random
from datetime import datetime

from consts import FINAL_GRAPH, GRAPH_FILE_NPY, SPECIES, REL_MAT, BASE_DIR
from Alg.CGC import CGC


def cluster():
    A_lst = map(FINAL_GRAPH, SPECIES)
    # Hs = load(BASE_DIR + '/clusters_2017-02-01T01:25:00.npz')
    # H_lst = map(Hs.__getitem__, SPECIES)
    # del Hs
    # A_lst = []
    # for sp in SPECIES:
    #     A_lst.append(load(GRAPH_FILE_NPY[sp]))
    S_dct = {}
    lambda_dct = {}
    for i in range(len(SPECIES)):
        for j in range(i+1, len(SPECIES)):
            t = (i,j)
            S_dct.update({t: REL_MAT(t).T})
            lambda_dct.update({t: random.uniform(0,1)+1})
    cgc = CGC(A_lst, S_dct, [100,100,100], lambda_dct, 'RSS')
    print cgc.iter()
    H_lst = cgc.H_lst
    data = dict(zip(SPECIES, H_lst))
    savez_compressed(BASE_DIR + '/clusters_{}'.format(datatime.now().isoformat()), **data)


if __name__ == '__main__':
    cluster()

# iter 27 err 4.81663557058e-05
