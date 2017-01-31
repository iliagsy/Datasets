# coding: utf-8
from numpy import *
import random

from consts import FINAL_GRAPH, SPECIES, REL_MAT, BASE_DIR
from Alg.CGC import CGC


def cluster():
    A_lst = map(FINAL_GRAPH, SPECIES)
    S_dct = {}
    lambda_dct = {}
    for i in range(len(SPECIES)):
        for j in range(i+1, len(SPECIES)):
            t = (i,j)
            S_dct.update({t: REL_MAT(t).T})
            lambda_dct.update({t: random.uniform(0,1)*2})
    cgc = CGC(A_lst, S_dct, [100,100,100], lambda_dct, 'RSS')
    print cgc.iter()
    H_lst = cgc.H_lst
    data = dict(zip(SPECIES, H_lst))
    savez_compressed(BASE_DIR + '/clusters', **data)


if __name__ == '__main__':
    cluster()
