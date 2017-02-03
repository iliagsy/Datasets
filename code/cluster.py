# coding: utf-8
from numpy import *
from scipy.sparse import csr_matrix
import random
import json
from datetime import datetime
import logging
import logging.config

from consts import FINAL_GRAPH, GRAPH_FILE_NPY, SPECIES, REL_MAT, BASE_DIR
from Alg.CGC import CGC


logging.config.fileConfig(BASE_DIR+'/code/logging.conf')
logger = logging.getLogger('cluster')


def cluster():
    A_lst = []
    for sp in SPECIES:
    # 取绝对值 http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000117
        A_lst.append(csr_matrix(abs(FINAL_GRAPH(sp))))
    # Hs = load(BASE_DIR + '/clusters_2017-02-01T01:25:00.npz')
    # H_lst = map(Hs.__getitem__, SPECIES)
    # del Hs
    S_dct = {}
    lambda_arr = zeros((len(SPECIES), len(SPECIES)))
    for i in range(len(SPECIES)):
        for j in range(i+1, len(SPECIES)):
            t = (i,j)
            S_dct.update({t: csr_matrix(REL_MAT(t)).transpose()})
            lambda_arr[i, j] = random.uniform(0,1)+1
    k_lst = [99,100,101]
    lossFunc = 'CD'
    cgc = CGC(A_lst, S_dct, k_lst, lambda_arr, lossFunc)
    itr, err, obj = cgc.iter(TOL=0.5*10**(-4))
    H_lst = cgc.H_lst
    data = dict(zip(SPECIES, H_lst))
    fn = BASE_DIR + '/clusters_{}'.format(datetime.now().isoformat())
    savez_compressed(fn, **data)
    params = {
        'itr': itr,
        'obj': obj,
        'err': err,
        'lambdas': lambda_arr.tolist(),
        'file': fn+'.npz',
        'lossFunc': lossFunc,
        'k_lst': k_lst,
    }
    logger.warn(json.dumps(params, indent=4))


if __name__ == '__main__':
    cluster()
