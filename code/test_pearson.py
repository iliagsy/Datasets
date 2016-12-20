# coding: utf-8

from Alg.pearson import Pearson
from consts import BASE_DIR, NPY_FILE
import json
A = load(NPY_FILE['human'])
errs = []
for i in range(A.shape[0]):
    for j in range(i+1, A.shape[0]):
        try:
            p = Pearson.Pearson(A[i], A[j])
        except:
            errs.append([i, j])

with open(BASE_DIR + '/code/err_log.txt', 'w') as fh:
    json.dumps(errs, indent=4)
