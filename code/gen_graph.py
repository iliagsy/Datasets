import json
from numpy import *
from consts import GRAPH_FILE, SPECIES, NPY_FILE, BASE_DIR

from Alg.pearson import Pearson as _P


def get_pearson_matrix(data_arr, res_arr, savefn=None, ret_err=False, test=False):
    errs = []
    assert res_arr.shape == (data_arr.shape[0],
                             data_arr.shape[0])
    for i in range(data_arr.shape[0]-1):
        for j in range(i+1, data_arr.shape[0]):
            try:
                r = _P.Pearson(data_arr[i], data_arr[j])
            except:
                errs.append([i, j])
            res_arr[i][j] = res_arr[j][i] = r
    if savefn:
        savetxt(savefn, res_arr)
    if ret_err:
        return res_arr, errs
    else:
        return res_arr


def gen_graph(species, test=False):
    assert species.lower() in SPECIES
    data_arr = load(NPY_FILE[species]); l = data_arr.shape[0]
    if test:
        data_arr = data_arr[:1000]; l = data_arr.shape[0]
    res_arr = zeros([l, l])
    test_str = '.test' if test else ''
    res_arr, errs = get_pearson_matrix(data_arr, res_arr,
                                       savefn=GRAPH_FILE[species] + test_str,
                                       ret_err=True)
    with open(BASE_DIR + '/code/err_log_%s%s.txt' % (species, test_str), 'w') as fh:
        fh.write(json.dumps(errs, indent=4))


# for sp in SPECIES:
gen_graph('human', test=False)
