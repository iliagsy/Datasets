import json
from numpy import *
from consts import GRAPH_FILE, SPECIES, NPY_FILE, BASE_DIR

from Alg.pearson import Pearson as _P


def get_pearson_matrix(data_arr, savefn=None, ret_err=False):
    t = _P.PearsonMat(data_arr, ret_err=ret_err)
    if savefn:
        savetxt(savefn, t if not isinstance(t, tuple) else t[0])
    return t


def gen_graph(species, test=False):
    assert species.lower() in SPECIES
    data_arr = load(NPY_FILE[species.lower()])
    if test:
        data_arr = data_arr[:1000]
    test_str = '.test' if test else ''
    res_arr, errs = get_pearson_matrix(data_arr,
                                       savefn=GRAPH_FILE[species] + test_str,
                                       ret_err=True)
    with open(BASE_DIR + '/code/err_log_%s%s.txt' % (species, test_str), 'w') as fh:
        fh.write(json.dumps(errs, indent=4))


# for sp in SPECIES:
gen_graph('human', test=False)
