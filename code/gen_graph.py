# coding: utf-8
from numpy import *
from consts import (GRAPH_FILE, SPECIES, SPECIES_REL, NPY_FILE, BASE_DIR,
                    ZIDX_FILE, REL_FILE, ORI_GENE_FILE, REL_MAT_FILE)


from Alg.pearson import Pearson as _P


def get_graph_matrix(data_arr, pearsonFunc, savefn=None):
# a wrapper
    arr = pearsonFunc(data_arr)
    for r in range(arr.shape[0]):
        arr[r, r] = 0.
    if savefn:
        save(savefn, arr)
    return arr


def gen_graph(species, pf=_P.PearsonMat0):
    assert species.lower() in SPECIES
    data_arr = remove_zero(species)
    res_arr = get_graph_matrix(data_arr,
                               pearsonFunc=pf,
                               savefn=GRAPH_FILE[species])


def gen_relation_mat(species, savefn=None):
    assert len(species) == 2
    cols = map(SPECIES_REL.index, species)
    cvt = {}
    for i in range(2):
        cvt[i] = lambda s: s.strip('\"')
    arr = loadtxt(REL_FILE, delimiter=',', skiprows=1, dtype=str,
                  usecols=cols, converters=cvt)

    names = []
    for sp in species:
        names.append(
            loadtxt(ORI_GENE_FILE[sp], delimiter=',', skiprows=1,
                    dtype=str, usecols=[0]
                    ).tolist()
        )

    rel_mat = zeros([len(names[0]), len(names[1])])
    for r in range(arr.shape[0]):
        try:
            i = names[0].index(arr[r][0])
            j = names[1].index(arr[r][1])
        except ValueError:
            continue
        rel_mat[i, j] += 1

    # normalize non-zero rows
    Sum = add.reduce(rel_mat, axis=1, keepdims=True)
    Sum[Sum == 0.] = inf
    rel_mat = rel_mat / Sum

    if savefn:
        savez_compressed(savefn, rel_mat)
    return rel_mat


def remove_zero(species):
    data_arr = load(NPY_FILE[species.lower()])
    zidx_ = add.reduce(data_arr != 0., axis=1) == 0.
    zidx = arange(data_arr.shape[0])[zidx_]
    save(ZIDX_FILE[species], zidx)
    return data_arr[logical_not(zidx_), :]


def test_zero(species):
    assert species.lower() in SPECIES
    data_arr = load(NPY_FILE[species.lower()])
    print sum(add.reduce(data_arr != 0., axis=1) == 0.) / float(data_arr.shape[0])
# 0.013 worm
# 0.003 human


if __name__ == '__main__':
    # for i in range(len(SPECIES)):
    #     for j in range(i+1, len(SPECIES)):
    #         sp1, sp2 = SPECIES[i], SPECIES[j]
    #         arr = gen_relation_mat([sp1, sp2],
    #                                savefn=REL_MAT_FILE(sp1, sp2))
    for sp in SPECIES:
        gen_graph(sp)
