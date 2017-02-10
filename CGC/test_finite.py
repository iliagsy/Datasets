# coding: utf-8
# pearson matrix是否都是有限数
# 对角元为0
from consts import GRAPH_FILE_NPY, SPECIES
from numpy import *
import numpy as np

skip = []  # ['fly']

for sp in SPECIES:
    if sp in skip: continue
    arr = load(GRAPH_FILE_NPY[sp])
    assert np.all(diag(arr) == 0.)  # yes
    print sp, sum(isnan(arr)) / float(arr.shape[0]**2)
    del arr
# fly 0.0
# worm 0.025
# human 0.007
