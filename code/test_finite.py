# coding: utf-8
# pearson matrix是否都是有限数
# 对角元为0
from consts import GRAPH_FILE, SPECIES
from numpy import *
import numpy as np

skip = []  # ['fly']

for sp in SPECIES:
    if sp in skip: continue
    arr = loadtxt(GRAPH_FILE[sp])
    assert np.all(diag(arr) == 0.)
    # print sp, sum(isnan(arr)), sum(isinf(arr))
    del arr
