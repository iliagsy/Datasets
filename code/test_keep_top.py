# coding: utf-8
import numpy as np
import json
from numpy import *

from consts import GRAPH_FILE_GZ, SPECIES, GRAPH_PRF
from keep_top import keep_top_d

A = random.rand(10, 10)
A = A * 100 - 50

for d in range(1, 6):
    arr = keep_top_d(A, d)
    print d, float(sum(arr != 0)) / A.shape[0]
