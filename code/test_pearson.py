# coding: utf-8

from numpy import *
from Alg.pearson import Pearson as _P
from consts import BASE_DIR, NPY_FILE
import json
A = load(NPY_FILE['human'])
A = A[:2, :]

print A
savetxt('p1.txt', _P.PearsonMat1(A))
savetxt('p0.txt', _P.PearsonMat0(A))
savetxt('p.txt', _P.PearsonMat(A))
