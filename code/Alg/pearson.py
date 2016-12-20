# coding: utf-8
from numpy import *

class Pearson(object):

    @classmethod
    def Pearson(cls, vec1, vec2):
        cov = cls.Cov(vec1, vec2)
        var1, var2 = map(cls.Var, [vec1, vec2])
        prsn = cov / sqrt(var1 * var2)
        return prsn

    @classmethod
    def Exp(cls, vec):
        return add.reduce(vec) / vec.shape[0]

    @classmethod
    def Var(cls, vec):
        E = cls.Exp(vec)
        return cls.Exp((vec - E) ** 2)  # 方差

    @classmethod
    def Cov(cls, vec1, vec2):
        assert vec1.shape[0] == vec2.shape[0]
        mult = vec1 * vec2
        E1, E2 = map(cls.Exp, [vec1, vec2])
        cov = cls.Exp(mult) - E1 * E2
        return cov
