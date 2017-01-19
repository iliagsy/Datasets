# coding: utf-8
from numpy import *

class Pearson(object):

    @classmethod
    def Exp(cls, arr):
        # 期望
        # arr: vec / 2darray
        d = arr.ndim - 1  # 最高一维
        return add.reduce(arr, d) / arr.shape[d]

    @classmethod
    def Var(cls, arr=None, E=None, E_2=None):
        # 方差
        # arr: vec / 2darray
        assert arr is not None or (E is not None
                                   and E_2 is not None)
        if E is None:
            E = cls.Exp(arr)
        if E_2 is None:
            E_2 = cls.Exp(arr ** 2)
        return E_2 - E ** 2

    @classmethod
    def Sd(cls, arr=None, E=None, E_2=None):
        # 标准差
        Var = cls.Var(arr, E, E_2)
        return sqrt(Var)

    @classmethod
    def PearsonMat(cls, arr, fill=True, ret_err=False):
        # fill: fill the matrix, making it symmetric
        # ret_err: return error log
        errs = []
        E = cls.Exp(arr)
        # Var = cls.Var(arr, E=E)
        Sd = cls.Sd(arr, E=E)
        res_arr = zeros([arr.shape[0], arr.shape[0]])
        for i in range(arr.shape[0]-1):
            for j in range(i+1, arr.shape[0]):
                try:
                    E_i_j = cls.Exp(arr[i] * arr[j])
                    # r = (E_i_j - E[i] * E[j]) / sqrt((Var[i] * Var[j]))
                    r = (E_i_j - E[i] * E[j]) / (Sd[i] * Sd[j])
                except:
                    errs.append([i, j])
                res_arr[i][j] = r
                if fill: res_arr[j][i] = r
        if ret_err:
            return res_arr, errs
        return res_arr
