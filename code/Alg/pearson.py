# coding: utf-8
import numpy as np
from numpy import *
from scipy.stats import pearsonr

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
    def PearsonMat0(cls, arr, fill=True, ret_err=False):
        # fill: fill the matrix, making it symmetric
        # ret_err: return error log
        errs = []
        E = cls.Exp(arr)
        # Var = cls.Var(arr, E=E)
        Sd = cls.Sd(arr, E=E)
        res_arr = zeros([arr.shape[0], arr.shape[0]])
        for i in range(arr.shape[0]):
            for j in range(i, arr.shape[0]):
                try:
                    E_i_j = cls.Exp(arr[i] * arr[j])
                    r = (E_i_j - E[i] * E[j]) / (Sd[i] * Sd[j])
                except:
                    errs.append([i, j])
                res_arr[i][j] = r
                if fill: res_arr[j][i] = r
        if ret_err:
            return res_arr, errs
        return res_arr


    @classmethod
    def H(cls, arr):
        SX = add.reduce(arr, axis=1)
        SX.shape = (1, SX.shape[0])  # 0-dim array -> 1-dim array
        S = np.dot(SX.T, SX)
        N = arr.shape[1]
        H = np.dot(arr, arr.T) - S / N
        return H


    @classmethod
    def PearsonMat1(cls, arr):
        H = cls.H(arr)
        D = np.diag(H)
        D.shape = (1, D.shape[0])
        P = H / sqrt(np.dot(D.T, D))
        return P


    @classmethod
    def PearsonMat(cls, arr, fill=True, ret_err=False, ret_p=False):
        res_arr = zeros([arr.shape[0], arr.shape[0]])
        p = zeros([arr.shape[0], arr.shape[0]])
        for i in range(arr.shape[0]):
            for j in range(i, arr.shape[0]):
                res_arr[i, j], p[i, j] = pearsonr(arr[i], arr[j])
                if fill: res_arr[j, i] = res_arr[i, j]
        if ret_p:
            return res_arr, p
        return res_arr
