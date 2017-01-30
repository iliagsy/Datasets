# coding: utf-8
from numpy import *
import numpy as np


class CGC(object):
    def __init__(self, A_lst, S_dct, k_lst, lambda_dct, lossFunc):
        d = self.d = len(A_lst)
        try:
            assert (len(k_lst) == d
                    and len(S_dct) == d*(d-1)/2
                    and len(lambda_dct) == d*(d-1)/2
                    and lossFunc in ('CD', 'RSS'))
        except AssertionError:
            print len(k_lst), len(S_dct), len(lambda_dct), lossFunc; exit()
        self.n_lst = map(lambda arr: arr.shape[0], A_lst)
        for t in S_dct.keys():
            try:
                assert S_dct[t].shape == (self.n_lst[t[1]], self.n_lst[t[0]])
            except AssertionError:
                print t, S_dct[t].shape
        self.A_lst = map(CGC._normalize, A_lst)
        self.S_dct = S_dct
        self.k_lst = k_lst
        self.lambda_dct = lambda_dct
        self.H_lst = []
        for i in range(d):
            H = -random.rand(self.n_lst[i], self.k_lst[i]) + 1
            self.H_lst.append(H)
        self.lossFunc = lossFunc

    def iter(self, itertimes=100):
        for itr in range(itertimes):
            err = 0.

            for i in range(self.d):
                H = self.H_lst[i]
                Psi = zeros_like(H)
                Xi = zeros_like(H)
                self._increParam(Xi, Psi, H, i)
                Psi += dot(self.A_lst[i], H)
                Xi += H.dot(H.T).dot(H)
                self.H_lst[i] = H * (Psi / Xi)**(1/4)
                del Psi, Xi
                err = max(err, np.max(abs(H - self.H_lst[i])))
            if err < 10**(-6):
                return itr
        return -1

    def _increParam(self, Xi, Psi, H, pi_):
        for j in range(self.d):
            if j == pi_: continue
            small, large = map(lambda f: f([pi_, j]), [min, max])
            Lambda = self.lambda_dct[(small, large)]
            S = self.S_dct[(small, large)]

            if self.lossFunc == 'RSS':
                if j < pi_:
                    Xi += Lambda/2 * H
                elif j > pi_:
                    Xi += Lambda/2 * S.T.dot(S).dot(H)
                    S = S.T
                Psi += Lambda / 2 * dot(S, self.H_lst[j])
            elif self.lossFunc == 'CD':
                # increment Xi
                if j < pi_:
                    F1 = H
                elif j > pi_:
                    F1 = S.T.dot(S).dot(H)
                    S = S.T
                Xi += Lambda * F1.dot(F1.T).dot(H)
                # increment Psi
                F2 = S.dot(self.H_lst[j])
                Psi += Lambda * F2.dot(F2.T).dot(H)

    @classmethod
    def _normalize(cls, A):
    # 取绝对值 http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000117
    # normalize by Frobenius norm
        return abs(A) / sum(A**2)
