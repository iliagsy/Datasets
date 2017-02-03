# coding: utf-8
from numpy import *
from scipy.sparse import csr_matrix
import scipy.sparse.linalg as scila
import numpy as np
from datetime import *
import logging
import logging.config
from consts import *

logging.config.fileConfig(BASE_DIR+'/code/logging.conf')
logger = logging.getLogger('CGC')


class CGC(object):
    def __init__(self, A_lst, S_dct, k_lst, lambda_arr, lossFunc, H_lst=None):
        d = self.d = len(A_lst)
        assert (len(k_lst) == d
                and len(S_dct) == d*(d-1)/2
                and lambda_arr.shape[0] == d
                and lossFunc in ('CD', 'RSS'))
        self.n_lst = map(lambda arr: arr.shape[0], A_lst)
        for t in S_dct.keys():
            assert S_dct[t].shape == (self.n_lst[t[1]], self.n_lst[t[0]])
        self.A_lst = map(CGC._normalize, A_lst)
        self.S_dct = S_dct
        self.k_lst = k_lst
        self.lambda_arr = lambda_arr
        if H_lst is None:
            self.H_lst = []
            for i in range(d):
                H = -random.rand(self.n_lst[i], self.k_lst[i]) + 1
                self.H_lst.append(matrix(H))
        else:
            self.H_lst = H_lst
        self.lossFunc = lossFunc

    def iter(self, itertimes=100, TOL=0.5 * 10**(-5)):
        ob = self.objective
        for itr in range(itertimes):
            for i in range(self.d):
                H = np.copy(self.H_lst[i])
                Psi = zeros_like(H)
                Xi = zeros_like(H)
                self._increParam(Xi, Psi, H, i)
                Psi += self.A_lst[i].dot(H)
                Xi += H.dot(H.T).dot(H)
                assert all(Psi/Xi >= 0)  # check
                self.H_lst[i] = H * (Psi / Xi)**(1./4)
                del Psi, Xi, H
            ob1 = self.objective; err = ob1-ob
            ob = ob1
            logger.warn('iter {} obj {} err {}'.format(itr, ob, err))  # 用来画图
            if abs(err) < TOL:
                break
        return itr, err, ob

    @property
    def objective(self):
        ob = 0.
        for i in range(self.d):
            ob += linalg.norm(self.A_lst[i] - self.H_lst[i].dot(self.H_lst[i].T))**2
        for i,j in self.S_dct.keys():
            if self.lossFunc == 'RSS':
                F = self.S_dct[(i,j)].dot(self.H_lst[i]) - self.H_lst[j]
            elif self.lossFunc == 'CD':
                F_ = self.S_dct[(i,j)].dot(self.H_lst[i])
                F = F_.dot(F_.transpose()) - self.H_lst[j].dot(self.H_lst[j].transpose())
                del F_
            ob += self.lambda_arr[i, j] * linalg.norm(F)**2
        return ob

    def _increParam(self, Xi, Psi, H, pi_):
        for j in range(self.d):
            if j == pi_: continue
            small, large = map(lambda f: f([pi_, j]), [min, max])
            Lambda = self.lambda_arr[small, large]
            S = self.S_dct[(small, large)]

            if self.lossFunc == 'RSS':
                if j < pi_:
                    Xi += Lambda/2 * H
                elif j > pi_:
                    Xi += Lambda/2 * S.transpose().dot(S).dot(H)
                    S = S.transpose()
                Psi += Lambda / 2 * S.dot(self.H_lst[j])
            elif self.lossFunc == 'CD':
                # increment Xi
                if j < pi_:
                    F1 = H
                elif j > pi_:
                    F1 = S.transpose().dot(S).dot(H)
                    S = S.transpose()
                Xi += Lambda * F1.dot(F1.transpose()).dot(H)
                # increment Psi
                F2 = S.dot(self.H_lst[j])
                Psi += Lambda * F2.dot(F2.transpose()).dot(H)

    @classmethod
    def _normalize(cls, A):
    # normalize by Frobenius norm
        return A / scila.norm(A)

    @classmethod
    def _normH(cls, H):
        return H / add.reduce(H, axis=1, keepdims=True)
