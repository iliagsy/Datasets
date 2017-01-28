from numpy import *
import numpy as np


class CGC(object):
    def __init__(self, A_lst, S_dct, k_lst, lambda_dct):
        d = self.d = len(A_lst)
        assert (len(k_lst) == d
                and len(S_dct) <= d*(d-1)/2
                and len(lambda_dct) == d*(d-1)/2)
        self.A_lst = map(cls._normalize, A_lst)
        self.n_lst = map(lambda arr: arr.shape[0], self.A_lst)
        self.S_dct = S_dct
        self.k_lst = k_lst
        self.lambda_dct = lambda_dct
        self.H_lst = []
        for i in range(d):
            H = -random.rand(n_lst[i], k_lst[i]) + 1
            self.H_lst.append(H)

    def iter_RSS(self, itertimes=100):
        for itr in range(itertimes):
            err = 0.
            for i in range(self.d):
                H = self.H_lst[i]
                Psi = dot(self.A_lst[i], H)
                Xi = H.dot(H.T).dot(H)
                for j in range(d):
                    if j == i: continue
                    small, large = map(lambda f: f([i, j]), [min, max])
                    Lambda = self.lambda_dct.get((small, large))
                    S = self.S_dct.get((small, large)), zeros_like(H))
                    if j < i:
                        Xi += Lambda/2 * H
                    elif j > i:
                        Xi += Lambda/2 * S.T.dot(S).dot(H)
                        S = S.T
                    Psi += Lambda / 2 * dot(S, self.H_lst[j])
                self.H_lst[i] = H * (Psi / Xi)**(1/4)
                err += mean(abs(H - self.H_lst[i]))
            if err < 10**(-6):
                return True
        return False

    def iter_CD(self, itertimes=100):
        for itr in range(itertimes):
            err = 0.
            for i in range(self.d):
                H = self.H_lst[i]
                Psi = dot(self.A_lst[i], H)
                Xi = H.dot(H.T).dot(H)
                for j in range(d):
                    if j == i: continue
                    small, large = map(lambda f: f([i, j]), [min, max])
                    Lambda = self.lambda_dct.get((small, large))
                    S = self.S_dct.get((small, large)), zeros_like(Xi))
                    # increment Xi
                    if j < i:
                        F1 = H
                    elif j > i:
                        F1 = S.T.dot(S).dot(H)
                        S = S.T
                    Xi += Lambda * F1.dot(F1.T).dot(H)
                    # increment Psi
                    F2 = S.dot(self.H_lst[j])
                    Psi += Lambda * F2.dot(F2.T).dot(H)
                self.H_lst[i] = H * (Psi / Xi)**(1/4)
                err += mean(abs(H - self.H_lst[i]))
            if err < 10**(-6):
                return True
        return False

    @classmethod
    def _normalize(cls, A):
    # normalize affinity matrix by Frobenius norm
        return A / sum(A**2)
