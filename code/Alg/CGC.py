from numpy import *


class CGC(object):
    def __init__(self, A_lst, S_lst, k_lst, lambda_lst):
        d = self.d = len(A_lst)
        assert (len(k_lst) == d
                and len(lambda_lst) == d*(d-1)/2
                and len(S_lst) == d*(d-1)/2)
        self.A_lst = map(cls._normalize, A_lst)
        self.n_lst = map(lambda arr: arr.shape[0], self.A_lst)
        self.S_lst = S_lst
        self.k_lst = k_lst
        self.lambda_lst = lambda_lst
        self.H_lst = []
        for i in range(d):
            H = -random.rand(n_lst[i], k_lst[i]) + 1
            self.H_lst.append(H)

    def iterate(self):
        pass
        
    @classmethod
    def _normalize(cls, A):
    # normalize affinity matrix by Frobenius norm
        return A / sum(A**2)
