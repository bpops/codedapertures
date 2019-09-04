import numpy as np
import matplotlib.pyplot as plt
import pyprimes

class mask():
    
    def __init__(self):
        A_ij = None
    
    def show(self, inverse=False):
        if inverse:
            plt.imshow(self.A_ij, cmap="binary")
        else:
            plt.imshow(self.A_ij, cmap="binary_r")
        plt.axis('off')
        plt.show()

    def get_pattern(self, inverse=False):
        if inverse:
            return 1-self.A_ij
        else:
            return self.A_ij
        
class URA(mask):
    
    def __init__(self, rank=4, mult=2):
        self.rank = rank
        self.mult = mult
        
        # get r, s
        r, s = self._get_prime_pairs(self.rank)
        self.r = r
        self.s = s
        
        # generate C_r(I)
        C_r_I = np.zeros(r) - 1
        C_s_J = np.zeros(s) - 1
        for x in range(1, r):
            C_r_I[x**2 % r] = 1
        for y in range(1, s):
            C_s_J[y**2 % s] = 1

        # generate A_IJ
        A_IJ = np.zeros([r,s])
        for I in range(r):
            for J in range(s):
                if I == 0:
                    A_IJ[I,J] = 0
                elif J == 0:
                    A_IJ[I,J] = 1
                elif C_r_I[I] * C_s_J[J] == 1:
                    A_IJ[I,J] = 1

        # generate A_ij
        m = self.mult
        A_ij = np.zeros([m*r,m*s])
        for i in range(m*r):
            for j in range(m*s):
                A_ij[i,j] = A_IJ[i%r,j%s]
        A_ij = np.roll(A_ij, int((r+1)/2), axis=0)
        A_ij = np.roll(A_ij, int((s+1)/2), axis=1)
        self.A_ij = A_ij
        
        self.info()
        
    def info(self):
        print("Uniformly Redundant Array")
        print("r, s: %i, %i" % (self.r, self.s))
        print("multiplier: %i" % self.mult)
        
    def _get_prime_pairs(self, rank):
        pit = pyprimes.primes()

        # intialize
        p1 = next(pit)
        this_rank = -1

        # loop through
        while True:
            p2 = next(pit)
            if (p2-p1) == 2:
                this_rank += 1
            else:
                p1 = p2
            if this_rank == rank:
                break

        return p1, p2