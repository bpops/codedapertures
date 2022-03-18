#                                                                                 
#            _|_|_|    _|_|    _|_|_|    _|_|_|    _|      _|                
#          _|        _|    _|  _|    _|  _|    _|    _|  _|                  
#          _|        _|_|_|_|  _|_|_|    _|_|_|        _|                    
#          _|        _|    _|  _|        _|            _|                    
#            _|_|_|  _|    _|  _|        _|            _|                    
#                                                                                 
#                Coded Aperture Production in PYthon                                                                              
#
#                           MIT license
#                  https://github.com/bpops/cappy
#

import numpy             as np
import matplotlib.pyplot as plt
import pyprimes
import random

class mask():
    """
    Mask

    Holds one of several types of coded aperture patterns.
    """

    def __init__(self):
        self.A_ij = None
    
    def show(self, inverse=False):
        """
        Plots the mask to the screen

        Parameters
        ----------
        inverse : bool
            if True, will invert the array before plotting
        """
        plt.rcParams['figure.figsize'] = [8,8]
        cmap = "binary_r" if inverse else "binary"
        plt.imshow(self.A_ij, cmap=cmap, aspect=1)
        plt.axis('off')
        plt.show()

    def get_pattern(self, inverse=False):
        """
        Returns the pattern as an array

        Parameters
        ----------
        inverse : bool
            if True, will invert the array before returning

        Returns
        -------
        A_ij : ndarrray
            the 2d boolean mask
        """
        mask = 1-self.A_ij if inverse else self.A_ij
        return mask

    def add_border(self, width, empty=False):
        """
        Adds a border with the given width

        Parameters
        ----------
        width : int
            width of border
        """
        
        new_width  = self.width  + width*2
        new_height = self.height + width*2
        new_mask   = np.zeros((new_width, new_height))
        if not empty: new_mask = new_mask +1
        new_mask[width:-width,width:-width] = self.A_ij
        self.A_ij = new_mask

            
class ura(mask):
    """
    Uniformly Redundant Array

    Parameters
    ----------
    rank : int
        the rank of prime pairs to use
    mult : int
        the number of times to tile the pattern in both dimensions
    quiet : bool
        if True, will print information about the array upon creation
    """

    
    def __init__(self, rank=4, mult=2, quiet=False):
        self.rank = rank
        self.mult = mult
        
        # get r, s
        r, s = self.__get_prime_pairs(self.rank)
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
        
        # get width/height
        self.width = self.A_ij.shape[0]
        self.height = self.A_ij.shape[1]

        if not quiet: self.report()
        
    def report(self):
        """
        Report the array info
        """
        print("Uniformly Redundant Array")
        print("r, s: %i, %i (rank %i)" % (self.r, self.s, self.rank))
        print("multiplier: %i" % self.mult)
        
    def __get_prime_pairs(self, rank):
        """
        Determine prime pairs at specified rank

        Parmeters
        ---------
        rank : int
            the rank of prime pairs to determine
        """
        pit = pyprimes.primes()

        # intialize
        p1 = next(pit)
        this_rank = -1

        # find primes
        while True:
            p2 = next(pit)
            if (p2-p1) == 2:
                this_rank += 1
            else:
                p1 = p2
            if this_rank == rank:
                break

        return p1, p2
        
class rand_array(mask):
    """
    Class to hold a randomly generate array

    Parameters
    ----------
    r : int
        number of 'x' elements in the array
    s : int
        number of 'y' elements in the array
    fill : float
        fill factor fraction
    quiet : bool
        if True, will print mask info upon creation
    """
    
    def __init__(self, r=10, s=10, fill=0.5, quiet=False):
        self.r = r
        self.s = s
        self.fill = fill
        
        # randomly fill
        A_ij = np.zeros([r, s])
        for i in range(r):
            for j in range(s):
                if random.random() < self.fill:
                    A_ij[i,j] = 1
        self.A_ij = A_ij
        self.actual_fill = np.sum(A_ij)/(self.r*self.s)
        
        # get width/height
        self.width = self.A_ij.shape[0]
        self.height = self.A_ij.shape[1]

        if not quiet: self.report()
        
    def report(self):
        """
        Report on the mask information
        """
        print("Random Array")
        print("r, s: %i, %i" % (self.r, self.s))
        print("desired fill factor: %.2f" % self.fill)
        print("actuall fill factor: %.2f" % self.actual_fill)
        
class mura(mask):
    """
    Modified Uniformly Redundant Array
    """
    
    def __init__(self, rank=5, quiet=False, mult=2):
        self.rank = rank
        self.L = self.__get_prime(rank)
        self.mult = mult
        
        # get r, s
        r = self.L
        s = self.L
        
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
        self.A_ij = A_ij

        # get width/height
        self.width = self.A_ij.shape[0]
        self.height = self.A_ij.shape[1]


        if not quiet: self.report()
        
    def report(self):
        """
        Report on the mask information
        """
        print("Modified Uniformly Redundant Array")
        print("L: %i (rank %i)" % (self.L, self.rank))
        
    def __get_prime(self, rank):
        """
        Determine prime of specified rank
        """
        m = 1
        this_rank = -1
        while True:
            L = 4*m + 1
            if pyprimes.isprime(L):
                this_rank += 1
            if this_rank == rank:
                break
            m += 1
        return L