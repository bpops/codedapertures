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

import numpy              as     np
import matplotlib.pyplot  as     plt
from   matplotlib.patches import RegularPolygon
import pyprimes
import random


class mask():
    """
    Mask

    Holds one of several types of coded aperture patterns.
    """

    def __init__(self):
        self.A_ij = None
    
    def show(self, inverse=False, size=8):
        """
        Plots the mask to the screen

        Parameters
        ----------
        inverse : bool
            if True, will invert the array before plotting
        size : int
            size of the plot (default 8)
        """
        plt.rcParams['figure.figsize'] = [size,size]
        cmap = "binary_r" if inverse else "binary"
        plt.imshow(np.transpose(self.A_ij), cmap=cmap, aspect=1)
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
            the 2D boolean mask
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
        self.A_ij  = new_mask


class rand_1d(mask):
    """
    Random 1-dimensional Array

    Parameters
    ----------
    x : int
        number of 'x' elements in the array
    y : int
        number of 'y' elements in the array
    fill : float
        fill factor fraction
    quiet : bool
        if True, will print mask info upon creation
    """
    
    def __init__(self, x=10, y=10, fill=0.5, quiet=False):
        self.r = x
        self.s = y
        self.fill = fill
        
        # randomly fill
        A_ij = np.zeros([self.r, self.s])
        for i in range(self.r):
            if random.random() < self.fill:
                A_ij[i,:] = 1
        self.A_ij = A_ij
        self.actual_fill = np.sum(A_ij[:,0])/(self.r)
        
        # get width/height
        self.width = self.A_ij.shape[0]
        self.height = self.A_ij.shape[1]

        if not quiet: self.report()
        
    def report(self):
        """
        Report on the mask information
        """
        print("Random 1D Array")
        print(f"x, y: {self.r}, {self.s}")
        print(f"desired fill factor: {self.fill:.2f}")
        print(f"actual  fill factor: {self.actual_fill:.2f}")

    def get_pattern(self, inverse=False, collapse=False):
        """
        Returns the pattern as an array

        Parameters
        ----------
        inverse : bool
            if True, will invert the array before returning
        collapse : bool
            if True, will collapse the array to one-dimension

        Returns
        -------
        A_ij : ndarrray
            the 2D (or 1D, if selected) boolean mask
        """
        mask = 1-self.A_ij if inverse else self.A_ij
        if collapse: mask = mask[:,0]
        return mask

class rand_2d(mask):
    """
    Random 2-dimensional Array

    Parameters
    ----------
    x : int
        number of 'x' elements in the array
    y : int
        number of 'y' elements in the array
    fill : float
        fill factor fraction
    quiet : bool
        if True, will print mask info upon creation
    """
    
    def __init__(self, x=10, y=10, fill=0.5, quiet=False):
        self.r = x
        self.s = y
        self.fill = fill
        
        # randomly fill
        A_ij = np.zeros([self.r, self.s])
        for i in range(self.r):
            for j in range(self.s):
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
        print("Random 2D Array")
        print(f"x, y: {self.r}, {self.s}")
        print(f"desired fill factor: {self.fill:.2f}")
        print(f"actual  fill factor: {self.actual_fill:.2f}")


class ura(mask):
    """
    Uniformly Redundant Array

    Parameters
    ----------
    rank : int
        the rank of prime pairs to use (0 -> (5,3) 1 -> (13,11) etc.)
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
        print(f"r, s: {self.r}, {self.s} (rank {self.rank})")
        print(f"multiplier: {self.mult}")
        
    def __get_prime_pairs(self, rank):
        """
        Determine prime pairs at specified rank

        Parmeters
        ---------
        rank : int
            the rank of prime pairs to determine (0 -> 5, 1 -> 13, etc.)
        """

        assert rank >= 0, f"rank must be great than or equal to zero, got {rank}"

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


class mura(mask):
    """
    Modified Uniformly Redundant Array

    Parameters
    ----------
    rank : int
        the rank of prime pairs to use (0 -> (5,3) 1 -> (13,11) etc.)
    mult : int
        the number of times to tile the pattern in both dimensions
    quiet : bool
        if True, will print information about the array upon creation
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
        print(f"L: {self.L} (rank {self.rank})")
        
    def __get_prime(self, rank):
        """
        Determine prime of specified rank

        Parameters
        ----------
        rank : int
            the rank of prime pairs (0 -> 5, 1 -> 13, etc.)
        """

        assert rank >= 0, f"rank must be great than or equal to zero, got {rank}"

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
    
class shura(mask):
    """
    Skew-Hadamard Uniformly Redundant Array

    Parameters
    ----------
    n : int
        determines the order, v, where v=4n-1 (default 6)
    r : int
        feeds into pattern (default 5)
    mult : int
        multiplier (default 10)
    quiet : bool
        if True, will print information about the array upon creation
    """

    def __init__(self, n=6, r=5, mult=4, quiet=False):
        self.n    = n
        self.r    = r
        self.mult = mult
        
        # calculate intermediates
        self.v   = 4*n-1
        self.k   = 2*n-1
        self.lam = n-1

        # construct cyclic difference set D
        self.D = np.zeros((int((self.v-1)/2)), dtype=np.int32)
        for i in range(len(self.D)):
            self.D[i] = ((i+1)**2) % self.v

        # determine labels
        rx = r*mult
        self.l = np.zeros((rx,rx))
        for i in range(rx):
            for j in range(rx):
                self.l[i,j] = (i + r*j) % self.v

        # calculate mask
        self.mask = np.zeros(self.l.shape)
        for i in range(self.mask.shape[0]):
            for j in range(self.mask.shape[1]):
                if self.l[i,j] in self.D:
                    self.mask[i,j] = 1

        if not quiet: self.report()

    def report(self):
        """
        Report the array info
        """
        print("Skew-Hadamard Uniformly Redundant Array")
        print(f"n: {self.n}")
        print(f"order (v): {self.v}")
        print(f"k: {self.k}")
        print(f"lambda: {self.lam}")
        print(f"r: {self.r}")
        print(f"multiplier: {self.mult}")

    def show(self):
        """
        Plot the mask
        """

        # determine open/closed pixels
        pix_close  = np.where(self.mask == 1)
        pix_open   = np.where(self.mask == 0)
        
        # determine open/closed pixels
        x_closed = pix_close[0]
        y_closed = pix_close[1]
        x_opened = pix_open[0]
        y_opened = pix_open[1]

        hex_vert = 1/(np.sqrt(3)/2)
        hex_radius = (hex_vert)/2.0

        fig, ax = plt.subplots(1)
        ax.set_aspect('equal')

        # closed pixels
        for x, y in zip (x_closed, y_closed):

            # determine patch origin
            x += y * 0.5
            y *= np.sqrt(3)/2

            # recenter
            x -= (self.mask.shape[0] + self.mask.shape[1]/2.0)/2.0
            y -= (self.mask.shape[1] * 1/hex_vert)/2.0

            # add hexagon
            hex = RegularPolygon((x, y), numVertices=6, radius=hex_radius, 
                                    orientation=np.radians(60), 
                                    facecolor='k', alpha=0.5, edgecolor='k')
            ax.add_patch(hex)

        # open pixels
        for x, y in zip (x_opened, y_opened):

            # determine patch origin
            x += y * 0.5
            y *= np.sqrt(3)/2

            # recenter
            x -= (self.mask.shape[0] + self.mask.shape[1]/2.0)/2.0
            y -= (self.mask.shape[1] * 1/hex_vert)/2.0

            # add hexagon
            hex = RegularPolygon((x, y), numVertices=6, radius=hex_radius, 
                                    orientation=np.radians(60), 
                                    facecolor='w', alpha=0.2, edgecolor='k')
            ax.add_patch(hex)

        plt.xlim(-self.mask.shape[0]/3.0,self.mask.shape[0]/3.0)
        plt.ylim(-self.mask.shape[0]/3.0,self.mask.shape[1]/3.0)
        plt.title(f"SHURA [o:{self.v}, r:{self.r}, x:{self.mult}]")
