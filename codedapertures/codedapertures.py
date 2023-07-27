#
#       _____       _       _ _____             _                            
#      |     |___ _| |___ _| |  _  |___ ___ ___| |_ _ _ ___ ___ ___          
#      |   --| . | . | -_| . |     | . | -_|  _|  _| | |  _| -_|_ -|         
#      |_____|___|___|___|___|__|__|  _|___|_| |_| |___|_| |___|___|         
#                                  |_|                                       
#
#            a python package for generating coded apertures                                                                        
#
#                             MIT license
#                    https://github.com/bpops/cappy
#

#import copy
import numpy              as     np
import matplotlib.pyplot  as     plt
from   matplotlib.patches import RegularPolygon
import pyprimes
import random


class mask_sq():
    """
    Mask

    Holds one of several types of square coded aperture patterns.
    """

    def __init__(self):
        self.A_ij = None

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


class rand_1d(mask_sq):
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
        plt.title("Random 1D")
        plt.show()

class rand_2d(mask_sq):
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
        plt.title("Random 2D")
        plt.show()

class ura(mask_sq):
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
        plt.title("URA")
        plt.show()

class mura(mask_sq):
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
        plt.title("MURA")
        plt.show()

class shura():
    """
    Skew-Hadamard Uniformly Redundant Array

    Parameters
    ----------
    n : int
        determines the order, v, where v=4n-1 (default 6)
    r : int
        feeds into pattern (default 5)
    quiet : bool
        if True, will print information about the array upon creation
    """

    def __init__(self, n=6, r=5, radius=5, quiet=False):
        self.n    = n
        self.r    = r

        # calculate mask size
        self.radius       = radius
        self.diam         = self.radius*2+1
        self.side_width   = radius + 1

        # create axial matrix
        self.axial_matrix = np.zeros((self.diam,self.diam))
        self.loc_matrix   = np.zeros((2,self.diam,self.diam))

        # calculate intermediates
        self.v   = 4*n-1
        self.k   = 2*n-1
        self.lam = n-1

        # construct cyclic difference set D
        self.D = np.zeros((int((self.v-1)/2)), dtype=np.int32)
        for i in range(len(self.D)):
            self.D[i] = ((i+1)**2) % self.v

        # determine labels
        self.rx = self.diam
        self.ry = self.diam
        self.l = np.zeros((self.rx,self.ry))
        for i in range(self.rx):
            for j in range(self.ry):
                self.l[i,j] = (i + r*j) % self.v

        # calculate mask
        self.mask = np.zeros(self.l.shape)
        for i in range(self.mask.shape[0]):
            for j in range(self.mask.shape[1]):
                if self.l[i,j] in self.D:
                    self.mask[i,j] = 1

        # map to axial matrix
        for i in range(self.diam):
            for j in range(self.diam):
                if (i+j > (self.radius-1)) and (i+j < (self.diam+self.radius)):
                    # this next line should not work correctly. WHY does it work?
                    self.axial_matrix[i,j] = self.mask[i,j]
                else:
                    self.axial_matrix[i,j] = np.nan

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
        print(f"side: {self.side_width}")
  
    def show_rhombus(self):
        """
        Plot the mask rhombus
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
            y *= 1/hex_vert

            # recenter
            x -= (self.mask.shape[0] + self.mask.shape[1]/2.0)/2.0 -1/hex_vert
            y -= (self.mask.shape[1] * 1/hex_vert)/2.0

            # add hexagon
            hex = RegularPolygon((x, y), numVertices=6, radius=hex_radius, 
                                    orientation=np.radians(60), 
                                    facecolor='k', alpha=0.6, edgecolor='k')
            ax.add_patch(hex)

        # open pixels
        for x, y in zip (x_opened, y_opened):

            # determine patch origin
            x += y * 0.5
            y *= np.sqrt(3)/2

            # recenter
            x -= (self.mask.shape[0] + self.mask.shape[1]/2.0)/2.0 -1/hex_vert
            y -= (self.mask.shape[1] * 1/hex_vert)/2.0

            # add hexagon
            hex = RegularPolygon((x, y), numVertices=6, radius=hex_radius, 
                                    orientation=np.radians(60), 
                                    facecolor='w', alpha=0.2, edgecolor='k')
            ax.add_patch(hex)

        plt.xlim(-self.rx,self.rx)
        plt.ylim(-self.ry,self.ry)
        #plt.xlim(-10,10)
        #plt.ylim(-10,10)
        plt.title(f"SHURA rhombus [o:{self.v}, r:{self.r}]")
        plt.show()

    def show(self):
        """
        Plot the mask
        """

        # set up plot parameters
        fig, ax = plt.subplots(1)
        ax.set_aspect('equal')
        hex_width = 1.0 # face-to-face distance
        hex_vert  = (hex_width)*(2.0/np.sqrt(3))

        # draw hexagon array
        for y in range(self.diam):
            row_width = self.diam - abs(self.radius-y)
            start_i   = np.max((self.radius-y,0))
            for x in range(row_width):
                facecolor = 'k' if self.axial_matrix[x+start_i,y] == 1 else 'w'
                alpha     = 0.6 if self.axial_matrix[x+start_i,y] == 1 else 0.3
                hex = RegularPolygon((x+0.5*abs(y-self.radius)-self.radius,
                                      ((y-self.radius)*((3/2)*hex_vert/2.0))),
                                     numVertices=6, radius=hex_vert/2.0, 
                                     orientation=np.radians(60), 
                                     facecolor=facecolor, alpha=alpha, edgecolor='k')
                ax.add_patch(hex)

        # set axis limits
        plt.xlim(-self.radius*hex_vert,self.radius*hex_vert)
        plt.ylim(-self.radius,self.radius)
        plt.title(f"SHURA [o:{self.v}, r:{self.r}]")
        plt.show()

class rand_hex():
    """
    Random Hexagonal Array

    Parameters
    ----------
    radius: int
        vertex-to-vertex radius of the array, minus half the pixel width
    fill: float
        fraction fill
    quiet : bool
        if True, will print information about the array upon creation
    """
    def __init__(self, radius=3, fill=0.5, quiet=False):

        # get/determine mask properties
        self.radius       = radius
        self.diam         = self.radius*2+1
        self.side_width   = radius + 1
        self.axial_matrix = np.zeros((self.diam,self.diam))
        self.loc_matrix   = np.zeros((2,self.diam,self.diam))
        self.fill         = fill

        # generate mask pattern
        for i in range(self.diam):
            for j in range(self.diam):
                if (i+j > (self.radius-1)) and (i+j < (self.diam+self.radius)):
                    if random.random() < self.fill:
                        self.axial_matrix[i,j] = 1
                else:
                    self.axial_matrix[i,j] = np.nan

        # determine actual fill factor
        self.actual_fill = np.sum(self.axial_matrix == 1) / np.sum(~np.isnan(self.axial_matrix))

        # generate locations
        for i in range(self.diam):
            for j in range(self.diam):
                if not np.isnan(self.axial_matrix[i,j]):
                    self.loc_matrix[0,i,j] = i*np.sqrt(3) - abs(j-self.radius)/2.0
                    self.loc_matrix[1,i,j] = -i + j

        if not quiet: self.report()

    def report(self):
        """
        Report the array info
        """
        print("Random Hexagonal Array")
        print(f"radius: {self.radius}")
        print(f"diameter: {self.diam}")
        print(f"side_width: {self.side_width}")
        print(f"desired fill factor: {self.fill:.2f}")
        print(f"actual  fill factor: {self.actual_fill:.2f}")

    def show(self):
        """
        Plot the mask
        """

        # set up plot parameters
        fig, ax = plt.subplots(1)
        ax.set_aspect('equal')
        hex_width = 1.0 # face-to-face distance
        hex_vert  = (hex_width)*(2.0/np.sqrt(3))

        # draw hexagon array
        for y in range(self.diam):
            row_width = self.diam - abs(self.radius-y)
            start_i   = np.max((self.radius-y,0))
            for x in range(row_width):
                facecolor = 'k' if self.axial_matrix[x+start_i,y] == 1 else 'w'
                alpha     = 0.6 if self.axial_matrix[x+start_i,y] == 1 else 0.3
                hex = RegularPolygon((x+0.5*abs(y-self.radius)-self.radius,
                                      ((y-self.radius)*((3/2)*hex_vert/2.0))),
                                     numVertices=6, radius=hex_vert/2.0, 
                                     orientation=np.radians(60), 
                                     facecolor=facecolor, alpha=alpha, edgecolor='k')
                ax.add_patch(hex)

        # set axis limits
        plt.xlim(-self.radius*hex_vert,self.radius*hex_vert)
        plt.ylim(-self.radius,self.radius)
        plt.title(f"Random Hex Array [diam: {self.diam}, fill: {self.actual_fill:.2f}]")
        plt.show()
