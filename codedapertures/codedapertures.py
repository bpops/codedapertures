#
#           ____                __             __
#          /\  _`\             /\ \           /\ \
#          \ \ \/\_\    ___    \_\ \     __   \_\ \
#           \ \ \/_/_  / __`\  /'_` \  /'__`\ /'_` \
#      ______\ \ \L\ \/\ \L\ \/\ \L\ \__  __//\ \L\ \
#     /\  _  \\ \____/\ \____/\ \___,/\ \____\ \___,_\
#     \ \ \L\ \\/_____ \/_____ \/____\/\/,_\_/______ /_ __    __    ____
#      \ \  __ \/\ '__`\  /'__`\/\`'__\ \ \/ /\ \/\ \/\`'__\/'__`\ /',__\
#       \ \ \/\ \ \ \L\ \/\  __/\ \ \/ \ \ \_\ \ \_\ \ \ \//\  __//\__, `\
#        \ \_\ \_\ \ ,__/\ \____\\ \_\  \ \__\\ \____/\ \_\\ \____\/\____/
#         \/_/\/_/\ \ \/  \/____/ \/_/   \/__/ \/___/  \/_/ \/____/\/___/
#                  \ \_\
#                   \/_/
#
#              a python package for generating coded apertures                                                                        
#
#                               MIT license
#                      https://github.com/bpops/cappy
#

from   commpy             import pnsequence
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
        self.L = get_prime(rank)
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
    rank : int
        determines the order, v, a prime of the form v=4n-1 
        (default 6)
    r : int
        feeds into pattern (default 5)
    quiet : bool
        if True, will not print information about the array upon creation
    """

    def __init__(self, rank=4, r=5, radius=5, quiet=False):
        self.rank = rank
        self.v    = self.get_order(self.rank)
        self.n    = int((self.v+1)/4)
        self.r    = r

        # calculate mask size
        self.radius       = radius
        self.diam         = self.radius*2+1
        self.side_width   = radius + 1

        # create axial matrix
        self.axial_matrix = np.zeros((self.diam,self.diam))
        self.loc_matrix   = np.zeros((2,self.diam,self.diam))

        # calculate intermediates
        self.v   = 4*self.n-1
        self.k   = 2*self.n-1
        self.lam = self.n-1

        # construct cyclic difference set D
        self.D = np.zeros((int((self.v-1)/2)), dtype=np.int32)
        for i in range(len(self.D)):
            self.D[i] = ((i+1)**2) % self.v

        # determine labels
        self.rx = self.diam
        self.ry = self.diam
        self.l = np.zeros((self.rx,self.ry), dtype=np.int16)
        for i in range(self.rx):
            for j in range(self.ry):
                self.l[i,j] = (i + r*j) % self.v

        # calculate mask
        self.mask = np.zeros(self.l.shape, dtype=np.int16)
        for i in range(self.mask.shape[0]):
            for j in range(self.mask.shape[1]):
                if self.l[i,j] in self.D:
                    self.mask[i,j] = 1

        # map to axial matrix
        for i in range(self.diam):
            for j in range(self.diam):
                if (i+j > (self.radius-1)) and (i+j < (self.diam+self.radius)):
                    self.axial_matrix[i,j] = self.mask[i,j]
                else:
                    self.axial_matrix[i,j] = np.nan

        if not quiet: self.report()

    def get_order(self,rank):
        """
        Determine order from the given rank, n. Order is defined
        as a prime satisfying the condition v=4n-1

        Parameters
        ----------
        rank : int
            rank; atleast 1

        Returns
        -------
        v : int
            prime order
        """
        
        n = 1
        this_rank = -1
        while True:
            v = 4*n - 1
            if pyprimes.isprime(v):
                this_rank += 1
            if this_rank == rank:
                break
            n += 1
        return v

    def report(self):
        """
        Report the array info
        """
        print("Skew-Hadamard Uniformly Redundant Array")
        print(f"rank:       {self.rank}")
        print(f"n:          {self.n}")
        print(f"order (v):  {self.v}")
        print(f"k:          {self.k}")
        print(f"lambda:     {self.lam}")
        print(f"r:          {self.r}")
        print(f"side width: {self.side_width}")
  
    def show_rhombus(self, labels=False, labelsize=10):
        """
        Plot the mask rhombus

        Parameters
        ----------
        labels : bool
            if True, will show the labels on top of each pixel
        labelsize: int
            fontsize for labels
        """

        # determine hex vertex-to-vertex and radius
        hex_vert = 1/(np.sqrt(3)/2)
        hex_radius = (hex_vert)/2.0

        # set up plotting
        fig, ax = plt.subplots(1)
        ax.set_aspect('equal')

        for x_i in range(self.mask.shape[0]):
            for y_i in range(self.mask.shape[1]):

                # determine patch origin
                x = x_i + y_i * 0.5
                y = y_i/hex_vert  

                # recenter
                x -= (self.mask.shape[0] + self.mask.shape[1]/2.0)/2.0 -1/hex_vert
                y -= (self.mask.shape[1] * 1/hex_vert)/2.0 - hex_vert/2.0

                # add hexagon
                if self.mask[x_i,y_i] == 1:
                    facecolor='k'
                else:
                    facecolor='w'
                hex = RegularPolygon((x, y), numVertices=6, radius=hex_radius, 
                              orientation=np.radians(60), 
                               facecolor=facecolor, alpha=0.6, edgecolor='k')
                ax.add_patch(hex)

                # add label
                if labels:
                    plt.annotate(self.l[x_i,y_i], (x,y),
                                 ha='center',va='center', fontsize=labelsize,
                                 transform=ax.transAxes)

        plt.xlim(-self.rx/1.2,self.rx/1.2)
        plt.ylim(-self.ry/2.0,self.ry/2.0)
        plt.title(f"SHURA rhombus [o:{self.v}, r:{self.r}]")
        plt.show()

    def show(self, labels=False, labelsize=10):
        """
        Plot the mask

        Parameters
        ----------
        labels : bool
            if True, will show the labels on top of each pixel
        labelsize: int
            fontsize for labels
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
                alpha     = 0.9 if self.axial_matrix[x+start_i,y] == 1 else 0.3
                label     = self.l[x+start_i,y]
                hex = RegularPolygon((x+0.5*abs(y-self.radius)-self.radius,
                                      ((y-self.radius)*((3/2)*hex_vert/2.0))),
                                     numVertices=6, radius=hex_vert/2.0, 
                                     orientation=np.radians(60), 
                                     facecolor=facecolor, alpha=alpha,
                                     edgecolor='k')
                ax.add_patch(hex)
                if labels:
                    plt.annotate(label, (x+0.5*abs(y-self.radius)-self.radius,
                                         (y-self.radius)*((3/2)*hex_vert/2.0)),
                                 ha='center',va='center', fontsize=labelsize,
                                 transform=ax.transAxes)

        # set axis limits
        plt.xlim(-self.radius*hex_vert,self.radius*hex_vert)
        plt.ylim(-self.radius,self.radius)
        plt.title(f"SHURA [o:{self.v}, r:{self.r}]")
        plt.show()

class hura():
    """
    Hexagonal Uniformly Redundant Array

    Parameters
    ----------
    rank : int
        determines the order, v, a prime of the form 3 or 12n+7
    quiet : bool
        if True, will not print information about the array upon creation
    """

    def __init__(self, rank=4, radius=5, quiet=False):
        self.rank = rank
        self.v    = self.get_order_from_rank(self.rank)
        self.n    = int((self.v-7)/12)
        self.r    = self.get_valid_r()

        # calculate mask size
        self.radius       = radius
        self.diam         = self.radius*2+1
        self.side_width   = radius + 1

        # create axial matrix
        self.axial_matrix = np.zeros((self.diam,self.diam))
        self.loc_matrix   = np.zeros((2,self.diam,self.diam))

        # calculate intermediates
        self.k   = 2*self.n-1
        self.lam = self.n-1

        # construct cyclic difference set D
        self.D = np.zeros((int((self.v-1)/2)), dtype=np.int32)
        for i in range(len(self.D)):
            self.D[i] = ((i+1)**2) % self.v

        # determine labels
        self.rx = self.diam
        self.ry = self.diam
        self.l = np.zeros((self.rx,self.ry), dtype=np.int16)
        for i in range(self.rx):
            for j in range(self.ry):
                self.l[i,j] = (i + self.r*j) % self.v

        # calculate mask
        self.mask = np.zeros(self.l.shape, dtype=np.int16)
        for i in range(self.mask.shape[0]):
            for j in range(self.mask.shape[1]):
                if self.l[i,j] in self.D:
                    self.mask[i,j] = 1

        # map to axial matrix
        for i in range(self.diam):
            for j in range(self.diam):
                if (i+j > (self.radius-1)) and (i+j < (self.diam+self.radius)):
                    self.axial_matrix[i,j] = self.mask[i,j]
                else:
                    self.axial_matrix[i,j] = np.nan

        if not quiet: self.report()

    def get_order_from_rank(self, rank):
        """
        Get the Order, v, from specified rank

        Parameters
        ----------
        rank : int
            rank, or nth order that satisfies 12n+7
        """
        if rank == 0:
            return 3
        else:
            n = 0
            this_rank = -1
            while True:
                v = 12*n +7
                if pyprimes.isprime(v):
                    this_rank += 1
                if this_rank == rank:
                    break
                n += 1
        return v
    
    def get_valid_r(self):
        """
        Determines the valid r from v

        Returns
        r : int
            the value that satifies r**2 % v == (r-1) % v
        """

        r_not_found = True;
        r = 0
        while r_not_found:
            if (r**2 % self.v) == ((r-1) % self.v):
                r_not_found = False
                break
            r += 1
        return r

    def report(self):
        """
        Report the array info
        """
        print("Hexagonal Uniformly Redundant Array")
        print(f"rank:       {self.rank}")
        print(f"n:          {self.n}")
        print(f"order (v):  {self.v}")
        print(f"k:          {self.k}")
        print(f"lambda:     {self.lam}")
        print(f"r:          {self.r}")
        print(f"side width: {self.side_width}")
  
    def show_rhombus(self, labels=False, labelsize=10):
        """
        Plot the mask rhombus

        Parameters
        ----------
        labels : bool
            if True, will show the labels on top of each pixel
        labelsize: int
            fontsize for labels
        """

        # determine hex vertex-to-vertex and radius
        hex_vert = 1/(np.sqrt(3)/2)
        hex_radius = (hex_vert)/2.0

        # set up plotting
        fig, ax = plt.subplots(1)
        ax.set_aspect('equal')

        for x_i in range(self.mask.shape[0]):
            for y_i in range(self.mask.shape[1]):

                # determine patch origin
                x = x_i + y_i * 0.5
                y = y_i/hex_vert  

                # recenter
                x -= (self.mask.shape[0] + self.mask.shape[1]/2.0)/2.0 -1/hex_vert
                y -= (self.mask.shape[1] * 1/hex_vert)/2.0 - hex_vert/2.0

                # add hexagon
                if self.mask[x_i,y_i] == 1:
                    facecolor='k'
                else:
                    facecolor='w'
                hex = RegularPolygon((x, y), numVertices=6, radius=hex_radius, 
                              orientation=np.radians(60), 
                               facecolor=facecolor, alpha=0.6, edgecolor='k')
                ax.add_patch(hex)

                # add label
                if labels:
                    plt.annotate(self.l[x_i,y_i], (x,y),
                                 ha='center',va='center', fontsize=labelsize,
                                 transform=ax.transAxes)

        plt.xlim(-self.rx/1.2,self.rx/1.2)
        plt.ylim(-self.ry/2.0,self.ry/2.0)
        plt.title(f"SHURA rhombus [o:{self.v}, r:{self.r}]")
        plt.show()

    def show(self, labels=False, labelsize=10):
        """
        Plot the mask

        Parameters
        ----------
        labels : bool
            if True, will show the labels on top of each pixel
        labelsize: int
            fontsize for labels
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
                alpha     = 0.9 if self.axial_matrix[x+start_i,y] == 1 else 0.3
                label     = self.l[x+start_i,y]
                hex = RegularPolygon((x+0.5*abs(y-self.radius)-self.radius,
                                      ((y-self.radius)*((3/2)*hex_vert/2.0))),
                                     numVertices=6, radius=hex_vert/2.0, 
                                     orientation=np.radians(60), 
                                     facecolor=facecolor, alpha=alpha,
                                     edgecolor='k')
                ax.add_patch(hex)
                if labels:
                    plt.annotate(label, (x+0.5*abs(y-self.radius)-self.radius,
                                         (y-self.radius)*((3/2)*hex_vert/2.0)),
                                 ha='center',va='center', fontsize=labelsize,
                                 transform=ax.transAxes)

        # set axis limits
        plt.xlim(-self.radius*hex_vert,self.radius*hex_vert)
        plt.ylim(-self.radius,self.radius)
        plt.title(f"HURA [o:{self.v}, r:{self.r}]")
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
        print(f"radius:       {self.radius}")
        print(f"diameter:     {self.diam}")
        print(f"side width:   {self.side_width}")
        print(f"desired fill: {self.fill:.2f}")
        print(f"actual  fill: {self.actual_fill:.2f}")

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
                                     facecolor=facecolor, alpha=alpha,
                                     edgecolor='k')
                ax.add_patch(hex)

        # set axis limits
        plt.xlim(-self.radius*hex_vert,self.radius*hex_vert)
        plt.ylim(-self.radius,self.radius)
        plt.title(f"Random Hex Array [diam: {self.diam}, fill: {self.actual_fill:.2f}]")
        plt.show()

class pnp():
    """
    Pseudo-Noise Product Array

    Paramters
    ---------
    m : int
        degree of a_i
    n : int
        degree of b_i
    mult : int
        how many times to tile the array (default: 2)
    """

    def __init__(self, m, n, mult=2):

        # generate primitive polynomials
        self.m = m
        self.n = n
        a = prim_poly(m)
        b = prim_poly(n)
        self.r = len(a)
        self.s = len(b)

        # generate mask
        self.mask = np.zeros((self.r,self.s))
        for i in range(self.r):
            for j in range(self.s):
                self.mask[i,j] = a[i]*b[j]
        self.mask = np.roll(self.mask,(-(m-2),-(n-2)), axis=(0,1))

        # tile array
        self.mask = np.tile(self.mask, (mult, mult))
        self.mask = np.roll(self.mask, (-int(self.mask.shape[0]/4.0)-1,
                                        -int(self.mask.shape[1]/4.0)-1),
                                        axis=(0,1))

        self.report()

    def report(self):
        """
        Report the array info
        """
        print("Pseudo-Noise Product Array")
        print(f"m: {self.m}")
        print(f"n: {self.n}")
        print(f"r (m width): {self.r}")
        print(f"s (n width): {self.s}")

    def show(self, inverse=False, size=5):
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
        plt.imshow(self.mask, cmap="gray", origin='lower', aspect=1)
        plt.axis('off')
        plt.title(f"Pseudo-Noise Product Array [m: {self.m}, n: {self.n}]")
        plt.show()
        

def prim_poly(m):
    """
    Primitive Polynomial

    Parameters
    ----------
    m : int
        degree (between 1 and 40); large numbers will take a very long time

    Returns
    -------
    pnsequence : ndarray
        a pseudo-random sequence satisfying the primitive polynomial of the
        degree specified
    """

    length = 2**m-1

    # define the first 40 primitive polynomial indices here
    h_x = {1:(0), 2:(1,0), 3:(1,0), 4:(1,0), 5:(2,0), 6:(1,0), 7:(1,0),
           8:(6,5,1,0), 9:(4,0), 10:(3,0), 11:(2,0), 12:(7,4,3,0),
           13:(4,3,1,0), 14:(12,11,1,0), 15:(1,0), 16:(5,3,2,0), 17:(3,0),
           18:(7,0), 19:(6,5,1,0), 20:(3,0), 21:(2,0), 22:(1,0), 23:(5,0),
           24:(4,3,1,0), 25:(3,0), 26:(8,7,1,0), 27:(8,7,1,0), 28:(3,0),
           29:(2,0), 30:(16,15,1,0), 31:(3,0), 32:(28,27,1,0), 33:(13,0), 
           34:(15,14,1,0), 35:(2,0), 36:(11,0), 37:(12,10,2,0), 38:(6,5,1,0),
           39:(4,0), 40:(21,19,2,0)}
    min_m = np.min(list(int(key) for key in h_x.keys()))
    max_m = np.max(list(int(key) for key in h_x.keys()))

    # check the degree exists
    if (m < min_m) or (m > max_m):
        raise ValueError(f"degree must be betweeen {min_m} and {max_m}")

    # generate mask for this degree
    mask = np.zeros(m)
    for i in h_x[m]:
        mask[m-i-1] = 1

    # initialize seed to match results from [MacWilliams 1976]
    seed = np.zeros(m)
    seed[0] = 1

    return pnsequence(m,seed,mask,2**m-1)

def get_prime(rank):
    """
    Determine prime of specified rank

    Parameters
    ----------
    rank : int
        the rank of the prime number
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