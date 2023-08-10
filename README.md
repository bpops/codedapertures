# CodedApertures
#### a python package for generating coded apertures  

[![Downloads](https://static.pepy.tech/personalized-badge/codedapertures?period=total&units=international_system&left_color=black&right_color=blue&left_text=Downloads)](https://pepy.tech/project/codedapertures)

CodedApertures is a python module that allows one to easily generate and display common coded aperture patterns. Coded apertures are a spatial encoding technique for straight-line optics, wherein traditional lensing (e.g., visible light) is not possible. Even wherein tradiational lensing is possible, there may be other advantages (infinite depth of field). Coded apertures may therefore be used for hard x-ray and gamma-ray imaging for astrophysics, medical imaging, and homeland security applications.

New to coded apertures? Here's a nice article: https://www.paulcarlisle.net/codedaperture/.

### Usage

Install from [PyPI](https://pypi.org/project/codedapertures/) with PIP:
```
pip install codedapertures
```

Coded patterns that this package can generate:

| class | coded aperture |
|---|---|
| rand1d | random array (1D) |
| rand2d | random array (2D) |
| ura | uniformly redundant array |
| mura | modified uniformly redundant array |
| pnp | pseudo-noise product array |
| randhex | random array (hexagonal) |
| shura | skew-hadamard uniformly redundant array |
| hura | hexagonal uniformly redundant array |

Note that for consistency with relevant decoding algorithms, we define the
binary meaning as such
| val | mask pixel |
|---|---|
| 0 | closed |
| 1 | open |

See [demo.ipynb](https://github.com/bpops/codedapertures/blob/master/demo.ipynb) for examples of use.

### References

URA pattern: E. E. Fenimore and T. M. Cannon, "Coded aperture imaging with uniformly redundant arrays," Appl. Opt. 17, 337-347 (1978).

MURA pattern:  E.E. Fenimore and S. R. Gottesman, "New family of binary arrays for coded aperture imaging" Appl. Opt. 28 (20): 4344-4352 (1989).

SHURA and HURA pattern: M.H. Finger and T.A. Prince, "Hexagonal Uniformly Redundant Arrays for Coded-Aperture Imaging," Proc. 19th Int. Cosmic Ray Conf., 3: 295-298 (1985).

PNP pattern: "PNP - A new class of coded aperture arrays," S. Gottesman and E. Schneid, IEEE Trans. Nucl. Sci., 33(1): 745-749 (1986).

Pseudo-Random Sequences and Primitive Polynomials: F.J. MacWilliams and N.J.A. Sloane, "Pseudo-Random Sequences and Arrays", Proc. of the IEEE, 64, 1715 (1976).

This package uses an axial coordinate system for the hexagonal grids. It follows the concept outlined excellently at https://www.redblobgames.com/grids/hexagons/#map-storage.