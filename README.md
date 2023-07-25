# CodedApertures
#### a python package for generating coded apertures  

[![Downloads](https://static.pepy.tech/personalized-badge/codedapertures?period=total&units=international_system&left_color=black&right_color=blue&left_text=Downloads)](https://pepy.tech/project/codedapertures)

CodedApertures is a python module that allows one to easily generate and display common coded aperture patterns. Coded apertures are a spatial encoding technique for straight-line optics, wherein traditional lensing (e.g., visible light) is not possible. Even wherein tradiational lensing is possible, there may be other advantages (infinite depth of field). Coded apertures may therefore be used for hard x-ray and gamma-ray imaging for astrophysics, medical imaging, and homeland security applications.

New to coded apertures? Here's a nice article: https://www.paulcarlisle.net/codedaperture/

PyPI URL: https://pypi.org/project/codedapertures/

### Usage

Install with PIP:
```
pip install codedapertures
```

See [demo.ipynb](https://github.com/bpops/codedapertures/blob/master/demo.ipynb) for examples of use.

### Credits

URA pattern: E. E. Fenimore and T. M. Cannon, "Coded aperture imaging with uniformly redundant arrays," Appl. Opt. 17, 337-347 (1978).

MURA pattern:  E.E. Fenimore and S. R. Gottesman, "New family of binary arrays for coded aperture imaging" Appl. Opt. 28 (20): 4344-4352 (1989).

SHURA and HURA pattern: M.H. Finger and T.A. Prince, "Hexagonal Uniformly Redundant Arrays for Coded-Aperture Imaging," Proc. 19th Int. Cosmic Ray Conf., 3: 295-298 (1985).

This package will use an axial coordinate system for the hexagonal grids. It will follow the concept outlined at https://www.redblobgames.com/grids/hexagons/#map-storage.
