import setuptools

setuptools.setup(name='codedapertures',
      version='0.2',
      description='Coded Aperture Production in PYthon (CAPPY)',
      long_description='CAPPY is a python module that allows one to easily generate and display common coded aperture patterns. Coded apertures are a spatial encoding technique for straight-line optics, wherein traditional lensing (e.g., visible light) is not possible. Even wherein tradiational lensing is possible, there may be other advantages (infinite depth of field). Coded apertures may therefore be used for hard x-ray and gamma-ray imaging for astrophysics, medical imaging, and homeland security applications.',
      url='https://github.com/bpops/cappy',
      author='bpops',
      license='MIT',
      install_requires=['pyprimes',
                'numpy',
                'matplotlib'],
      zip_safe=False)
