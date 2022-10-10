from setuptools import setup, find_packages
from setuptools.command.install import install
import os

setup(
    name='sister-sbg',
    description= 'SISTER: Space-based Imaging Spectroscopy and Thermal PathfindER',
    version='1.1.0',
    url='https://github.com/EnSpec/sister',
    author = 'Adam Chlus',
    packages=find_packages(),
    install_requires=['earthengine-api',
                      'gdal',
                      'h5py',
                      'hy-tools',
                      'numba',
                      'numpy',
                      'ray',
                      'rtree',
                      'pandas',
                      'pyproj',
                      'pysolar',
                      'scikit-image',
                      'scipy',
                      'statsmodels'],
    python_requires='>=3.6, !=3.9.*',
    include_package_data=True,
    package_data={'': ['data/*']})
