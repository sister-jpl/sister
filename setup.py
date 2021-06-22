from setuptools import setup, find_packages

setup(
    name='sister-sbg',
    description= 'SISTER: Space-based Imaging Spectroscopy and Thermal PathfindER',
    version='0.1',
    license='GPL_v3',
    url='https://github.com/EnSpec/sister',
    author = 'Adam Chlus',
    packages=find_packages(),
    install_requires=['ee',
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
    python_requires='>=3.6, !=3.9.*'
    )
