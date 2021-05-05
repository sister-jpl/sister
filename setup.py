from setuptools import setup, find_packages

setup(
    name='sister-sbg',
    description= 'SISTER: Space-based Imaging Spectroscopy and Thermal PathfindER',
    version='0.1',
    license='MIT',
    url='https://github.com/EnSpec/hytools',
    author = 'Adam Chlus',
    packages=find_packages(),
    install_requires=['h5py',
                      'pyproj',
                      'numpy',
                      'pandas',
                      'ray',
                      'scikit-learn',
                      'scipy',
                      'hy-tools'],
    python_requires='>=3.6, !=3.9.*'
    )


