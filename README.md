# SISTER
## Space-based Imaging Spectroscopy and Thermal pathfindER

This repository contains code for implementing prototype algorithm workflows
for the Space-based Imaging Spectroscopy and Thermal pathfindER (SISTER).

This repository is under active development and currently contains
code for preprocessing imaging spectroscopy data from airborne and spaceborne
sensors for input into higher level algorithms including atmospheric, topographic
and BRDF correction algorithms.

### Installation
To install, clone
```bash
git clone https://github.com/EnSpec/sister.git
```

and install with setuptools
```bash
python setup.py install
```

### Examples

#### PRISMA HDF to ENVI

The following code takes as input a PRISMA L1 radiance image along with ESA Copernicus DEM tiles and exports
three ENVI formated files:
1. Merged VNIR+SWIR radiance datacube
2. Location datacube (longitude, latitude, altitude)
3. Observables datacube (sensor, solar geometry......)

[PRISMA Algorithm Workflow](https://raw.githubusercontent.com/EnSpec/sister/master/figures/prisma_workflow.svg)

```python

import os
from sister.sensors import prisma

base_dir = '/'

# Copernicus DEM directory
elev_dir = '/data/sister/data/cop_dsm/'
base_name = '20200621003500_20200621003505_0001'
l1_zip  = '/data/prisma/PRS_L1_STD_OFFL_%s.zip'% base_name
out_dir = '%s/data/prisma/rad/PRISMA_%s/' % (base_dir,base_name)
temp_dir =  '%s/temp/' % base_dir

#Wavelength shift surface
shift = '%s/sister/data/prisma/wavelength_shift/PRISMA_20200721104249_20200721104253_0001_wavelength_shift_surface' % base_dir

#Perform Landsat image matching (recommended)
match = True

#Project image to UTM
proj = True

#Output resolution in meters
res = 90

if not os.path.isdir(out_dir):
    os.mkdir(out_dir)

if not os.path.isdir(temp_dir):
    os.mkdir(temp_dir)

prisma.he5_to_envi(l1_zip,out_dir,temp_dir,elev_dir,
            shift = shift,match=match,proj = proj, res = res)

```
