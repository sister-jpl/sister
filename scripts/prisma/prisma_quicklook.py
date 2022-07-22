#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SISTER
Space-based Imaging Spectroscopy and Thermal PathfindER
Author: Adam Chlus

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, version 3 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

import argparse
from collections import OrderedDict
import glob
import numpy as np
import os
import shutil
import zipfile
import h5py
from PIL import Image
import PIL
import matplotlib.pyplot as plt

home = os.path.expanduser('~')
zip_dir = '/data2/prisma/raw/'
out_dir = '%s/Dropbox/rs/sister/figures/prisma_quicklooks/' % home
temp_dir = '/data2/temp/'

for l1_zip in glob.glob('%s*L1*.zip' % zip_dir):

    base_name = os.path.basename(l1_zip).replace('.zip','')
    image_file ="%s/%s.jpeg" % (out_dir,base_name)
    if os.path.isfile(image_file):
        continue

    with zipfile.ZipFile(l1_zip,'r') as zipped:
        zipped.extractall(temp_dir)
    l1  = '%s/%s.he5' % (temp_dir,base_name)

    l1_obj = h5py.File(l1,'r')
    vnir_data =  l1_obj['HDFEOS']["SWATHS"]['PRS_L1_HCO']['Data Fields']['VNIR_Cube']
    vnir_waves = l1_obj.attrs.get('List_Cw_Vnir')
    vnir_fwhm = l1_obj.attrs.get('List_Fwhm_Vnir')

    band3 =   vnir_data[:,np.argmin(np.abs(vnir_waves-550)),:]
    band2 = vnir_data[:,np.argmin(np.abs(vnir_waves-660)),:]
    band1 = vnir_data[:,np.argmin(np.abs(vnir_waves-850)),:]
    rgb=  np.stack([band1,band2,band3])

    rgb = np.moveaxis(rgb,0,-1).astype(float)
    bottom = np.nanpercentile(rgb,5,axis = (0,1))
    top = np.nanpercentile(rgb,95,axis = (0,1))
    rgb = np.clip(rgb,bottom,top)
    rgb = (rgb-np.nanmin(rgb,axis=(0,1)))/(np.nanmax(rgb,axis= (0,1))-np.nanmin(rgb,axis= (0,1)))
    rgb =(rgb*255).astype(np.uint8)

    print(base_name)
    print(l1_obj.attrs['Processor_Version'])
    plt.imshow(rgb)
    plt.show()
    plt.close()

    im = Image.fromarray(rgb)
    im.save(image_file)

    os.remove(l1)
