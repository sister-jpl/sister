#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SISTER
@author: Adam Chlus

This script uncompresses Copernicus DSM tiles, extracts the DEM raster and recompresses it.
"""

import os
import tarfile
import glob

tar_dir = '/data2/cop_dem/2019_1/'
output_dir = '/data2/cop_dsm/'
temp_dir = '/data1/temp/'

for tar_file in glob.glob(tar_dir + '*.tar'):
    print(os.path.basename(tar_file))
    file = tarfile.open(tar_file)
    # Find dt2 file
    for member in file.getmembers():
        if '.dt2' in member.name:
            file_name = os.path.basename(member.name)
            extracted = file.extractfile(member)
            output_filename = "%s%s.tar.gz" % (output_dir,file_name[:-4])
            member.name = file_name
            with tarfile.open(output_filename, "w:gz") as tar:
                tar.addfile(member,extracted)
    file.close()
