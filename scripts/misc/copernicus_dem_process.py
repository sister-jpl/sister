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

import os
import tarfile
import glob

tar_dir = '/data2/cop_dem/2019_1/'
output_dir = '/data2/cop_dsm/'
temp_dir = '/data1/temp/'

for tar_file in glob.glob(tar_dir + '*.tar'):
    try:
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

    except:
        print("ERROR: %s " % os.path.basename(tar_file))
