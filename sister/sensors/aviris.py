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
import shutil
import tarfile
import hytools as ht
from hytools.io.envi import WriteENVI
import numpy as np
from ..utils.geometry import resample



def create_loc_ort(loc_file,glt_file):
    loc = ht.HyTools()
    loc.read_file(loc_file,'envi')

    glt = ht.HyTools()
    glt.read_file(glt_file,'envi')

    samples = np.abs(glt.get_band(0))
    lines = np.abs(glt.get_band(1))

    lon = np.copy(loc.get_band(0))
    lat= np.copy(loc.get_band(1))
    elv = np.copy(loc.get_band(2))

    lon_proj = lon[lines.flatten()-1,samples.flatten()-1]
    lon_proj = lon_proj.reshape(lines.shape)

    lat_proj = lat[lines.flatten()-1,samples.flatten()-1]
    lat_proj = lat_proj.reshape(lines.shape)

    elv_proj = elv[lines.flatten()-1,samples.flatten()-1]
    elv_proj = elv_proj.reshape(lines.shape)

    lon_proj[samples ==0] = -9999
    lat_proj[samples ==0] = -9999
    elv_proj[samples ==0] = -9999

    # Create output file
    header_dict =glt.get_header()
    header_dict['bands'] = 3
    header_dict['data type'] = 4
    header_dict['band names'] = loc.get_header()['band names']
    header_dict['description'] = loc.get_header()['description']
    writer = WriteENVI(loc_file +'_ort',header_dict)
    writer.write_band(lon_proj,0)
    writer.write_band(lat_proj,1)
    writer.write_band(elv_proj,2)
    writer.close()

def preprocess(input_tar,out_dir,temp_dir,res = 0):
    '''
    input_tar = '/data2/avcl/raw/f080709t01p00r15.tar.gz'
    input_tar = '/data2/avng/rdn/ang20170901t195659.tar.gz'
    out_dir ='/data1/temp/ang_pre/output/'
    temp_dir ='/data1/temp/ang_pre/temp/'

    '''

    base_name = os.path.basename(input_tar)

    if base_name.startswith('ang'):
        base_name = base_name[:18]
    elif base_name.startswith('f'):
        base_name = base_name[:16]
    else:
        raise ValueError('Unrecognized sensor')

    out_dir = "%s/%s/" % (out_dir,base_name)

    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)

    file = tarfile.open(input_tar)
    tar_contents = [temp_dir+c for c in file.getnames()]
    file.extractall(temp_dir)
    file.close()

    #AVIRIS NG
    if base_name.startswith('ang'):
        obs_ort_file = [x for x in tar_contents if x.endswith('obs_ort')][0]
        rdn_file = [x for x in tar_contents if x.endswith('img')][0]
        loc_file = [x for x in tar_contents if x.endswith('loc')][0]
        glt_file = [x for x in tar_contents if x.endswith('glt')][0]
        create_loc_ort(loc_file,glt_file)

    #AVIRIS Classic
    elif base_name.startswith('f'):
        obs_ort_file =  [x for x in tar_contents if x.endswith('obs_ort')][0]
        rdn_file = [x for x in tar_contents if x.endswith('ort_img')][0]
        igm_file = [x for x in tar_contents if x.endswith('igm')][0]
        glt_file = [x for x in tar_contents if x.endswith('ort_glt')][0]

        #Rename files for consistency
        loc_file = igm_file.replace('igm','loc')
        os.rename(igm_file,loc_file)
        os.rename(igm_file+'.hdr',loc_file+'.hdr')

        #Scale radiance
        os.rename(rdn_file,rdn_file+'_unscale')
        os.rename(rdn_file + '.hdr',rdn_file+'_unscale.hdr')

        rdn = ht.HyTools()
        rdn.read_file(rdn_file+'_unscale','envi')

        prefix = rdn.base_name[:3]
        if prefix in ['f95', 'f96', 'f97', 'f98', 'f99', 'f00',
                      'f01', 'f02', 'f03', 'f04', 'f05']:
            gains = np.r_[50.0 * np.ones(160), 100.0 * np.ones(64)]
        elif prefix in ['f06', 'f07', 'f08', 'f09', 'f10', 'f11',
                        'f12', 'f13', 'f14', 'f15', 'f16', 'f17',
                        'f18', 'f19', 'f20', 'f21']:
            gains = np.r_[300.0 * np.ones(110), 600.0 * np.ones(50), 1200.0 * np.ones(64)]

        rdn_header = rdn.get_header()
        rdn_header['byte order'] = 0
        rdn_header['no data value'] = -9999
        rdn_header['data type'] = 4
        writer = WriteENVI(rdn.file_name.replace('_unscale',''),rdn_header)
        iterator =rdn.iterate(by = 'band')

        while not iterator.complete:
            band = np.copy(iterator.read_next()).astype(float)
            band /= gains[iterator.current_band]
            band[~rdn.mask['no_data']] = -9999
            writer.write_band(band,iterator.current_band)

        create_loc_ort(loc_file,glt_file)

    loc_ort_file = loc_file+'_ort'

    for file in [obs_ort_file,rdn_file,loc_ort_file]:
        if res==0:
            shutil.move(file,out_dir)
            shutil.move(file + '.hdr',out_dir)
        else:
            resample(file,out_dir,res)

    shutil.rmtree(tar_contents[0])
