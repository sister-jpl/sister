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
import requests
import numpy as np
import hytools as ht
from hytools.io.envi import WriteENVI
from hytools.io.envi import envi_header_dict
import h5py
from skimage.util import view_as_blocks
import pyproj

def get_neon_radiance(site,date,line,out_dir):
    '''Given a site, date and line number this scripts retrieves the line via
    the NEON API data portal

    Args:
        site (str): Four letter site code.
        date (str): Date in format YYYMMDD.
        line (int): Line number.
        out_dir (str): Output directory path.

    Returns:
        filename (str): Pathname of downloaded HDF file.

    '''

    year = date[:4]
    month =  date[4:6]

    product_request = requests.get('https://data.neonscience.org/api/v0/data/DP1.30008.001/%s/%s-%s' % (site,year,month)).json()
    files= product_request['data']['files']

    filename = None

    # Cycle until matching file is found
    for file in files:
        if ('L%03d' % line in file['name']) & (date in file['name']):
            # Download image to disk
            base_name = file['name'].replace('_radiance.h5','')
            out_dir = out_dir+  base_name + '/'
            if not os.path.isdir(out_dir):
                os.mkdir(out_dir)
            url = file['url']
            filename = '%s/%s' % (out_dir,file['name'])

            if not os.path.isfile(filename):
                with requests.get(url, stream=True) as r:
                    print("Downloading %s" % file['name'])
                    r.raise_for_status()
                    with open(filename, 'wb') as f:
                        for chunk in r.iter_content(chunk_size=int(1E8)):
                            f.write(chunk)
    if not filename:
        print('%s %s %s not found' % (site,date,line))
    return filename

def h5_radiance_to_envi(filename,resolution = 1):
    '''Convert a NEON HDF radiance file to ENVI formated
    image along with observables and location data cubes

    TODO: Recalculate terrain azimuth, provided product may be
    incorrect

    Args:
        filename (str): Path to HDF file.
        resolution (int, optional): Output image resolution. Defaults to 1.

    Returns:
        None.

    '''

    # Load HDF file
    hdf_obj = h5py.File(filename,'r')

    key = [key for key in hdf_obj.keys()][0]
    rad_dec = hdf_obj[key]['Radiance']['RadianceDecimalPart']
    rad_int =hdf_obj[key]['Radiance']['RadianceIntegerPart']
    obs = hdf_obj[key]['Radiance']['Metadata']['Ancillary_Rasters']['OBS_Data']
    igm = hdf_obj[key]['Radiance']['Metadata']['Ancillary_Rasters']['IGM_Data']

    wavelengths =hdf_obj[key]['Radiance']['Metadata']['Spectral_Data']['Wavelength'][:].tolist()
    fwhm =hdf_obj[key]['Radiance']['Metadata']['Spectral_Data']['FWHM'][:].tolist()

    map_info = hdf_obj[key]['Radiance']['Metadata']['Coordinate_System']['Map_Info'][()].decode("utf-8").split(',')
    epsg = hdf_obj[key]['Radiance']['Metadata']['Coordinate_System']['EPSG Code'][()].decode("utf-8")

    new_lines = rad_dec.shape[0]//resolution
    new_cols = rad_dec.shape[1]//resolution

    map_info[5] =  resolution
    map_info[6] =  resolution

    map_info = [str(info).strip() for info in map_info]

    # Export integer and decimal radiance components
    # to temporary ENVI files
    rad_dict = envi_header_dict()
    rad_dict['lines']= rad_dec.shape[0]
    rad_dict['samples']= rad_dec.shape[1]
    rad_dict['bands']=  rad_dec.shape[2]
    rad_dict['interleave']= 'bsq'
    rad_dict['data type'] = 12
    rad_dict['byte order'] = 0
    dec_temp = filename.replace('radiance.h5','rad_dec')
    writer = WriteENVI(dec_temp,rad_dict )
    writer.write_chunk(rad_dec, 0,0)

    int_temp = filename.replace('radiance.h5','rad_int')
    writer = WriteENVI(int_temp,rad_dict )
    writer.write_chunk(rad_int, 0,0)

    int_obj = ht.HyTools()
    int_obj.read_file(int_temp, 'envi')

    dec_obj = ht.HyTools()
    dec_obj.read_file(dec_temp, 'envi')

    # Export radiance
    ##################
    rad_dict = envi_header_dict()
    rad_dict['lines']= new_lines
    rad_dict['samples']= new_cols
    rad_dict['bands']= rad_dec.shape[2]
    rad_dict['wavelength']= wavelengths
    rad_dict['fwhm']= fwhm
    rad_dict['interleave']= 'bil'
    rad_dict['data type'] = 4
    rad_dict['wavelength units'] = "nanometers"
    rad_dict['byte order'] = 0
    rad_dict['data ignore value']= -9999
    rad_dict['map info'] =map_info

    output_name = filename.replace('radiance.h5','rad')
    writer = WriteENVI(output_name,rad_dict)

    for band_num in range(rad_dict['bands']):
        print(band_num)
        band_int = int_obj.get_band(band_num).astype(float)
        band_dec = dec_obj.get_band(band_num)/50000
        band = band_int + band_dec
        band[band_int==255] = np.nan
        band = band[:new_lines*resolution,:new_cols*resolution]
        band  = view_as_blocks(band, (resolution,resolution)).mean(axis=(2,3))
        band[np.isnan(band)] = -9999
        writer.write_band(band,band_num)

    os.remove(dec_temp)
    os.remove(int_temp)

    # Export observables
    ####################
    obs_dict = envi_header_dict()
    obs_dict['band_names'] = ['path length', 'to-sensor azimuth',
                                'to-sensor zenith','to-sun azimuth',
                                  'to-sun zenith','phase', 'slope',
                                  'aspect', 'cosine i','UTC time']
    obs_dict['data type'] = 4
    obs_dict['lines']= new_lines
    obs_dict['samples']= new_cols
    obs_dict['bands']= 10
    obs_dict['fwhm']= fwhm
    obs_dict['interleave']= 'bil'
    obs_dict['data type'] = 4
    obs_dict['byte order'] = 0
    obs_dict['data ignore value']= -9999
    obs_dict['map info'] =map_info

    output_name = filename.replace('radiance.h5','obs')
    writer = WriteENVI(output_name,obs_dict)

    for band_num in range(obs_dict['bands']):
        print(band_num)
        band = obs[:,:,band_num]
        band[band==-9999] = np.nan
        band = band[:new_lines*resolution,:new_cols*resolution]
        band  = view_as_blocks(band, (resolution,resolution)).mean(axis=(2,3))
        band[np.isnan(band)] = -9999
        writer.write_band(band,band_num)

    # Export location datacube (lon,lat,elevation)
    ##############################################
    loc_dict = envi_header_dict()
    loc_dict['band_names'] = ['longitude','latitude','elevation']
    loc_dict['data type'] = 4
    loc_dict['lines']= new_lines
    loc_dict['samples']= new_cols
    loc_dict['bands']= 3
    loc_dict['fwhm']= fwhm
    loc_dict['interleave']= 'bil'
    loc_dict['data type'] = 4
    loc_dict['byte order'] = 0
    loc_dict['data ignore value']= -9999
    loc_dict['map info'] =map_info

    output_name = filename.replace('radiance.h5','loc')
    writer = WriteENVI(output_name,loc_dict)

    in_proj = pyproj.Proj("+init=EPSG:%s" % epsg)
    out_proj= pyproj.Proj("+init=EPSG:4326")

    longitude,latitude = pyproj.transform(in_proj, out_proj,
                                          igm[:,:,0],
                                          igm[:,:,1])

    elevation = igm[:,:,2]
    mask = elevation==-9999

    for band_num,band in enumerate([longitude,latitude,elevation]):
        print(band_num)
        band[mask] = np.nan
        band = band[:new_lines*resolution,:new_cols*resolution]
        band  = view_as_blocks(band, (resolution,resolution)).mean(axis=(2,3))
        band[np.isnan(band)] = -9999
        writer.write_band(band,band_num)

    os.remove(filename)
