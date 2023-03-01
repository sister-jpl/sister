#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SISTER
Space-based Imaging Spectroscopy and Thermal PathfindER
Author: Adam Chlus

This module uses code from the EMIT-SDS

https://github.com/emit-sds/emit-utils/blob/develop/emit_utils/reformat.py

"""

import os
import glob
import netCDF4 as nc
import hytools as ht
import numpy as np
import datetime as dt
from ..utils.geometry import utm_zone


NO_DATA_VALUE = -9999.0

def nc_to_envi(img_file,out_dir,temp_dir,obs_file=None,export_loc=False,crid='000'):
    '''
    This function exports up to three UTM projected EMIT datasets:
        *_RDN/*_RFL : Radiance or reflectance

        Optional:
        *_OBS : Observation file in the format of JPL obs files:
                1. Pathlength (m)
                2. To-sensor view azimuth angle (degrees)
                3. To-sensor view zenith angle (degrees)
                4. To-sun azimuth angle (degrees)
                5. To-sun zenith angle (degrees)
                6. Phase
                7. Slope (Degrees)
                8. Aspect (Degrees)
                9. Cosine i
                10. UTC decimal hours
        Optional:
        *LOC : Location file in the following format:
                1. Longitude (decimal degrees)
                2. Longitude (decimal degrees)
                3. Elevation (m)

    img_file(str): EMIT L1B_RDN or L2A_RFL netCDF
    out_dir(str): Output directory of ENVI datasets
    temp_dir(str): Temporary directory
    obs_file(str): EMIT L1B_OBS netCDF, will export if provided. Default: None.
    export_loc(bool): Export location datacube; Default: False

    '''
    img_nc = nc.Dataset(img_file)

    if 'radiance' in img_nc.variables.keys():
        data = img_nc['radiance']
        product = "L1B_RDN"
    elif 'reflectance' in img_nc.variables.keys():
        data = img_nc['reflectance']
        product = "L2A_RFL"
    else:
        print("Unrecognized input image dataset")
        return

    fwhm = img_nc['sensor_band_parameters']['fwhm'][:].data
    waves = img_nc['sensor_band_parameters']['wavelengths'][:].data

    # Get geospatial data (from EMIT-SDS)
    gt = np.array(img_nc.__dict__["geotransform"])
    glt = np.zeros(list(img_nc.groups['location']['glt_x'].shape) + [2], dtype=np.int32)
    glt[...,0] = np.array(img_nc.groups['location']['glt_x'])
    glt[...,1] = np.array(img_nc.groups['location']['glt_y'])
    valid_glt = np.all(glt != 0, axis=-1)
    glt[valid_glt] -= 1
    map_info =  f'{{Geographic Lat/Lon, 1, 1, {gt[0]}, {gt[3]}, {gt[1]}, {gt[5]*-1},WGS-84}}'

    # Get spatial and temporal extents
    latitude = np.array(img_nc.groups['location']['lat'])
    longitude = np.array(img_nc.groups['location']['lon'])

    corner_1 = [longitude[0,0],  latitude[0,0]]
    corner_2 = [longitude[0,-1], latitude[0,-1]]
    corner_3 = [longitude[-1,-1],latitude[-1,-1]]
    corner_4 = [longitude[-1,0], latitude[-1,0]]

    datetime = img_file.split('_')[4]
    date = dt.datetime.strptime(datetime[:8],'%Y%m%d')

    if obs_file:
        obs_nc = nc.Dataset(obs_file)
        utc_time = np.copy(obs_nc['obs'][:,:,9])
        start_time = utc_time.min()
        start_hour = int(start_time)
        start_minute = (start_time-start_hour)*60
        start_second = round((start_minute - int(start_minute))*60)
        start_minute = int(start_minute)
        start_delta = dt.timedelta(hours = start_hour,
                                    minutes = start_minute,
                                    seconds = start_second)

        end_time = utc_time.max()
        end_hour = int(end_time)
        end_minute = (end_time-end_hour)*60
        end_second = round((end_minute - int(end_minute))*60)
        end_minute = int(end_minute)
        end_delta = dt.timedelta(hours = end_hour,
                                    minutes = end_minute,
                                    seconds = end_second)
        start_time = date+start_delta
        end_time = date+end_delta
    else:
        start_time = date
        end_time = date

    # Process radiance/reflectace datacubes
    data_header = ht.io.envi_header_dict()
    data_header['lines']= glt.shape[0]
    data_header['samples']= glt.shape[1]
    data_header['bands']= data.shape[2]
    data_header['byte order']= 0
    data_header['data ignore value']= NO_DATA_VALUE
    data_header['data type']= 4
    data_header['interleave']= 'bil'
    data_header['map info'] = map_info

    data_gcs = f'{temp_dir}/data_gcs'

    writer = ht.io.WriteENVI(data_gcs,data_header)
    data_prj = np.full((glt.shape[0], glt.shape[1]),NO_DATA_VALUE)

    print(f"Exporting EMIT {product} dataset")

    for band_num in range(data.shape[2]):
        data_prj[valid_glt] = data[:,:,band_num][glt[valid_glt, 1], glt[valid_glt,0]]
        writer.write_band(data_prj, band_num)
        data_prj[:] = NO_DATA_VALUE

    zone,direction  = utm_zone(longitude, latitude)
    if direction == 'North':
        epsg_dir = 6
    else:
        epsg_dir = 7

    out_crs = "EPSG:32%s%02d" % (epsg_dir,zone)

    # Reproject to UTM
    print(f'Projecting data to UTM Zone {zone} {direction} at 60m resolution')
    data_utm = f'{out_dir}/SISTER_EMIT_{product}_{datetime}_{crid}.bin'
    warp_command = f'gdalwarp -overwrite -t_srs {out_crs} -tr 60 60 -r near -of ENVI {data_gcs} {data_utm}'
    os.system(warp_command)

    #Update header file
    data_header_file = data_utm.replace('.bin','.hdr')
    data_header = ht.io.envi.parse_envi_header(data_header_file)
    data_header['description']= 'Radiance micro-watts/cm^2/nm/sr'
    data_header['band names']=[]
    data_header['fwhm']=fwhm.tolist()
    data_header['wavelength']= waves.tolist()
    data_header['wavelength units'] = 'nanometers'
    data_header['start acquisition time'] = start_time.strftime('%Y-%m-%dT%H:%M:%SZ')
    data_header['end acquisition time'] = end_time.strftime('%Y-%m-%dT%H:%M:%SZ')
    data_header['bounding box'] =[corner_1,corner_2,corner_3,corner_4]
    data_header['coordinate system string'] = []
    data_header['sensor type'] ='EMIT'
    ht.io.envi.write_envi_header(data_header_file,data_header)

    if export_loc:
        print(f"Exporting EMIT location dataset")
        # Process location datacube
        loc_header = ht.io.envi_header_dict()
        loc_header['lines']= glt.shape[0]
        loc_header['samples']= glt.shape[1]
        loc_header['bands']= 3
        loc_header['byte order']= 0
        loc_header['data ignore value']= NO_DATA_VALUE
        loc_header['data type']= 4
        loc_header['interleave']= 'bil'
        loc_header['map info'] = map_info

        loc_gcs = f'{temp_dir}/loc_gcs'
        writer = ht.io.WriteENVI(loc_gcs,loc_header)

        loc_band = np.full((glt.shape[0], glt.shape[1]),NO_DATA_VALUE)

        for band_num,band_name in enumerate(['lon','lat','elev']):
            band = np.copy(img_nc.groups['location'][band_name])
            loc_band[valid_glt] = band[glt[valid_glt, 1], glt[valid_glt,0]]
            writer.write_band(loc_band, band_num)
            loc_band[:] = NO_DATA_VALUE

        loc_utm =   f'{out_dir}/SISTER_EMIT_{product}_{datetime}_{crid}_LOC.bin'

        print(f'Projecting location datacube to UTM Zone {zone} {direction} at 60m resolution')
        warp_command = f'gdalwarp -overwrite -t_srs {out_crs} -tr 60 60 -r near -of ENVI {loc_gcs} {loc_utm}'
        os.system(warp_command)

        #Update header file
        loc_header_file = loc_utm.replace('.bin','.hdr')
        loc_header = ht.io.envi.parse_envi_header(loc_header_file)
        loc_header['band names'] = ['longitude', 'latitude','elevation']
        loc_header['description']= 'Location datacube'
        loc_header['coordinate system string'] = []
        loc_header['start acquisition time'] = start_time.strftime('%Y-%m-%dT%H:%M:%SZ')
        loc_header['end acquisition time'] = end_time.strftime('%Y-%m-%dT%H:%M:%SZ')
        loc_header['bounding box'] =[corner_1,corner_2,corner_3,corner_4]
        loc_header['sensor type'] ='EMIT'

        ht.io.envi.write_envi_header(loc_header_file,loc_header)

    if obs_file:
        # Process observation datacube
        print(f"Exporting EMIT observation dataset")

        obs_header = ht.io.envi_header_dict()
        obs_header['lines']= glt.shape[0]
        obs_header['samples']= glt.shape[1]
        obs_header['bands']= 11
        obs_header['byte order']= 0
        obs_header['data ignore value']= NO_DATA_VALUE
        obs_header['data type']= 4
        obs_header['interleave']= 'bil'
        obs_header['map info'] = map_info

        obs_gcs = f'{temp_dir}/obs_gcs'
        writer = ht.io.WriteENVI(obs_gcs,obs_header)

        obs_band = np.full((glt.shape[0], glt.shape[1]),NO_DATA_VALUE)

        for band_num in range(11):
            band = np.copy(obs_nc['obs'][:,:,band_num])
            obs_band[valid_glt] = band[glt[valid_glt, 1], glt[valid_glt,0]]
            writer.write_band(obs_band, band_num)
            loc_band[:] = NO_DATA_VALUE

        obs_utm = f'{out_dir}/SISTER_EMIT_{product}_{datetime}_{crid}_OBS.bin'

        print(f'Projecting observation datacube to UTM Zone {zone} {direction} at 60m resolution')
        warp_command = f'gdalwarp -overwrite -t_srs {out_crs} -tr 60 60 -r near -of ENVI {obs_gcs} {obs_utm}'
        os.system(warp_command)

        #Update header file
        obs_header_file = obs_utm.replace('.bin','.hdr')
        obs_header = ht.io.envi.parse_envi_header(obs_header_file)
        obs_header['band names'] = ['path length', 'to-sensor azimuth',
                                    'to-sensor zenith','to-sun azimuth',
                                      'to-sun zenith','phase', 'slope',
                                      'aspect', 'cosine i','UTC time','earth-sun distance']
        obs_header['description'] = 'Observation datacube'
        obs_header['coordinate system string'] = []
        obs_header['start acquisition time'] = start_time.strftime('%Y-%m-%dT%H:%M:%SZ')
        obs_header['end acquisition time'] = end_time.strftime('%Y-%m-%dT%H:%M:%SZ')
        obs_header['bounding box'] =[corner_1,corner_2,corner_3,corner_4]
        obs_header['sensor type'] ='EMIT'

        ht.io.envi.write_envi_header(obs_header_file,obs_header)


    #Delete GDAL generated XMLs
    for xml in glob.glob(f'{out_dir}/*.xml'):
        os.remove(xml)
