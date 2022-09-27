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
import datetime as dt
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

    # Get image bounds coordinates
    corner_1 = [lon[0,0],  lat[0,0]]
    corner_2 = [lon[0,-1], lat[0,-1]]
    corner_3 = [lon[-1,-1],lat,-1]]
    corner_4 = [lon[-1,0], latitude[-1,0]]


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

    return corner_1,corner_2,corner_3,corner_4

def time_correct(obs_ort_file):
    obs = ht.HyTools()
    obs.read_file(obs_ort_file,'envi')

    utc_time = obs.get_band(9)[obs.mask['no_data']]
    if (utc_time.max() > 24) | (utc_time.min() < 0):
        hour = np.ones((obs.lines,obs.columns))
        hour += int(obs.base_name[12:14])
        hour += int(obs.base_name[14:16])/60
        hour += int(obs.base_name[16:18])/3600
        hour[~obs.mask['no_data']] = obs.no_data

        obs.load_data(mode = 'r+')

        if obs.interleave == 'bip':
            obs.data[:,:,9] =hour
        elif obs.interleave == 'bil':
            obs.data[:,9,:] =hour
        elif obs.interleave == 'bsq':
            obs.data[9,:,:] =hour

        obs.close_data()

def get_spatiotemporal_extents(loc_file,obs_ort_file):
    ''' Get image acquiition start and end time and bounding box of image
    '''

    loc = ht.HyTools()
    loc.read_file(loc_file,'envi')

    bounds = [loc.get_band(0).min(),
             loc.get_band(0).max(),
             loc.get_band(1).min(),
             loc.get_band(1).max()]

    obs = ht.HyTools()
    obs.read_file(obs_ort_file,'envi')

    utm_time = obs.get_band(9)
    start_time = utm_time[utm_time != obs.no_data].min()
    start_hour = int(start_time)
    start_minute = (start_time-start_hour)*60
    start_second = round((start_minute - int(start_minute))*60)
    start_minute = int(start_minute)

    start_delta = dt.timedelta(hours = start_hour,
                                minutes = start_minute,
                                seconds = start_second)


    utm_time = obs.get_band(9)
    end_time = utm_time[utm_time != obs.no_data].max()
    end_hour = int(end_time)
    end_minute = (end_time-end_hour)*60
    end_second = round((end_minute - int(end_minute))*60)
    end_minute = int(end_minute)

    end_delta = dt.timedelta(hours = end_hour,
                                minutes = end_minute,
                                seconds = end_second)

    return bounds,start_delta,end_delta


def preprocess(input_tar,out_dir,temp_dir,res = 0):
    '''
    '''

    base_name = os.path.basename(input_tar)

    if base_name.startswith('ang'):
        base_name = base_name[:18]
        date = dt.datetime.strptime(base_name[3:10],'%Y%m%d')
    elif base_name.startswith('f'):
        base_name = base_name[:16]
        date = dt.datetime.strptime("20%s" % base_name[1:7],'%Y%m%d')
    else:
        raise ValueError('Unrecognized sensor')

    file = tarfile.open(input_tar)
    tar_contents = [temp_dir+c for c in file.getnames()]
    file.extractall(temp_dir)
    file.close()

    #AVIRIS NG
    if base_name.startswith('ang'):
        instrument = 'AVNG'
        obs_ort_file = [x for x in tar_contents if x.endswith('obs_ort')][0]
        rdn_file = [x for x in tar_contents if x.endswith('img')][0]
        loc_file = [x for x in tar_contents if x.endswith('loc')][0]
        glt_file = [x for x in tar_contents if x.endswith('glt')][0]
        corner_1,corner_2,corner_3,corner_4 = create_loc_ort(loc_file,glt_file)
        time_correct(obs_ort_file)

    #AVIRIS Classic
    elif base_name.startswith('f'):
        instrument = 'AVCL'
        obs_ort_file = [x for x in tar_contents if x.endswith('obs_ort')][0]
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

        corner_1,corner_2,corner_3,corner_4 = create_loc_ort(loc_file,glt_file)

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

    loc_ort_file = loc_file+'_ort'

    #Get spatial and temporal extents of datsaet
    bounds,start_delta,end_delta = get_spatiotemporal_extents(loc_file,obs_ort_file)
    start_time =date+start_delta
    end_time =date+end_delta
    datetime = start_time.strftime('%Y%m%dT%H%M%S')
    lon_min,lon_max,lat_min,lat_max = bounds

    out_dir = "%s/SISTER_%s_%s_L1B_RDN_000/" % (out_dir,instrument,datetime)

    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)

    for file in [obs_ort_file,rdn_file,loc_ort_file]:

        if 'loc' in file:
            product = 'LOC'
        elif 'obs' in file:
            product = 'OBS'
        else:
            product = 'RDN'

        new_file = "%s/SISTER_%s_%s_L1B_%s_000" % (os.path.dirname(file),instrument,datetime,product)
        new_file_hdr = new_file+ '.hdr'

        os.rename(file,new_file)
        os.rename(file+ '.hdr',new_file_hdr)

        #Add spatial and temporal extents to header file
        header = ht.io.envi.parse_envi_header(new_file_hdr)

        # Create new header
        clean_header = {}
        clean_header['bands'] = header['bands']
        clean_header['samples'] = header['samples']
        clean_header['lines'] =  header['lines']
        clean_header['interleave'] =  header['interleave']
        clean_header['data type'] =  header['data type']
        clean_header['map info'] =  header['map info']
        clean_header['sensor type'] =instrument

        if product == 'RDN':
            clean_header['description'] = 'Radiance micro-watts/cm^2/nm/sr'
            clean_header['wavelength'] =  header['wavelength']
            clean_header['fwhm'] =  header['fwhm']
            clean_header['wavelength units'] =  'nanometers'
        elif product == 'OBS':
            clean_header['description'] = 'Observation datacube'
            clean_header['band names'] = ['path length','to-sensor azimuth',
                                          'to-sensor zenith','to-sun azimuth',
                                          'to-sun zenith','phase','slope','aspect',
                                          'cosine i','UTC time']
        elif product == 'LOC':
            clean_header['description'] = 'Location datacube'
            clean_header['band names'] = ['longitude','latitude','elevation']

        clean_header['start acquisition time'] = start_time.strftime('%Y-%m-%dT%H:%M:%SZ')
        clean_header['end acquisition time'] = end_time.strftime('%Y-%m-%dT%H:%M:%SZ')
        clean_header['bounding box'] =[corner_1,corner_2,corner_3,corner_4]

        ht.io.envi.write_envi_header(new_file,clean_header,mode = 'w')

        if res==0:
            shutil.move(new_file,out_dir)
            shutil.move(new_file_hdr,out_dir)
        else:
            resample(new_file,out_dir,res)

    shutil.rmtree(tar_contents[0])
