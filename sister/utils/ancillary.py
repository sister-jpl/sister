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

import logging
from hytools.io.envi import WriteENVI,envi_header_dict
import numpy as np

def loc_export(loc_file,longitude,latitude,elevation):
    '''Export location datasets to disk
    '''
    loc_header = envi_header_dict()
    loc_header['description']= 'Location datacube'
    loc_header['lines']= longitude.shape[0]
    loc_header['samples']= longitude.shape[1]
    loc_header['bands']= 3
    loc_header['interleave']= 'bil'
    loc_header['data type'] = 4
    loc_header['band names'] = ['longitude', 'latitude','elevation']
    loc_header['byte order'] = 0

    loc = np.array([longitude,latitude,elevation])
    loc = np.moveaxis(loc,0,-1)

    WriteENVI(loc_file,loc_header)
    outfile = open(loc_file, 'rb+')

    line_num=0
    for line in loc:
        outfile.seek(line_num * loc_header['samples'] *  loc_header['bands'] * np.dtype('float32').itemsize)
        outfile.write(line.T.astype('float32').tobytes())
        line_num+=1
    outfile.close()


def obs_export(obs_file,pathlength,sensor_az,sensor_zn,solar_az,solar_zn,phase,slope,aspect,cosine_i,utc_time):
    '''Export observables datasets to disk
    '''
    obs_header = envi_header_dict()
    obs_header['description'] = 'Observation datacube'
    obs_header['lines']= sensor_az.shape[0]
    obs_header['samples']= sensor_az.shape[1]
    obs_header['bands']= 10
    obs_header['interleave']= 'bil'
    obs_header['data type'] = 4
    obs_header['byte order'] = 0
    obs_header['band names'] = ['path length', 'to-sensor azimuth',
                                'to-sensor zenith','to-sun azimuth',
                                  'to-sun zenith','phase', 'slope',
                                  'aspect', 'cosine i','UTC time']

    WriteENVI(obs_file,obs_header)

    obs = np.array([pathlength,sensor_az,sensor_zn,
                    solar_az,solar_zn,phase,
                    slope,aspect,cosine_i,utc_time])
    obs = np.moveaxis(obs,0,-1)

    outfile = open(obs_file, 'rb+')

    line_num=0
    for line in obs:
        outfile.seek(line_num * obs_header['samples'] *  obs_header['bands'] * np.dtype('float32').itemsize)
        outfile.write(line.T.astype('float32').tobytes())
        line_num+=1
    outfile.close()
