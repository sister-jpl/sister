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


def loc_export(loc_file,longitude,latitude,elevation):
    '''Export location datasets to disk
    '''
    loc_header = envi_header_dict()
    loc_header['lines']= longitude.shape[0]
    loc_header['samples']= longitude.shape[1]
    loc_header['bands']= 3
    loc_header['interleave']= 'bil'
    loc_header['data type'] = 4
    loc_header['band_names'] = ['Longitude', 'Latitude','Elevation']
    loc_header['byte order'] = 0

    writer = WriteENVI(loc_file,loc_header)
    writer.write_band(longitude,0)
    writer.write_band(latitude,1)
    writer.write_band(elevation,2)


def obs_export(obs_file,pathlength,sensor_az,sensor_zn,solar_az,solar_zn,phase,slope,aspect,cosine_i,utc_time):
    '''Export observables datasets to disk
    '''
    obs_header = envi_header_dict()
    obs_header['lines']= pathlength.shape[0]
    obs_header['samples']= pathlength.shape[1]
    obs_header['bands']= 10
    obs_header['interleave']= 'bil'
    obs_header['data type'] = 4
    obs_header['byte order'] = 0
    obs_header['band_names'] = ['path length', 'to-sensor azimuth',
                                'to-sensor zenith','to-sun azimuth',
                                  'to-sun zenith','phase', 'slope',
                                  'aspect', 'cosine i','UTC time']

    writer = WriteENVI(obs_file,obs_header)
    writer.write_band(pathlength,0)
    writer.write_band(sensor_az,1)
    writer.write_band(sensor_zn,2)
    writer.write_band(solar_az,3)
    writer.write_band(solar_zn,4)
    writer.write_band(phase,5)
    writer.write_band(slope,6)
    writer.write_band(aspect,7)
    writer.write_band(cosine_i,8)
    writer.write_band(utc_time,9)










