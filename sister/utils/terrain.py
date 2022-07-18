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


https://copernicus-dem-30m.s3.amazonaws.com/

"""

import os
import glob
import tarfile
import logging
import numpy as np
import hytools as ht
from hytools.io.envi import WriteENVI,envi_header_dict
import pandas as pd
from rtree import index
from scipy.spatial import cKDTree
from .misc import download_file


def dem_generate(longitude,latitude,elev_dir,temp_dir):
    '''
    Args:
        longitude (float): Longitude array
        latitude (float): Latitude array
        elev_dir (str): Directory of zipped elevation tiles
        temp_dir (str): Temporary output directory

    Returns:
        dem (np.array): Elevation array.

    '''
    # Get extents of image
    lon_min = longitude.min()
    lon_max = longitude.max()
    lat_min = latitude.min()
    lat_max = latitude.max()

    if 'aws' in elev_dir:
        tiles = pd.read_csv(elev_dir + 'tileList.txt',header = None).values.flatten()
    else:
        tiles = glob.glob(elev_dir + '*.tar.gz')

    idx = index.Index(properties=index.Property())

    #Get list of intersecting tiles
    for i, tile in enumerate(tiles):
        lat,ign,lon = os.path.basename(tile).replace('_COG','').split('_')[3:6]
        if 'W' in lon:
            lon = -1*float(lon[1:])
        else:
            lon = float(lon[1:])
        if 'S' in lat:
            lat = -1*float(lat[1:])
        else:
            lat = float(lat[1:])
        idx.insert(i,(lon,lat,lon+1,lat+1))
    tiles_intersect = [tiles[n] for n in idx.intersection((lon_min, lat_min, lon_max, lat_max))]

    if len(tiles_intersect) == 0:
        constant_elev = float(input("No overlapping tiles found, enter constant elevation for scene (m): "))
        elevation = np.ones(longitude.shape) * constant_elev
    else:
        tile_string = "Found %s intersecting elevation tiles:" % len(tiles_intersect)
        for tile in tiles_intersect:
            tile_string+= '\n\t%s' % tile

            if 'aws' in elev_dir:
                tile_url = "%s%s/%s.tif" % (elev_dir,tile,tile)
                tile_file = "%s%s.tif" % (temp_dir,tile)
                download_file(tile_file,tile_url)
            else:
                with tarfile.open(tile, 'r') as tar_ref:
                    tar_ref.extractall(temp_dir)
        logging.info(tile_string)

    logging.info('Merging DEM tiles')
    dem_file  = '%stemp_dem' % temp_dir
    os.system('gdal_merge.py -o %s -of ENVI %sCopernicus_DSM*' % (dem_file,temp_dir))

    dem_obj = ht.HyTools()
    dem_obj.read_file(dem_file, 'envi')

    ulx = float(dem_obj.map_info[3])
    uly = float(dem_obj.map_info[4])
    pix_x = float(dem_obj.map_info[5])
    pix_y = float(dem_obj.map_info[6])

    dem_lat,dem_lon = np.indices((dem_obj.lines,dem_obj.columns))

    dem_xl = int((lon_min-ulx)//pix_x)
    dem_xr = int((lon_max-ulx)//pix_x)
    dem_yu = int((uly-lat_max)//pix_y)
    dem_yd = int((uly-lat_min)//pix_y)

    dem_subset = dem_obj.get_chunk(dem_xl,dem_xr,dem_yu,dem_yd)
    dem_lat,dem_lon = np.indices(dem_subset.shape[:2])
    dem_lat = (lat_max- dem_lat*pix_y).flatten()
    dem_lon = (lon_min+ dem_lon*pix_x).flatten()

    #Create spatial index and nearest neighbor sample
    src_points =np.concatenate([np.expand_dims(dem_lon,axis=1),
                                np.expand_dims(dem_lat,axis=1)],axis=1)
    tree = cKDTree(src_points,balanced_tree= False)

    dst_points = np.concatenate([longitude.flatten()[:,np.newaxis],
                                 latitude.flatten()[:,np.newaxis]],
                                 axis=1)

    indexes = tree.query(dst_points,k=1)[1]
    indices_int = np.unravel_index(indexes,(dem_subset.shape[0],
                                            dem_subset.shape[1]))
    elevation = dem_subset[indices_int[0],indices_int[1]].reshape(longitude.shape)

    #Set negative elevations to 0
    if np.sum(elevation<0) > 0:
        logging.warning('Elevations below sea level found, setting to 0m')
        elevation[elevation<0] =0

    return elevation

def slope_aspect(elevation,temp_dir):
    ''' Use GDAL to calculte slope and aspect
    '''

    dem_dict  = envi_header_dict()
    dem_dict ['lines']= elevation.shape[0]
    dem_dict ['samples']= elevation.shape[1]
    dem_dict ['bands']= 1
    dem_dict ['interleave']= 'bsq'
    dem_dict ['data type'] = 4

    dem_file = '%stemp_dem_clip' % temp_dir
    writer = WriteENVI(dem_file,dem_dict)
    writer.write_band(elevation,0)

    slope_file =  '%s_slope' % temp_dir
    aspect_file =  '%s_aspect' % temp_dir

    logging.info('Calculating slope')
    os.system('gdaldem slope -of ENVI %s %s'% (dem_file,slope_file))

    logging.info('Calculating aspect')
    os.system('gdaldem aspect -f ENVI %s %s' % (dem_file,aspect_file))

    asp_obj = ht.HyTools()
    asp_obj.read_file(aspect_file, 'envi')
    aspect =asp_obj.get_band(0)

    slp_obj = ht.HyTools()
    slp_obj.read_file(slope_file, 'envi')
    slope =slp_obj.get_band(0)

    return slope,aspect
