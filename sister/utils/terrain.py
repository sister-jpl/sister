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


def terrain_generate(longitude,latitude,elev_dir,temp_dir):
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

        # Retry reading tile list if fails
        retry = 0
        while (retry < 10) & (len(tiles) != 26450):
            tiles = pd.read_csv(elev_dir + 'tileList.txt',header = None).values.flatten()
            retry+=1

        if len(tiles) != 26450:
            raise ValueError('Failed to download tile list.')

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
        raise ValueError('No overlapping Copernicus DEM tiles found.')

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

    slope_file =  '%stemp_slope' % temp_dir
    aspect_file =  '%stemp_aspect' % temp_dir

    logging.info('Calculating slope')
    os.system('gdaldem slope -compute_edges -of ENVI %s %s'% (dem_file,slope_file))

    logging.info('Calculating aspect')
    os.system('gdaldem aspect -compute_edges -of ENVI %s %s' % (dem_file,aspect_file))

    terrain_arrs = []

    for file in [dem_file,slope_file,aspect_file]:

        trr_obj = ht.HyTools()
        trr_obj.read_file(file, 'envi')

        ulx = float(trr_obj.map_info[3])
        uly = float(trr_obj.map_info[4])
        pix_x = float(trr_obj.map_info[5])
        pix_y = float(trr_obj.map_info[6])

        trr_lat,trr_lon = np.indices((trr_obj.lines,trr_obj.columns))

        trr_xl = int((lon_min-ulx)//pix_x)
        trr_xr = int((lon_max-ulx)//pix_x)
        trr_yu = int((uly-lat_max)//pix_y)
        trr_yd = int((uly-lat_min)//pix_y)

        trr_subset = trr_obj.get_chunk(trr_xl,trr_xr,trr_yu,trr_yd)
        trr_lat,trr_lon = np.indices(trr_subset.shape[:2])
        trr_lat = (lat_max- trr_lat*pix_y).flatten()
        trr_lon = (lon_min+ trr_lon*pix_x).flatten()

        #Create spatial index and nearest neighbor sample
        src_points =np.concatenate([np.expand_dims(trr_lon,axis=1),
                                    np.expand_dims(trr_lat,axis=1)],axis=1)
        tree = cKDTree(src_points,balanced_tree= False)

        dst_points = np.concatenate([longitude.flatten()[:,np.newaxis],
                                     latitude.flatten()[:,np.newaxis]],
                                     axis=1)

        indexes = tree.query(dst_points,k=1)[1]
        indices_int = np.unravel_index(indexes,(trr_subset.shape[0],
                                                trr_subset.shape[1]))
        terrain = trr_subset[indices_int[0],indices_int[1]].reshape(longitude.shape)

        #Set negative elevations to 0
        if (np.sum(terrain<0) > 0) and 'dem' in file:
            logging.warning('Elevations below sea level found, setting to 0m')
            terrain[terrain<0] =0

        terrain_arrs.append(terrain)

    return terrain_arrs
