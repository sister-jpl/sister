#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SISTER
Space-based Imaging Spectroscopy and Thermal PathfindER
Author: Adam Chlus
"""
import os
import glob
import tarfile
import numpy as np
import hytools as ht
from hytools.io.envi import WriteENVI,envi_header_dict
from rtree import index
from scipy.spatial import cKDTree


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

    # Create a simple spatial index to find intersecting tiles
    idx = index.Index(properties=index.Property())
    tiles = glob.glob(elev_dir + '*.tar.gz')
    for i, tile in enumerate(tiles):
        lat,ign,lon = os.path.basename(tile).split('_')[3:6]
        if 'W' in lon:
            lon = -1*float(lon[1:])
        else:
            lon = float(lon[1:])
        if 'S' in lat:
            lat = -1*float(lat[1:])
        else:
            lat = float(lat[1:])
        idx.insert(i,(lon,lat,lon+1,lat+1))

    tiles_inter = [tiles[n] for n in idx.intersection((lon_min, lat_min, lon_max, lat_max))]

    if len(tiles_inter) > 0:
        print("\nIntersecting elevation tiles:")
        for tile in tiles_inter:
            print('\t' + os.path.basename(tile))
            with tarfile.open(tile, 'r') as tar_ref:
                tar_ref.extractall(temp_dir)

        print('Merging DEM tiles')
        dem_file  = '%stemp_dem' % temp_dir
        os.system('gdal_merge.py -o %s -of ENVI %s*.dt2' % (dem_file,temp_dir))

        dem_obj = ht.HyTools()
        dem_obj.read_file(dem_file, 'envi')

        ulx = float(dem_obj.map_info[3])
        uly = float(dem_obj.map_info[4])
        pix = float(dem_obj.map_info[5])

        dem_lat,dem_lon = np.indices((dem_obj.lines,dem_obj.columns))

        dem_xl = int((lon_min-ulx)//pix)
        dem_xr = int((lon_max-ulx)//pix)
        dem_yu = int((uly-lat_max)//pix)
        dem_yd = int((uly-lat_min)//pix)

        dem_subset = dem_obj.get_chunk(dem_xl,dem_xr,dem_yu,dem_yd)
        dem_lat,dem_lon = np.indices(dem_subset.shape[:2])
        dem_lat = (lat_max- dem_lat*pix).flatten()
        dem_lon = (lon_min+ dem_lon*pix).flatten()

        #Create spatial index and nearest neighbor sample
        src_points =np.concatenate([np.expand_dims(dem_lon,axis=1),np.expand_dims(dem_lat,axis=1)],axis=1)
        tree = cKDTree(src_points,balanced_tree= False)

        dst_points = np.concatenate([longitude.flatten()[:,np.newaxis],
                                     latitude.flatten()[:,np.newaxis]],
                                     axis=1)

        dists, indexes = tree.query(dst_points,k=1)
        indices_int = np.unravel_index(indexes,(dem_subset.shape[0],
                                                dem_subset.shape[1]))
        elevation = dem_subset[indices_int[0],indices_int[1]].reshape(longitude.shape)

    else:
        constant_elev = float(input("\nNo overlapping tiles found, enter constant elevation for scene (m): "))
        elevation = np.ones(longitude.shape.shape) * constant_elev

    #Set negative elevations to 0
    if np.sum(elevation<0) > 0:
        print('Elevations below sea level found, setting to 0m')
        elevation[elevation<0] =0

    return elevation

def slope_aspect(elevation,temp_dir):

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

    print('Calculating slope')
    os.system('gdaldem slope -of ENVI %s %s'% (dem_file,slope_file))

    print('Calculating aspect')
    os.system('gdaldem aspect -f ENVI %s %s' % (dem_file,aspect_file))

    asp_obj = ht.HyTools()
    asp_obj.read_file(aspect_file, 'envi')
    aspect =asp_obj.get_band(0)

    slp_obj = ht.HyTools()
    slp_obj.read_file(slope_file, 'envi')
    slope =slp_obj.get_band(0)

    return slope,aspect










