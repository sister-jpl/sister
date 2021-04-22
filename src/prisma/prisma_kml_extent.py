#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SISTER
Space-based Imaging Spectroscopy and Thermal PathfindER
Author: Adam Chlus
"""

import argparse
import datetime as dt
import glob
import os
import shutil
import zipfile
import h5py
from collections import OrderedDict
import fiona
from fiona.crs import from_epsg

def main():
    '''
    This function exports the extents of the scene a KML to be used for Landsat Scene selection

    base_name: Scene basename (ex: 20200205001044_20200205001049_0001)
    zip_dir: Directory of zipped L2c data
    out_dir: Output directory
    temp_dir: Temporary directory

    '''
    parser = argparse.ArgumentParser(description = "Export PRISMA scene extent to KML")
    parser.add_argument('base_name',help="Image basename", type = str)
    parser.add_argument('zip_dir',help="Output directory", type = str)
    parser.add_argument('out_dir',help="Output directory", type = str)
    parser.add_argument('temp_dir',help="Temporary data directory", type = str)

    args = parser.parse_args()
    zip_dir  = args.zip_dir
    base_name = args.base_name
    temp_dir =  "%s/temp_%s/" % (args.temp_dir,base_name)
    out_dir = '%s/PRS_%s/' %  (args.out_dir,base_name)

    ##Testing
    # parser = argparse.ArgumentParser(description = "Export PRISMA scene extent to KML")
    # args = parser.parse_args()
    # zip_dir  = '/data2/prisma/zip/'
    # base_name = "20200712223611_20200712223615_0001"
    # out_dir = '/data2/prisma/envi/PRS_%s/' % base_name
    # temp_dir = '/data2/temp/temp_%s/' % base_name

    for file in glob.glob("%s*L2*%s*" % (zip_dir,base_name)):
        zip_base  =os.path.basename(file)
        print('Unzipping %s' % zip_base)
        with zipfile.ZipFile(file,'r') as zipped:
            zipped.extractall(temp_dir)
    l2  = '%sPRS_L2C_STD_%s.he5' % (temp_dir,base_name)

    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)
    if not os.path.isdir(temp_dir):
        os.mkdir(temp_dir)

    l_obj = h5py.File(l2,'r')
    subdir = [x for x in l_obj['HDFEOS']["SWATHS"].keys() if 'HCO' in x][0]
    geo =  l_obj['HDFEOS']["SWATHS"][subdir]['Geolocation Fields']

    # Get extents of image
    lon1 = geo['Longitude'][0,0]
    lon2 = geo['Longitude'][0,-1]
    lon3 = geo['Longitude'][-1,-1]
    lon4 = geo['Longitude'][-1,0]

    lat1 = geo['Latitude'][0,0]
    lat2 = geo['Latitude'][0,-1]
    lat3 = geo['Latitude'][-1,-1]
    lat4 = geo['Latitude'][-1,0]

    extent = {'geometry' : {'coordinates': [[(lon1,lat1),
                            (lon2,lat2),
                            (lon3,lat3),
                            (lon4,lat4),
                            (lon1,lat1)]],
                            'type': 'Polygon'},
              'properties' : OrderedDict([])}

    schema = {'geometry': 'Polygon',
              'properties': OrderedDict([])}

    temp_geojson = '%s/temp.geojson' % temp_dir
    with fiona.open(temp_geojson,
                    'w',
                    driver="GeoJSON",
                    crs=from_epsg(4326),
                    schema=schema) as c:
        c.write(extent)

    #Convert to KML
    print('Converting to KML')
    os.system('ogr2ogr %sPRS_%s_extent.kml %s' % (out_dir,base_name,temp_geojson))
    print('Deleting temporary files')
    shutil.rmtree(temp_dir)

if __name__== "__main__":
    main()


















