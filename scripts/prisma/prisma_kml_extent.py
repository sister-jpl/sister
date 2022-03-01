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

import argparse
from collections import OrderedDict
import os
import shutil
import zipfile
import h5py
import fiona
from fiona.crs import from_epsg

def main():
    '''
    This function exports the extents of a PRISMA scene to a KML file.

    l1_zip: Path to zipped L1 data
    out_dir: Output directory
    temp_dir: Temporary directory

    '''
    parser = argparse.ArgumentParser(description = "Export PRISMA scene extent to KML")
    parser.add_argument('l1_zip',help="Input zipped L1 HDF path", type = str)
    parser.add_argument('out_dir',help="Output directory", type = str)
    parser.add_argument('temp_dir',help="Temporary data directory", type = str)

    args = parser.parse_args()
    l1_zip  = args.l1_zip
    base_name = os.path.basename(l1_zip).replace('.zip','')
    temp_dir =  "%s/temp_%s/" % (args.temp_dir,base_name)
    out_dir = args.out_dir
    kml_file =  "%s/%s" % (out_dir,
                            os.path.basename(l1_zip).replace('zip','kml'))
    if os.path.isfile(kml_file):
        print('%s exists' % kml_file)
        return

    print('Unzipping %s' % base_name)
    with zipfile.ZipFile(l1_zip,'r') as zipped:
        zipped.extractall(temp_dir)
    l1  = '%s/%s.he5' % (temp_dir,base_name)

    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)
    if not os.path.isdir(temp_dir):
        os.mkdir(temp_dir)

    l_obj = h5py.File(l1,'r')
    geo =  l_obj['HDFEOS']["SWATHS"]['PRS_L1_HCO']['Geolocation Fields']

    # Get extents of image
    lon1 = geo['Longitude_VNIR'][0,0]
    lon2 = geo['Longitude_VNIR'][0,-1]
    lon3 = geo['Longitude_VNIR'][-1,-1]
    lon4 = geo['Longitude_VNIR'][-1,0]

    lat1 = geo['Latitude_VNIR'][0,0]
    lat2 = geo['Latitude_VNIR'][0,-1]
    lat3 = geo['Latitude_VNIR'][-1,-1]
    lat4 = geo['Latitude_VNIR'][-1,0]

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
    print('Converting to KML')
    os.system('ogr2ogr %s %s >/dev/null 2>&1' % (kml_file,temp_geojson))
    print('Deleting temporary files')
    shutil.rmtree(temp_dir)

if __name__== "__main__":
    main()
