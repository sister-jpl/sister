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
from sister.sensors import prisma


def main():
    '''
    This function is a comman line wrapper around the PRISMA hdf to envi converter
    function he5_to_envi and generates the following files:

        *_rad* : Merged and optionally shift corrected radiance cube
        *_obs* : Observables file in the format of JPL obs files:
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
        *_loc* : Location file in the following format:
                1. Longitude (decimal degrees)
                2. Longitude (decimal degrees)
                3. Elevation (m)

    l1(str): L1 zipped radiance data product path
    out_dir(str): Output directory of ENVI radiance datasets
    temp_dir(str): Temporary directory for intermediate files
    elev_dir (str): Directory path to zipped Copernicus DSM tiles
    shift (str) : Pathname of wavelength shift correction surface file
    match (str or list) : Perform landsat image matching (recommended)
    proj (bool) : Project image to UTM grid
    res (int) : Resolution of projected image, 30 should be one of its factors (90,120,150.....)

    '''

    parser = argparse.ArgumentParser(description = "Convert PRISMA HDF to ENVI format")
    parser.add_argument('l1_zip',help="Path to zipped L1 HDF radiance file", type = str)
    parser.add_argument('out_dir',help="Output directory", type = str)
    parser.add_argument('temp_dir',help="Temporary directory", type = str)
    parser.add_argument('elev_dir',help="Directory path to zipped Copernicus DSM tiles", type = str)
    parser.add_argument("-shift", help="Path to wavelength shift surface", type = str,required=False,default=None)
    parser.add_argument("-match", help="Perform landsat image matching", action='store_true',required=False)
    parser.add_argument("-proj", help="Project image to UTM grid",action='store_true', required=False)
    parser.add_argument("-res", help="Projected image resolution in meters",type = int, required=False,default= 30)

    args = parser.parse_args()

    prisma.he5_to_envi(args.l1_zip,args.out_dir,args.temp_dir,args.elev_dir,
                       shift = args.shift,match=args.match,proj = args.proj,
                       res = args.res)


if __name__=='__main__':

    main()
