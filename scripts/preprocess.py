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
import os
from sister.sensors import prisma,aviris


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
    elev_dir (str): Directory zipped Copernicus elevation tiles or url to AWS Copernicus data
                    ex : 'https://copernicus-dem-30m.s3.amazonaws.com/'
    shift (str) : URL of wavelength shift correction surface file
    match (str or list) : Perform landsat image matching (recommended)
    proj (bool) : Project image to UTM grid
    res (int) : Resolution of projected image, 30 should be one of its factors (90,120,150.....)

    '''

    parser = argparse.ArgumentParser(description = "Convert PRISMA HDF to ENVI format")
    parser.add_argument('input',help="Path to compressed input file", type = str)
    parser.add_argument('out_dir',help="Output directory", type = str)
    parser.add_argument('temp_dir',help="Temporary directory", type = str)
    #parser.add_argument('resolution',help="Output resolution", type = int)

    args = parser.parse_args()

    base_name = os.path.basename(args.input)

    if base_name.startswith('PRS'):
        aws_cop_url='https://copernicus-dem-30m.s3.amazonaws.com/'
        shift_surface='https://github.com/EnSpec/sister/raw/master/data/prisma/wavelength_shift/PRISMA_20200721104249_20200721104253_0001_wavelength_shift_surface'
        prisma.he5_to_envi(args.input,args.out_dir,args.temp_dir,
                           aws_cop_url,
                           shift = shift_surface,match=False,
                           proj = True,res = 30)
    elif base_name.startswith('ang') or base_name.startswith('f'):
        aviris.preprocess(args.input,args.out_dir,args.temp_dir,
                          resolution = 30)
    elif base_name.startswith('DESIS'):
        print('DESIS')
    else:
        print('Unrecognized input sensor')


if __name__=='__main__':

    main()
