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
import sister
from sister.sensors import prisma,aviris,desis

def main():
    parser = argparse.ArgumentParser(description = "Convert PRISMA HDF to ENVI format")
    parser.add_argument('input',help="Path to compressed input file", type = str)
    parser.add_argument('out_dir',help="Output directory", type = str)
    parser.add_argument('temp_dir',help="Temporary directory", type = str)
    parser.add_argument('resolution',help="Output resample resolution",type=int, default = 0)
    parser.add_argument('smile', nargs='?',help="Path to smile wavelengths", default = False)
    parser.add_argument('rad_coeff', nargs='?',help="Path to radiometric coeffs",default = False)
    parser.add_argument('landsat', nargs='?',help="Landsat reference file",default = False)

    args = parser.parse_args()
    base_name = os.path.basename(args.input)
    aws_cop_url='https://copernicus-dem-30m.s3.amazonaws.com/'

    if base_name.startswith('PRS'):
        prisma.he5_to_envi(args.input,args.out_dir,args.temp_dir,
                           aws_cop_url,
                           shift = args.smile,
                           rad_coeff =args.rad_coeff,
                           proj = True,
                           match=args.landsat)
    elif base_name.startswith('ang') or base_name.startswith('f'):
        aviris.preprocess(args.input,args.out_dir,args.temp_dir,
                          res = args.resolution)
    elif base_name.startswith('DESIS'):
        desis.l1c_process(args.input,args.out_dir,args.temp_dir,
                           aws_cop_url)
    else:
        print('Unrecognized input sensor')


if __name__=='__main__':

    main()
