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
import yaml
from sister import Sister

def main():

    parser = argparse.ArgumentParser(description = "Convert PRISMA HDF to ENVI format")
    parser.add_argument('base_name',help="File base name", type = str, required=True)
    parser.add_argument('config_file',help="Configuration YAML file", type = str, required=True)
    parser.add_argument('-steps',help="Workflow steps", nargs='+', required=True)
    args = parser.parse_args()

    with open(args.config_file,'r') as file:
        configs = yaml.load(file)

    manager = Sister(args.base_name,configs)

    workflow =args.steps

    if 'rdn' in workflow:
        manager.radiance()
    if 'rfl' in workflow:
        manager.reflectance()


if __name__=='__main__':
    main()
