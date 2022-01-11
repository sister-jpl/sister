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
import json
import os
import hytools as ht
import pandas as pd
import numpy as np

home = os.path.expanduser("~")

name,url = 'filtered_other','https://data.ecosis.org/dataset/dea65562-994e-47d8-be37-2a8a3aaffdc9/resource/eb2b6493-7bde-42c7-9fa4-75ba036c671f/download/filtered_other.csv'

surface_params = {'moderate:':   [{ "interval":[300,420],
                              "regularizer":10,
                              "correlation":"decorrelated"
                            },
                            {
                              "interval":[420,785],
                              "regularizer":10,
                              "correlation":"decorrelated"
                            },
                            {
                              "interval":[785,1050],
                              "regularizer":1e-6,
                              "correlation":"EM"
                            },
                            {
                              "interval":[1050,1050],
                              "regularizer":10,
                              "correlation":"decorrelated"
                            },
                            {
                              "interval":[1050,1250],
                              "regularizer":1e-6,
                              "correlation":"EM"
                            },
                            {
                              "interval":[1250,2500],
                              "regularizer":10,
                              "correlation":"decorrelated"
                            }],

            'weak':    [ {"interval":[300,2500],
                          "regularizer":10,
                          "correlation":"decorrelated"
                        }],

            'strong':  [ {"interval":[300,2500],
                            "regularizer":1e-6,
                            "correlation":"EM"
                          }]}


def get_surface_spectra(output_file,url):
    '''Download surface spectra from ECOSIS and save
    as ENVI file
    '''
    spectra = pd.read_csv(url)

    header= ht.io.envi_header_dict()
    header['samples'] =1
    header['bands'] = len(spectra.columns)
    header['lines'] =len(spectra.index)
    header['interleave'] = 'bip'
    header['data type'] = 4
    header['byte order'] = 0
    header['wavelength'] = spectra.columns.values.tolist()
    writer = ht.io.WriteENVI(output_file, header)

    line= 0
    for r,row in spectra.iterrows():
        writer.write_line(row.values,line)
        line+=1

def gen_wavelength_file(rad_file,out_dir):
    '''Export radiance image wavelengths and FWHM to file
    '''

    radiance = ht.HyTools()
    radiance.read_file(rad_file,'envi')
    wave_fit = pd.DataFrame(index = [x for x in range(radiance.bands)])
    wave_fit[0] =np.round(radiance.wavelengths/1000,5)
    wave_fit[1] = np.round(radiance.fwhm/1000,5)
    wavelength_file = '%s/wavelength_fit.txt' % out_dir
    wave_fit.to_csv(wavelength_file , sep = '\t',header=False)
    return wavelength_file


def surface_config_gen(surface_data_path,windows,wavelength_file,surface_file,out_config,prior = 'weak'):

    surfaces= ['filtered_other','filtered_veg','filtered_ocean','surface_Liquids']

    surface_config = {}
    surface_config["output_model_file"]= surface_file
    surface_config["wavelength_file"]= wavelength_file
    surface_config["normalize"]="Euclidean"
    surface_config["reference_windows"]= windows
    surface_config["sources"] = []

    for s,surface in enumerate(surfaces):

        surface_config["sources"][s] = {}
        surface_config["sources"][s]["input_spectrum_files"] = ["%s/%s" % (surface_data_path,surface)]
        surface_config["sources"][s]["n_components"] = 1
        surface_config["sources"][s]["windows"] = surface_params[prior]

    with open(out_config, 'w') as outfile:
        json.dump(surface_config,outfile,indent=3)
