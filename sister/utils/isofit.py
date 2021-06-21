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
import glob
import json
import os
import shutil
import subprocess
import hytools as ht
import isofit
from isofit.utils import surface_model
import pandas as pd
import numpy as np

home = os.path.expanduser("~")
isofit_surface_dir =  '%s/isotest/data/neon/surface/surface_20210215_emit/' % home
isofit_surface_dir =  '%s/Dropbox/rs/sister/data/isofit/surface_models/surface_20210215_emit/' % home
windows = [[380.0, 1300.0], [1450, 1780.0], [1950.0, 2450.0]]

def surface_config_gen(surface_data_path,surface_type,windows,wavelength_file,surface_file,out_config):
    ''' Export surface model config file
    '''
    surface_config = {}
    surface_config["output_model_file"]= surface_file
    surface_config["wavelength_file"]= wavelength_file
    surface_config["normalize"]="Euclidean"
    surface_config["reference_windows"]= windows
    if surface_type == 'terrestrial':

        surface_config["sources"] = [[] for x in range(4)]
        surface_config["sources"][0] = {}
        surface_config["sources"][0]["input_spectrum_files"] = ["%s/filtered_other" % surface_data_path]
        surface_config["sources"][0]["n_components"] = 1
        surface_config["sources"][0]["windows"] = [
                                                    {
                                                      "interval":[300,420],
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
                                                    }
                                                  ]

        surface_config["sources"][1] = {}
        surface_config["sources"][1]["input_spectrum_files"] = ["%s/filtered_veg" % surface_data_path]
        surface_config["sources"][1]["n_components"] = 1
        surface_config["sources"][1]["windows"]= [
                                                  {
                                                    "interval":[300,420],
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
                                                  }
                                                ]


        surface_config["sources"][2] = {}
        surface_config["sources"][2]["input_spectrum_files"] = ["%s/filtered_ocean" % surface_data_path]
        surface_config["sources"][2]["n_components"] = 1
        surface_config["sources"][2]["windows"] = [
                                                  {
                                                    "interval":[300,420],
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
                                                  }
                                                ]
        surface_config["sources"][3] = {}
        surface_config["sources"][3]["input_spectrum_files"] = ["%s/surface_Liquids" % surface_data_path]
        surface_config["sources"][3]["n_components"] = 2
        surface_config["sources"][3]["windows"] =  [
                                                      {
                                                        "interval":[300,420],
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
                                                      }
                                                    ]

    elif surface_type == 'decorr':

        surface_config["sources"][0] = {}
        surface_config["sources"][0]["input_spectrum_files"] = ["%s/filtered_other" % surface_data_path]
        surface_config["sources"][0]["n_components"] = 1
        surface_config["sources"][0]["windows"] = [
                                                    {
                                                      "interval":[300,2500],
                                                      "regularizer":10,
                                                      "correlation":"decorrelated"
                                                    }
                                                  ]

        surface_config["sources"][1] = {}
        surface_config["sources"][1]["input_spectrum_files"] = ["%s/filtered_veg" % surface_data_path]
        surface_config["sources"][1]["n_components"] = 1
        surface_config["sources"][1]["windows"]= [
                                                     {
                                                      "interval":[300,2500],
                                                      "regularizer":10,
                                                      "correlation":"decorrelated"
                                                    }
                                                ]


        surface_config["sources"][2] = {}
        surface_config["sources"][2]["input_spectrum_files"] = ["%s/filtered_ocean" % surface_data_path]
        surface_config["sources"][2]["n_components"] = 1
        surface_config["sources"][2]["windows"] = [
                                                    {
                                                      "interval":[300,2500],
                                                      "regularizer":10,
                                                      "correlation":"decorrelated"
                                                    }
                                                ]
        surface_config["sources"][3] = {}
        surface_config["sources"][3]["input_spectrum_files"] = ["%s/surface_Liquids" % surface_data_path]
        surface_config["sources"][3]["n_components"] = 2
        surface_config["sources"][3]["windows"] =  [
                                                    {
                                                      "interval":[300,2500],
                                                      "regularizer":10,
                                                      "correlation":"decorrelated"
                                                    }
                                                    ]

    with open(out_config, 'w') as outfile:
        json.dump(surface_config,outfile,indent=3)


def gen_wavelength_file(rad_file,out_dir):
    '''
    Export radiance image wavelengths and FWHM to file
    '''

    radiance = ht.HyTools()
    radiance.read_file(rad_file,'envi')
    wave_fit = pd.DataFrame(index = [x for x in range(radiance.bands)])
    wave_fit[0] =np.round(radiance.wavelengths/1000,5)
    wave_fit[1] = np.round(radiance.fwhm/1000,5)
    wavelength_file = '%s/wavelength_fit.txt' % out_dir
    wave_fit.to_csv(wavelength_file , sep = '\t',header=False)
    return wavelength_file


def run_isofit(in_dir,out_dir,temp_dir,segment_size=15,surface_type = 'terrestrial'):

    # Testing
    # in_dir =  '/data2/prisma/rad/PRS_20210616165230_20210616165235_0001_flat/'
    # out_dir = '/data2/prisma/rfl/PRS_20210616165230_20210616165235_0001_flat/'
    # temp_dir = '/data2/temp/PRS_20210616165230_20210616165235_0001_flat/'

    #Don't like how this is implemented...works for now
    rad_files = glob.glob("%s/*_rad_*" % in_dir)
    rad_files.sort()
    obs_files = glob.glob("%s/*_obs_*" % in_dir)
    obs_files.sort()
    loc_files = glob.glob("%s/*_loc_*" % in_dir)
    loc_files.sort()

    base_name = os.path.basename(rad_files[0])

    if base_name.startswith('PRS'):
        date  = base_name.split('_')[1][:8]
        base_name = '_'.join(base_name.split('_')[:4])
    elif base_name.startswith('DESIS'):
        date  = base_name.split('_')[1][:8]
        base_name = '_'.join(base_name.split('_')[:4])

    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)
    else:
        shutil.rmtree(out_dir)
        os.mkdir(out_dir)

    if not os.path.isdir(temp_dir):
        os.mkdir(temp_dir)
    else:
        shutil.rmtree(temp_dir)
        os.mkdir(temp_dir)

    surface_config = "%s/surface_config.json" % temp_dir
    surface_file = '%s/surface_filtered.mat'% temp_dir
    wavelength_file  = gen_wavelength_file(rad_files[0],temp_dir)

    surface_config_gen(isofit_surface_dir,
                       surface_type,
                       windows,
                       wavelength_file,
                       surface_file,
                       surface_config)

    surface_model(surface_config)

    apply_oe  = ['python']
    apply_oe.append('%s/Dropbox/rs/sister/algorithms/atc/isofit/isofit/utils//apply_oe.py' % home)
    apply_oe.append(rad_files[0])
    apply_oe.append(loc_files[0])
    apply_oe.append(obs_files[0])
    apply_oe.append(temp_dir)
    apply_oe.append('NA-%s' % date)
    apply_oe += ['--surface_path', surface_file]
    apply_oe += ['--n_cores','8']
    apply_oe += ['--empirical_line','1']
    apply_oe += ['--presolve','1']
    apply_oe += ['--segment_size',str(segment_size)]
    apply_oe += ['--ray_temp_dir',temp_dir]
    apply_oe += ['--log_file','%s/%s_logfile' % (out_dir,base_name)]
    apply_oe += ['--rdn_factors_path',"%s/Dropbox/rs/sister/sister/data/prisma/radiance_factors/PRISMA_20200721104249_20200721104253_0001_radiance_factors.txt" % (home)]
    apply_oe += ['--emulator_base','/data1/isofit/sRTMnet/sRTMnet_v100']

    apply = subprocess.Popen(apply_oe)
    apply.wait()

    for suffix in ['rfl','rfl.hdr','uncert','uncert.hdr']:
        shutil.move('%s/output/%s_rad_geo_%s'%(temp_dir,base_name,suffix),
                    '%s/%s_rad_geo_%s'%(out_dir,base_name,suffix))

    shutil.rmtree(temp_dir)






