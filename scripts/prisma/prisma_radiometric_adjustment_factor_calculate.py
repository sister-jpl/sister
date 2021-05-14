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

import os
import hytools as ht
import pandas as pd
import isofit
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import glob

from isofit.core.forward import ForwardModel
from isofit.core.geometry import Geometry
from isofit.inversion.inverse import Inversion
from isofit.configs.configs import Config
from scipy.interpolate import interp1d
from isofit.utils import surface_model
import datetime as dt
import ray
from skimage.filters import threshold_otsu
from scipy.interpolate import interp1d
from scipy.signal import savgol_filter
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import minimize
from hytools.io.envi import WriteENVI,envi_header_dict
from scipy.signal import gaussian
import shutil

home = os.path.expanduser('~')

# Potentential reference scenes
base_names = ['20200715103934_20200715103938_0001',
              '20210209104043_20210209104047_0001',
              '20200721104249_20200721104253_0001']
export_figures = True

for base_name in base_names:

    iso_base = '%s/isotest/isofit/' % home
    figure_dir = '%s/sister/figures/prisma/radiance_factors/' % home
    doy = dt.datetime.strptime(base_name[:8],'%Y%m%d').timetuple().tm_yday
    base_file =  '%s/data/prisma/rad/PRISMA_%s_rad_factor/PRISMA_%s_' % (home, base_name,base_name)
    cal_dir = '%s/temp/%s' % (home,base_name)

    if os.path.isdir(cal_dir):
        shutil.rmtree(cal_dir)

    os.mkdir(cal_dir)
    for subdir in ['/configs','/data','/output','/lut_full/']:
        os.mkdir(cal_dir+subdir)

    #Load datasets
    radiance = ht.HyTools()
    radiance.read_file(base_file + 'rad_unprj', 'envi')
    radiance.create_bad_bands([[300,380.0],[1300.0,1500],[1780.0,1975.0],[2450.0,2600]])

    observables = ht.HyTools()
    observables.read_file(base_file + 'obs_unprj', 'envi')

    location = ht.HyTools()
    location.read_file(base_file + 'loc_unprj', 'envi')

    #Get center of image
    obs_mean  =  observables.get_chunk(0,-1,0,-1).mean(axis =(0,1))
    loc_mean  =  location.get_chunk(0,-1,0,-1).mean(axis =(0,1))

    # Export wavelength file
    wavelength_file = '%s/data/wavelengths.txt' % cal_dir
    wave_fit = pd.DataFrame(index = [x for x in range(radiance.bands)])
    wave_fit[0] =np.round(radiance.wavelengths/1000,5)
    wave_fit[1] = np.round(radiance.fwhm/1000,5)
    wave_fit.to_csv(wavelength_file, sep = '\t',header=False)

    # Create decorrelated surface model
    windows= [[380.0, 1300.0], [1500, 1780.0], [1975.0, 2450.0]]
    surface_data_path = "%s/isotest/data/surface/surface_20210215_emit" % home
    surface_file = "%s/data/surface_emit.mat" % cal_dir
    surface_config = "%s/configs/surface.json" % cal_dir
    surface_config_gen_decorr(surface_data_path,
                       windows,
                       wavelength_file,
                       surface_file,
                       surface_config)
    surface_model(surface_config)

    # Create modtran template file
    engine = 'sRTMnet'
    template_file = '%s/configs/radtran_template.json' % cal_dir
    modtran_template_create(obs_mean.astype(np.float64),
                            loc_mean.astype(np.float64),
                            doy,
                            template_file)

    # Create forward model configs
    rtm_config,surface_config,instrument_config =fw_configs(iso_base,
                                                            cal_dir,
                                                            template_file,
                                                            wavelength_file,
                                                            surface_file,
                                                            1,
                                                            "multicomponent_surface",
                                                            engine,
                                                            observables)

    full_config = Config({"forward_model":{"instrument":instrument_config,
                                           "surface":surface_config,
                                           "radiative_transfer": rtm_config}})
    #Start up ray
    if ray.is_initialized():
        ray.shutdown()
    ray.init(num_cpus = 48)

    # Run intial forward model
    fm = ForwardModel(full_config)

    inversion_settings = {"implementation": {"mode": "inversion",
                                             "inversion": {
                                              "windows": windows}}}
    inverse_config = Config(inversion_settings)
    iv = Inversion(inverse_config, fm)

    # Refine water vapor LUT
    water = []
    window  =5
    for line in range(10,990,100):
        print(line)
        obs_mean  =  observables.get_chunk(2,990,
                                          line-window,line+window).mean(axis =(0,1))
        loc_mean  =  location.get_chunk(2,990,
                                          line-window,line+window).mean(axis =(0,1))
        rad_mean  = radiance.get_chunk(2,990,
                                     line-window,line+window).mean(axis =(0,1))
        geom = Geometry(obs = obs_mean,loc = loc_mean)
        state_trajectory = iv.invert(rad_mean, geom)
        water.append(state_trajectory[-1][-1])

    #Create new water table
    water_min = np.min(water)*.9
    water_max = np.max(water)*1.1
    water_vals = np.linspace(water_min,water_max,5).tolist()

    # Recreate configs with refined water vapor
    for file in glob.glob(cal_dir+"/lut_full/*"):
        os.remove(file)

    rtm_config,surface_config,instrument_config =fw_configs(iso_base,
                                                            cal_dir,
                                                            template_file,
                                                            wavelength_file,
                                                            surface_file,
                                                            1,
                                                            "multicomponent_surface",
                                                            engine,
                                                            observables,
                                                            water_vals)

    full_config = Config({"forward_model":{"instrument":instrument_config,
                                           "surface":surface_config,
                                           "radiative_transfer": rtm_config}})

    fm = ForwardModel(full_config)
    inverse_config = Config(inversion_settings)
    iv = Inversion(inverse_config, fm)

    # Use the middle of the scene for calibration, least impacted by residual
    # wavelength shift
    line  =500

    # Get mean line data
    obs_mean  =  observables.get_chunk(0,-1,
                                      line,line+1).mean(axis =(0,1))
    loc_mean  =  location.get_chunk(0,-1,
                                      line,line+1).mean(axis =(0,1))
    rad_mean  = radiance.get_chunk(0,-1,
                                 line,line+1).mean(axis =(0,1))

    geom = Geometry(obs = obs_mean,loc = loc_mean)
    state_trajectory = iv.invert(rad_mean, geom)
    state_est = state_trajectory[-1]
    rfl_est1 = state_est[:radiance.bands]


    # Generate smooth reflectance spectrum
    rfl_smooth = np.copy(rfl_est1)
    rfl_smooth[:30]= savgol_filter(rfl_smooth[:30],7,3)
    d_start,d_end = 30,94
    rfl_smooth[d_start:d_end] = savgol_filter(rfl_smooth[d_start:d_end],31,2)
    d_start,d_end = 112,140
    rfl_smooth[d_start:d_end] = savgol_filter(rfl_smooth[d_start:d_end],25,1)
    d_start,d_end = 161,223
    rfl_smooth[d_start:d_end] = savgol_filter(rfl_smooth[d_start:d_end],17,2)

    #Foward model to at sensor radiance
    new_state = state_est.copy()
    new_state[:radiance.bands] =rfl_smooth
    smooth_rad = fm.calc_rdn(new_state,geom)

    rfl_est1[radiance.bad_bands] = np.nan
    rfl_smooth[radiance.bad_bands] = np.nan
    smooth_rad[radiance.bad_bands] = np.nan

    #Plot initial inversion
    fig = plt.figure(figsize = (6,4))
    ax = fig.add_subplot(111)
    ax.plot(radiance.wavelengths,rfl_est1,c='k',lw=2)
    ax.set_xlabel("Wavelength (nm)",fontsize = 15)
    ax.set_ylabel("Reflectance",fontsize = 15)
    ax.set_xlim(375,2500)
    for direction in ['left','right','top','bottom']:
        ax.spines[direction].set_linewidth(1.5)

    if export_figures:
        plt.savefig("%s/PRISMA_%s_initial_rfl.png" % (figure_dir,base_name),
                    bbox_inches = 'tight', dpi = 500)
        plt.show()
        plt.close()

    #Plot initial inversion and smooth reflectance
    fig = plt.figure(figsize = (6,4))
    ax = fig.add_subplot(111)
    ax.plot(radiance.wavelengths,rfl_est1,c='k',lw=2,alpha=.5)
    ax.plot(radiance.wavelengths,rfl_smooth,c='r',lw=2)
    ax.set_xlabel("Wavelength (nm)",fontsize = 15)
    ax.set_ylabel("Reflectance",fontsize = 15)
    ax.set_xlim(375,2500)
    for direction in ['left','right','top','bottom']:
        ax.spines[direction].set_linewidth(1.5)
    if export_figures:
        plt.savefig("%s/PRISMA_%s_smoothed_rfl.png" % (figure_dir,base_name),
        plt.show()
                    bbox_inches = 'tight', dpi = 500)
        plt.close()

    #Calculate radiance factors, set water absorption bands to 1
    ratio = smooth_rad/rad_mean
    ratio[radiance.bad_bands] = 1
    rad_mean[radiance.bad_bands] = np.nan

    fig = plt.figure(figsize = (10,4),facecolor='w')
    ax1 = fig.add_subplot(211)
    ax1.plot(radiance.wavelengths, rad_mean,label= 'ISOFIT', c = 'r')
    ax1.plot(radiance.wavelengths, smooth_rad,label= 'ISOFIT smooth', c = 'k')
    ax1.legend(frameon = False)
    ax1.set_xlabel('Wavelength (nm)')
    ax1.set_ylabel('Reflectance')
    ax1.legend(bbox_to_anchor=(.9,-.17), frameon= False, ncol= 2, fontsize = 12)

    ax2 = fig.add_subplot(212)
    ax2.plot(radiance.wavelengths, ratio, c = 'gray')
    ax2.legend(frameon = False)
    ax2.set_xlabel('Wavelength (nm)')
    ax2.set_ylabel('Correction factor')
    ax2.legend(bbox_to_anchor=(.9,-.17), frameon= False, ncol= 2, fontsize = 12)

    if export_figures:
        plt.savefig("%s/PRISMA_%s_smoothed_atsensor_radiance_factors.png" % (figure_dir,base_name),
                    bbox_inches = 'tight', dpi = 500)
        plt.show()
        plt.close()


    rad_corr_file = "%s/sister/data/prisma/radiance_factors/PRISMA_%s_radiance_factors.txt" % (home,base_name)
    rad_corr = pd.DataFrame(index = [x for x in range(radiance.bands)])
    rad_corr[0] = ratio
    rad_corr.to_csv(rad_corr_file , sep = '\t',header=False)







