#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SISTER
Space-based Imaging Spectroscopy and Thermal PathfindER
Author: Adam Chlus
"""
import os
import logging
from multiprocessing import cpu_count
import shutil
import subprocess
import hytools as ht
from isofit.utils import surface_model
import numpy as np
from .sensors import prisma # ,desis
from .utils.isofit import gen_wavelength_file,surface_config_gen,get_surface_spectra
from .utils.misc import download_file
import yaml

class Sister:
    def __init__(self,base_name,configs):

        self.base_name = base_name
        self.project = True

        if base_name.startswith('PRS'):
            self.sensor = 'prisma'
            self.date = base_name[4:12]
        elif base_name.startswith('DESIS'):
            self.sensor = 'desis'
            self.date = base_name[4:12]

        self.rdn_cfg = configs[self.sensor]
        self.rdn_dir = configs[self.sensor]['rdn']
        self.rfl_dir = configs[self.sensor]['rfl']
        self.raw_dir = configs[self.sensor]['raw']
        self.tmp_dir = configs['tmp_dir']
        self.dsm_dir = configs['dsm_dir']

        if self.project:
            suffix = '_prj'
        else:
            suffix = ''

        #Assign file path names
        self.rdn_file = self.rdn_dir + self.base_name + '/' + self.base_name + '_rdn' + suffix
        self.loc_file = self.rdn_dir + self.base_name + '/' + self.base_name + '_loc' + suffix
        self.obs_file = self.rdn_dir + self.base_name + '/' + self.base_name + '_obs' + suffix
        self.rfl_file = self.rdn_file.replace('rdn','rfl')
        self.unc_file = self.rfl_file.replace('_rfl','_unc')
        self.isofit = configs['isofit']

    def exists(self,file):
        return os.path.isfile(file)

    def radiance(self):
        #First convert to radiance
        if self.sensor == 'prisma':
            l1_zip = '%s%s.zip' % (self.raw_dir,self.base_name.replace('PRS','PRS_L1_STD_OFFL'))
            prisma.he5_to_envi(l1_zip,self.rdn_dir,
                               self.tmp_dir,self.dsm_dir,
                               shift = self.rdn_cfg['shift'],
                               match = self.rdn_cfg['match'],
                               proj = self.rdn_cfg['project'],
                               res = self.rdn_cfg['resolution'])
        elif self.sensor == 'desis':
            #l1_zip = '%s%s.zip' % (self.raw_dir,self.base_name.replace('PRS','PRS_L1_STD_OFFL'))
            desis.geotiff_to_envi(self.base_name,l1b_zip,self.rdn_dir,
                                 self.tmp_dir,self.dsm_dir,
                                 match = self.rdn_cfg['match'],
                                 proj = self.rdn_cfg['project'],
                                 res = self.rdn_cfg['resolution'])

    def reflectance(self):
        '''Run ISOFIT reflectance inversion
        '''

        # Check for input files
        if not (self.exists(self.rdn_file) & self.exists(self.loc_file) & self.exists(self.obs_file)):
            print('Radiance data not found, generating radiance datasets')
            self.radiance()

        rfl_temp = self.tmp_dir + '%s_isofit/' % self.base_name
        if not os.path.isdir(rfl_temp):
            os.mkdir(rfl_temp)
        else:
            shutil.rmtree(rfl_temp)
            os.mkdir(rfl_temp)

        for key in self.isofit['surface']:
            get_surface_spectra(rfl_temp + key,
                                self.isofit['surface'][key])

        surface_config = "%s/surface_config.json" % rfl_temp
        surface_file = '%s/surface_filtered.mat'% rfl_temp
        wavelength_file = gen_wavelength_file(self.rdn_file,
                                               rfl_temp)
        surface_config_gen(rfl_temp,
                           self.isofit['surface_type'],
                           self.isofit['windows'],
                           wavelength_file,
                           surface_file,
                           surface_config)

        surface_model(surface_config)

        apply_oe  = ['python']
        apply_oe.append('%s/isofit/utils/apply_oe.py' % self.isofit['base'])
        apply_oe.append(self.rdn_file)
        apply_oe.append(self.loc_file)
        apply_oe.append(self.obs_file)
        apply_oe.append(rfl_temp)

        if self.sensor in ['prisma','desis']:
            apply_oe.append('NA-%s' % self.date)
        apply_oe += ['--surface_path', surface_file]
        apply_oe += ['--n_cores',str(cpu_count())]
        apply_oe += ['--empirical_line','1']
        apply_oe += ['--presolve','1']
        apply_oe += ['--segment_size', str(self.isofit['segment_size'])]
        apply_oe += ['--ray_temp_dir',rfl_temp]
        apply_oe += ['--log_file','%s/%s_logfile' % (rfl_temp,self.base_name)]
        apply_oe += ['--emulator_base',self.isofit['emulator']]

        if self.isofit['radiance_factors']:
            rad_factors = "%s/radiance_factors.txt" % rfl_temp
            download_file(rad_factors,self.isofit['radiance_factors'])
            apply_oe += ['--rdn_factors_path',rad_factors]

        apply = subprocess.Popen(apply_oe)
        apply.wait()

        if not os.path.isdir("%s/%s" % (self.rfl_dir,self.base_name)):
            os.mkdir("%s/%s" % (self.rfl_dir,self.base_name))

        # Mask windows and rename ISOFIT output files
        rfl = '%s/output/%s_rfl_prj_rfl'% (rfl_temp,self.base_name)
        uncert = '%s/output/%s_uncert_prj_uncert'% (rfl_temp,self.base_name)

        for new,old in [(self.rfl_file,rfl),(self.unc_file,uncert)]:
            hy_obj = ht.HyTools()
            hy_obj.read_file(old,'envi')
            mask = []
            for wave in hy_obj.wavelengths:
                window = False
                for start,end in test.isofit['windows']:
                    if (wave>start) and (wave<end):
                        window |= True
                mask.append(window)
            header_dict = hy_obj.get_header()
            writer = ht.io.WriteENVI(new, header_dict)
            for band_num in range(hy_obj.bands):
                band  = hy_obj.get_band(band_num)
                if not mask[band_num]:
                    band = np.zeros(band.shape)
                writer.write_band(band,band_num)

        shutil.rmtree(rfl_tmp)






