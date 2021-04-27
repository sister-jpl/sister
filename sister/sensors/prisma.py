#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SISTER
Space-based Imaging Spectroscopy and Thermal PathfindER
Author: Adam Chlus
"""

import datetime as dt
import os
import zipfile
import h5py
import hytools as ht
from hytools.io.envi import WriteENVI,envi_header_dict
from hytools.topo.topo import calc_cosine_i
import numpy as np
from rtree import index
from scipy.interpolate import interp1d
from scipy.spatial import cKDTree
import pyproj
from pysolar import solar
from skimage.util import view_as_blocks
from scipy.ndimage import uniform_filter
from ..utils.terrain import *
from ..utils.geometry import *
from ..utils.ancillary import *


home = os.path.expanduser("~")

def he5_to_envi(base_name,l1_zip,l2c_zip,out_dir,elev_dir,temp_dir,
                smile = None,match=None,rfl = False, project = True, res = 30):
    '''
    This function exports three unprojected files:
        *_rad_unpj : Merged and optionaly smile corrected radiance cube
        *_obs_unprj : Observables file in the format of JPL obs files:
                1. Pathlength (m)
                2. Sensor view azimuth angle (degrees)
                3. Sensor view zenith angle (degrees)
                4. Solar azimuth angle (degrees)
                5. Solar zenith angle (degrees)
                6. Sensor view azimuth angle in degrees
                7. Sensor view azimuth angle in degrees
        *_loc_unprj : Location file in the following format:
                1. Longitude (decimal degrees)
                2. Longitude (decimal degrees)
                3. Elevation (m)

    l1(str): L1 zipped data product path
    l2(str): L2C zipped data product path
    base_name :
    out_dir(str): Output directory of ENVI datasets
    elev_dir (str): Directory zipped elevation tiles
    smile (str) : Pathname of smile correction surface file
    match (str or list) : Pathname to Landsat image(s) for image re-registration (recommended)
    rfl (bool) : Export ASI L2C surface reflectance
    '''

    # home = os.path.expanduser("~")
    # elev_dir = '/data2/cop_dsm/'
    # base_name = "20200703022145_20200703022149_0001"
    # l1_zip  = '/data2/prisma/zip/PRS_L1_STD_OFFL_%s.zip'% base_name
    # l2c_zip  = '/data2/prisma/zip/PRS_L2C_STD_%s.zip' % base_name
    # out_dir = '/data1/temp//PRISMA/'
    # temp_dir =  '/data1/temp//temp_PRISMA/'
    # smile = '%s/Dropbox/rs/sister/data/prisma/PRS_20200721104249_20200721104253_0001_smile' % home
    # match = '/data2/landsat/LC08_L2SP_115029_20190708_20200827_02_T1_SR_B5.TIF'
    # rfl = False
    # project = True
    # res = 30

    for file in [l1_zip,l2c_zip]:
        zip_base  =os.path.basename(file)
        print('Unzipping %s' % zip_base)
        with zipfile.ZipFile(file,'r') as zipped:
            zipped.extractall(temp_dir)

    l1  = '%sPRS_L1_STD_OFFL_%s.he5' % (temp_dir,base_name)
    l2c  = '%sPRS_L2C_STD_%s.he5' % (temp_dir,base_name)

    file_suffixes =  ['rad','loc','obs']

    if smile:
        smile_obj = ht.HyTools()
        smile_obj.read_file(smile, 'envi')
        shift_surf_smooth = smile_obj.get_band(0)
        smile_correct = True

    for product in [l1,l2c]:
        l_obj = h5py.File(product,'r')
        subdir = [x for x in l_obj['HDFEOS']["SWATHS"].keys() if 'HCO' in x][0]

        if 'L2C' in subdir:
            if not rfl:
                continue
            print('Exporting reflectance data')
            measurement = 'rfl'
            file_suffixes.append(measurement)

        else:
            measurement = 'rad'
            print('Exporting radiance data')

        # Export VNIR to temporary ENVI
        vnir_data =  l_obj['HDFEOS']["SWATHS"][subdir]['Data Fields']['VNIR_Cube']
        vnir_waves = l_obj.attrs.get('List_Cw_Vnir')
        vnir_fwhm = l_obj.attrs.get('List_Fwhm_Vnir')

        rad_dict = envi_header_dict ()
        rad_dict['lines']= vnir_data.shape[0]
        rad_dict['samples']= vnir_data.shape[2]
        rad_dict['bands']=  vnir_data.shape[1]
        rad_dict['wavelength']= vnir_waves
        rad_dict['fwhm']= vnir_fwhm
        rad_dict['interleave']= 'bsq'
        rad_dict['data type'] = 12
        rad_dict['wavelength units'] = "nanometers"
        rad_dict['byte order'] = 0
        vnir_temp = '%sPRISMA_-%s_%s_vnir' % (temp_dir,base_name,measurement)

        writer = WriteENVI(vnir_temp,rad_dict )
        writer.write_chunk(np.moveaxis(vnir_data[:,:,:],1,2), 0,0)

        # Export SWIR to temporary ENVI
        swir_data =  l_obj['HDFEOS']["SWATHS"][subdir]['Data Fields']['SWIR_Cube']
        swir_waves = l_obj.attrs.get('List_Cw_Swir')
        swir_fwhm = l_obj.attrs.get('List_Fwhm_Swir')

        rad_dict = envi_header_dict ()
        rad_dict['lines']= swir_data.shape[0]
        rad_dict['samples']= swir_data.shape[2]
        rad_dict['bands']=  swir_data.shape[1]
        rad_dict['wavelength']= swir_waves
        rad_dict['fwhm']= swir_fwhm
        rad_dict['interleave']= 'bil'
        rad_dict['data type'] = 12
        rad_dict['wavelength units'] = "nanometers"
        rad_dict['byte order'] = 0
        swir_temp = '%sPRISMA_%s_%s_swir' % (temp_dir,base_name,measurement)

        writer = WriteENVI(swir_temp,rad_dict )
        writer.write_chunk(np.moveaxis(swir_data[:,:,:],1,2), 0,0)

        vnir_waves = np.flip(vnir_waves[6:])
        swir_waves = np.flip(swir_waves[:-3])

        vnir_fwhm = np.flip(vnir_fwhm[6:])
        swir_fwhm = np.flip(swir_fwhm[:-3])

        vnir_obj = ht.HyTools()
        vnir_obj.read_file(vnir_temp, 'envi')

        swir_obj = ht.HyTools()
        swir_obj.read_file(swir_temp, 'envi')

        if project:
            output_name = '%sPRISMA_%s_%s_unprj' % (temp_dir,base_name,measurement)
        else:
            output_name = '%sPRISMA_%s_%s_unprj' % (out_dir,base_name,measurement)

        rad_dict  = envi_header_dict()
        rad_dict ['lines']= vnir_obj.lines-4 #Clip edges of array
        rad_dict ['samples']=vnir_obj.columns-4  #Clip edges of array
        rad_dict ['bands']= len(vnir_waves.tolist() + swir_waves.tolist())
        rad_dict ['wavelength']= vnir_waves.tolist() + swir_waves.tolist()
        rad_dict ['fwhm']= vnir_fwhm.tolist() + swir_fwhm.tolist()
        rad_dict ['interleave']= 'bil'
        rad_dict ['data type'] = 4
        rad_dict ['wavelength units'] = "nanometers"
        rad_dict ['byte order'] = 0

        writer = WriteENVI(output_name,rad_dict)
        iterator_v =vnir_obj.iterate(by = 'line')
        iterator_s =swir_obj.iterate(by = 'line')

        while not iterator_v.complete:
            chunk_v = iterator_v.read_next()[:,6:]
            chunk_v =np.flip(chunk_v,axis=1)
            chunk_s = iterator_s.read_next()[:,:-3]
            chunk_s =np.flip(chunk_s,axis=1)

            if (iterator_v.current_line >=2) and (iterator_v.current_line <= 997):
                if (measurement == 'rad') & smile_correct:
                    vnir_interpolator = interp1d(shift_surf_smooth[iterator_v.current_line,:60],
                                                   chunk_v,fill_value = "extrapolate",kind='cubic')
                    chunk_v = vnir_interpolator(vnir_waves)
                    swir_interpolator = interp1d(shift_surf_smooth[iterator_v.current_line,60:],
                                                   chunk_s,fill_value = "extrapolate",kind='cubic')
                    chunk_s = swir_interpolator(swir_waves)

                line = np.concatenate([chunk_v,chunk_s],axis=1)/1000.
                writer.write_line(line[2:-2,:], iterator_v.current_line-2)

    #Load ancillary datasets
    geo =  l_obj['HDFEOS']["SWATHS"][subdir]['Geolocation Fields']
    geom = l_obj['HDFEOS']['SWATHS'][subdir]['Geometric Fields']
    pvs =  l_obj['Info']["Ancillary"]['PVSdata']

    # Time
    '''1. Convert from MJD2000 to UTC hours
       2. Fit line to estimate continous time.
    '''

    def dhour(day):
        epoch = dt.datetime(2000,1, 1,)
        epoch = epoch.replace(tzinfo=dt.timezone.utc)

        hour =  (day-day//1)*24
        minute =  (hour-hour//1)*60
        second= (minute-minute//1)*60
        microsecond= (second-second//1)*1000000
        time = epoch + dt.timedelta(days=day//1,hours=hour//1,
                                    minutes=minute//1,seconds=second,
                                    microseconds =microsecond)
        return time.hour + time.minute/60. + time.second/3600.

    v_dhour = np.vectorize(dhour)
    utc_time = v_dhour(np.array(geo['Time'][:]))

    utc_time = np.array(utc_time)
    X = np.concatenate([np.arange(1000)[:,np.newaxis], np.ones(utc_time.shape)],axis=1)
    slope, intercept = np.linalg.lstsq(X,utc_time,rcond=-1)[0].flatten()
    utc_t_linear = slope*np.arange(1000)+ intercept
    utc_time = np.ones(geo['Longitude'][:,:].shape[0]) *utc_t_linear[:,np.newaxis]
    utc_time = utc_time[2:-2,2:-2]

    # Solar geometries
    '''Solar geometry is calculated based on the mean scene acquisition time
    which varies by less than 5 seconds from start to end of the scene and is
    computationally more efficient.
    '''
    epoch = dt.datetime(2000,1, 1,)
    epoch = epoch.replace(tzinfo=dt.timezone.utc)
    time = epoch + dt.timedelta(days=np.array(geo['Time'][:]).mean())

    solar_az = solar.get_azimuth(geo['Latitude'][:,:],geo['Longitude'][:,:],time)[2:-2,2:-2]
    solar_zn = 90-solar.get_altitude(geo['Latitude'][:,:],geo['Longitude'][:,:],time)[2:-2,2:-2]
    sensor_az =geom['Rel_Azimuth_Angle'][2:-2,2:-2] +solar_az # Not exactly right...need to fix minor...
    sensor_zn = geom['Observing_Angle'][2:-2,2:-2]

    longitude= geo['Longitude'][2:-2,2:-2]
    latitude= geo['Latitude'][2:-2,2:-2]

    #Create initial elevation raster
    elevation= dem_generate(longitude,latitude,elev_dir,temp_dir)

    zone,direction = utm_zone(longitude,latitude)

    if match:
        easting,northing =  dd2utm(longitude,latitude)
        coords =np.concatenate([np.expand_dims(easting.flatten(),axis=1),
                                np.expand_dims(northing.flatten(),axis=1)],axis=1)

        project = Projector()
        project.create_tree(coords,easting.shape)
        project.query_tree(easting.min()-100,northing.max()+100,30)

        sensor_az_prj = project.project_band(sensor_az,-9999)
        sensor_zn_prj = project.project_band(sensor_zn,-9999)
        elevation_prj = project.project_band(elevation.astype(np.float),-9999)

        rad_file = '%sPRISMA_%s_rad_unprj' % (temp_dir,base_name)
        radiance = ht.HyTools()
        radiance.read_file(rad_file, 'envi')

        #Average over Landsat 8 Band 5 bandwidth
        warp_band = np.zeros(longitude.shape)
        for wave in range(850,890,10):
            warp_band += radiance.get_wave(wave)/7.
        warp_band = project.project_band(warp_band,-9999)
        warp_band = 16000*(warp_band-warp_band.min())/warp_band.max()

        #Calculate optimal shift
        y_model,x_model = image_match(match,warp_band,
                                      easting.min()-100,northing.max()+100,
                                      sensor_zn_prj,sensor_az_prj,elevation_prj)

        #Apply uniform filter
        smooth_elevation = uniform_filter(elevation,25)
        smooth_az = uniform_filter(sensor_az,25)
        smooth_zn = uniform_filter(sensor_zn,25)

        i,a,b,c = y_model
        y_offset = i + a*smooth_zn +b*smooth_az + c*smooth_elevation

        i,a,b,c= x_model
        x_offset = i + a*smooth_zn +b*smooth_az + c*smooth_elevation

        new_easting = easting+  30*x_offset
        new_northing = northing- 30*y_offset

        zone,direction = utm_zone(longitude,latitude)
        longitude,latitude = utm2dd(new_easting,new_northing,zone,direction)

        #Recalculate elevation with new coordinates
        elevation= dem_generate(longitude,latitude,elev_dir,temp_dir)

    if project:
        loc_file = '%sPRISMA_%s_loc_unprj' % (temp_dir,base_name)
    else:
        loc_file = '%sPRISMA_%s_loc_unprj' % (out_dir,base_name)

    loc_export(loc_file,longitude,latitude,elevation)

    # Path length
    ''' 1. Convert Lat, Lon and Altitude images to X,Y,Z (ECEF) coordinates
        2. Calculate satellite position, original data are sampled at 1Hz
        resulting in steps in the position data, a line is
        fit to each dimension to estimate continuous position.
        3. Calculate satellite-to-pixel distance.
    '''

    ecef = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
    lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')
    x, y, z = pyproj.transform(lla, ecef,
                               longitude,
                               latitude,
                               elevation, radians=False)

    gps_seconds = []
    for second,week in  zip(pvs['GPS_Time_of_Last_Position'][:].flatten(),pvs['Week_Number'][:].flatten()):
        gps_second = (week*7*24*60*60) + second
        time = dt.datetime(1980, 1, 6) + dt.timedelta(seconds=gps_second - (37 - 19))
        gps_seconds.append(time.hour*3600 + time.minute*60. + time.second)

    line_seconds = []
    for day in geo['Time'][:].flatten():
        epoch = dt.datetime(2000,1, 1,)
        epoch = epoch.replace(tzinfo=dt.timezone.utc)
        time = epoch + dt.timedelta(days=day)
        line_seconds.append(time.hour*3600 + time.minute*60. + time.second)

    measurements = np.arange(len(pvs['Wgs84_pos_x'][:]))

    sat_t = np.array(gps_seconds)[:,np.newaxis]
    X = np.concatenate([measurements[:,np.newaxis], np.ones(sat_t.shape)],axis=1)
    slope, intercept = np.linalg.lstsq(X,sat_t,rcond=-1)[0].flatten()
    sat_t_linear = slope*measurements+ intercept

    line_t = np.array(line_seconds)[:,np.newaxis]
    X = np.concatenate([np.arange(1000)[:,np.newaxis], np.ones(line_t.shape)],axis=1)
    slope, intercept = np.linalg.lstsq(X,line_t,rcond=-1)[0].flatten()
    line_t_linear = slope*np.arange(1000)+ intercept

    pathlength = np.zeros(x.shape)
    satelite_xyz = []

    for sat_pos,grd_pos in zip(['x','y','z'],[x,y,z]):
        sat_p = np.array(pvs['Wgs84_pos_%s' % sat_pos][:])
        X = np.concatenate([measurements[:,np.newaxis], np.ones(sat_p.shape)],axis=1)
        slope, intercept = np.linalg.lstsq(X,sat_p,rcond=-1)[0].flatten()
        sat_p_linear = slope*measurements+ intercept
        interpolator = interp1d(sat_t_linear,sat_p_linear,
                                fill_value="extrapolate",kind = 'linear')
        sat_interp = interpolator(line_t_linear)
        pathlength += (sat_interp[2:-2]-grd_pos)**2
        satelite_xyz.append(sat_interp)
    pathlength = np.sqrt(pathlength)

    slope,aspect = slope_aspect(elevation,temp_dir)

    cosine_i = calc_cosine_i(np.radians(solar_zn),
                             np.radians(solar_az),
                             np.radians(slope),
                             np.radians(aspect))

    phase =  np.arccos(np.cos(np.radians(solar_zn)))*np.cos(np.radians(solar_zn))
    phase += np.sin(np.radians(solar_zn))*np.sin(np.radians(solar_zn))*np.cos(np.radians(geom['Rel_Azimuth_Angle'][2:-2,2:-2]))

    # Export observable datacube
    if project:
        obs_file = '%sPRISMA_%s_obs_unprj' % (temp_dir,base_name)
    else:
        obs_file = '%sPRISMA_%s_obs_unprj' % (out_dir,base_name)

    obs_export(obs_file,pathlength,sensor_az,sensor_zn,
               solar_az,solar_zn,phase,slope,aspect,
               cosine_i,utc_time)

    if project:

        #Create new projector with corrected coordinates
        new_coords =np.concatenate([np.expand_dims(new_easting.flatten(),axis=1),
                        np.expand_dims(new_northing.flatten(),axis=1)],axis=1)

        project = Projector()
        project.create_tree(new_coords,new_easting.shape)
        project.query_tree(new_easting.min()-100,new_northing.max()+100,30)

        blocksize = int(res/30)
        map_info = ['UTM', 1, 1, new_easting.min()-100, new_northing.max()+100,res,
                           res,zone,direction, 'WGS-84' , 'units=Meters']
        out_cols = int(blocksize* (project.output_shape[1]//blocksize))
        out_lines = int(blocksize* (project.output_shape[0]//blocksize))

        print('Georeferencing datasets')
        for file in file_suffixes:
            print(file)
            input_name = '%sPRISMA_%s_%s_unprj' % (temp_dir,base_name,file)
            hy_obj = ht.HyTools()
            hy_obj.read_file(input_name, 'envi')
            iterator =hy_obj.iterate(by = 'band')

            out_header = hy_obj.get_header()
            out_header['lines']= project.output_shape[0]//blocksize
            out_header['samples']=project.output_shape[1]//blocksize
            out_header['data ignore value'] = -9999
            out_header['map info'] = map_info

            output_name = '%sPRISMA_%s_%s_geo' % (out_dir,base_name,file)
            writer = WriteENVI(output_name,out_header)

            while not iterator.complete:
                band = project.project_band(iterator.read_next(),-9999)
                band[band == -9999] = np.nan
                band = np.nanmean(view_as_blocks(band[:out_lines,:out_cols], (blocksize,blocksize)),axis=(2,3))
                band[band<0] = 0
                band[np.isnan(band)] = -9999
                writer.write_band(band,iterator.current_band)












