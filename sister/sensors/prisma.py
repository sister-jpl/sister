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

import datetime as dt
import logging
import os
import zipfile
import shutil
import h5py
import hytools as ht
import pandas as pd
from hytools.io.envi import WriteENVI,envi_header_dict
from hytools.topo.topo import calc_cosine_i
import numpy as np
from scipy.interpolate import interp1d
from pysolar import solar
from scipy.ndimage import uniform_filter
from ..utils.terrain import *
from ..utils.geometry import *
from ..utils.ancillary import *
from ..utils.misc import download_file

home = os.path.expanduser("~")

def he5_to_envi(l1_zip,out_dir,temp_dir,elev_dir,shift = False, rad_coeff = False,
                match=False,proj = False):
    '''
    This function exports three files:
        *_rdn* : Merged and optionally shift corrected radiance cube
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
    out_dir(str): Output directory of ENVI datasets
    temp_dir(str): Temporary directory for intermediate
    elev_dir (str): Directory zipped Copernicus elevation tiles or url to AWS Copernicus data
                    ex : 'https://copernicus-dem-30m.s3.amazonaws.com/'
    shift (bool) : Apply wavelength shift correction surface file
    rad_coeff (bool) : Apply radiometric correction coefficients file
    match (bool or string) : Perform landsat image matching, if string path to reference file
    proj (bool) : Project image to UTM grid
    '''

    base_name = os.path.basename(l1_zip)[16:-4]
    out_dir = f'{out_dir}/PRS_{base_name}/'

    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)

    logging.basicConfig(filename= f'{out_dir}/PRS_{base_name}.log',
            format='%(asctime)s: %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S',
            level=logging.NOTSET)

    temp_dir = f'{temp_dir}/tmpPRS_{base_name}/'
    if not os.path.isdir(temp_dir):
        os.mkdir(temp_dir)

    zip_base  =os.path.basename(l1_zip)
    logging.info('Unzipping %s' % zip_base)
    with zipfile.ZipFile(l1_zip,'r') as zipped:
        zipped.extractall(temp_dir)

    l1_obj = h5py.File(f'{temp_dir}/PRS_L1_STD_OFFL_{base_name}.he5','r')
    version = l1_obj.attrs['Processor_Version'].decode('UTF-8')
    version_str = version.replace('.','').replace('-','')

    logging.info('PRISMA processor: %s' % version)

    apply_shift = False
    if shift:
        apply_shift = True
        shift_obj = np.load(shift)
        if f'shifts_v{version_str}' in shift_obj.keys():
            shift_surface = shift_obj[f'shifts_v{version_str}']
            interp_kind =str(shift_obj['interp_kind'])

        else:
            print(f'Smile: Processor version {version_str} not found.')
            apply_shift = False

    apply_coeff = False
    if rad_coeff:
        apply_coeff = True
        coeff_obj = np.load(rad_coeff)
        if f'coeffs_v{version_str}' in coeff_obj.keys():
            coeff_arr = coeff_obj[f'coeffs_v{version_str}']
        else:
            print(f'Rad coefficients: Processor version {version_str} not found.')
            apply_coeff = False


    #Define output paths
    if proj:
        rdn_file = f'{temp_dir}/PRS_{base_name}_rdn'
        loc_file = f'{temp_dir}/PRS_{base_name}_loc'
        obs_file = f'{temp_dir}/PRS_{base_name}_obs'
    else:
        rdn_file = f'{out_dir}/PRS_{base_name}_rdn'
        loc_file = f'{out_dir}/PRS_{base_name}_loc'
        obs_file = f'{out_dir}/PRS_{base_name}_obs'

    measurement = 'rdn'
    logging.info('Exporting radiance data')

    # Export VNIR to temporary ENVI
    vnir_data =  l1_obj['HDFEOS']["SWATHS"]['PRS_L1_HCO']['Data Fields']['VNIR_Cube']
    vnir_waves = l1_obj.attrs.get('List_Cw_Vnir')
    vnir_fwhm = l1_obj.attrs.get('List_Fwhm_Vnir')

    rdn_dict = envi_header_dict ()
    rdn_dict['lines']= vnir_data.shape[0]
    rdn_dict['samples']= vnir_data.shape[2]
    rdn_dict['bands']=  vnir_data.shape[1]
    rdn_dict['wavelength']= vnir_waves
    rdn_dict['fwhm']= vnir_fwhm
    rdn_dict['interleave']= 'bsq'
    rdn_dict['data type'] = 12
    rdn_dict['wavelength units'] = "nanometers"
    rdn_dict['byte order'] = 0
    vnir_temp = f'{temp_dir}/PRS_{base_name}_{measurement}_vnir'

    writer = WriteENVI(vnir_temp,rdn_dict )
    writer.write_chunk(np.moveaxis(vnir_data[:,:,:],1,2), 0,0)

    # Export SWIR to temporary ENVI
    swir_data =  l1_obj['HDFEOS']["SWATHS"]['PRS_L1_HCO']['Data Fields']['SWIR_Cube']
    swir_waves = l1_obj.attrs.get('List_Cw_Swir')
    swir_fwhm = l1_obj.attrs.get('List_Fwhm_Swir')

    rdn_dict['lines']= swir_data.shape[0]
    rdn_dict['samples']= swir_data.shape[2]
    rdn_dict['bands']=  swir_data.shape[1]
    rdn_dict['wavelength']= swir_waves
    rdn_dict['fwhm']= swir_fwhm
    swir_temp = f'{temp_dir}/PRS_{base_name}_{measurement}_swir'

    writer = WriteENVI(swir_temp,rdn_dict )
    writer.write_chunk(np.moveaxis(swir_data[:,:,:],1,2), 0,0)

    vnir_waves = np.flip(vnir_waves[3:]) #6
    swir_waves = np.flip(swir_waves[:-6]) #-3

    vnir_fwhm = np.flip(vnir_fwhm[3:])
    swir_fwhm = np.flip(swir_fwhm[:-6])

    vnir_obj = ht.HyTools()
    vnir_obj.read_file(vnir_temp, 'envi')

    swir_obj = ht.HyTools()
    swir_obj.read_file(swir_temp, 'envi')

    rdn_dict  = envi_header_dict()
    rdn_dict ['description'] = f"Radiance v{version} micro-watts/cm^2/nm/sr"
    rdn_dict ['lines']= vnir_obj.lines-4 #Clip edges of array
    rdn_dict ['samples']=vnir_obj.columns-4  #Clip edges of array
    rdn_dict ['bands']= len(vnir_waves.tolist() + swir_waves.tolist())
    rdn_dict ['wavelength']= vnir_waves.tolist() + swir_waves.tolist()
    rdn_dict ['fwhm']= vnir_fwhm.tolist() + swir_fwhm.tolist()
    rdn_dict ['interleave']= 'bil'
    rdn_dict ['data type'] = 4
    rdn_dict ['wavelength units'] = "nanometers"
    rdn_dict ['byte order'] = 0
    rdn_dict ['default bands'] = [int(vnir_obj.wave_to_band(850)),
                                  int(vnir_obj.wave_to_band(660)),
                                  int(vnir_obj.wave_to_band(560))]

    writer = WriteENVI(rdn_file,rdn_dict)
    iterator_v =vnir_obj.iterate(by = 'line')
    iterator_s =swir_obj.iterate(by = 'line')

    while not iterator_v.complete:
        chunk_v = iterator_v.read_next()[:,3:]
        chunk_v =np.flip(chunk_v,axis=1)
        chunk_s = iterator_s.read_next()[:,:-6]
        chunk_s =np.flip(chunk_s,axis=1)

        if (iterator_v.current_line >=2) and (iterator_v.current_line <= 997):
            if (measurement == 'rdn') & apply_shift:
                vnir_interpolator = interp1d(vnir_waves+shift_surface[iterator_v.current_line-2,:63],
                                               chunk_v[2:-2,:],fill_value = "extrapolate",
                                               kind=interp_kind)
                chunk_v = vnir_interpolator(vnir_waves)
                swir_interpolator = interp1d(swir_waves+shift_surface[iterator_v.current_line-2,63:],
                                               chunk_s[2:-2,:],fill_value = "extrapolate",
                                               kind=interp_kind)
                chunk_s = swir_interpolator(swir_waves)

                line = np.concatenate([chunk_v,chunk_s],axis=1)/1000.

            else:
                line = np.concatenate([chunk_v,chunk_s],axis=1)[2:-2,:]/1000.

            #Apply rad coeffs
            if apply_coeff:
                line*=coeff_arr[iterator_v.current_line-2,:]

            writer.write_line(line, iterator_v.current_line-2)

    #Load ancillary datasets
    geo =  l1_obj['HDFEOS']["SWATHS"]['PRS_L1_HCO']['Geolocation Fields']
    pvs =  l1_obj['Info']["Ancillary"]['PVSdata']

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
    utc_time = np.ones(geo['Longitude_VNIR'][:,:].shape[0]) *utc_time[:,np.newaxis]
    utc_time = np.rot90(utc_time[2:-2,2:-2] ,k=1)

    mjd2000_epoch = dt.datetime(2000,1, 1,)
    mjd2000_epoch = mjd2000_epoch.replace(tzinfo=dt.timezone.utc)
    mean_time = mjd2000_epoch + dt.timedelta(days=np.array(geo['Time'][:]).mean())

    longitude= geo['Longitude_VNIR'][2:-2,2:-2]
    latitude= geo['Latitude_VNIR'][2:-2,2:-2]

    #Create initial terrain datsets
    elevation,slope,aspect = terrain_generate(longitude,latitude,elev_dir,temp_dir)
    zone,direction = utm_zone(longitude,latitude)

    # Calculate satellite X,Y,Z position for each line
    ''' GPS data are sampled at 1Hz resulting in steps in the
        position data, a line is fit to each dimension to estimate
        continuous position.

        There are more GPS samples than there are lines, to allign
        the GPS signal with the line, we use the provided 'Time'
        information for each line to match with the GPS data.

        When converting GPS time to UTC we use 17 sec difference
        instead of 18 sec because it matches the time provided in
        the time array.
    '''

    # Convert satellite GPS position time to UTC
    sat_t = []
    for second,week in zip(pvs['GPS_Time_of_Last_Position'][:].flatten(),pvs['Week_Number'][:].flatten()):
        gps_second = week*7*24*60*60 + second
        gps_epoch = dt.datetime(1980, 1, 6)
        gps_time  = gps_epoch+ dt.timedelta(seconds=gps_second - 17)
        sat_t.append(gps_time.hour*3600 + gps_time.minute*60. + gps_time.second)
    sat_t = np.array(sat_t)[:,np.newaxis]

    # Convert line MJD2000 to UTC
    grd_t = []
    for day in geo['Time'][:].flatten():
        time = mjd2000_epoch + dt.timedelta(days=day)
        grd_t.append(time.hour*3600 + time.minute*60. + time.second)
    grd_t = np.array(grd_t)[:,np.newaxis]

    #Fit a line to ground time
    X = np.concatenate([np.arange(1000)[:,np.newaxis], np.ones(grd_t.shape)],axis=1)
    time_slope, intercept = np.linalg.lstsq(X,grd_t,rcond=-1)[0].flatten()
    line_t_linear = time_slope*np.arange(1000)+ intercept

    #Fit a line to satellite time
    measurements = np.arange(len(sat_t))
    X = np.concatenate([measurements[:,np.newaxis], np.ones(sat_t.shape)],axis=1)
    time_slope, intercept = np.linalg.lstsq(X,sat_t,rcond=-1)[0].flatten()
    sat_t_linear = time_slope*measurements+ intercept

    # Interpolate x,y,z satelite positions
    sat_xyz = []
    for sat_pos in ['x','y','z']:
        sat_p = np.array(pvs[f'Wgs84_pos_{sat_pos}'][:])
        sat_slope, intercept = np.linalg.lstsq(X,sat_p,rcond=-1)[0].flatten()
        sat_p_linear = sat_slope*measurements+ intercept
        interpolator = interp1d(sat_t_linear,sat_p_linear,
                                fill_value="extrapolate",kind = 'linear')
        sat_interp = interpolator(line_t_linear)
        sat_xyz.append(sat_interp[2:-2])
    sat_xyz = np.array(sat_xyz)

    # Calculate sensor to ground pathlength
    grd_xyz = np.array(dda2ecef(longitude,latitude,elevation))
    path = pathlength(sat_xyz,grd_xyz)

    # Export satellite position to csv
    sat_lon,sat_lat,sat_alt = ecef2dda(sat_xyz[0],sat_xyz[1],sat_xyz[2])
    satellite_df = pd.DataFrame()
    satellite_df['lat'] = sat_lat
    satellite_df['lon'] = sat_lon
    satellite_df['alt'] = sat_alt
    satellite_df.to_csv(f'{out_dir}/PRS_{base_name}_satellite_loc.csv')

    # Convert satellite coords to local ENU
    sat_enu  = np.array(dda2utm(sat_lon,sat_lat,sat_alt,
                       utm_zone(longitude,latitude)))
    # Convert ground coords to local ENU
    easting,northing,up  =dda2utm(longitude,latitude,
                                elevation)

    # Calculate sensor geometry
    sensor_zn,sensor_az = sensor_view_angles(sat_enu,
                                             np.array([easting,northing,up]))

    # Perform image matching
    if match:
        coords =np.concatenate([np.expand_dims(easting.flatten(),axis=1),
                                np.expand_dims(northing.flatten(),axis=1)],axis=1)
        warp_east = easting.min()-100
        warp_north =northing.max()+100
        pixel_size = 30

        project = Projector()
        project.create_tree(coords,easting.shape)
        project.query_tree(warp_east,warp_north,pixel_size)

        # Project independent variables
        sensor_az_prj = project.project_band(sensor_az,-9999)
        sensor_zn_prj = project.project_band(sensor_zn,-9999)
        elevation_prj = project.project_band(elevation.astype(np.float),-9999)

        radiance = ht.HyTools()
        radiance.read_file(rdn_file, 'envi')

        #Average over Landsat 8 Band 5 bandwidth and warp
        unwarp_band = np.zeros(longitude.shape)
        for wave in range(850,890,10):
            unwarp_band += radiance.get_wave(wave)/7.
        warp_band = project.project_band(unwarp_band,-9999)
        warp_band = 16000*(warp_band-warp_band.min())/warp_band.max()

        if isinstance(match,bool):
            landsat,land_east,land_north = get_landsat_image(longitude,latitude,
                                                             mean_time.month,
                                                             max_cloud = 5)
        else:
            lst = ht.HyTools()
            lst.read_file(match,'envi')
            landsat = lst.get_band(0)
            land_east = float(lst.map_info[3])
            land_north = float(lst.map_info[4])

        #Calculate offsets between reference and input images
        offset_x = int((warp_east-land_east)//pixel_size)
        offset_y = int((land_north-warp_north)//pixel_size)

        #Calculate optimal shift
        y_model,x_model = image_match(landsat,warp_band,
                                      offset_x,offset_y,
                                      sensor_zn_prj,sensor_az_prj,elevation_prj)

        #Apply uniform filter
        smooth_elevation = uniform_filter(elevation,25)
        smooth_az = uniform_filter(sensor_az,25)
        smooth_zn = uniform_filter(sensor_zn,25)

        # Generate y and x offset surfaces
        i,a,b,c = y_model
        y_offset = i + a*smooth_zn +b*smooth_az + c*smooth_elevation

        i,a,b,c= x_model
        x_offset = i + a*smooth_zn +b*smooth_az + c*smooth_elevation

        # Calculate updated coordinates
        easting = easting+  30*x_offset
        northing = northing- 30*y_offset

        zone,direction = utm_zone(longitude,latitude)
        longitude,latitude = utm2dd(easting,northing,zone,direction)

        #Recalculate elevation with new coordinates
        logging.info('Rebuilding DEM')
        elevation,slope,aspect= terrain_generate(longitude,latitude,elev_dir,temp_dir)
    # Solar geometries
    '''Solar geometry is calculated based on the mean scene acquisition time
    which varies by less than 5 seconds from start to end of the scene and is
    computationally more efficient.
    '''

    #Calculate solar geometry
    solar_az = solar.get_azimuth(geo['Latitude_VNIR'][:,:],geo['Longitude_VNIR'][:,:],mean_time)[2:-2,2:-2]
    solar_zn = 90-solar.get_altitude(geo['Latitude_VNIR'][:,:],geo['Longitude_VNIR'][:,:],mean_time)[2:-2,2:-2]

    # Recalculate sensor geometry with updated coordinates
    sensor_zn,sensor_az = sensor_view_angles(sat_enu,
                                             np.array([easting,northing,up]))


    # Export location datacube
    loc_export(loc_file,longitude,latitude,elevation)


    cosine_i = calc_cosine_i(np.radians(solar_zn),
                             np.radians(solar_az),
                             np.radians(slope),
                             np.radians(aspect))
    rel_az = np.radians(solar_az-sensor_az)
    phase =  np.arccos(np.cos(np.radians(solar_zn)))*np.cos(np.radians(solar_zn))
    phase += np.sin(np.radians(solar_zn))*np.sin(np.radians(solar_zn))*np.cos(rel_az)

    # Export observables datacube
    obs_export(obs_file,path,sensor_az,sensor_zn,
               solar_az,solar_zn,phase,slope,aspect,
               cosine_i,utc_time)

    # Get image bounds coordinates
    corner_1 = [longitude[0,0],  latitude[0,0]]
    corner_2 = [longitude[0,-1], latitude[0,-1]]
    corner_3 = [longitude[-1,-1],latitude[-1,-1]]
    corner_4 = [longitude[-1,0], latitude[-1,0]]

    start_time = dt.datetime.strptime(base_name.split('_')[0],'%Y%m%d%H%M%S')
    end_time = dt.datetime.strptime(base_name.split('_')[1],'%Y%m%d%H%M%S')

    if proj:
        #Create new projector with corrected coordinates
        new_coords =np.concatenate([np.expand_dims(easting.flatten(),axis=1),
                        np.expand_dims(northing.flatten(),axis=1)],axis=1)

        res = 30

        project = Projector()
        project.create_tree(new_coords,easting.shape)
        project.query_tree(easting.min()-600,northing.max()+600,res)

        no_data_mask = project.project_band(np.zeros(project.output_shape),-9999) == -9999
        start_col = np.argwhere(no_data_mask.sum(axis=0) != project.output_shape[0])[0][0]
        end_col = start_col + len(np.argwhere(no_data_mask.sum(axis=0) != project.output_shape[0]))

        start_line = np.argwhere(no_data_mask.sum(axis=1) != project.output_shape[1])[0][0]
        end_line = start_line + len(np.argwhere(no_data_mask.sum(axis=1) != project.output_shape[1]))

        map_info = ['UTM',
                    1,
                    1,
                    easting.min()-600 - (res/2) + res*start_col,
                    northing.max()+600 + (res/2) - res*start_line,
                    res,
                    res,
                    zone,
                    direction,
                    'WGS-84' ,
                    'units=Meters']

        logging.info('Projecting datasets to WGS84 UTM at %sm resolution' % res)
        for file in ['rdn','loc','obs']:
            input_name = f'{temp_dir}/PRS_{base_name}_{file}'
            hy_obj = ht.HyTools()
            hy_obj.read_file(input_name, 'envi')
            iterator =hy_obj.iterate(by = 'band')

            out_header = hy_obj.get_header()
            out_header['lines']=  end_line - start_line
            out_header['samples']= end_col - start_col
            out_header['data ignore value'] = -9999
            out_header['map info'] = map_info
            out_header['start acquisition time'] = start_time.strftime('%Y-%m-%dT%H:%M:%SZ')
            out_header['end acquisition time'] = end_time.strftime('%Y-%m-%dT%H:%M:%SZ')
            out_header['bounding box'] =[corner_1,corner_2,corner_3,corner_4]
            out_header['sensor type'] ='PRISMA'

            output_name = f'{out_dir}/PRS_{base_name}_{file}_prj'
            writer = WriteENVI(output_name,out_header)

            while not iterator.complete:
                band = project.project_band(iterator.read_next(),-9999)
                if file == 'rdn':
                    band[(band<0) & (band>-9999)] = 0
                band[np.isnan(band)] = -9999
                writer.write_band(band[start_line:end_line,start_col:end_col],iterator.current_band)

    logging.info('Deleting temporary files')
    shutil.rmtree(temp_dir)
