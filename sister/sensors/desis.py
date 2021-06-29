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
import os
import logging
import shutil
import zipfile
import xml.etree.ElementTree as ET
try:
    import gdal
except:
    from osgeo import gdal
import hytools as ht
from hytools.io.envi import WriteENVI,envi_header_dict
from hytools.topo.topo import calc_cosine_i
import numpy as np
import pandas as pd
import pyproj
from pysolar import solar
from scipy.interpolate import griddata,interp1d
from scipy.ndimage import uniform_filter
from scipy.optimize import curve_fit
from skimage.util import view_as_blocks
from ..utils.terrain import *
from ..utils.geometry import *
from ..utils.ancillary import *

def gaussian(x, mu, fwhm):
    sig = fwhm/(2* np.sqrt(2*np.log(2)))
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

def geotiff_to_envi(l1b_zip,out_dir,temp_dir,elev_dir,
                    match=None,proj = True, res = 30):

    '''
     This function exports three files:
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

     l1(str): L1B zipped radiance data product path
     out_dir(str): Output directory of ENVI datasets
     temp_dir(str): Temporary directory for intermediate
     elev_dir (str): Directory zipped Copernicus elevation tiles
     shift (str) : Pathname of wavelength shift correction surface file
     match (str or list) : Pathname to Landsat image for image re-registration (recommended)
     proj (bool) : Project image to UTM grid
     res (int) : Resolution of projected image, 30 should be one of its factors (90,120,150.....)

    This functions assumes that the L1B and L1C zipped radiance products are located in the
    same directory
    '''

    base_name = os.path.basename(l1b_zip)[14:-4]

    out_dir = '%s/DESIS_%s/'% (out_dir,base_name)
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)

    temp_dir = '%s/DESIS_%s/'% (temp_dir,base_name)
    if not os.path.isdir(temp_dir):
        os.mkdir(temp_dir)

    zip_base  =os.path.basename(l1b_zip)
    logging.info('Unzipping %s' % zip_base)
    with zipfile.ZipFile(l1b_zip,'r') as zipped:
        zipped.extractall(temp_dir)

    l1b_file = gdal.Open('%s/DESIS-HSI-L1B-%s-SPECTRAL_IMAGE.tif' % (temp_dir,base_name))

    # Parse relevant metadata from XML file, assume metadata are in same directory as iamges
    tree = ET.parse('%s/DESIS-HSI-L1B-%s-METADATA.xml' % (temp_dir,base_name))
    root = tree.getroot()
    specific =  root[3]
    band_meta = {}
    for item in specific.findall('bandCharacterisation'):
        for band in item:
            for meta in band:
                if meta.tag not in band_meta.keys():
                    band_meta[meta.tag] = []
                string =  str(meta.text.encode('utf8'))
                string = string.replace('\\n',' ')
                string = string.replace('b',' ')
                string = string.replace("'",' ')
                values= string.split(',')
                values = [float(x) for x in values]
                if len(values) == 1:
                    values= values[0]
                band_meta[meta.tag].append(values)
    offset = np.array(band_meta['offsetOfBand'])
    gain = np.array(band_meta['gainOfBand'])
    waves = np.array(band_meta['wavelengthCenterOfBand'])
    fwhm = np.array(band_meta['wavelengthWidthOfBand'])
    response = np.array(band_meta['response'])
    response_waves = np.array(band_meta['wavelengths'])

    # Fit a gaussian to the reponse to determine center wavelength and fhwm
    opt_waves = []
    opt_fwhm = []
    for i,wave in enumerate(waves):
        popt, pcov = curve_fit(gaussian, response_waves[i],np.array(response[i])/max(response[i]),[waves[i],fwhm[i]])
        opt_waves.append(popt[0])
        opt_fwhm.append(popt[1])

    scene_az = float(specific.findall('sceneAzimuthAngle')[0].text)
    scene_zn = float(specific.findall('sceneIncidenceAngle')[0].text)

    # Get acquisition start and end time
    base =  root[2]
    time_str = base.findall('temporalCoverage')[0].findall('startTime')[0].text
    time_str=time_str.replace('T',' ').replace('Z','')
    start_time = dt.datetime.strptime(time_str,"%Y-%m-%d %H:%M:%S.%f")
    start_time = start_time.replace(tzinfo=dt.timezone.utc)

    time_str = base.findall('temporalCoverage')[0].findall('endTime')[0].text
    time_str=time_str.replace('T',' ').replace('Z','')
    end_time = dt.datetime.strptime(time_str,"%Y-%m-%d %H:%M:%S.%f")
    end_time = end_time.replace(tzinfo=dt.timezone.utc)
    date = dt.datetime.strftime(start_time , "%Y%m%d")

    # Get orbital data
    orbit_data = []
    for item in specific.findall('orbit'):
        for line in item:
            for point in line.findall('point'):
                time_str = line.findall('timeUTC')[0].text
                time_str=time_str.replace('T',' ').replace('Z','')
                orbit_time = dt.datetime.strptime(time_str,"%Y-%m-%d %H:%M:%S.%f")
                orbit_time = orbit_time.replace(tzinfo=dt.timezone.utc)

                #The metadata contains orbital info beyond the collection window
                # we use the acquisition start and end times to filter points
                if (orbit_time>=start_time) & (orbit_time<=end_time):
                    for location in point.findall('location'):
                        x = float(location.findall('X')[0].text)
                        y = float(location.findall('Y')[0].text)
                        z = float(location.findall('Z')[0].text)
                        orbit_data.append([x,y,z])

    orbit_data = np.array(orbit_data)

    l1b_band = l1b_file.ReadAsArray().mean(axis=0)

    # Get bounding coordinates of scene
    coord_dict= {}
    polygon = base.findall('spatialCoverage')[0].findall('boundingPolygon')[0]
    for point in polygon:
        name= point.findall('frame')[0].text
        lat= float(point.findall('latitude')[0].text)
        lon= float(point.findall('longitude')[0].text)
        coord_dict[name] = [lat,lon]

    # Get ISS altitude
    altitude_m = float(base.findall('altitudeCoverage')[0].text)

    raster = l1b_file.ReadAsArray()
    mask = raster[1].astype(float)
    mask = mask==mask[0][0]

    rad_dict  = envi_header_dict()
    rad_dict ['lines']= l1b_file.RasterYSize
    rad_dict ['samples']= l1b_file.RasterXSize-85
    rad_dict ['bands']= len(waves)-1
    rad_dict ['wavelength']= opt_waves[1:]
    rad_dict ['fwhm']= opt_fwhm[1:]
    rad_dict ['interleave']= 'bil'
    rad_dict ['data type'] = 4
    rad_dict ['wavelength units'] = "nanometers"
    rad_dict ['byte order'] = 0
    rad_dict ['data ignore value'] = -9999
    rad_dict ['default bands'] = [np.argmin(np.abs(waves-660)),
                                  np.argmin(np.abs(waves-560)),
                                  np.argmin(np.abs(waves-460))]

    #Define output paths
    if proj:
        rad_file = '%sDESIS_%s_rdn' % (temp_dir,base_name)
        loc_file = '%sDESIS_%s_loc' % (temp_dir,base_name)
        obs_file = '%sDESIS_%s_obs' % (temp_dir,base_name)
    else:
        rad_file = '%sDESIS_%s_rdn' % (out_dir,base_name)
        loc_file = '%sDESIS_%s_loc' % (out_dir,base_name)
        obs_file = '%sDESIS_%s_obs' % (out_dir,base_name)

    writer = WriteENVI(rad_file,rad_dict )

    #Write VNIR cube
    logging.info('Exporting radiance data')
    for line_num in range(l1b_file.RasterYSize):
        line = raster[:,line_num,:].astype(float)
        line = line*gain[:,np.newaxis] + offset[:,np.newaxis]
        line = line[1:,85:].T
        writer.write_line(line,line_num)

    del raster

    # Location datacube
    ###########################################################################
    lines,columns = np.indices((l1b_file.RasterYSize,l1b_file.RasterXSize))

    lat_vals= []
    lon_vals = []
    points = [[0,0],
              [l1b_file.RasterYSize,0],
              [l1b_file.RasterYSize,l1b_file.RasterXSize],
              [0,l1b_file.RasterXSize]]

    for point in [1,2,3,4]:
        lat_vals.append(coord_dict['point_%s' % point][0])
        lon_vals.append(coord_dict['point_%s' % point][1])

    longitude = griddata(points, lon_vals, (lines,columns), method='linear')[:,85:]
    latitude = griddata(points, lat_vals, (lines,columns), method='linear')[:,85:]

    #Create initial elevation raster
    elevation= dem_generate(longitude,latitude,elev_dir,temp_dir)
    zone,direction = utm_zone(longitude,latitude)

    solar_az = solar.get_azimuth(latitude,longitude,start_time)
    solar_zn = 90-solar.get_altitude(latitude,longitude,start_time)

    ecef = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
    lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')

    grd_xyz = np.array(pyproj.transform(lla, ecef,
                               longitude,
                               latitude,
                               elevation, radians=False))

    # Calculate satellite XYZ position
    sat_xyz = []
    line_grid = np.linspace(0,1,orbit_data.shape[0])*longitude.shape[0]
    for sat_coord in orbit_data.T:
        interpolator = interp1d(line_grid,sat_coord)
        sat_interp = interpolator(np.arange(longitude.shape[0]))
        sat_xyz.append(sat_interp)
    sat_xyz = np.array(sat_xyz)
    path =np.linalg.norm(sat_xyz[:,:,np.newaxis]-grd_xyz,axis=0)
    # Export satellite position to csv
    sat_lon,sat_lat,sat_alt = ecef2dda(sat_xyz[0],sat_xyz[1],sat_xyz[2])
    satellite_df = pd.DataFrame()
    satellite_df['lat'] = sat_lat
    satellite_df['lon'] = sat_lon
    satellite_df['alt'] = sat_alt
    satellite_df.to_csv('%sDESIS_%s_satellite_loc.csv' % (out_dir,base_name))


    # Convert satellite coords to local ENU
    sat_enu  = np.array(dda2utm(sat_lon,sat_lat,sat_alt,
                       utm_zone(longitude,latitude)))
    # Convert ground coords to local ENU
    easting,northing,up  =dda2utm(longitude,latitude,
                                elevation)

    # Calculate sensor geometry
    sensor_zn,sensor_az = sensor_view_angles(sat_enu,
                                             np.array([easting,northing,up]))


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
        radiance.read_file(rad_file, 'envi')

        #Average over Landsat 8 Band 5 bandwidth and warp
        warp_band = np.zeros(longitude.shape)
        for wave in range(850,890,10):
            warp_band += radiance.get_wave(wave)/7.
        warp_band = project.project_band(warp_band,-9999)
        warp_band = 16000*(warp_band-warp_band.min())/warp_band.max()

        landsat,land_east,land_north = get_landsat_image(longitude,latitude,
                                                         end_time.month,max_cloud = 5)

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
        elevation= dem_generate(longitude,latitude,elev_dir,temp_dir)

    loc_export(loc_file,longitude,latitude,elevation)

    # Generate remaining observable layers
    slope,aspect = slope_aspect(elevation,temp_dir)
    cosine_i = calc_cosine_i(np.radians(solar_zn),
                             np.radians(solar_az),
                             np.radians(slope),
                             np.radians(aspect))
    rel_az = np.radians(solar_az-sensor_az)
    phase =  np.arccos(np.cos(np.radians(solar_zn)))*np.cos(np.radians(solar_zn))
    phase += np.sin(np.radians(solar_zn))*np.sin(np.radians(solar_zn))*np.cos(rel_az)

    utc_time = (lines/(l1b_file.RasterYSize) * (end_time-start_time).seconds)/60/60
    utc_time+= start_time.hour + start_time.minute/60
    utc_time = utc_time[:,85:]

    obs_export(obs_file,path,sensor_az,sensor_zn,
               solar_az,solar_zn,phase,slope,aspect,
               cosine_i,utc_time)

    if proj:
        #Create new projector with corrected coordinates
        new_coords =np.concatenate([np.expand_dims(easting.flatten(),axis=1),
                        np.expand_dims(northing.flatten(),axis=1)],axis=1)

        project = Projector()
        project.create_tree(new_coords,easting.shape)
        project.query_tree(easting.min()-100,northing.max()+100,30)

        blocksize = int(res/30)
        map_info = ['UTM', 1, 1, easting.min()-100, northing.max()+100,res,
                           res,zone,direction, 'WGS-84' , 'units=Meters']
        out_cols = int(blocksize* (project.output_shape[1]//blocksize))
        out_lines = int(blocksize* (project.output_shape[0]//blocksize))

        logging.info('Georeferencing datasets')
        for file in ['rdn','loc','obs']:
            logging.info(file)
            input_name = '%sDESIS_%s_%s' % (temp_dir,base_name,file)
            hy_obj = ht.HyTools()
            hy_obj.read_file(input_name, 'envi')
            iterator =hy_obj.iterate(by = 'band')

            out_header = hy_obj.get_header()
            out_header['lines']= project.output_shape[0]//blocksize
            out_header['samples']=project.output_shape[1]//blocksize
            out_header['data ignore value'] = -9999
            out_header['map info'] = map_info

            output_name = '%sDESIS_%s_%s_prj' % (out_dir,base_name,file)
            writer = WriteENVI(output_name,out_header)

            while not iterator.complete:
                band = project.project_band(iterator.read_next(),-9999)
                band[band == -9999] = np.nan
                band = np.nanmean(view_as_blocks(band[:out_lines,:out_cols], (blocksize,blocksize)),axis=(2,3))
                if file == 'rdn':
                    band[band<0] = 0
                band[np.isnan(band)] = -9999
                writer.write_band(band,iterator.current_band)
    logging.info('Deleting temporary files')
    shutil.rmtree(temp_dir)
