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

import xml.etree.ElementTree as ET
import numpy as np
import gdal
import hytools as ht
import pyproj
import datetime as dt
import os
import zipfile
import shutil
import pandas as pd
import hytools as ht
from hytools.io.envi import WriteENVI,envi_header_dict
from hytools.topo.topo import calc_cosine_i
from scipy.interpolate import interp1d
from pysolar import solar
from scipy.interpolate import griddata
from scipy.optimize import curve_fit
from skimage.util import view_as_blocks
from scipy.ndimage import uniform_filter
from ..utils.terrain import *
from ..utils.geometry import *
from ..utils.ancillary import *


def gaussian(x, mu, fwhm):
    sig = fwhm/(2* np.sqrt(2*np.log(2)))
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

def desis_to_envi(base_name,l1b_zip,l1c_zip,out_dir,elev_dir,temp_dir,
                  match=None,project = True, res = 30):
    '''
    This function exports a PRISMA L1 radiance products to an ENVI formatted
    binary file along with ancillary datasets including location data and geometry.
    This scripts requires L1 radiance, L2c reflectance and SRTM elevation products for a scene.
    '''

    #Testing
    # elev_dir = '/data2/cop_dsm/'
    # base_name = "DT0488344520_005-20200818T141910"
    # l1b_zip = '/data2/desis/zip/DESIS-HSI-L1B-%s-V0210.zip' %  base_name
    # l1c_zip = '/data2/desis/zip/DESIS-HSI-L1C-%s-V0210.zip' % base_name
    # out_dir = '/data2/desis/envi/DESIS_%s/' % base_name
    # temp_dir = '/data2/temp/'
    # match = '/data2/landsat/LC08_L2SP_025028_20200617_20200824_02_T1_SR_B5.TIF'
    # rfl = False
    # project = True
    # res = 30

    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)

    temp_dir = '%s/DESIS_%s/'% (temp_dir,base_name)
    if not os.path.isdir(temp_dir):
        os.mkdir(temp_dir)

    for file in [l1b_zip,l1c_zip]:
        zip_base  =os.path.basename(file)
        print('Unzipping %s' % zip_base)
        with zipfile.ZipFile(file,'r') as zipped:
            zipped.extractall(temp_dir)

    l1  = '%sPRS_L1_STD_OFFL_%s.he5' % (temp_dir,base_name)
    l2c  = '%sPRS_L2C_STD_%s.he5' % (temp_dir,base_name)


    l1b_file = gdal.Open('%s/DESIS-HSI-L1B-%s-V0210-SPECTRAL_IMAGE.tif' % (temp_dir,base_name))
    l1c_file = gdal.Open('%s/DESIS-HSI-L1C-%s-V0210-SPECTRAL_IMAGE.tif' % (temp_dir,base_name))

    # Parse relevant metadata from XML file, assume metadata are in same directory as iamges
    tree = ET.parse('%s/DESIS-HSI-L1B-%s-V0210-METADATA.xml' % (temp_dir,base_name))
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
    mirror_angle = float(specific.findall('pointingMirrorAngle')[0].text)

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
    l1c_band = l1c_file.ReadAsArray().mean(axis=0)
    rfl_mask = (l1c_band !=l1c_band[0][0]).astype(int)

    # Get bounding coordinates of scene
    coord_dict= {}
    polygon = base.findall('spatialCoverage')[0].findall('boundingPolygon')[0]
    for point in polygon:
        name= point.findall('frame')[0].text
        lat= float(point.findall('latitude')[0].text)
        lon= float(point.findall('longitude')[0].text)
        coord_dict[name] = [lat,lon]

    # Get corner coordinates of projects L2 scene
    top = np.argmax(np.cumsum(rfl_mask[0]))
    bottom = np.argmax(np.cumsum(rfl_mask[-2]))
    left = np.argmax(np.cumsum(rfl_mask[:,0]))
    right = np.argmax(np.cumsum(rfl_mask[:,-1]))

    ulx,pixel,a,uly,b,c =l1c_file.GetGeoTransform()

    coords = np.array([[top,0],[l1c_file.RasterXSize,right],[bottom,l1c_file.RasterYSize],[0,left]]).T
    coords[0] = ulx + coords[0]*pixel
    coords[1] = uly - coords[1]*pixel

    map_proj = pyproj.Proj(l1c_file.GetProjection())
    lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')
    corner_lon,corner_lat = pyproj.transform(map_proj, lla,
                                coords[0],
                                coords[1],
                                radians=False)

    #Update coord_dict with close matching corner coordinate
    for point in coord_dict:
        if 'point' in point:
            lat,lon = coord_dict[point]
            dist = np.sqrt((corner_lon-lon)**2 + (corner_lat-lat)**2)
            closest = np.argmin(dist)
            print(point)
            print(lat,lon)
            print(corner_lat[closest],corner_lon[closest])
            coord_dict[point]= [corner_lat[closest],corner_lon[closest]]
            print(corner_lat[closest]-lat,corner_lon[closest]-lon)

    # Get ISS altitude
    altitude_m = float(base.findall('altitudeCoverage')[0].text)

    raster = l1b_file.ReadAsArray()
    mask = raster[1].astype(float)
    mask = mask==mask[0][0]

    rad_dict  = envi_header_dict ()
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

    if project:
        output_name = '%sDESIS_%s_rad_unprj' % (temp_dir,base_name)
    else:
        output_name = '%sDESIS_%s_rad_unprj' % (out_dir,base_name)

    writer = WriteENVI(output_name,rad_dict )

    #Write VNIR cube
    print('Exporting radiance data')
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

    # Export observable datacube
    ################################################################
    solar_az = solar.get_azimuth(latitude,longitude,start_time)
    solar_zn = 90-solar.get_altitude(latitude,longitude,start_time)

    ecef = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
    lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')

    x, y, z = pyproj.transform(lla, ecef,
                               longitude,
                               latitude,
                               elevation, radians=False)

    # Calculate sensor-to-ground pathlength
    satelite_xyz = []
    pathlength = np.zeros(x.shape)
    for sat_coord,grd_coord in zip(orbit_data.T,[x,y,z]):
        line_grid = np.linspace(0,1,len(sat_coord))*longitude.shape[0]
        interpolator = interp1d(line_grid,sat_coord)
        sat_interp = interpolator(np.arange(longitude.shape[0]))
        pathlength += (sat_interp[:,np.newaxis]-grd_coord)**2
        satelite_xyz.append(sat_interp)
    pathlength = np.sqrt(pathlength)

    # Export satellite position to csv
    sv_lon,sv_lat,sv_alt = pyproj.transform(ecef,lla,
                                satelite_xyz[0],
                                satelite_xyz[1],
                                satelite_xyz[2],
                                radians=False)
    satellite_df = pd.DataFrame()
    satellite_df['lat'] = sv_lat
    satellite_df['lon'] = sv_lon
    satellite_df['alt'] = sv_alt
    satellite_df.to_csv('%sDESIS_%s_satellite_loc.csv' % (out_dir,base_name))

    #Calculate sensor zenith angle
    sv_x = (x-satelite_xyz[0][:,np.newaxis])[:,:,np.newaxis]
    sv_y = (y-satelite_xyz[1][:,np.newaxis])[:,:,np.newaxis]
    sv_z = (z-satelite_xyz[2][:,np.newaxis])[:,:,np.newaxis]
    sv_xyz = np.concatenate([sv_x,sv_y,sv_z],axis=2)
    grnd_xyz = np.concatenate([x[:,:,np.newaxis],
                                y[:,:,np.newaxis],
                                z[:,:,np.newaxis],],axis=2)
    dot = np.einsum('ijk,ijk->ij', sv_xyz,grnd_xyz)
    denom = np.linalg.norm(sv_xyz,axis=2)* np.linalg.norm(grnd_xyz,axis=2)
    sensor_zn = 180-np.degrees(np.arccos(dot/denom))


    # Calculate sensor angles.
    DX = x-np.expand_dims(satelite_xyz[0], axis=1)
    DY = y-np.expand_dims(satelite_xyz[1], axis=1)
    DZ = z-np.expand_dims(satelite_xyz[2], axis=1)

    view_azimuth = np.arcsin(np.abs(DX)/np.sqrt(DX**2+DY**2))
    ind = (DX>0)&(DY<0)
    view_azimuth[ind]=np.pi-view_azimuth[ind]
    ind = (DX<0)&(DY<0)
    view_azimuth[ind]=np.pi+view_azimuth[ind]
    ind = (DX<0)&(DY>0)
    view_azimuth[ind]=2*np.pi-view_azimuth[ind]
    sensor_az = np.degrees(view_azimuth)

    print("####################")
    print("Scene zenith:",scene_zn)
    print("Calc zenith:",sensor_zn[511,511])
    print("Zenith diff:",scene_zn-sensor_zn[511,511])
    print("####################")
    print("Scene azimuth:",scene_az)
    print("Calc azimuth:",sensor_az[511,511])
    print("Azimuth diff:",scene_az-sensor_az[511,511])

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

        rad_file = '%sDESIS_%s_rad_unprj' % (temp_dir,base_name)
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

        easting = easting+  30*x_offset
        northing = northing- 30*y_offset

        zone,direction = utm_zone(longitude,latitude)
        longitude,latitude = utm2dd(easting,northing,zone,direction)

        #Recalculate elevation with new coordinates
        elevation= dem_generate(longitude,latitude,elev_dir,temp_dir)

    if project:
        loc_file = '%sDESIS_%s_loc_unprj' % (temp_dir,base_name)
    else:
        loc_file = '%sDESIS_%s_loc_unprj' % (out_dir,base_name)

    loc_export(loc_file,longitude,latitude,elevation)

    slope,aspect = slope_aspect(elevation,temp_dir)

    # Calculate illumination datasets
    cosine_i = calc_cosine_i(np.radians(solar_zn),
                             np.radians(solar_az),
                             np.radians(slope),
                             np.radians(aspect))

    phase =  np.arccos(np.cos(np.radians(solar_zn)))*np.cos(np.radians(sensor_zn))
    phase += np.sin(np.radians(solar_zn))*np.sin(np.radians(sensor_zn))*np.cos(np.radians(sensor_az-solar_az))

    utc_time = (lines/(l1b_file.RasterYSize) * (end_time-start_time).seconds)/60/60
    utc_time+= start_time.hour + start_time.minute/60
    utc_time = utc_time[:,85:]

    # Export observable datacube
    if project:
        obs_file = '%sDESIS_%s_obs_unprj' % (temp_dir,base_name)
    else:
        obs_file = '%sDESIS_%s_obs_unprj' % (out_dir,base_name)

    obs_export(obs_file,pathlength,sensor_az,sensor_zn,
               solar_az,solar_zn,phase,slope,aspect,
               cosine_i,utc_time)

    if project:

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

        print('Georeferencing datasets')
        for file in ['rad','loc','obs']:
            print(file)
            input_name = '%sDESIS_%s_%s_unprj' % (temp_dir,base_name,file)
            hy_obj = ht.HyTools()
            hy_obj.read_file(input_name, 'envi')
            iterator =hy_obj.iterate(by = 'band')

            out_header = hy_obj.get_header()
            out_header['lines']= project.output_shape[0]//blocksize
            out_header['samples']=project.output_shape[1]//blocksize
            out_header['data ignore value'] = -9999
            out_header['map info'] = map_info

            output_name = '%sDESIS_%s_%s_geo' % (out_dir,base_name,file)
            writer = WriteENVI(output_name,out_header)

            while not iterator.complete:
                band = project.project_band(iterator.read_next(),-9999)
                band[band == -9999] = np.nan
                band = np.nanmean(view_as_blocks(band[:out_lines,:out_cols], (blocksize,blocksize)),axis=(2,3))
                band[band<0] = 0
                band[np.isnan(band)] = -9999
                writer.write_band(band,iterator.current_band)

    shutil.rmtree(temp_dir)