#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SISTER
Space-based Imaging Spectroscopy and Thermal PathfindER
Author: Adam Chlus

"""
import os
from itertools import tee
import logging
import numpy as np
from rtree import index
from scipy.spatial import cKDTree
from numba import jit
import statsmodels.api as sm
import ray
import pyproj
from skimage.util import view_as_blocks
import ee

try:
    import gdal
except:
    from osgeo import gdal


class Projector():
    """Projector class"""

    def __init__(self):
        self.tree = None
        self.indices = None
        self.pixel_size = None
        self.input_shape = None
        self.inter_shape = None
        self.output_shape = None

        self.mask = None

    def create_tree(self,coords,input_shape):
        self.input_shape = input_shape
        self.tree = cKDTree(coords,balanced_tree= False)

    def query_tree(self,ulx,uly,pixel_size):
        half_pix = pixel_size/2
        self.pixel_size = pixel_size

        lines = (uly-self.tree.mins[1])//half_pix
        if lines%2 !=0:
            lines+=1

        columns = (self.tree.maxes[0]-ulx)//half_pix
        if columns%2 !=0:
            columns+=1

        self.inter_shape = (int(lines),int(columns))
        int_north,int_east = np.indices(self.inter_shape)
        int_east = (int_east*half_pix + ulx).flatten()
        int_north = (uly-int_north*half_pix).flatten()
        int_north= np.expand_dims(int_north,axis=1)
        int_east= np.expand_dims(int_east,axis=1)
        dest_points = np.concatenate([int_east,int_north],axis=1)

        distances,self.indices =  self.tree.query(dest_points,k=1)
        self.indices = np.unravel_index(self.indices,self.input_shape)
        distances =  distances.reshape(self.inter_shape)
        self.mask = distances > 2*pixel_size
        self.output_shape = (int(lines/2),int(columns/2))

    def project_band(self,band,no_data):
        band = np.copy(band[self.indices[0],self.indices[1]].reshape(self.inter_shape))
        band[self.mask] = np.nan
        band = np.nanmean(view_as_blocks(band, (2,2)),axis=(2,3))
        band[np.isnan(band)] = no_data
        return band


def utm_zone(longitude, latitude):
    ''' Returns UTM zone and direction
    '''
    zone = int(np.ceil((longitude.min()+180)/6))
    if latitude.mean() >0:
        direction = 'N'
    else:
        direction = 'S'
    return zone,direction


def dda2utm(longitude,latitude,altitude,zn_dir =None):
    '''Convert coordinates in decimal degrees into
    UTM eastings and northings
    '''

    if zn_dir:
        zone,direction =zn_dir
    else:
        zone,direction = utm_zone(longitude, latitude)

    if direction == 'N':
        epsg_dir = 6
    else:
        epsg_dir = 7

    out_crs = pyproj.Proj("+init=EPSG:32%s%02d" % (epsg_dir,zone))
    in_crs= pyproj.Proj("+init=EPSG:4326")

    easting,northing,up = pyproj.transform(in_crs,out_crs,
                                           longitude,
                                           latitude,
                                           altitude)
    return easting,northing,up


def dda2ecef(longitude,latitude,altitude):
    '''Convert longitude,latitude and altitude to Earth Centered Earth Fixed
    coordinates (x,y,x)

    '''
    in_crs = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')
    out_crs = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')

    x, y, z = pyproj.transform(in_crs, out_crs,
                               longitude,
                               latitude,
                               altitude, radians=False)
    return x,y,z

def ecef2dda(x,y,z):
    '''Convert Earth Centered Earth Fixed
    coordinates (x,y,z) to longitude,latitude and altitude

    '''
    in_crs = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
    out_crs = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')

    lon, lat, alt = pyproj.transform(in_crs, out_crs,
                               x,y,z, radians=False)
    return lon, lat, alt



def utm2dd(easting,northing,zone,direction):
    '''Convert coordinates in UTM eastings and northings into
        decimal degrees
    '''
    if direction == 'N':
        epsg_dir = 6
    else:
        epsg_dir = 7

    in_crs = pyproj.Proj("+init=EPSG:32%s%02d" % (epsg_dir,zone))
    out_crs= pyproj.Proj("+init=EPSG:4326")
    longitude,latitude= pyproj.transform(in_crs,out_crs,easting,northing)
    return longitude,latitude


@jit(nopython=True)
def optimal_shift(indices,offset_x,offset_y,warp_band,ref_band,shift_max):
    '''Calculate optimal X and Y shift based on correlation between
    reference and warp bands. Use numba to speed up processing.
    '''

    sx,sy,px1,px2,py1,py2 = indices
    warp_sub = warp_band[py1:py2,px1:px2]
    mask = warp_sub.flatten() != -9999
    x = warp_sub.flatten()[mask]
    x_error = x-x.mean()
    x_error_sqr_sum = np.sum(x_error**2)

    opt_y,opt_x = 0,0
    opt_corr = -100

    for y_shift in range(-shift_max,shift_max):
        for x_shift in range(-shift_max,shift_max):
            lx1,lx2,ly1,ly2 = px1+offset_x+x_shift,px2+offset_x+x_shift,py1+offset_y+y_shift,py2+offset_y+y_shift
            ref_sub = ref_band[ly1:ly2,lx1:lx2]
            y = ref_sub.flatten()[mask]
            num  = np.sum(x_error*(y-y.mean()))
            denom = np.sqrt(x_error_sqr_sum*np.sum((y-y.mean())**2))
            if denom !=0:
                corr_coeff= num/denom
            else:
                corr_coeff= 0
            if corr_coeff > opt_corr:
                opt_corr = corr_coeff
                opt_y,opt_x = y_shift,x_shift
    return [sy,sx,opt_y,opt_x,opt_corr]

@ray.remote
def ray_optimal_shift(indices,offset_x,offset_y,warp_band,ref_band,shift_max):
    '''Wrapper around numba shift optimizer
    '''
    return optimal_shift(indices,offset_x,offset_y,warp_band,ref_band,shift_max)


def image_match(ref_band,warp_band,offset_x,offset_y,sensor_zn_prj,sensor_az_prj,elevation_prj,shift_max=15):
    ''' This function takes as input a single landsat band 5 (850 nm) image
    and an image to be warped and calculates the offset surface in the X and Y
    direction as a function of the sensor zenith, sensor azimuth and ground elevation.

    Args:
        ref_file (str): Landsat B5 (850nm) pathanem.
        warp_band (np.array): Projected band to be warped.
        ulx (float): Upper left easting coordinate
        uly (float): Upper left northing coordinate
        sensor_zn_prj (np.array): Projected sensor zennith array.
        sensor_az_prj (np.array): Projected sensor azimuth array.
        elevation_prj_prj (np.array): Projected elevation array.
        shift_max (int): Maximum pixel shift to test.

    Returns:
        Y and X offset model coefficients.

    '''

    logging.info('Image matching')

    if ray.is_initialized():
        ray.shutdown()
    ray.init()

    # Share arrays
    warp_band_r = ray.put(warp_band)
    ref_band_r = ray.put(ref_band)

    window = 25
    step = 5

    covar_surf = np.zeros((1+warp_band.shape[0]//step,1+warp_band.shape[1]//step,3))

    # Run correlation matching, use numba and ray to speed up
    results = []
    for sy,py1 in enumerate(range(0,warp_band.shape[0],step)):
        py2 = min(py1+window,warp_band.shape[0])
        for sx,px1 in enumerate(range(0,warp_band.shape[1],step)):
            px2 = min(px1+window,warp_band.shape[1])
            warp_sub = warp_band[py1:py2,px1:px2]
            mask = warp_sub.flatten() != -9999
            elev_sub =  elevation_prj[py1:py2,px1:px2].flatten()
            zn_sub =  sensor_zn_prj[py1:py2,px1:px2].flatten()
            az_sub =  sensor_az_prj[py1:py2,px1:px2].flatten()
            if np.sum(mask)> 100:
                covar_surf[sy,sx,:] = [zn_sub[mask].mean(),az_sub[mask].mean(),elev_sub[mask].mean()]
                indices = [sx,sy,px1,px2,py1,py2]
                results.append(ray_optimal_shift.remote(indices,
                                                        offset_x,offset_y,
                                                        warp_band_r,ref_band_r,
                                                        shift_max))
    results = ray.get(results)

    # Put results into arrays
    shift_surf = np.zeros((1+warp_band.shape[0]//step,1+warp_band.shape[1]//step,2))
    corr_surf = np.zeros((1+warp_band.shape[0]//step,1+warp_band.shape[1]//step))
    for sy,sx,opt_y,opt_x,opt_corr in results:
        shift_surf[int(sy),int(sx),:] = opt_y,opt_x
        corr_surf[int(sy),int(sx)] = opt_corr

    ray.shutdown()

    # Calculate shift slopes
    px, py = np.gradient(shift_surf[:,:,0])
    slope_y = np.sqrt(px ** 2 + py ** 2)

    px, py = np.gradient(shift_surf[:,:,1])
    slope_x = np.sqrt(px ** 2 + py ** 2)

    # Mask
    corr_thres = .4
    mask = (slope_y==0) & (corr_surf  > corr_thres) & (slope_x==0)

    # Fit X and Y offset linear models
    X = covar_surf[mask]
    y = shift_surf[:,:,1][mask]
    model_x = sm.OLS(y,sm.add_constant(X)).fit()
    logging.info('X offset model')
    logging.info(model_x.summary())

    y = shift_surf[:,:,0][mask]
    model_y = sm.OLS(y,sm.add_constant(X)).fit()
    logging.info('Y offset model')
    logging.info(model_y.summary())

    return (model_y.params,model_x.params)

def pathlength(sat_xyz,grd_xyz):
    '''Calculate pathlength from satellite
    to ground

    '''
    return np.linalg.norm(sat_xyz[:,np.newaxis]-grd_xyz,axis=0)


def sensor_view_angles(sat_enu,grd_enu):
    '''Calculates sensor zenith and azimuth angle
    in degrees

    TODO: Confirm correct results in all
    '''


    p = (sat_enu[:,np.newaxis]-grd_enu)/pathlength(sat_enu,grd_enu)
    sensor_zn = 90-np.degrees(np.arcsin(p[2]))
    sensor_az = np.degrees(np.arctan(p[0]/p[1]))

    DX,DY,DZ= grd_enu - sat_enu[:,np.newaxis]
    sensor_az[DY>0]= -sensor_az[DY>0]

    return sensor_zn,sensor_az

def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)


def get_landsat_image(longitude,latitude,month,max_cloud = 5):
    '''Given a set of coordinates and a month this function uses
    Google Earth Engine to generate a landsat scene using scenes from
    +/- 1 month of the input month
    '''

    # Get extents of image, extend to allow for shifting
    lat1 = float(np.max(latitude)) + .02
    lat2= float(np.min(latitude)) - .02

    lon1 = float(np.max(longitude)) + .02
    lon2= float(np.min(longitude)) - .02

    ee.Initialize()

    bounds  = ee.Geometry.Polygon(list([(lon1,lat1),
                                        (lon2,lat1),
                                        (lon2,lat2),
                                        (lon1,lat2),
                                        (lon1,lat1)]))

    #Retrieve Landsat 8 collection and average
    landsat8 = ee.ImageCollection("LANDSAT/LC08/C01/T1")
    landsat8_bounds = landsat8.filterBounds(bounds)
    landsat8_month = landsat8_bounds.filter(ee.Filter.calendarRange(month-1,month+1,'month'))
    landsat8_cloud = landsat8_month.filterMetadata('CLOUD_COVER','less_than',max_cloud).sort('CLOUDY_PIXEL_PERCENTAGE')
    landsat_mean = landsat8_cloud.select('B5').mean()
    latlon = ee.Image.pixelLonLat().addBands(landsat_mean)

    lats = []
    lons = []
    values = []

    #Split up area into smaller chunks to prevent exceeding max pixels
    mini_lats = np.linspace(lat2,lat1,5)
    mini_lons = np.linspace(lon2,lon1,5)

    for mini_lat1,mini_lat2 in pairwise(mini_lats):
        for mini_lon1,mini_lon2 in pairwise(mini_lons):

            mini_bounds =  ee.Geometry.Polygon(list([(mini_lon1,mini_lat1),
                                                       (mini_lon2,mini_lat1),
                                                       (mini_lon2,mini_lat2),
                                                       (mini_lon1,mini_lat2),
                                                       (mini_lon1,mini_lat1)]))
            latlon_reducer = latlon.reduceRegion(
                              reducer=ee.Reducer.toList(),
                              geometry=mini_bounds,
                              scale=30)
            values+= np.array((ee.Array(latlon_reducer.get("B5")).getInfo())).tolist()
            lats += np.array((ee.Array(latlon_reducer.get("latitude")).getInfo())).tolist()
            lons+= np.array((ee.Array(latlon_reducer.get("longitude")).getInfo())).tolist()


    lats = np.array(lats)[:,np.newaxis]
    lons= np.array(lons)[:,np.newaxis]
    values= np.array(values)[:,np.newaxis]

    easting,northing,up = dda2utm(lons,lats,[0 for x in lats],zn_dir =None)

    coords =np.concatenate([easting,northing],axis=1)

    project = Projector()
    project.create_tree(coords,np.expand_dims(easting.flatten(),axis=1).shape)

    ulx = easting.min()-100
    uly = northing.max()+100
    pixel_size = 30

    project.query_tree(ulx,uly,pixel_size)

    values_prj = project.project_band(np.expand_dims(values.flatten(),axis=1),-9999)

    return values_prj,ulx,uly



