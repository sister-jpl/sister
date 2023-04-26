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
from itertools import tee
import logging
import numpy as np
from scipy.spatial import cKDTree
from numba import jit
import statsmodels.api as sm
import ray
import pyproj
from skimage.util import view_as_blocks
import hytools as ht
from hytools.io.envi import WriteENVI
from scipy.stats import circmean
from scipy.ndimage import binary_erosion


#Temporary fix.....
try:
    import ee
except:
    print('Unable to import Google Earth Engine API')

class Projector():
    """Nearest Neighbor Projector class"""

    def __init__(self):
        self.tree = None
        self.indices = None
        self.pixel_size = None
        self.input_shape = None
        self.output_shape = None

        self.mask = None

    def create_tree(self,coords,input_shape):
        self.input_shape = input_shape
        self.tree = cKDTree(coords,balanced_tree= False)

    def query_tree(self,ulx,uly,pixel_size):
        self.pixel_size = pixel_size

        lines = (uly-self.tree.mins[1])//pixel_size
        if lines%2 !=0:
            lines+=20
        else:
            lines+=19

        columns = (self.tree.maxes[0]-ulx)//pixel_size
        if columns%2 !=0:
            columns+=20
        else:
            columns+=19

        self.output_shape = (int(lines),int(columns))
        north,east = np.indices(self.output_shape)
        east = (east*pixel_size + ulx).flatten()
        north = (uly-north*pixel_size).flatten()
        north= np.expand_dims(north,axis=1)
        east= np.expand_dims(east,axis=1)
        dest_points = np.concatenate([east,north],axis=1)

        distances,self.indices =  self.tree.query(dest_points,k=1)
        self.indices = np.unravel_index(self.indices,self.input_shape)
        self.distances =  distances.reshape(self.output_shape)
        self.mask = ~binary_erosion(self.distances < pixel_size,
                                   structure=np.ones((5,5)),
                                   border_value=1)

    def project_band(self,band,no_data):
        band = np.copy(band[self.indices[0],self.indices[1]].reshape(self.output_shape))
        band[self.mask] = no_data
        return band


def utm_zone(longitude, latitude):
    ''' Returns UTM zone and direction
    '''
    zone = int(np.ceil((longitude.min()+180)/6))
    if latitude.mean() >0:
        direction = 'North'
    else:
        direction = 'South'
    return zone,direction


def dda2utm(longitude,latitude,altitude,zn_dir =None):
    '''Convert coordinates in decimal degrees into
    UTM eastings and northings
    '''

    if zn_dir:
        zone,direction =zn_dir
    else:
        zone,direction = utm_zone(longitude, latitude)

    if direction == 'North':
        epsg_dir = 6
    else:
        epsg_dir = 7

    out_crs = pyproj.Proj("epsg:32%s%02d" % (epsg_dir,zone))
    in_crs= pyproj.Proj("epsg:4326")

    easting,northing,up = pyproj.transform(in_crs,out_crs,
                                           latitude,
                                           longitude,
                                           altitude)
    return easting,northing,up


def dda2ecef(longitude,latitude,altitude):
    '''Convert longitude,latitude and altitude to Earth Centered Earth Fixed
    coordinates (x,y,x)

    '''
    in_crs = pyproj.Proj(proj='longlat', ellps='WGS84', datum='WGS84')
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
    if direction.startswith('N'):
        epsg_dir = 6
    else:
        epsg_dir = 7

    in_crs = pyproj.Proj("epsg:32%s%02d" % (epsg_dir,zone))
    out_crs= pyproj.Proj("epsg:4326")
    latitude,longitude = pyproj.transform(in_crs,out_crs,easting,northing)
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
    ray.init(num_cpus=13)
    logging.info('Ray intitialized')

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

    logging.info('Correlation matching complete')
    results = ray.get(results)
    logging.info('Ray results retrieved')

    # Put results into arrays
    shift_surf = np.zeros((1+warp_band.shape[0]//step,1+warp_band.shape[1]//step,2))
    corr_surf = np.zeros((1+warp_band.shape[0]//step,1+warp_band.shape[1]//step))
    for sy,sx,opt_y,opt_x,opt_corr in results:
        shift_surf[int(sy),int(sx),:] = opt_y,opt_x
        corr_surf[int(sy),int(sx)] = opt_corr

    ray.shutdown()
    logging.info('Ray shutdown')

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
    return np.linalg.norm(sat_xyz[:,:,np.newaxis]-grd_xyz,axis=0)


def sensor_view_angles(sat_enu,grd_enu):
    '''Calculates sensor zenith and azimuth angle
    in degrees
    '''
    grd_enu = np.rot90(grd_enu,k=3,axes=(1,2))

    p = (sat_enu[:,:,np.newaxis]-grd_enu)/pathlength(sat_enu,grd_enu)
    sensor_zn = 90-np.degrees(np.arcsin(p[2]))
    sensor_az = np.degrees(np.arctan(p[0]/p[1]))
    sensor_az[sensor_az<0]+=360

    DX,DY,DZ= grd_enu - sat_enu[:,:,np.newaxis]
    sensor_az[(DX>0) & (DY>0)]= 180+sensor_az[(DX>0) & (DY>0)]
    sensor_az[(DX<0) & (DY>=0)] = sensor_az[(DX<0) & (DY>=0)] - 180

    sensor_zn =  np.rot90(sensor_zn)
    sensor_az =  np.rot90(sensor_az)

    return sensor_zn,sensor_az

def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)

def get_landsat_image(longitude,latitude,month,max_cloud = 5,band=5,project = True,month_window = 1):
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
    month_min = max(month-month_window,1)
    month_max = min(month+month_window,12)

    landsat8_month = landsat8_bounds.filter(ee.Filter.calendarRange(month_min,month_max,'month'))
    landsat8_cloud = landsat8_month.filterMetadata('CLOUD_COVER','less_than',max_cloud).sort('CLOUDY_PIXEL_PERCENTAGE')
    landsat_mean = landsat8_cloud.select('B%s' % band).mean()
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
            values+= np.array((ee.Array(latlon_reducer.get("B%s" % band)).getInfo())).tolist()
            lats += np.array((ee.Array(latlon_reducer.get("latitude")).getInfo())).tolist()
            lons+= np.array((ee.Array(latlon_reducer.get("longitude")).getInfo())).tolist()


    lats = np.array(lats)[:,np.newaxis]
    lons= np.array(lons)[:,np.newaxis]
    values= np.array(values)[:,np.newaxis]


    if project:
        project = Projector()
        easting,northing,up = dda2utm(lons,lats,[0 for x in lats],zn_dir =None)
        coords =np.concatenate([easting,northing],axis=1)
        project.create_tree(coords,np.expand_dims(easting.flatten(),axis=1).shape)

        ulx = easting.min()-100
        uly = northing.max()+100
        pixel_size = 30

        project.query_tree(ulx,uly,pixel_size)

        values_prj = project.project_band(np.expand_dims(values.flatten(),axis=1),-9999)
        return values_prj,ulx,uly


    return values,lons,lats


def rotate_coords(x_i,y_i,x_p,y_p,theta):
    '''Rotate coordinates (xi,yi) theta degrees around point (xp,yp)
    '''
    theta = np.radians(theta)
    x_r = x_p + ((x_i-x_p)*np.cos(theta)) - ((y_i-y_p)*np.sin(theta))
    y_r = y_p + ((x_i-x_p)*np.sin(theta)) + ((y_i-y_p)*np.cos(theta))
    return x_r,y_r



def resample(in_file,out_dir,resolution,kind = 'closest', verbose = True, unrotate = False):
    ''' Perform a two-step spatial resampling to . First, pixels are aggregated and
    averaged, next a nearest neighbor algorithm is used to resample images to resolution.
    '''

    out_image = out_dir + '/' + os.path.basename(in_file)
    image = ht.HyTools()
    image.read_file(in_file,'envi')

    x_ul = float(image.map_info[3])
    y_ul = float(image.map_info[4])
    pixel_res = float(image.map_info[5])

    if 'rotation' in image.map_info[-1]:
        rotation = 360 + float(image.map_info[-1].split('=')[-1])
    else:
        rotation = 0

    # Calculate rotated and unrotated coords
    y_ind,x_ind = np.indices((image.lines,image.columns))
    y_rcoord = y_ul - y_ind*pixel_res
    x_rcoord = x_ul + x_ind*pixel_res

    if unrotate:
        x_coord,y_coord = rotate_coords(x_rcoord,y_rcoord,x_ul,y_ul,rotation)
    else:
        x_coord,y_coord = x_rcoord,y_rcoord

    bin_size = int(np.round(resolution/pixel_res))
    if verbose:
        print("Aggregating every %s pixels" % bin_size)

    lines =bin_size*(image.lines//bin_size)
    columns =bin_size*(image.columns//bin_size)

    y_coord_bin = np.nanmean(view_as_blocks(y_coord[:lines,:columns],
                                     (bin_size,bin_size)),axis=(2,3))
    x_coord_bin= np.nanmean(view_as_blocks(x_coord[:lines,:columns],
                                     (bin_size,bin_size)),axis=(2,3))


    # Get extent of output array
    xmax = int(resolution * (x_coord_bin.max()//resolution)) + resolution
    ymax = int(resolution * (y_coord_bin.max()//resolution)) + resolution
    xmin = int(resolution * (x_coord_bin.min()//resolution)) - resolution
    ymin = int(resolution * (y_coord_bin.min()//resolution)) - resolution

    out_columns = int((xmax-xmin)/resolution)
    out_lines =  int((ymax-ymin)/resolution)

    #Calculate coordinates of output array
    image_shape = (out_lines,out_columns)
    y_coord_out,x_coord_out = np.indices(image_shape)*resolution
    y_coord_out = ymax -  y_coord_out
    x_coord_out = xmin + x_coord_out

    #Create tree to convert pixels to geolocated pixels
    src_points =np.concatenate([np.expand_dims(x_coord_bin.flatten(),axis=1),
                                np.expand_dims(y_coord_bin.flatten(),axis=1)],axis=1)
    tree = cKDTree(src_points,balanced_tree= True)

    dst_points = np.concatenate([np.expand_dims(x_coord_out.flatten(),axis=1),
                                 np.expand_dims(y_coord_out.flatten(),axis=1)],axis=1)

    dists, indexes = tree.query(dst_points,k=1)
    indices_int = np.unravel_index(indexes,x_coord_bin.shape)
    mask = dists.reshape(image_shape) >resolution

    out_header = image.get_header()
    out_header['lines'] = out_lines
    out_header['samples'] = out_columns
    out_header['map info'][3] = str(xmin)
    out_header['map info'][4] = str(ymax)
    out_header['map info'][5:7] = resolution,resolution
    if unrotate:
        out_header['map info'][-1] ='rotation=0.0000'
    out_header['byte order'] = 0
    out_header['data ignore value'] = image.no_data

    writer = WriteENVI(out_image,out_header)
    iterator =image.iterate(by = 'band')

    while not iterator.complete:
        if verbose &  (iterator.current_band%10 == 0):
            print("%s/%s" % (iterator.current_band,image.bands))
        band = np.copy(iterator.read_next()).astype(float)
        band[~image.mask['no_data']] = np.nan
        bins  = view_as_blocks(band[:lines,:columns],(bin_size,bin_size))

        if (iterator.current_band in [1,2,3,4,7]) and ('OBS' in image.base_name):
            bins = np.radians(bins)
            band = circmean(bins,axis=2,nan_policy = 'omit')
            band = circmean(band,axis=2,nan_policy = 'omit')
            band = np.degrees(band)
        else:
            band = np.nanmean(bins,axis=(2,3))

        band = band[indices_int[0],indices_int[1]].reshape(image_shape)
        band[mask]= image.no_data
        band[np.isnan(band)]= image.no_data
        writer.write_band(band,iterator.current_band)
