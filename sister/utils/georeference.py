#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SISTER
Space-based Imaging Spectroscopy and Thermal PathfindER
Author: Adam Chlus
"""
import os
import glob
import tarfile
import numpy as np
import hytools as ht
from hytools.io.envi import WriteENVI,envi_header_dict
from rtree import index
from scipy.spatial import cKDTree


def utm_zone(longitude, latitude):
    ''' Returns UTM zone and direction
    '''
    zone = int(np.ceil((longitude.min()+180)/6))
    if latitude.mean() >0:
        direction = 'N'
    else:
        direction = 'S'
    return zone,direction


def dd2utm(longitude,latitude):
    '''Convert coordinates in decimal degrees int
    UTM eastings and northings
    '''
    zone,direction = utm_zone(longitude, latitude)

    if direction == 'N':
        epsg_dir = 6
    else:
        epsg_dir = 7

    outPCS = pyproj.Proj("+init=EPSG:32%s%02d" % (epsg_dir,zone))
    inGCS= pyproj.Proj("+init=EPSG:4326")

    # Convert to easting and northing,
    easting,northing = pyproj.transform(inGCS,outPCS,longitude,latitude)
    return easting,northing


def utm2dd(easting,northing,zone,direction):
    '''Convert coordinates in UTM eastings and northings into
        decimal degrees
    '''
    if direction == 'N':
        epsg_dir = 6
    else:
        epsg_dir = 7

    inPCS = pyproj.Proj("+init=EPSG:32%s%02d" % (epsg_dir,zone))
    outGCS= pyproj.Proj("+init=EPSG:4326")

    # Convert to easting and northing,
    longitude,latitude= pyproj.transform(inPCS,outGCS,easting,northing)
    return longitude,latitude



class Projector():
    """Projector class"""

    def __init__(self):
        """Constructor method
        """
        self.tree = None
        self.indices = None
        self.pixel_size = None
        self.input_shape = None
        self.ouput_shape = None
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

        self.output_shape = (int(lines),int(columns))
        int_north,int_east = np.indices(self.output_shape)
        int_east = (int_east*half_pix + ulx).flatten()
        int_north = (uly-int_north*half_pix).flatten()
        int_north= np.expand_dims(int_north,axis=1)
        int_east= np.expand_dims(int_east,axis=1)
        dest_points = np.concatenate([int_east,int_north],axis=1)

        distances,self.indices =  self.tree.query(dest_points,k=1)
        self.indices = np.unravel_index(self.indices,self.input_shape)
        distances =  distances.reshape(self.output_shape)
        self.mask = distances > 2*pixel_size

    def project(self,band,no_data):
        band = np.copy(band[self.indices[0],self.indices[1]].reshape(self.output_shape))
        band[self.mask] = np.nan
        band = np.nanmean(view_as_blocks(band, (2,2)),axis=(2,3))
        band[np.isnan(band)] = no_data
        return band





def

















