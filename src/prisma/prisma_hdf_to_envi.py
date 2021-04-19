# -*- coding: utf-8 -*-
''' prisma_radiance_export.py
'''
import argparse
import datetime as dt
import glob
import os
import shutil
import zipfile
import h5py
import hytools as ht
from hytools.io.envi import WriteENVI,envi_header_dict
from hytools.misc import progbar
from hytools.topo.topo import calc_cosine_i
import numpy as np
from rtree import index
from scipy.interpolate import interp1d
from scipy.spatial import cKDTree
import pyproj
from pysolar import solar
from skimage.util import view_as_blocks
import pandas as pd

def main():
    '''
    This function exports a PRISMA L1 radiance products to an ENVI formatted
    binary file along with ancillary datasets including location data and geometry.
    This scripts requires L1 radiance, L2c reflectance and SRTM elevation products for a scene.

    The VNIR and SWIR dectector cubes are stitched together at 950nm.

    base_name: Scene basename (ex: 20200205001044_20200205001049_0001)
    zip_dir: Directory of zipped L1 and L2c data products
    srtm_dir: Directory of zipped SRTM tiles
    out_dir: Output directory
    temp_dir: Temporary file directory
    sml (optional): Path name of smile correction surface
    sml (optional): Georeference resolution, if not specifiec data will not be geolocated,
                    30 must be a factor of output resolution.
    rfl (optional): Export ASI reflectance product
    '''
    parser = argparse.ArgumentParser(description = "Converts PRISMA data ENVI format")
    parser.add_argument('base_name',help="Image basename", type = str)
    parser.add_argument('zip_dir',help="Directory path of L1 radiance data", type = str)
    parser.add_argument('srtm_dir',help="Directory path of SRTM tiles", type = str)
    parser.add_argument('out_dir',help="Output directory", type = str)
    parser.add_argument('temp_dir',help="Temporary data directory", type = str)
    parser.add_argument("-sml", help="Path to smile surface", required=False,type = str)
    parser.add_argument("-geo", help="Georeference resolution, if not specifiec data will not be geolocated", required=False, type = int)
    parser.add_argument("-rfl", help="Export ASI reflectance product", required=False, action='store_true')

    args = parser.parse_args()
    srtm_dir = args.srtm_dir
    zip_dir  = args.zip_dir
    base_name = args.base_name
    temp_dir =  "%s/temp_%s/" % (args.temp_dir,base_name)
    out_dir = '%s/PRS_%s/' %  (args.out_dir,base_name)
    smile_correct = False

    # # ##Testing
    # parser = argparse.ArgumentParser(description = "Converts PRISMA data ENVI format")
    # args = parser.parse_args()
    # srtm_dir = '/data2/srtm/'
    # l1_dir  = '/data2/prisma/l1/'
    # l2_dir = '/data2/prisma/l2c/'
    # base_name = "20200621005743_20200621005747_0001"
    # # out_dir = '/data2/prisma/envi/PRS_%s/' % base_name
    # temp_dir = '/data2/temp/prisma/'
    # args.smile =True
    # args.geo =True

    for file in glob.glob("%s*%s*" % (zip_dir,base_name)):
        zip_base  =os.path.basename(file)
        print('Unzipping %s' % zip_base)
        with zipfile.ZipFile(file,'r') as zipped:
            zipped.extractall(temp_dir)
    l1  = '%sPRS_L1_STD_OFFL_%s.he5' % (temp_dir,base_name)
    l2  = '%sPRS_L2C_STD_%s.he5' % (temp_dir,base_name)

    file_suffixes =  ['rad','loc','obs']

    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)

    if not os.path.isdir(temp_dir):
        os.mkdir(temp_dir)

    if args.sml:
        smile = ht.HyTools()
        smile.read_file(args.sml, 'envi')
        shift_surf_smooth = smile.get_band(0)
        smile_correct = True

    for product in [l1,l2]:
        l_obj = h5py.File(product,'r')
        subdir = [x for x in l_obj['HDFEOS']["SWATHS"].keys() if 'HCO' in x][0]

        if 'L2C' in subdir:
            if not args.rfl:
                continue
            print('\nExporting reflectance data')
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
        vnir_temp = '%sNA-%s_%s_vnir' % (temp_dir,base_name,measurement)

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
        swir_temp = '%sNA-%s_%s_swir' % (temp_dir,base_name,measurement)

        writer = WriteENVI(swir_temp,rad_dict )
        writer.write_chunk(np.moveaxis(swir_data[:,:,:],1,2), 0,0)

        #Export merged radiance/reflectance datacube
        vnir_waves = np.flip(vnir_waves[6:])
        swir_waves = np.flip(swir_waves[:-3])

        vnir_fwhm = np.flip(vnir_fwhm[6:])
        swir_fwhm = np.flip(swir_fwhm[:-3])

        vnir_obj = ht.HyTools()
        vnir_obj.read_file(vnir_temp, 'envi')

        swir_obj = ht.HyTools()
        swir_obj.read_file(swir_temp, 'envi')

        if args.geo:
            output_name = '%sNA-%s_%s' % (temp_dir,base_name,measurement)
        else:
            output_name = '%sNA-%s_%s' % (out_dir,base_name,measurement)

        rad_dict  = envi_header_dict()
        rad_dict ['lines']= vnir_obj.lines
        rad_dict ['samples']=vnir_obj.columns
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

            if (measurement == 'rad') & smile_correct:
                vnir_interpolator = interp1d(shift_surf_smooth[iterator_v.current_line,:60],
                                               chunk_v,fill_value = "extrapolate",kind='cubic')
                chunk_v = vnir_interpolator(vnir_waves)
                swir_interpolator = interp1d(shift_surf_smooth[iterator_v.current_line,60:],
                                               chunk_s,fill_value = "extrapolate",kind='cubic')
                chunk_s = swir_interpolator(swir_waves)
            line = np.concatenate([chunk_v,chunk_s],axis=1)/1000.
            writer.write_line(line, iterator_v.current_line)
            progbar(iterator_v.current_line,1000)

    # Location datacube
    #################################
    geo =  l_obj['HDFEOS']["SWATHS"][subdir]['Geolocation Fields']

    # Get extents of image
    lon_min = geo['Longitude'][:,:].min()
    lon_max = geo['Longitude'][:,:].max()
    lat_min = geo['Latitude'][:,:].min()
    lat_max = geo['Latitude'][:,:].max()

    # Create a simple spatial index to find intersecting tiles
    idx = index.Index(properties=index.Property())
    tiles = glob.glob(srtm_dir + '*.zip')
    for i, tile in enumerate(tiles):
        lat,lon = os.path.basename(tile).split('_')[:2]
        if 'w' in lon:
            lon = -1*float(lon[1:])
        else:
            lon = float(lon[1:])
        if 's' in lat:
            lat = -1*float(lat[1:])
        else:
            lat = float(lat[1:])
        idx.insert(i,(lon,lat,lon+1,lat+1))

    tiles_inter = [tiles[n] for n in idx.intersection((lon_min, lat_min, lon_max, lat_max))]

    if len(tiles_inter) > 0:
        print("\nIntersecting SRTM tiles:")
        for tile in tiles_inter:
            print('\t' + os.path.basename(tile))
            with zipfile.ZipFile(tile, 'r') as zip_ref:
                zip_ref.extractall(temp_dir)

        print('Merging SRTM tiles')
        dem  = '%stemp_dem' % temp_dir
        os.system('gdal_merge.py -o %s -of ENVI %s*bil' % (dem,temp_dir))

        dem_obj = ht.HyTools()
        dem_obj.read_file(dem, 'envi')

        ulx = float(dem_obj.map_info[3])
        uly = float(dem_obj.map_info[4])
        pix = float(dem_obj.map_info[5])

        dem_lat,dem_lon = np.indices((dem_obj.lines,dem_obj.columns))

        dem_xl = int((lon_min-ulx)//pix)
        dem_xr = int((lon_max-ulx)//pix)
        dem_yu = int((uly-lat_max)//pix)
        dem_yd = int((uly-lat_min)//pix)

        dem_subset = dem_obj.get_chunk(dem_xl,dem_xr,dem_yu,dem_yd)
        dem_lat,dem_lon = np.indices(dem_subset.shape[:2])
        dem_lat = (lat_max- dem_lat*pix).flatten()
        dem_lon = (lon_min+ dem_lon*pix).flatten()

        #Create spatial index and sample
        src_points =np.concatenate([np.expand_dims(dem_lon,axis=1),np.expand_dims(dem_lat,axis=1)],axis=1)
        tree = cKDTree(src_points,balanced_tree= False)

        dst_points = np.concatenate([geo['Longitude'][:,:].flatten()[:,np.newaxis],
                                     geo['Latitude'][:,:].flatten()[:,np.newaxis]],
                                     axis=1)

        dists, indexes = tree.query(dst_points,k=1)
        indices_int = np.unravel_index(indexes,(dem_subset.shape[0],
                                                dem_subset.shape[1]))
        prisma_dem = dem_subset[indices_int[0],indices_int[1]].reshape(geo['Longitude'][:,:].shape)

    else:
        elevation = float(input("\nNo overlapping SRTM tiles found, enter constant elevation for scene (m): "))
        prisma_dem = np.ones(geo['Longitude'][:,:].shape) * elevation
        prisma_dem =  prisma_dem + (np.arange(prisma_dem.shape[0])*0.0001)

    #Set negative elevations to 0
    if np.sum(prisma_dem<0) > 0:
        print('Elevations below sea level found, setting to 0m')
        prisma_dem[prisma_dem<0] =0

    loc_dict  = envi_header_dict()
    loc_dict['lines']= vnir_data.shape[0]
    loc_dict['samples']= vnir_data.shape[2]
    loc_dict['bands']= 3
    loc_dict['interleave']= 'bsq'
    loc_dict['data type'] = 4
    loc_dict['band_names'] = ['Longitude', 'Latitude','Elevation']
    loc_dict['byte order'] = 0

    if args.geo:
        output_name = '%sNA-%s_loc' % (temp_dir,base_name)
    else:
        output_name = '%sNA-%s_loc' % (out_dir,base_name)

    print('Exporting location datacube')
    writer = WriteENVI(output_name,loc_dict)
    writer.write_band(geo['Longitude'][:,:],0)
    writer.write_band(geo['Latitude'][:,:],1)
    writer.write_band(prisma_dem,2)

    dem_dict  = envi_header_dict()
    dem_dict ['lines']= vnir_data.shape[0]
    dem_dict ['samples']= vnir_data.shape[2]
    dem_dict ['bands']= 1
    dem_dict ['interleave']= 'bsq'
    dem_dict ['data type'] = 4

    output_name = '%sPRS_dem' % temp_dir
    writer = WriteENVI(output_name,dem_dict)
    writer.write_band(prisma_dem,0)

    slope_name =  '%sPRS_slp' % temp_dir
    aspect_name =  '%sPRS_asp' % temp_dir

    print('Calculating slope')
    os.system('gdaldem slope -of ENVI %s %s'% (output_name,slope_name))
    print('Calculating aspect')
    os.system('gdaldem aspect -f ENVI %s %s' % (output_name,aspect_name))

    asp_obj = ht.HyTools()
    asp_obj.read_file(aspect_name, 'envi')

    slp_obj = ht.HyTools()
    slp_obj.read_file(slope_name, 'envi')

    # Observables datacube
    #################################
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

    # Solar geometries
    '''Solar geometry is calculated based on the mean scene acquisition time
    which varies by less than 5 seconds from start to end of the scene and is
    computationally more efficient.
    '''
    epoch = dt.datetime(2000,1, 1,)
    epoch = epoch.replace(tzinfo=dt.timezone.utc)
    time = epoch + dt.timedelta(days=np.array(geo['Time'][:]).mean())

    solar_az = solar.get_azimuth(geo['Latitude'][:,:],geo['Longitude'][:,:],time)
    solar_zn = 90-solar.get_altitude(geo['Latitude'][:,:],geo['Longitude'][:,:],time)
    sensor_az =geom['Rel_Azimuth_Angle'][:,:] +solar_az # No exactly right...need to fix minor...

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
                               geo['Longitude'][:,:],
                               geo['Latitude'][:,:],
                               prisma_dem, radians=False)

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
        pathlength += (sat_interp-grd_pos)**2
        satelite_xyz.append(sat_interp)
    pathlength = np.sqrt(pathlength)

    cosine_i = calc_cosine_i(np.radians(solar_zn),
                             np.radians(solar_az),
                             np.radians(slp_obj.get_band(0)),
                             np.radians(asp_obj.get_band(0)))

    phase =  np.arccos(np.cos(np.radians(solar_zn)))*np.cos(np.radians(geom['Observing_Angle'][:,:]))
    phase += np.sin(np.radians(solar_zn))*np.sin(np.radians(geom['Observing_Angle'][:,:]))*np.cos(np.radians(geom['Rel_Azimuth_Angle'][:,:]))

    # Export satellite position to csv
    lon,lat,alt = pyproj.transform(ecef,lla,
                                satelite_xyz[0],
                                satelite_xyz[1],
                                satelite_xyz[2],
                                radians=False)
    satellite_df = pd.DataFrame()
    satellite_df['lat'] = lat
    satellite_df['lon'] = lon
    satellite_df['alt'] = alt
    satellite_df.to_csv('%sNA-%s_satellite_loc.csv' % (out_dir,base_name))

    # Export observable datacube
    print('Exporting observables datacube')
    anc_header = envi_header_dict()
    anc_header['lines']= vnir_data.shape[0]
    anc_header['samples']= vnir_data.shape[2]
    anc_header['bands']= 10
    anc_header['interleave']= 'bil'
    anc_header['data type'] = 4
    anc_header['band_names'] = ['path length', 'to-sensor azimuth',
                                'to-sensor zenith','to-sun azimuth',
                                  'to-sun zenith','phase', 'slope',
                                  'aspect', 'cosine i','UTC time']
    anc_header['byte order'] = 0

    if args.geo:
        output_name = '%sNA-%s_obs' % (temp_dir,base_name)
    else:
        output_name = '%sNA-%s_obs' % (out_dir,base_name)

    writer = WriteENVI(output_name,anc_header)
    writer.write_band(pathlength,0)
    writer.write_band(sensor_az,1)
    writer.write_band(geom['Observing_Angle'][:,:],2)
    writer.write_band(solar_az,3)
    writer.write_band(solar_zn,4)
    writer.write_band(phase,5)
    writer.write_band(slp_obj.get_band(0),6)
    writer.write_band(asp_obj.get_band(0),7)
    writer.write_band(cosine_i,8)
    writer.write_band(utc_time,9)

    #Export georeferenced image
    if args.geo:
        location = ht.HyTools()
        location.read_file('%sNA-%s_loc' % (temp_dir,base_name), 'envi')
        #Clip bada lines at edges of detector
        lon = location.get_band(0)[2:-2,2:-2]
        lat = location.get_band(1)[2:-2,2:-2]

        # Determine UTM zone
        zone = int(np.ceil((lon.min()+180)/6))
        if lat.mean() >0:
            direction = 'N'
        else:
            direction = 'S'

        outPCS = pyproj.Proj("+init=EPSG:326%02d" % zone)
        inGCS= pyproj.Proj("+init=EPSG:4326")

        # Convert to easting and northing,
        easting,northing = pyproj.transform(inGCS,outPCS,lon,lat)
        easting = easting.flatten()
        northing = northing.flatten()

        resolution = args.geo

        if resolution%30 == 0:
            blocksize = int(resolution/30)
        else:
            print('Resolution error: 30 is not a factor of %s, setting resolution to 30m' % resolution)
            blocksize = 1

        image_shape = (int((northing.max()-northing.min())//30),int((easting.max()-easting.min())//30))
        int_north,int_east = np.indices(image_shape)
        int_east = (int_east*30 + easting.min()).flatten()
        int_north = (northing.max()-int_north*30).flatten()
        int_north= np.expand_dims(int_north,axis=1)
        int_east= np.expand_dims(int_east,axis=1)

        #Create spatial index
        src_points =np.concatenate([np.expand_dims(easting,axis=1),np.expand_dims(northing,axis=1)],axis=1)
        tree = cKDTree(src_points,balanced_tree= False)
        dst_points = np.concatenate([int_east,int_north],axis=1)
        dists, indexes = tree.query(dst_points,k=1)
        dists = dists.reshape(image_shape)
        indices_int = np.unravel_index(indexes,lon.shape)
        mask = dists > 45

        map_info_string = ['UTM', 1, 1, easting.min(), northing.max(),resolution,
                           resolution,zone,direction, 'WGS-84' , 'units=Meters']

        out_cols = blocksize* (mask.shape[1]//blocksize)
        out_lines = blocksize* (mask.shape[0]//blocksize)

        print('Georeferencing datasets')
        total_bands = rad_dict['bands'] + anc_header['bands'] + loc_dict['bands']
        curr_band = 0
        for file in file_suffixes:
            input_name = '%sNA-%s_%s' % (temp_dir,base_name,file)
            hy_obj = ht.HyTools()
            hy_obj.read_file(input_name, 'envi')
            iterator =hy_obj.iterate(by = 'band')

            out_header = hy_obj.get_header()
            out_header ['lines']= int(out_lines/blocksize)
            out_header ['samples']= int(out_cols/blocksize)
            out_header['data ignore value'] = -9999
            out_header['map info'] = map_info_string

            output_name = '%sNA-%s_%s_geo' % (out_dir,base_name,file)
            writer = WriteENVI(output_name,out_header)

            while not iterator.complete:
                band = iterator.read_next()[2:-2,2:-2]
                band = np.copy(band[indices_int[0],indices_int[1]].reshape(mask.shape))
                band[mask] = np.nan
                band = np.nanmean(view_as_blocks(band[:out_lines,:out_cols], (blocksize,blocksize)),axis=(2,3))
                band[np.isnan(band)] = -9999
                writer.write_band(band,iterator.current_band)
                curr_band+=1
                progbar(curr_band,total_bands)

    print('\nDeleting temporary files.')
    shutil.rmtree(temp_dir)

if __name__== "__main__":
    main()


















