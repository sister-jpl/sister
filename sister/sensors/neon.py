import requests
import os
import tarfile
import numpy as np
from hytools.io.envi import WriteENVI
from hytools.io.envi import envi_header_dict
import h5py
from skimage.util import view_as_blocks
import pyproj

def get_neon_radiance(site,date,line,out_dir):
    '''Given a site, date and line number this scripts retreives the line via
    the NEON API data portal'''

    product_request = requests.get('https://data.neonscience.org/api/v0/data/DP1.30008.001/%s/%s-%s' % (site,year,month)).json()
    files= product_request['data']['files']

    # Cycle until matching file is found
    for file in files:
        if 'L%03d' % line in file['name']:
            # Download image to disk
            url = file['url']
            filename = '%s/%s' % (out_dir,file['name'])
            if not os.path.isfile(filename):
                with requests.get(url, stream=True) as r:
                    print("Downloading %s" % file['name'])
                    r.raise_for_status()
                    with open(filename, 'wb') as f:
                        for chunk in r.iter_content(chunk_size=int(1E8)):
                            f.write(chunk)
    return filename

def neon_to_envi(filename,resolution = 1,compress=True):
    '''Convert a NEON HDF radiance file to ENVI formated
    image along with observables and location data cube
    '''

    # Load HDF file
    hdf_obj = h5py.File(filename,'r')

    key = [key for key in hdf_obj.keys()][0]

    rad_dec = hdf_obj[key]['Radiance']['RadianceDecimalPart']
    rad_int =hdf_obj[key]['Radiance']['RadianceIntegerPart']
    obs = hdf_obj[key]['Radiance']['Metadata']['Ancillary_Rasters']['OBS_Data']
    igm = hdf_obj[key]['Radiance']['Metadata']['Ancillary_Rasters']['IGM_Data']

    wavelengths =hdf_obj[key]['Radiance']['Metadata']['Spectral_Data']['Wavelength'][:].tolist()
    fwhm =hdf_obj[key]['Radiance']['Metadata']['Spectral_Data']['FWHM'][:].tolist()

    map_info = hdf_obj[key]['Radiance']['Metadata']['Coordinate_System']['Map_Info'][()].decode("utf-8").split(',')
    epsg = hdf_obj[key]['Radiance']['Metadata']['Coordinate_System']['EPSG Code'][()].decode("utf-8")

    new_lines = rad_dec.shape[0]//resolution
    new_cols = rad_dec.shape[1]//resolution

    map_info[5] =  resolution
    map_info[6] =  resolution

    map_info = [str(info).strip() for info in map_info]

    # Export radiance
    ####################################################
    rad_dict = envi_header_dict()
    rad_dict['lines']= new_lines
    rad_dict['samples']= new_cols
    rad_dict['bands']= rad_dec.shape[2]
    rad_dict['wavelength']= wavelengths
    rad_dict['fwhm']= fwhm
    rad_dict['interleave']= 'bil'
    rad_dict['data type'] = 4
    rad_dict['wavelength units'] = "nanometers"
    rad_dict['byte order'] = 0
    rad_dict['data ignore value']= -9999
    rad_dict['map info'] =map_info

    output_name = filename.replace('radiance.h5','rad')
    writer = WriteENVI(output_name,rad_dict)

    for band_num in range(rad_dict['bands']):
        print(band_num)
        band_int = rad_int[:,:,band_num].astype(float)
        band_dec = rad_dec[:,:,band_num]/50000
        band = band_int + band_dec
        band[band_int==255] = np.nan
        band = band[:new_lines*resolution,:new_cols*resolution]
        band  = view_as_blocks(band, (resolution,resolution)).mean(axis=(2,3))
        band[np.isnan(band)] = -9999
        writer.write_band(band,band_num)

    # Export observables
    ####################################################
    obs_dict = envi_header_dict()
    obs_dict['band_names'] = ['path length', 'to-sensor azimuth',
                                'to-sensor zenith','to-sun azimuth',
                                  'to-sun zenith','phase', 'slope',
                                  'aspect', 'cosine i','UTC time']
    obs_dict['data type'] = 4
    obs_dict['lines']= new_lines
    obs_dict['samples']= new_cols
    obs_dict['bands']= 10
    obs_dict['fwhm']= fwhm
    obs_dict['interleave']= 'bil'
    obs_dict['data type'] = 4
    obs_dict['byte order'] = 0
    obs_dict['data ignore value']= -9999
    obs_dict['map info'] =map_info

    output_name = filename.replace('radiance.h5','obs')
    writer = WriteENVI(output_name,obs_dict)

    for band_num in range(obs_dict['bands']):
        print(band_num)
        band = obs[:,:,band_num]
        band[band==-9999] = np.nan
        band = band[:new_lines*resolution,:new_cols*resolution]
        band  = view_as_blocks(band, (resolution,resolution)).mean(axis=(2,3))
        band[np.isnan(band)] = -9999
        writer.write_band(band,band_num)

    # Export location datacube (lon,lat,elevation)
    ####################################################
    loc_dict = envi_header_dict()
    loc_dict['band_names'] = ['Longitude','Latitude','Elevation']
    loc_dict['data type'] = 4
    loc_dict['lines']= new_lines
    loc_dict['samples']= new_cols
    loc_dict['bands']= 3
    loc_dict['fwhm']= fwhm
    loc_dict['interleave']= 'bil'
    loc_dict['data type'] = 4
    loc_dict['byte order'] = 0
    loc_dict['data ignore value']= -9999
    loc_dict['map info'] =map_info

    output_name = filename.replace('radiance.h5','loc')
    writer = WriteENVI(output_name,loc_dict)

    in_proj = pyproj.Proj("+init=EPSG:%s" % epsg)
    out_proj= pyproj.Proj("+init=EPSG:4326")

    longitude,latitude = pyproj.transform(in_proj, out_proj,
                                          igm[:,:,0],
                                          igm[:,:,1])

    elevation = igm[:,:,2]
    mask = elevation==-9999

    for band_num,band in enumerate([longitude,latitude,elevation]):
        print(band_num)
        band[mask] = np.nan
        band = band[:new_lines*resolution,:new_cols*resolution]
        band  = view_as_blocks(band, (resolution,resolution)).mean(axis=(2,3))
        band[np.isnan(band)] = -9999
        writer.write_band(band,band_num)

    # if compress:
    #     with tarfile.open(filename.replace('_radiance.h5','.tar.gz'), "w:gz") as tar:
    #         tar.add(filename.replace('radiance.h5','%s.hdr' % suffix),
    #                 arcname=os.path.basename(filename.replace('radiance.h5','%s.hdr' % suffix)))

    #     for suffix in ['rad','loc','obs']:















