# SISTER L1 Preprocess PGE Documentation

## Description

The L1 preprocess PGE takes as input imaging spectroscopy datasets in their native formats and converts them to a common set of file types for downstream processing. Currently, preprocessing is supported for four sensors:
 
- AVIRIS Classic
- AVIRIS Next Generation
- DESIS 
- PRISMA

Processing steps vary by sensor and are descibed below.
 
### AVIRIS Classic (AVCL)

Input data for AVCL is a tar'ed and gzipped archived containing multiple ENVI formated datasets, including radiance, geolocation and geometry data.

File example:

	f080709t01p00r15.tar.gz

Preprocessing of AVCL data includes application of radiance gains to generate radiance data in the units of microwatts per centimeter squared per steradian (μW/cm<sup>2</sup>/sr). Preprocessing of AVCL also includes optional spatial resampling. Spatial resampling is performed by aggregationg and averaging pixels to the closest resolution to the target resolution. For example, for a target pixel size of 30m and source pixel size of 16m, pixel will be averaged in 2x2 blocks of pixels for an output resolution of 30m.
	
	
### AVIRIS Next Generation (AVNG)

Input data for AVNG is a tar'ed and gzipped archived containing multiple ENVI formated datasets, including radiance, geolocation and geometry data.

File example:

	ang20191027t204454.tar.gz
	
Preprocessing of AVNG also includes optional spatial resampling. Spatial resampling is performed by aggregationg and averaging pixels to the closest resolution to the target resolution. For example, for a target pixel size of 30m and source pixel size of 5.6m, pixel will be averaged in 5x5 blocks of pixels for an output resolution of 28m.

### DESIS (DESIS)
DESIS L1C radiance data are provided by the German Space Agency (DLR) and Teledyne as a zipped archive containing radiance data in a GeoTIFF file along with metadata in an XML file. 

Example file:

	DESIS-HSI-L1C-DT0700655132_004-20220317T070333-V0215.zip
			
Provided band gains and offsets are using to convert radiance data to physical units of μW/cm<sup>2</sup>/sr. Per-pixel sensor geometry data is not provided, instead a scene mean value is included in the metadata, this value is assigned to all pixels in the image. Per-pixel solar geometry is calculated at the start time of image collection. An elevation dataset is not provided with DESIS imagery and is generated during runtime. Copernicus DEM tiles that overlap the DESIS image extent are downloaded from AWS servers (['https://copernicus-dem-30m.s3.amazonaws.com/']('https://copernicus-dem-30m.s3.amazonaws.com/')), mosiacked and clipped to the geographic extent of the input dataset.

### PRISMA (PRISMA)

PRISMA L1 radiance data are provided by the Italian Space Agency (ASI) as a zipped HDF file containing unprojected radiance, pixel geolocation data and sensor and solar geometry data. 

Example file:

	PRS_L1_STD_OFFL_20210204092812_20210204092816_0001.zip

Prior to data export a series of correction routines are applied to the dataset to improve geometric registration and radiometry. First a smile correction is applied by resampling the radiance data using a precalculated wavelength center array, next a pseudo flat field correction is applied to the radiance data using a precalculated array of radiometric adjustment coefficients. Using the input Landsat image as a reference image pixel coordinates are then adjusted using an image matching algorithm. An elevation dataset is not provided with PRISMA imagery and is generated during runtime. Copernicus DEM tiles that overlap the PRISMA image extent are downloaded from AWS servers (['https://copernicus-dem-30m.s3.amazonaws.com/']('https://copernicus-dem-30m.s3.amazonaws.com/')), mosiacked and clipped to the geographic extent of the input dataset. Finally all datasets are projected in the appropriate WGS84 UTM zone at a spatial resolution of 30m.


	
## PGE Arguments

In addition to required MAAP job submission arguments the L1 preprocess PGE also takes the following arguments:


|Argument| Type |  Description | Default|
|---|---|---|---|
| l1_granule| string |URL to input L1 dataset granule| -|
| landsat | string |URL to composite Landsat reference image, required only for PRISMA datasets| 'None'|
| resolution| integer |Ouput pixel resolution in meters. Applies only to AVIRIS sensors, PRISMA and DESIS datasets are output at native resolution of 30m | 30|


## Outputs

The L1 preprocess PGE exports 3 ENVI formatted datacubes along with their associated header files. Output files follow the AVIRIS Next Generation formatting structure and include the following: 


|Product name| Description |  Units |
|---|---|---|---|
| \*RDN\*| ENVI Radiance |μW/cm<sup>2</sup>/sr|
| \*RDN\*.hdr| ENVI radiance header file  | - |
| \*LOC\*| ENVI Location datacube |-|
|| 1. WGS-84 longitude |decimal degrees|
|| 2. WGS-84 latitude |decimal degrees|
|| 3. Ground elevation |meters|
| \*LOC\*.hdr| ENVI location header file  | - |
| \*OBS\*| ENVI Observation datacube |-|
|| 1. path length |meters|
|| 2. to-sensor-azimuth |0 to 360 degrees clockwise from N|
|| 3. to-sensor-zenith |0 to 90 degrees from zenith|
|| 4. to-sun-azimuth |0 to 360 degrees clockwise from N|
|| 5. to-sun-zenith |0 to 90 degrees from zenith|
|| 6. solar phase |degrees between to-sensor and to-sun vectors in principal plane|
|| 7. slope |decimal degrees|
|| 8. aspect |0 to 360 degrees clockwise from N|
|| 9. cosine i |unitless|
|| 10. UTC time |decimal hours|
|| 11. Earth-sun distance |astronomical unit|
| \*OBS*.hdr| ENVI observable header file  | - |


File and band descriptions taken directly from [AVIRIS NG Data Product Readme]
(https://avirisng.jpl.nasa.gov/dataportal/ANG_L1B_L2_Data_Product_Readme_v02.txt)



## Examples

**PRISMA**	
	
	l1p_job_response = maap.submitJob(algo_id="sister-l1_preprocess",
										    version="1.0.0",
										    l1_granule= '../PRS_L1_STD_OFFL_20200917091806_20200917091810_0001.zip',
										    landsat='.../PRS_20200917091806_20200917091810_0001_landsat.tar.gz',
										    publish_to_cmr=False,
										    cmr_metadata={},
										    queue="sister-job_worker-32gb",
										    identifier="l1_preprocess_PRISMA_20200917T091806")
		 
**AVCL, AVNG, DESIS** 

Landsat argument not required, will default to 'None' 
 
 	l1p_job_response = maap.submitJob(algo_id="sister-l1_preprocess",
										    version="1.0.0",
										    l1_granule= '../ang20170827t175432.tar',
										    resolution=30,
										    publish_to_cmr=False,
										    cmr_metadata={},
										    queue="sister-job_worker-32gb",
										    identifier="l1_preprocess_AVNG_20170827T175432")



 
 





