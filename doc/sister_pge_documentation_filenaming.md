# SISTER Filenaming convention

## Basename

All output files use a common naming convention beginning with the basename which identifies the instrument, image collection date and time, product processing level, product name and version:

	INSTRUMENT_YYYYMMDDTHHMMSS_LEVEL_PRODUCT_VERSION

for example:

	AVNG_20220502T180901_L1B_RDN_001

### **INSTRUMENT**
Instrument identifier code

| Instrument | Identifier code | 
| ---|---|
| AVIRIS Classic | AVCL |
| AVIRIS Next Generation| AVNG |
| DESIS | DESIS |
| PRISMA| PRISMA |
 
### **YYYYMMDDTHHMMSS** 

UTC time at start of image collection: year, month, day, hour, minute and second.
 
### LEVEL and PRODUCT
 
Data processing level and product identifier
 
| PGE Name | LEVEL | Product | 
| ---|---|---|
| L1 preprocess| L1B | RDN |
| ISOFIT atmospheric correction | L2A | RFL |
| Spectral resampling | L2A | RSRFL |
| Reflectance correction| L2A | CORFL |
 
### Version
PGE version derived from repository release version or workflow identifier
 
 
## L1B Preprocess PGE
All outputs of the L1 PGE processor are compressed into a single tar.gz file using the following naming structure:
 
 		INSTRUMENT_YYYYMMDDTHHMMSS_L1B_RDN_VERSION.tar.gz
 	
for example:

		AVNG_20220502T180901_L1B_RDN_100.tar.gz

 	
The outputs of the PGE use the following naming convention: 

		INSTRUMENT_YYYYMMDDTHHMMSS_L1B_SUBPRODUCT_VERSION
	
| Subproduct code | Description | Example |
| ---|---|---|
| RDN | Radiance datacube | AVNG\_20220502T180901\_L1B\_RDN\_100 | 
| OBS | Observables datacube | AVNG\_20220502T180901\_L1B\_OBS\_100 |
| LOC | Location datacube | AVNG\_20220502T180901\_L1B\_LOC\_100 | 


## L2A ISOFIT atmospheric correction
All outputs of the L2 atmospheric correction are compressed into a single tar.gz file using the following naming structure:
 
 	 	INSTRUMENT_YYYYMMDDTHHMMSS_L2A_RFL_VERSION.tar.gz
for example:

		AVNG_20220502T180901_L2A_RFL_100.tar.gz

The outputs of the PGE use the following naming convention: 

		INSTRUMENT_YYYYMMDDTHHMMSS_L2A_SUBPRODUCT_VERSION
	 	 
| Subproduct code | Description | Example | 
| ---|---|---|
| RFL | Reflectance datacube | AVNG\_20220502T180901\_L2A\_RFL\_100 |
| UNC | Uncertainty datacube | AVNG\_20220502T180901\_L2A\_UNC\_100 | 
| STATE | State datacube | AVNG\_20220502T180901\_L2A\_STATE\_100 | 
| SEG | Segment datacube | AVNG\_20220502T180901\_L2A\_SEG\_100 |


## L2A Spectral resampling
All outputs of the L2 atmospheric correction are compressed into a single tar.gz file using the following naming structure:
 
 	 	INSTRUMENT_YYYYMMDDTHHMMSS_L2A_RSRFL_VERSION.tar.gz
 	 	
for example:

		AVNG_20220502T180901_L2A_RSRFL_100.tar.gz

The outputs of the PGE use the following naming convention: 
 
		INSTRUMENT_YYYYMMDDTHHMMSS_L2A_SUBPRODUCT_VERSION
		
| Subproduct code | Description |  Example | 
| ---|---|---|
| RSRFL | Resampled reflectance datacube | AVNG\_20220502T180901\_L2A\_RSRFL\_100 |
| RSUNC | Resampled uncertainty datacube | AVNG\_20220502T180901\_L2A\_RSUNC\_100 |


## L2A Reflectance correction
All outputs of the L2a reflectance correction are compressed into a single tar.gz file using the following naming structure:
 
 	 	INSTRUMENT_YYYYMMDDTHHMMSS_L2A_CORFL_VERSION.tar.gz
 	 	
example:

		AVNG_20220502T180901_L2A_CORFL_100.tar.gz	

The outputs of the PGE use the following naming convention: 

	INSTRUMENT_YYYYMMDDTHHMMSS_L2A_SUBPRODUCT_VERSION
	

| Subproduct code | Description | Example | 
| ---|---|---|
| CORFL | Corrected reflectance datacube | AVNG\_20220502T180901\_L1B\_CORFL\_100 |








 
 
 