# SISTER Filenaming convention

## Basename

All output files use a common naming convention beginning with the basename which identifies the instrument, image collection date and time, product processing level, product name and version:

	SISTER_INSTRUMENT_YYYYMMDDTHHMMSS_LEVEL_PRODUCT_CRID

for example:

	SISTER_AVNG_20220502T180901_L1B_RDN_001

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
 
### CRID
Unique workflow identification code 
 
## L1B Preprocess PGE
All outputs of the L1 PGE processor are compressed into a single tar.gz file using the following naming structure:
 
 		SISTER_INSTRUMENT_YYYYMMDDTHHMMSS_L1B_RDN_CRID.tar.gz
 	
for example:

		SISTER_AVNG_20220502T180901_L1B_RDN_001.tar.gz

 	
The outputs of the PGE use the following naming convention: 

		SISTER_INSTRUMENT_YYYYMMDDTHHMMSS_L1B_SUBPRODUCT_CRID
	
| Subproduct code | Description | Example |
| ---|---|---|
| RDN | Radiance datacube | SISTER_AVNG\_20220502T180901\_L1B\_RDN\_001 | 
| OBS | Observables datacube | SISTER_AVNG\_20220502T180901\_L1B\_OBS\_001 |
| LOC | Location datacube | SISTER_AVNG\_20220502T180901\_L1B\_LOC\_001 | 


## L2A ISOFIT atmospheric correction
All outputs of the L2 atmospheric correction are compressed into a single tar.gz file using the following naming structure:
 
 	 	SISTER_INSTRUMENT_YYYYMMDDTHHMMSS_L2A_RFL_CRID.tar.gz
for example:

		SISTER_AVNG_20220502T180901_L2A_RFL_001.tar.gz

The outputs of the PGE use the following naming convention: 

		SISTER_INSTRUMENT_YYYYMMDDTHHMMSS_L2A_SUBPRODUCT_CRID
	 	 
| Subproduct code | Description | Example | 
| ---|---|---|
| RFL | Reflectance datacube | SISTER_AVNG\_20220502T180901\_L2A\_RFL\_001 |
| UNC | Uncertainty datacube | SISTER_AVNG\_20220502T180901\_L2A\_UNC\_001 | 
| STATE | State datacube | SISTER_AVNG\_20220502T180901\_L2A\_STATE\_001 | 
| SEG | Segment datacube | SISTER_AVNG\_20220502T180901\_L2A\_SEG\_001 |


## L2A Spectral resampling
All outputs of the L2 atmospheric correction are compressed into a single tar.gz file using the following naming structure:
 
 	 	SISTER_INSTRUMENT_YYYYMMDDTHHMMSS_L2A_RSRFL_CRID.tar.gz
 	 	
for example:

		SISTER_AVNG_20220502T180901_L2A_RSRFL_001.tar.gz

The outputs of the PGE use the following naming convention: 
 
		SISTER_INSTRUMENT_YYYYMMDDTHHMMSS_L2A_SUBPRODUCT_CRID
		
| Subproduct code | Description |  Example | 
| ---|---|---|
| RSRFL | Resampled reflectance datacube | SISTER_AVNG\_20220502T180901\_L2A\_RSRFL\_001 |
| RSUNC | Resampled uncertainty datacube | SISTER_AVNG\_20220502T180901\_L2A\_RSUNC\_001 |


## L2A Reflectance correction
All outputs of the L2a reflectance correction are compressed into a single tar.gz file using the following naming structure:
 
 	 	SISTER_INSTRUMENT_YYYYMMDDTHHMMSS_L2A_CORFL_CRID.tar.gz
 	 	
example:

		SISTER_AVNG_20220502T180901_L2A_CORFL_001.tar.gz	

The outputs of the PGE use the following naming convention: 

	SISTER_INSTRUMENT_YYYYMMDDTHHMMSS_L2A_SUBPRODUCT_CRID
	

| Subproduct code | Description | Example | 
| ---|---|---|
| CORFL | Corrected reflectance datacube | SISTER_AVNG\_20220502T180901\_L2A_CORFL\_001 |








 
 
 