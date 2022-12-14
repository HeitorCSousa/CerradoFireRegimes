ESA AVHRR-LTDR Fire_cci Burned Area data: Pixel product
-----------------------------------------------
-----------------------------------------------

The AVHRR-LTDR Fire_cci Burned Area (BA) product v1.1 (also called FireCCILT11 for short) was obtained based on spectral information from the AVHRR Land Long Term Data Record (LTDR) product version 5 dataset (https://ltdr.modaps.eosdis.nasa.gov/cgi-bin/ltdr/ltdrPage.cgi)

FireCCILT11 is available for years 1982 to 2018 at pixel resolution (0.05 degrees) and grid resolution (0.25 degrees). Year 1994 is not included in the dataset because there was not enough input data from the LTDR Archive to process the dataset for that year.


Details of the Pixel product
---------------------------
The pixel product is composed of 5 files:

 - ...JD.tif: Day of first detection (Julian Day) of the burned area 

 - ...CL.tif: Confidence level of burned area detection

 - ...BA.tif: Burned area corresponding to the proporcion of the pixel that is calculated as burned.

 - ...OB.tif: Number of observations, i.e. times that the pixel has been observed in the month.

 - ....xml: Metadata of the product

The Spatial resolution of this BA product is 0.05 degrees, which is the resolution of the AVHRR-LTDR input data. 

The Coordinate Reference System (CRS) used is a geographic coordinate system (GCS) based on the World Geodetic System 84 (WGS84) reference ellipsoid and using a Plate Carrée projection with geographical coordinates of equal pixel size.

This product is distributed in global monthly files, grouped by year.


File format
-----------
The product is delivered in GeoTIFF format, and compressed into tar.gz files to reduce downloading file sizes.


Further information
-------------------
For further information see the Product User Guide at https://climate.esa.int/en/projects/fire/key-documents/.


Example citation:
----------------  

Chuvieco, E.; Pettinari, M.L.; Otón, G. (2020): ESA Fire Climate Change Initiative (Fire_cci): AVHRR-LTDR Fire_cci Burned Area Pixel product, version 1.1. Centre for Environmental Data Analysis, 28 December 2020. doi:10.5285/b1bd715112ca43ab948226d11d72b85e.
