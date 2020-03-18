# -*- coding: utf-8 -*-
"""
Example script demonstrating batch radiometric calibration using vegetation
LUTs.

Calls functions contained in radiocal module to handle the batch processing.

Actual correction is performed by calling C programs.

Created on Tue May 26 11:32:00 2015

@author: mdenbina
"""
import radiocal




# Path to the UAVSAR data files:
datapath = '/Users/mdenbina/UAVSAR/Louisiana/'

# Path to the folder containing the radiometric calibration programs
# (e.g., uavsar_calib_veg_v2 and geocode_uavsar)
programpath = '/Users/mdenbina/UAVSARRadioCal/'

# Calibration program:
calibprog = programpath+'uavsar_calib'

# Geocoding program:
geocodeprog = programpath+'uavsar_geocode'

# root name of the calibration table to use:
caltblroot = programpath+'caltbl_LA_GulfCo_Wetlands'

# calname='area_only' # name to append to the filename of the calibrated files


# min and max look angles, if post processing is enabled...
# look angles outside these bounds will be set to zero:
minlook = 24
maxlook = 64

# Polarizations to correct:
pol = [0, 1, 2]


# Root names pointing to the UAVSAR data to use for LUT creation, excluding the polarization and correction type (which get appended to this string to produce the full filename).
sardata = ['gulfco_14011_15058_109_150509','gulfco_32017_15058_008_150509','gulfco_14018_15058_007_150509','gulfco_27010_15058_011_150509','gulfco_09010_15058_012_150509','gulfco_32010_15058_110_150509']


# Subpaths pointing to a land cover or mask image to use for each UAVSAR scene.
# len() of maskdata needs to be the same as the len() of sardata.
maskdata = ['CCAP Clipped/CCAP_14011.grd','CCAP Clipped/CCAP_32017.grd','CCAP Clipped/CCAP_14018.grd','CCAP Clipped/CCAP_27010.grd','CCAP Clipped/CCAP_09010.grd','CCAP Clipped/CCAP_32010.grd']


# Path to save the LUT:
LUTpath = '/Users/mdenbina/UAVSARRadioCal/'

# A name to describe the LUT:
LUTname = 'LA_GulfCo_Wetlands'


# The SAR image and the mask should have the same extents, pixel size, etc.
# Trim the mask to the correct size (e.g., in QGIS) before running this.

# Array of allowed values in the mask -- values other than these will be
# excluded from the process.  For a boolean mask this should equal True, to
# include all pixels where the mask is True.  For a land cover image, list
# the class ID numbers of interest.
# Note: For Louisiana data using CCAP land cover, classes 15 and 18 are both
# emergent wetland (18: Estuarine Emergent Wetland, and 15: Palustrine
# Emergent Wetland).
allowed = [15,18]


# String corresponding to the correction we wish to load for LUT creation.  Presumably, generate LUT using area only corrected images, since that's the data we are applying the LUT to.
corrstr = 'area_only'


# These settings determine which pixels we use to generate the LUT, and which pixels are excluded, based on backscatter.
# Note, these cutoff values should be based on HV backscatter.  To be consistent between the polarizations, we always mask out
# the same pixels for each polarization.  The pixels excluded based on backscatter use the HVHV.
max_cutoff = 0.03 # pixels above this value will be excluded
min_cutoff = 0 # pixels below this will be excluded


# Set to true to assume range slope is zero, false otherwise:
flatdemflag = True

# Constant height value for the created flat DEM:
hgtval = -26.5883

# Note that this example script was created for wetland areas along the Gulf
# Coast in Louisiana, United States.  Here we have set flatdemflag to True
# since there is minimal topography in this area.  However, in areas with
# significant topography, one should set flatdemflag to False in order
# to calculate the terrain slope angle from the DEM.  Here we create a
# perfectly flat DEM which is approximately equal to the mean sea level
# height in this area.


# Savitzky-Golay filter to smooth LUT?
sgfilterflag = True # set to True to filter, False to leave alone
sgfilterwindow = 51 # filter window size--larger windows yield more smoothing



# Area Correction (in order to make the data to generate the LUT)
print ('DOING AREA CORRECTION...')
for num in range(0,len(sardata)):
    radiocal.batchcal(datapath, programpath, calibprog, geocodeprog, None,
                      calname='area_only', docorrectionflag=True, zerodemflag=True, 
                      createmaskflag=False, createlookflag=True, createslopeflag=False, 
                      overwriteflag=True, postprocessflag=False, pol=pol, hgtval=hgtval,
                      scene=sardata[num])


# LUT Creation
print ('CREATING LUT...')
radiocal.createlut(datapath, sardata, maskdata, LUTpath, LUTname, allowed,
              pol=pol, corrstr='area_only', min_cutoff=0,
              max_cutoff=0.03, flatdemflag=True, sgfilterflag=True, 
              sgfilterwindow=51, min_look=20, max_look=70, min_samples=100000)



# LUT Correction
print('DOING LUT CORRECTION...')

radiocal.batchcal(datapath, programpath, calibprog, geocodeprog, caltblroot,
             calname='area_veg', docorrectionflag=True, zerodemflag=True, 
             createmaskflag=True, createlookflag=True, createslopeflag=False, 
             overwriteflag=True, postprocessflag=True, minlook=minlook, 
             maxlook=maxlook, pol=pol, hgtval=hgtval)