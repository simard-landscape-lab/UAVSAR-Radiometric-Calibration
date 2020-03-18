Radiometric Calibration of UAVSAR Images Over Forested Terrain With Topography
Software Package
Last Updated: August 8, 2019
==============================================================================


This software implements radiometric calibration of UAVSAR MLC (multi-look complex intensity) images using the algorithms described in the following paper:

M. Simard, B. V. Riel, M. Denbina and S. Hensley, "Radiometric Correction of Airborne Radar Images Over Forested Terrain With Topography," in IEEE Transactions on Geoscience and Remote Sensing, vol. 54, no. 8, pp. 4488-4500, Aug. 2016.
doi: 10.1109/TGRS.2016.2543142
URL: http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=7447752&isnumber=7482858

The software first performs radiometric calibration using area normalization.  The illuminated area as seen by the SAR is calculated for each pixel using a digital elevation model (DEM) and a facet model of the ground terrain.  In addition, the software can correct for changes in radar backscatter in forested and vegetated areas as a function of the ground topography.  This is implemented using empirical look-up tables of vegetation backscatter in the HH, HV, and VV polarizations as a function of the viewing angle and the terrain slope angle in the range direction of the radar.

For more details on the algorithm, see the above paper.



Installation
============

From within the source code directory, the software can be compiled using:

1. make clean
2. make
3. make install

Note that make install simply copies the compiled executables into the /bin/ folder within the current directory.  If you wish to copy the compiled executables to a different folder (e.g., one that is on your PATH), change the BIN variable in the Makefile to the desired folder.  By default the Makefile uses g++ as the compiler.  Make sure that it is installed before attempting to compile the software.  On Mac OS X, g++ can be easily installed by execute "xcode-select --install" in a Terminal window.  This installs the Xcode command line tools, which include the g++ compiler (among other programs).  To use a different compiler, the CC variable in the Makefile can be changed to the desired compiler executable.


Usage
=====

UAVSAR data is freely available for download on the web at http://uavsar.jpl.nasa.gov/.

In order to perform radiometric calibration and geocode the resulting imagery, you will need to download the .ann file (metadata) and the .hgt file (DEM) for the scene of interest.  You will also need to download the .mlc (multi-look intensity) files for the polarizations you wish to calibrate (generally HHHH, HVHV, or VVVV).

There are two programs in this software package.  The first, uavsar_calib, performs the calibration on the original .mlc files.  The second program, uavsar_geocode, geocodes the calibrated .mlc file into a calibrated .grd file with constant spacing in geographic coordinates (latitude/longitude).  We will now discuss how to use these programs.  Note that there is also a helper Python script which streamlines the process a bit compared to using the programs directly.  This helper script is located in python/uavsar_radiocal_helper.py.  To see help and options for that script, run "python uavsar_radiocal_helper.py -h" in that folder.  We now proceed to a discussion of the uavsar_calib and uavsar_geocode programs directly.

An example usage of the uavsar_calib program is as follows:

uavsar_calib -c caltbl_NewHampshire_WhiteMountain_HH.flt -m Brtlet_07101_09061_001.mask -l Brtlet_07101_09061_001.look -s Brtlet_07101_09061_001.slope -u geomap.trans testdata/Brtlet_07101_09061_001_090814_L090_CX_01.ann HHHH Brtlet_07101_09061_001_HHHH_Calibrated.mlc

The -c option specifies the vegetation LUT to use for calibration.  This is an optional flag, and if not specified, the software will perform calibration using the area normalization, then stop, without performing the additional calibration.  In order to create a vegetation LUT file, see the createLUT() function in python/radiocal.py.  The LUT is essentially a file containing the average backscatter of a particular land cover class (e.g., forest, or wetland) as a function of both terrain slope and the SAR viewing angle.  In the vegetation_lut/ folder, there are previously created calibration files for a forested area of the US state of New Hampshire, and a wetland area of the US state of Louisiana.  In order to create a vegetation LUT, it is necessary to know the vegetation type of interest, and also have masks identifying which pixels of the UAVSAR imagery have that vegetation type and should be used in the LUT creation process.  For more details, see the createlut() function in python/radiocal.py, and example usage in python/radiocal_example_script.py.

Returning to the uavsar_calib command line options, the -m option specifies an output filename for a validity mask.  The validity mask will contain a value of 0 for pixels where the correction was performed.  For pixels where the correction could not be performed (e.g., incidence angle out of allowed range, negligible illuminated area), the value will be 1.  Note that the mask file is saved as a 4-byte float flat binary file with the same dimensions as the .mlc files.

The -l option flag specifies an output file to save the look angle for each pixel, measured between the SAR look vector and the nadir.  Similarly, the -s option flag specifies an output file to save the range-facing terrain slope angle for each pixel.  Note that the look angle and slope angle files are saved as flat binary files with the same dimensions and datatype (4-byte float) as the UAVSAR .hgt file containing the DEM used in the SAR processing.

The -u option flag specifies an output filename to store a transformation look up table which is then used in the geocoding process.  If you wish to geocode the calibrated results, this flag must be specified, an the filename given here must be given to the uavsar_geocode program in the later steps.

The -d option flag specifies an output filename to store a calibration difference factor layer which shows the change from the uncalibrated to calibrated files. The difference is calculated as: difference = calibrated/uncalibrated.

Finally, the uavsar_calib program has three required arguments: the annotation file of the data, the 4-letter polarization string (HHHH, HVHV, VVVV) you wish to calibrate, and a filename for the output calibrated .mlc file.  Note that the program assumes that the .mlc and .hgt files are in the same folder as the .ann file.  If they are not, the program will return an error.

Note that typing uavsar_calib without any arguments, or uavsar_calib -h, will give a more detailed list of all option flags.

Once calibration has finished, a new .mlc file containing the calibrated intensity will be created by uavsar_calib with the specified filename.  The next step is to geocode the calibrated .mlc file into a .grd file in geographic coordinates.  This is done using the uavsar_geocode program.  Running the program without any arguments will show the required arguments.  An example run of the program is:

uavsar_geocode Brtlet_07101_09061_001_HHHH_Calibrated.mlc 3300 geomap.trans Brtlet_07101_09061_001_HHHH_Calibrated.grd 24164 9293

The first argument is the input .mlc file to geocode.  The second argument is the number of columns in this file.  This number can be found in the annotation file, under “mlc_pwr.set_cols”.

The third argument is the transformation look up table created by uavsar_calib.  The fourth argument is the desired filename for the newly created geocoded .grd file.

The last two arguments represent the desired number of columns and rows in the output file, respectively.  Often, we will want the output .grd file to have the same dimensions as the .hgt file included with the UAVSAR data (and the look angle and slope angle files optionally created by uavsar_calib).  To find these dimensions, we can look in the annotation file for “hgt.set_cols” and “hgt.set_rows” respectively.

Note that the .grd file created by uavsar_geocode is a flat binary file in regular geographic coordinates which can be loaded into ENVI or other GIS software provided that a suitable header file is created containing the latitude and longitude of the upper left corner of the raster, as well as the latitude and longitude pixel spacing.  In the python/ directory there is a script called buildUAVSARhdr.py, which was created by Nathan Thomas for this purpose.  This script requires three arguments: the input annotation file, the geocoded .grd file, and the 4-letter string describing the polarization (e.g., HHHH, HVHV, or VVVV).  The script will create a .hdr file allowing the .grd file to be loaded.

In addition to buildUAVSARhdr.py, there are three other files in the Python/ directory.  The first is uavsar_radiocal_helper.py.  This is a Python script, executable from the command line, which streamlines and simplifies the usage of uavsar_calib and uavsar_geocode.  It performs radiometric calibration and geocoding, and also calls the Python ENVI .hdr creation, all with a single command line call.  For information on the options and usage of this script, you can execute "python uavsar_radiocal_helper.py -h".  Note that if this script is given a single annotation file as input, it will process that scene.  If it is given a directory as input, it will batch process all UAVSAR data within the given folder.

The other two files are radiocal.py and radiocal_example_script.py.  radiocal.py contains functions for batch processing as well as look-up table creation.  Its options are more in depth than the helper script.  radiocal_example_script.py is an example is an example script showing implementation of batch processing using the functions in radiocal.py.  The example script was used to perform radiometric calibration and geocoding of UAVSAR data located along the Gulf Coast in the US state of Louisiana.