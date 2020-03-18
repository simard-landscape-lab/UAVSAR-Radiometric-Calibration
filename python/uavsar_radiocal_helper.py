#! /usr/bin/env python
"""Contains a helper function allowing streamlined radiometric calibration
and geocoding of UAVSAR data with a single command line call."""

import os
import os.path
import subprocess
import argparse
from glob import glob

from buildUAVSARhdr import genHDRfromTXT


def runcal(annfile, name=None, caltbl=None, look=None, slope=None,
           mask=None, diff=None):
    """Performs radiometric calibration on a given UAVSAR dataset, and
        geocodes the result.
        
        Note that the actual calibration/geocoding is performed by the
        compiled C programs "uavsar_calib" and "uavsar_geocode".  These
        programs are called by this Python script with the necessary arguments.
        This function will first check to see if "uavsar_calib" and
        "uavsar_geocode" are on the system PATH environment variable and are
        callable from the terminal.  If not, the function will also check to
        see if they are in '../bin/' (relative to the location of this script).
        Please ensure that the programs are accessible through one of these
        methods.
        
        This function will perform the calibration for all three polarizations
        (HHHH, HVHV, and VVVV), or for whichever of those polarizations are
        present if not all of the MLC files can be found.
        
        Both .mlc files (in radar coordinates) and .grd files (in geographic
        coordinates, with ENVI .hdr file) will be created.  Both files will
        contain the calibrated backscattered power (not amplitude) in linear
        units.
        
        Arguments:
            annfile (str): Path and filename of the input annotation file.
            caltbl (str): Path and filename of vegetation look up table file.
                Default: None (perform area correction only).
            look (bool): Boolean flag that sets whether to create look angle
                file.
            slope (bool): Boolean flag that sets whether to save slope angle
                file.
            mask (bool): Boolean flag that sets whether to save mask file.
            diff (bool): Boolean flag that sets whether to create difference
                file.
        
    """
    # Find the programs to call.
    uavsar_calib_prog = subprocess.getoutput('which uavsar_calib')
    uavsar_geocode_prog = subprocess.getoutput('which uavsar_geocode')

    if uavsar_calib_prog == '':
        if os.path.isfile(os.path.dirname(os.path.realpath(__file__))+'/../bin/uavsar_calib'):
            uavsar_calib_prog = os.path.dirname(os.path.realpath(__file__))+'/../bin/uavsar_calib'
        else:
            print('uavsar_radiocal_helper.py -- Cannot find program uavsar_calib.  Aborting.')
            return
            
    if uavsar_geocode_prog == '':
        if os.path.isfile(os.path.dirname(os.path.realpath(__file__))+'/../bin/uavsar_geocode'):
            uavsar_geocode_prog = os.path.dirname(os.path.realpath(__file__))+'/../bin/uavsar_geocode'
        else:
            print('uavsar_radiocal_helper.py -- Cannot find program uavsar_geocode.  Aborting.')
            return
            
    if not os.path.isfile(annfile):
        print('uavsar_radiocal_helper.py -- Cannot find annotation file: "'+annfile+'".  Aborting.')
        return
        
    if name is None:
        name = 'Cal'
        
    # Load dimensions and filenames from annotation file.  
    with open(annfile, 'r') as ann:
        for line in ann:
            if 'mlc_mag.set_cols' in line:
                mlc_cols = line.split()[3]
            elif 'grd_mag.set_cols' in line:
                grd_cols = line.split()[3]
            elif 'grd_mag.set_rows' in line:
                grd_rows = line.split()[3]
            elif 'mlcHHHH' in line:
                mlcfileHH = line.split()[2]
            elif 'mlcHVHV' in line:
                mlcfileHV = line.split()[2]
            elif 'mlcVVVV' in line:
                mlcfileVV = line.split()[2]
        

    # Perform the processing for each polarization.
    datapath = os.path.dirname(annfile)
    
    for pol in range(3):
        if pol == 0:
            mlcfile = mlcfileHH
            polstr = 'HHHH'
            shortpol = 'HH'
        elif pol == 1:
            mlcfile = mlcfileHV
            polstr = 'HVHV'
            shortpol = 'HV'
        else:
            mlcfile = mlcfileVV
            polstr = 'VVVV'
            shortpol = 'VV'
        
        if os.path.isfile(datapath+'/'+mlcfile):
            basefile, ext = os.path.splitext(datapath+'/'+mlcfile)           
            basefile_nopol = basefile.split(polstr)[0] + basefile.split(polstr)[1]
            
            
            # String for call to calibration program.
            calib_exec = uavsar_calib_prog + ' '
            
            if caltbl is not None:
                calfiles = glob(caltbl+'*'+shortpol+'*.flt')

                if len(calfiles) == 1: 
                    calib_exec += '-c '+calfiles[0]+' '
                elif len(calfiles) > 1:
                    print('uavsar_radiocal_helper.py -- Too many calibration table files matching pattern: "'+caltbl+'*'+shortpol+'*.flt".  Aborting.')
                    os._exit(1)
                else:
                    print('uavsar_radiocal_helper.py -- Cannot find calibration table file matching pattern: "'+caltbl+'*'+shortpol+'*.flt".  Aborting.')
                    os._exit(1)

            if look:
                calib_exec += '-l '+basefile_nopol+'_look.grd '

            if slope:
                calib_exec += '-s '+basefile_nopol+'_slope.grd '
                
            if mask:
                calib_exec += '-m '+basefile_nopol+'_'+name+'_mask.mlc '

            if diff:
                calib_exec += '-d '+basefile+'_'+name+'_diff.mlc '
                
            calib_exec += '-u '+basefile+'_geomap.trans '       
            calib_exec += annfile+' '+polstr+' '+basefile+'_'+name+'.mlc'
            
            
            # String for call to geocoding program.
            geocode_exec = uavsar_geocode_prog + ' '
            geocode_exec += basefile+'_'+name+'.mlc '
            geocode_exec += str(mlc_cols) + ' '
            geocode_exec += basefile+'_geomap.trans '
            geocode_exec += basefile+'_'+name+'.grd '
            geocode_exec += str(grd_cols) + ' ' + str(grd_rows)
            
            
            print('uavsar_radiocal_helper.py -- Calibrating file: '+mlcfile)
            print(subprocess.getoutput(calib_exec))
            print('uavsar_radiocal_helper.py -- Geocoding file: '+mlcfile)
            print(subprocess.getoutput(geocode_exec))
            genHDRfromTXT(annfile, basefile+'_'+name+'.grd', polstr)
            
            if look:
                genHDRfromTXT(annfile, basefile_nopol+'_look.grd', polstr)
                look = False # no need to do this for more than one polarization
                
            if slope:
                genHDRfromTXT(annfile, basefile_nopol+'_slope.grd', polstr)
                slope = False # no need to do this for more than one polarization
            
            # Geocode mask file, if we created one.
            if mask:
                geocode_exec = uavsar_geocode_prog + ' '
                geocode_exec += basefile_nopol+'_'+name+'_mask.mlc '
                geocode_exec += str(mlc_cols) + ' '
                geocode_exec += basefile+'_geomap.trans '
                geocode_exec += basefile_nopol+'_'+name+'_mask.grd '
                geocode_exec += str(grd_cols) + ' ' + str(grd_rows)
                
                print('uavsar_radiocal_helper.py -- Geocoding mask file for: '+mlcfile)
                print(subprocess.getoutput(geocode_exec))
                genHDRfromTXT(annfile, basefile_nopol+'_'+name+'_mask.grd', polstr)
                mask = False # no need to do this for more than one polarization
            
            # Geocode difference file, if we created one.
            if diff:
               geocode_exec = uavsar_geocode_prog + ' '
               geocode_exec += basefile+'_'+name+'_diff.mlc '
               geocode_exec += str(mlc_cols) + ' '
               geocode_exec += basefile+'_geomap.trans '
               geocode_exec += basefile+'_'+name+'_diff.grd '
               geocode_exec += str(grd_cols) + ' ' + str(grd_rows)
               print('uavsar_radiocal_helper.py -- Geocoding difference file for: '+mlcfile)
               print(geocode_exec)
               print(subprocess.getoutput(geocode_exec))
               genHDRfromTXT(annfile, basefile+'_'+name+'_diff.grd', polstr)
               
            # Remove temporary geomap.trans file used for geocoding.
            rm_temp = 'rm '+basefile+'_geomap.trans'
            print(subprocess.getoutput(rm_temp))
            
    return



def main():
    print('uavsar_radiocal_helper.py -- Welcome.  Use "-h" flag to see help and options.')
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, help='Specify the input UAVSAR annotation file, or a directory containing multiple UAVSAR annotation files.')
    parser.add_argument('-c', '--cal', type=str, help='Specify the base filename for a vegetation calibration look up table, if desired.  Note that the HH, HV, and VV polarizations will have different files.  For this argument, please provide the part of the filename before the polarization identifier (e.g., "caltbl_test" instead of "caltbl_test_HH.flt").  This program will assume the calibration table filenames contain a two-letter polarization identifier (either "HH", "HV", or "VV) and have the extension ".flt".  If this argument is not specified, only area normalization correction will be performed.')
    parser.add_argument('-n', '--name', type=str, help='An optional string which will be at the end of the filenames of the calibrated output products, in order to identify them (e.g., "CalArea" or "CalVegLUT").  Default: "Cal".')
    parser.add_argument('-l', '--look', action='store_true', help='Toggle to save look angle file.')
    parser.add_argument('-s', '--slope', action='store_true', help='Toggle to save slope angle file.')
    parser.add_argument('-m', '--mask', action='store_true', help='Toggle to save validity mask file.')
    parser.add_argument('-d', '--diff', action='store_true', help='Toggle to save difference file.')
    args = parser.parse_args()
    
    if args.input == None:
        print('uavsar_radiocal_helper.py -- Input UAVSAR annotation file or data path not specified.  Use "-h" to see help and options.')
        os._exit(1)
    elif os.path.isfile(args.input):
        runcal(args.input, caltbl=args.cal, name=args.name, look=args.look,
               slope=args.slope, mask=args.mask, diff=args.diff)
    elif os.path.isdir(args.input):
        print('uavsar_radiocal_helper.py -- Input directory specified.  Batch processing all annotation files found in directory.')
        infiles = [file for file in os.listdir(args.input) if (file.endswith('.ann.txt') or file.endswith('.ann'))]
        
        for annfile in infiles:
            print('uavsar_radiocal_helper.py -- Processing "'+annfile+'"...')
            runcal(os.path.dirname(args.input)+'/'+annfile, caltbl=args.cal, name=args.name, look=args.look,
                   slope=args.slope, mask=args.mask, diff=args.diff)
    else:
        print("uavsar_radiocal_helper.py -- Input UAVSAR annotation file or data path does not exist.  Aborting.")
        os._exit(1)
        
    return


if __name__ == "__main__":
    main()