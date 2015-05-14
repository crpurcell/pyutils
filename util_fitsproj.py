#!/usr/bin/env python
#=============================================================================#
#                                                                             #
# NAME:     util_fitsproj.py                                                  #
#                                                                             #
# PURPOSE:  Utility functions to reproject or convert FITS data.              #
#                                                                             #
# MODIFIED: 14-May-2015 by C. Purcell                                         #
#                                                                             #
# CONTENTS:                                                                   #
#                                                                             #
#  regrid_equJ2000_swarp ... call SWARP to reproject an image into J2000      #
#  mosaic_fits_montage   ... call Montage to mosaic multiple fits files       #
#                                                                             #
#=============================================================================#
import os
import sys
import shutil
import math as m
import commands
import numpy as np
import montage as mo
try:
    import pyfits as pf
    import pywcs as pw
except ImportError:
    import astropy.io.fits as pf
    import astropy.wcs.wcs as pw


#-----------------------------------------------------------------------------#
def regrid_equJ2000_swarp(inFitsFile, outFitsFile):
    """Call SWARP to reproject an image into Equatorial J2000 coords."""
    
    fileRoot, fileExt = os.path.splitext(outFitsFile)
    outWeightFits = fileRoot +"_wei" + fileExt
    configFile = fileRoot + "_swarp.cfg"
    if os.path.exists(configFile): 
	    os.remove(configFile)
    if os.path.exists(outFitsFile): 
	    os.remove(outFitsFile)
    
    # Write a swarp config file
    CFILE = open(configFile,'w')
    CFILE.write( "IMAGEOUT_NAME          %s\n" % outFitsFile)
    CFILE.write( "WEIGHTOUT_NAME          %s\n" % outWeightFits)
    CFILE.write( "CELESTIAL_TYPE         EQUATORIAL\n")
    CFILE.write( "PROJECTION_TYPE        SIN\n")
    CFILE.write( "CENTER_TYPE            ALL\n")
    CFILE.write( "SUBTRACT_BACK          N\n")
    CFILE.write( "BLANK_BADPIXELS        N\n")
    CFILE.write( "WRITE_XML              N\n")
    CFILE.write( "VERBOSE_TYPE           QUIET\n")
    CFILE.close()
    
    # Construct the SWARP command and run
    try:
        command = 'swarp ' + inFitsFile + " -c " + configFile
	run_command(command)
    except Exception:            
        os.remove(configFile)
	raise Exception

    # Clean up
    os.remove(configFile)
    os.remove(outWeightFits)

    
#-----------------------------------------------------------------------------#
def mosaic_fits_montage(inFiles, outFile):
    """Call Montage to mosaic multiple fits files together"""

    # Split the paths and extensions
    outFilePath, outFileName = os.path.split(outFile)
    outFileRoot, outFileExt = os.path.splitext(outFileName)
    inFilePath, dummy = os.path.split(inFiles[0])
    if inFilePath == '':
        inFilePath = '.'        

    # Output a list of images
    if os.path.exists(outFileRoot + '.imglst'):
        os.remove(outFileRoot + '.imglst')
    LFILE = open(outFileRoot + '.imglst','w')
    LFILE.write( "|                                             fname|\n")
    LFILE.write( "|                                              char|\n")
    for fitsFile in inFiles:
        inFilePath,fileName = os.path.split(fitsFile)
        LFILE.write("%s\n" % fileName)
    LFILE.close()
    
    # Assemble the files for mosaicing
    if os.path.exists(outFilePath + '/' + outFileRoot):
        shutil.rmtree(outFilePath + '/' + outFileRoot, True)
    if os.path.exists(outFile):
        shutil.rmtree(outFile,True)
    try:
        if os.path.exists(outFileRoot):
            shutil.rmtree(outFileRoot,True)
        mo.mosaic(inFilePath, outFileRoot, imglist=outFileRoot+'.imglst',
                  exact_size=True)
        shutil.move(outFileRoot+'/mosaic.fits', outFile)
    except Exception:
        print "MONTAGE mosaic task failed!"
        if os.path.exists(outFileRoot + '.imglst'):
            os.remove(outFileRoot + '.imglst')
        if os.path.exists(outFileRoot):
            shutil.rmtree(outFileRoot, True)
        return 1

    # Clean up
    if os.path.exists(outFileRoot + '.imglst'):
        os.remove(outFileRoot + '.imglst')
    if os.path.exists(outFileRoot):
        shutil.rmtree(outFileRoot, True)
