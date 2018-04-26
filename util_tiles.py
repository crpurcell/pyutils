#!/usr/bin/env python
#=============================================================================#
#                                                                             #
# NAME:     util_tiles.py                                                     #
#                                                                             #
# PURPOSE:  Common functions for using mosaiced data.                         #
#                                                                             #
# MODIFIED: 14-Feb-2018 by C. Purcell                                         #
#                                                                             #
# CONTENTS:                                                                   #
#                                                                             #
# which_wise_tile     ... return WISE tilenames containing a coordinate       #
# which_mipsgal_tile  ... return MIPSGAL tilenames containing a coordinate    #
# which_glimpse_tile  ... return the GLIMPSE tilenames containing a coord     #
# which_glimpse_tiles ... which GLIMPSE tiles overlap a square patch of sky   #
# which_ukidss_cutout ... translate between CORNISH and UKIDSS names          #
# which_grs_cube      ... return the GRS filename containing a coordinate     #
#                                                                             #
#=============================================================================#

# Import standard python modules
import os
import sys
import re
import commands
import math as m
#import pyfits as pf
from astropy.io import fits as pf

from util_misc import *


#-----------------------------------------------------------------------------#
def which_wise_tile(l_deg, b_deg, coordFile):
    """
    Compare a coordinate to the bounds of a list of fits files and return 
    the root name of the WISE FITS tile on disk.
    """

    if not os.path.exists(coordFile):
        print "Coverage description file does not exist."
        return []

    fileLst = []; min_dist = []
    CONFIGFILE = open(coordFile, "r")
    
    # Compile a few useful regular expressions
    spaces = re.compile('\s+')
    comment = re.compile('#.*')
    
    for line in CONFIGFILE:
        line = line.rstrip("\n\r")

        # Filter for comments and blank lines
        if not comment.match(line):
            words = line.split()
            if float(words[1])>180.0:
                words[1]=float(words[1])-360.0
            if float(words[2])>180.0:
                words[2]=float(words[2])-360.0
            lRange=[words[1],words[2]]
            lRange.sort(numeric_compare)
            if l_deg>float(lRange[0]) and l_deg<float(lRange[1]):
                fileLst.append(words[0])
                if abs(l_deg-float(lRange[0]))<abs(l_deg-float(lRange[1])):
                    min_dist.append(abs(l_deg-float(lRange[0])))
                else:
                    min_dist.append(abs(l_deg-float(lRange[1])))

    CONFIGFILE.close()
    if len(fileLst)>0:
        return fileLst[0]
    else:
        return ""


#-----------------------------------------------------------------------------#
def which_mipsgal_tile(l_deg, b_deg, coordFile):
    """
    Compare a coordinate to the bounds of a list of fits files and return 
    the root name of the MIPSGAL FITS tile on disk.
    """

    if not os.path.exists(coordFile):
        print "Coverage description file does not exist."
        return []

    fileLst = []; min_dist = []
    CONFIGFILE = open(coordFile, "r")
    
    # Compile a few useful regular expressions
    spaces = re.compile('\s+')
    comment = re.compile('#.*')
    
    for line in CONFIGFILE:
        line = line.rstrip("\n\r")

        # Filter for comments and blank lines
        if not comment.match(line):
            words = line.split()
            if float(words[1])>180.0:
                words[1]=float(words[1])-360.0
            if float(words[2])>180.0:
                words[2]=float(words[2])-360.0
            lRange=[words[1],words[2]]
            lRange.sort(numeric_compare)
            if l_deg>float(lRange[0]) and l_deg<float(lRange[1]):
                fileLst.append(words[0])
                if abs(l_deg-float(lRange[0]))<abs(l_deg-float(lRange[1])):
                    min_dist.append(abs(l_deg-float(lRange[0])))
                else:
                    min_dist.append(abs(l_deg-float(lRange[1])))

    CONFIGFILE.close()
    if len(fileLst)>0:
        return fileLst[0]
    else:
        return ""


#-----------------------------------------------------------------------------#
def which_glimpse_tile(l_deg,b_deg):
    """Form the root name of the GLIMPSE IRAC tile given a Galctic coord"""

    lCent="%05.0f" % ((float("%05.0f" % (l_deg-0.5))+0.5)*100.0)

    if abs(float(lCent)) > 600:
        if abs(b_deg) <= 0.37625: 
            bCent="+0000"
        else: 
            if b_deg/abs(b_deg)<0: bCent = "%05.0f" % (-75.0)
            else: bCent = "+%04.0f" % (75.0)
    elif abs(float(lCent)) > 300:
        if abs(b_deg) <= 0.62825: 
            bCent="+0000"
        else: 
            if b_deg/abs(b_deg)<0: bCent = "%05.0f" % (-115.0)
            else: bCent = "+%04.0f" % (115.0)
    else:
        if abs(b_deg) <= 0.52175: 
            bCent="+0000"
        else: 
            if b_deg/abs(b_deg)<0: bCent = "%05.0f" % (-155.0)
            else: bCent = "+%04.0f" % (155.0)

    return "GLM_%s%s_mosaic" % (lCent,bCent)


#-----------------------------------------------------------------------------#
def which_glimpse_tiles(l_deg, b_deg, width_deg, height_deg, boundsFile):
    """Determine which GLIMPSE tiles overlap a square patch of sky"""

    fitsRootLst = []
    
    BF = open(boundsFile, 'r')
    for line in BF:
        e = line.strip('\r\n').split()
        e[1:] = [float(x) for x in e[1:]]

        # Tile rectangle
        height2_deg = e[4] - e[3]
        b2_deg = height2_deg / 2.0 + e[3]
        cosf = m.cos( m.radians( max(e[4], e[3]) ) )
        width2_deg = (e[2] - e[1]) / cosf
        l2_deg = (e[2] - e[1]) / 2.0 + e[1]
        
        # Check for overlap
        if check_overlap(l_deg, b_deg, width_deg, height_deg,
                      l2_deg, b2_deg, width2_deg, height2_deg):
            fitsRootLst.append(e[0])
    BF.close()
    
    return fitsRootLst


#-----------------------------------------------------------------------------#
def which_ukidss_cutout(cornishName, lookupFile):
    """Query a translation table between CORNISH names and UKIDSS names"""

    if not os.path.exists(lookupFile):
        print "Lookup file does not exist."
        return []

    lookupTable = {}
    LOOKUPFILE = open(lookupFile, "r")
    
    for line in LOOKUPFILE:        
        line = line.rstrip("\n\r")
        line = line.split()
        lookupTable[line[0]] = line[1]

    if lookupTable.has_key(cornishName):
        return lookupTable[cornishName]
    else:
        return ""


#-----------------------------------------------------------------------------#
def which_grs_cube(l_deg, b_deg):
    """Return the GRS filename containing a coordinate"""

    cubeFitsFile = None; minDist = None;
    
    coverage = [['grs_18-21_cube.fits',18.0,21.0],
                ['grs_19-23_cube.fits',19.0,23.0],
                ['grs_21-25_cube.fits',21.0,25.0],
                ['grs_23-27_cube.fits',23.0,27.0],
                ['grs_25-29_cube.fits',25.0,29.0],
                ['grs_27-31_cube.fits',27.0,31.0],
                ['grs_29-33_cube.fits',29.0,33.0],
                ['grs_31-35_cube.fits',31.0,35.0],
                ['grs_33-37_cube.fits',33.0,37.0],
                ['grs_35-39_cube.fits',35.0,39.0],
                ['grs_37-41_cube.fits',37.0,41.0],
                ['grs_39-43_cube.fits',39.0,43.0],
                ['grs_41-45_cube.fits',41.0,45.0],
                ['grs_43-47_cube.fits',43.0,47.0],
                ['grs_45-49_cube.fits',45.0,49.0],
                ['grs_47-51_cube.fits',47.0,51.0],
                ['grs_49-53_cube.fits',49.0,53.0],
                ['grs_51-55_cube.fits',51.0,55.0],
                ['grs_53-56_cube.fits',53.0,56.0]]

    for entry in coverage:
        
        if float(entry[1])>180.0:
            entry[1]=float(entry[1])-360.0
        if float(entry[2])>180.0:
            entry[2]=float(entry[2])-360.0
            
        lRange=[entry[1],entry[2]]
        lRange.sort(numeric_compare)
        l_mid = entry[2]-entry[1]
        
        if l_deg>float(lRange[0]) and l_deg<float(lRange[1]):

            if cubeFitsFile == None:
                cubeFitsFile = entry[0]
            else:
                if abs(l_deg-l_mid) < minDist:
                    cubeFitsFile = entry[0]
                    
    return cubeFitsFile


#-----------------------------------------------------------------------------#
def which_BOLOCAM_mosaic(l_deg, b_deg):
    """Return the BOLOCAM filenames containing a coordinate"""
    
    mosFitsFile = None
    minDist = None;
    coverage = [['G11.500+0.000_BOLOCAM_mosaic.fits',9.05,14.0],
                ['G15.500+0.000_BOLOCAM_mosaic.fits',13.0,18.0],
                ['G19.500+0.000_BOLOCAM_mosaic.fits',17.0,22.0],
                ['G23.500+0.000_BOLOCAM_mosaic.fits',21.0,26.0],
                ['G27.500+0.000_BOLOCAM_mosaic.fits',25.0,30.0],
                ['G31.500+0.000_BOLOCAM_mosaic.fits',29.0,34.0],
                ['G35.500+0.000_BOLOCAM_mosaic.fits',33.0,38.0],
                ['G39.500+0.000_BOLOCAM_mosaic.fits',37.0,42.0],
                ['G43.500+0.000_BOLOCAM_mosaic.fits',41.0,46.0],
                ['G47.500+0.000_BOLOCAM_mosaic.fits',45.0,50.0],
                ['G51.500+0.000_BOLOCAM_mosaic.fits',49.0,54.0],
                ['G55.500+0.000_BOLOCAM_mosaic.fits',53.0,58.0],
                ['G59.500+0.000_BOLOCAM_mosaic.fits',57.0,62.0],
                ['G63.500+0.000_BOLOCAM_mosaic.fits',61.0,66.0]]

    # Return None if outside the relevant b range
    if l_deg > 14.47 and l_deg < 15.47:
        if abs(b_deg) > 0.996:
            return mosFitsFile
    elif l_deg > 29.47 and l_deg < 31.55:
        if abs(b_deg) > 0.996:
            return mosFitsFile
    elif l_deg > 51.5 and l_deg < 58.4:
        if abs(b_deg) > 0.52:
            return mosFitsFile
    else:
        if abs(b_deg) > 0.6:
            return mosFitsFile
    
    for entry in coverage:
        
        if float(entry[1])>180.0:
            entry[1]=float(entry[1])-360.0
        if float(entry[2])>180.0:
            entry[2]=float(entry[2])-360.0
            
        l_range=[entry[1],entry[2]]
        l_range.sort(numeric_compare)
        l_mid = entry[2]-entry[1]
        
        if l_deg>float(l_range[0]) and l_deg<float(l_range[1]):

            if mosFitsFile == None:
                mosFitsFile = entry[0]
            else:
                if abs(l_deg-l_mid) < minDist:
                    mosFitsFile = entry[0]
                    
    return mosFitsFile

