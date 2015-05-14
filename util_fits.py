#!/usr/bin/env python
#=============================================================================#
#                                                                             #
# NAME:     util_fits.py                                                      #
#                                                                             #
# PURPOSE:  Utility functions to operate on FITS data.                        #
#                                                                             #
# MODIFIED: 14-May-2015 by C. Purcell                                         #
#                                                                             #
# CONTENTS:                                                                   #
#                                                                             #
#  mkWCSDict         ... parse fits header using wcslib                       #
#  downsize_wcsdict  ... scale the values in a WCS dict by a factor           #
#  get_wcs           ... older, non WCSLIB version of mkWCSDict               #
#  world2pix_lin     ... simple linear conversion world->pix                  #
#  pix2world_lin     ... simple linear conversion pix->world                  #
#  chan2world_lin    ... simple linear conversion chan->world                 #
#  world2chan_lin    ... simple linear conversion world->chan                 #
#  strip_fits_dims   ... strip header and / or data dimensions                #
#  fits_nans_to_blank .. convert NaNs to magic number in a FITS file          #
#  fits_blank_to_nans .. convert magic number to NaNs in a FITS file          #
#  fits_blank_zeros  ... set values of 0.0 to NaNs in a FITS file             #
#  calc_beam_area    ... calculate the beam area in pixels                    #
#  get_beam_from_history ... search the AIPS history for the clean beam       #
#  save_fits_1D_spec ... save a simple 1D spectrum to a FITS file             #
#  fits_make_lin_axis .. create an array with absica values (linear)          #
#  get_subfits       ... cut a sub portion from a FITS cube (new 2012 code!)  #
#  get_cube_spec     ... extract a spectrum from a FITS cube (needs updating  #
#  biscuit_cutter    ... cut a sub portion from a FITS cube (OLD)             #
#  biscuit_cutter_swarp  cut a sub portion from a FITS map using SWARP        #
#  query_fits_map    ... query the value at a position in a FITS image        #
#                                                                             #
# TODO:                                                                       #
#                                                                             #
#  * Update older wcs.sky2pix to wcs.world2pix                                #
#                                                                             #
#                                                                             #
#=============================================================================#
import os
import sys
import shutil
import copy
import re
import math as m
import commands
import numpy as np
try:
    import pyfits as pf
    import pywcs as pw
except ImportError:
    import astropy.io.fits as pf
    import astropy.wcs.wcs as pw

from util_misc import calc_stats


#-----------------------------------------------------------------------------#
def mkWCSDict(header, wrapX180=False, forceCheckDims=5):
    """Return a dictionary of key WCS parameters using WCSLIB"""

    w = {}    
    header2D = strip_fits_dims(header=header, minDim=2,
                               forceCheckDims=forceCheckDims)
    wcs2D = pw.WCS(header2D)
    w['header2D'] = header2D
    if header['NAXIS']>=3:
        header3D = strip_fits_dims(header=header, minDim=3,
                                   forceCheckDims=forceCheckDims)
        wcs3D = pw.WCS(header3D,)
        w['header3D'] = header3D

    # Header values
    w['xrpix']  = header['CRPIX1']
    w['xrval']  = header['CRVAL1']
    w['xnaxis'] = header['NAXIS1']
    w['yrpix']  = header['CRPIX2']
    w['yrval']  = header['CRVAL2']
    w['ynaxis'] = header['NAXIS2']
    if header['NAXIS']>=3:
        w['zrpix']  = header['CRPIX3']
        w['zrval']  = header['CRVAL3']
        w['znaxis'] = header['NAXIS3']

    # Calculate the centre coordinates
    xCent_pix = float(header['NAXIS1']) / 2.0 + 0.5
    yCent_pix = float(header['NAXIS2']) / 2.0 + 0.5
    if header['NAXIS']>=3:
        zCent_pix = float(header['NAXIS3']) /2.0 + 0.5
        [[w['xcent'], w['ycent'], w['zcent']]] = \
                       wcs3D.wcs_pix2sky([(xCent_pix, yCent_pix, zCent_pix)],1)
    else:
        [[w['xcent'], w['ycent']]] = \
                       wcs2D.wcs_pix2sky([(xCent_pix,yCent_pix)],1)
    if wrapX180:
        if w['xcent']>180.0:
            w['xcent'] = w['xcent'] - 360.0
    
    # Calculate the image bounds in world coords
    try:
        a = wcs2D.calc_footprint()
    except Exception:
        a = wcs2D.calcFootprint()
    w['xmin'] = np.hsplit(a, 2)[0].min()
    w['xmax'] = np.hsplit(a, 2)[0].max()
    w['ymin'] = np.hsplit(a, 2)[1].min()
    w['ymax'] = np.hsplit(a, 2)[1].max()
    if header['NAXIS']>=3:
        [[dummy1, dummy2, z1]] = \
               wcs3D.wcs_pix2sky([(xCent_pix, yCent_pix, 1.0)], 1)
        [[dummy1, dummy2, z2]] = \
               wcs3D.wcs_pix2sky([(xCent_pix, yCent_pix, header['NAXIS3'])], 1)
        w['zmin'] = min([z1, z2])
        w['zmax'] = max([z1, z2])
    if wrapX180:
        if w['xmin']>180.0:
            w['xmin'] = w['xmin'] - 360.0
        if w['xmax']>180.0:
            w['xmax'] = w['xmax'] - 360.0
    
    # Set the type of position coordinate axes
    if wcs2D.wcs.lngtyp == 'RA':
        w['coord_type'] = 'EQU'
    elif wcs2D.wcs.lngtyp == 'GLON':
        w['coord_type'] = 'GAL'
    else:
        w['coord_type'] = ''
        
    # Determine the pixel scale at the refpix
    if header['NAXIS']>=3:
        crpix = np.array(wcs3D.wcs.crpix)
        [[x1, y1, z1]] = wcs3D.wcs_pix2sky([crpix], 1)
        crpix += 1.0
        [[x2, y2, z2]] = wcs3D.wcs_pix2sky([crpix], 1)
        w['zdelt'] = z2 - z1
    else:
        crpix = np.array(wcs2D.wcs.crpix)
        [[x1, y1]] = wcs2D.wcs_pix2sky([crpix], 1)
        crpix -= 1.0
        [[x2, y2]] = wcs2D.wcs_pix2sky([crpix], 1)
    cosy = m.cos( m.radians((y1 + y2) / 2.0) )
    w['xdelt'] = (x2 - x1) * cosy
    w['ydelt'] = y2 - y1
    w['pixscale'] =  (abs(w['xdelt'])*abs(w['xdelt']))**0.5

    return w
    

#-----------------------------------------------------------------------------#
def downsize_wcsdict(w, downsize):
    """Scale the values in a WCS dict to match a resampled data array"""

    w['xrpix']  = float(w['xrpix']) / downsize
    w['xnaxis'] = float(w['xnaxis']) / downsize
    w['yrpix']  = float(w['yrpix']) / downsize
    w['ynaxis'] = float(w['ynaxis']) / downsize
    w['xdelt']  = float(w['xdelt']) * downsize
    w['ydelt']  = float(w['ydelt']) * downsize
    w['pixscale'] = abs( (w['xdelt'] + w['xdelt']) / 2.0 )
    
    return w


#-----------------------------------------------------------------------------#
def get_wcs(header, verbose=False):
    """Parse the simple WCS data from the fits header (assumes linear)"""
    
    w = {}; x_indx=None; y_indx=None; z_indx=None
    
    # Determine the coordinate type and the corresponding keyword index
    ra_regx = re.compile('^RA')
    dec_regx = re.compile('^DEC')
    glon_regx = re.compile('^GLON')
    glat_regx = re.compile('^GLAT')
    velo_regx = re.compile('^VELO')
    freq_regx = re.compile('^FREQ')
    for i in range(int(header['NAXIS'])):
        keyword = "CTYPE"+str(i+1)
        if ra_regx.match(header[keyword]):
            w['coord_type']="EQU";
            x_indx=i+1
        if dec_regx.match(header[keyword]):
            y_indx=i+1
        if glon_regx.match(header[keyword]):
            w['coord_type']="GAL";
            x_indx=i+1
        if glat_regx.match(header[keyword]):
            y_indx=i+1
        if velo_regx.match(header[keyword]):
            w['spec_type']="VELO";
            z_indx=i+1
        if velo_regx.match(header[keyword]):
            w['spec_type']="FREQ";
            z_indx=i+1
    if not x_indx or not y_indx:
        if verbose:
            print "Failed to find determine coordinate system."
        return 1
    if not z_indx and int(header['NAXIS'])>2:
        if verbose:
            print "WARNING: Failed to identify units on the spectral axis."
        z_indx = 0
    
    # Read the axis scaling keywords
    w['x_indx'] = x_indx
    w['y_indx'] = y_indx
    w['xrpix']  = header['CRPIX'+str(x_indx)]
    w['xrval']  = header['CRVAL'+str(x_indx)]
    w['xnaxis'] = header['NAXIS'+str(x_indx)]
    w['yrpix']  = header['CRPIX'+str(y_indx)]
    w['yrval']  = header['CRVAL'+str(y_indx)]
    w['ynaxis'] = header['NAXIS'+str(y_indx)]
    if z_indx:
        w['z_indx'] = z_indx
        w['znaxis'] = header['NAXIS'+str(z_indx)]
        w['zrpix']  = header['CRPIX'+str(z_indx)]
        w['zrval']  = header['CRVAL'+str(z_indx)]
        w['zdelt']  = header['CDELT'+str(z_indx)]
    try:
        xdelt=(-1)*(header['PIXSCAL'+str(x_indx)])/3600.0
        ydelt=(header['PIXSCAL'+str(y_indx)])/3600.0
    except Exception:
        xdelt = None; ydelt=None
    if not xdelt or not ydelt:
        try:
            xdelt=header['CDELT'+str(x_indx)]
            ydelt=header['CDELT'+str(y_indx)]
        except Exception:
            xdelt = None; ydelt=None
    if not xdelt or not ydelt:        
        try:
            xdelt=header['CD1_'+str(x_indx)] 
            ydelt=header['CD2_'+str(y_indx)] 
        except Exception:
            xdelt = None; ydelt=None        
    if not xdelt or not ydelt:
        return 1
    else:
        w['xdelt'] = xdelt
        w['ydelt'] = ydelt
    w['pixscale'] = abs((w['xdelt']+w['xdelt'])/2.0)

    # Calculate the image bounds in world coordinates (center of edge pixels)
    xminmax = []; zminmax = []
    w['ymin'] = (w['yrval']-((w['yrpix']-1)*w['ydelt']))
    w['ymax'] = (w['yrval']+((w['ynaxis']-w['yrpix'])*w['ydelt']))
    xminmax.append((w['xrval']-((w['xrpix']-1)*w['xdelt']
                                /m.cos(m.radians(w['yrval'])))))
    xminmax.append((w['xrval']+((w['xnaxis']-w['xrpix'])*w['xdelt']
                                /m.cos(m.radians(w['yrval'])))))
    w['xmin'] = min(xminmax)
    w['xmax'] = max(xminmax)
    
    if z_indx:
        zminmax.append(w['zrval']-((w['zrpix']-1)*w['zdelt']))        
        zminmax.append(w['zrval']+((w['znaxis']-w['zrpix'])*w['zdelt']))
        w['zmin'] = min(zminmax)
        w['zmax'] = max(zminmax)

    # Calculate the centre coordinates
    w['ycent'] = (w['yrval']+((w['ynaxis']/2.0-w['yrpix']+1)*w['ydelt']))
    cos_factor= m.cos(m.radians(w['ycent']))    
    w['xcent'] = (w['xrval']+((w['xnaxis']/2.0-w['xrpix']+1)
                              *w['xdelt']/cos_factor))

    return w


#-----------------------------------------------------------------------------#
def world2pix_lin(w, x_deg, y_deg):
    """Simple linear conversion from deg to pixel (1-based pix numbering)"""
    
    x_deg = float(x_deg)
    y_deg = float(y_deg)    
    cosF = m.cos( m.radians(y_deg) )
    x_pix = w['xrpix'] + (x_deg - w['xrval']) / (w['xdelt'] / cosF )
    y_pix = w['yrpix'] + (y_deg - w['yrval']) / w['ydelt'] 
    
    return x_pix, y_pix


#-----------------------------------------------------------------------------#
def pix2world_lin(w, x_pix, y_pix):
    """Simple linear conversion from pixel to deg (1-based pix numbering)"""
    
    x_pix = float(x_pix)
    y_pix = float(y_pix)   
    y_deg = w['yrval'] + ( (y_pix - w['yrpix']) * w['ydelt'] )
    cosF  = m.cos( m.radians(y_deg) )
    x_deg = w['xrval'] + ( (x_pix - w['xrpix']) * w['xdelt'] / cosF)

    return x_deg, y_deg


#-----------------------------------------------------------------------------#
def chan2world_lin(w, z_pix):    
    """Simple linear conversion from channel to world coordinates (1-based)"""
    
    z_deg = w['zrval'] + ( (z_pix - w['zrpix']) * w['zdelt'] )

    return z_deg



#-----------------------------------------------------------------------------#
def world2chan_lin(w, z_deg):
    """Simple linear conversion from world to channel coordinates (1-based)"""

    z_pix = w['zrpix'] + ( (z_deg - w['zrval']) / w['zdelt'] )
    
    return z_pix


#-----------------------------------------------------------------------------#
def strip_fits_dims(data=None, header=None, minDim=2, forceCheckDims=5):
    """
    Strip array and / or header dimensions from a FITS data-array or header.
    """

    xydata = None

    # Strip unused dimensions from the header
    if not data is None:
        
        naxis = len(data.shape)
        extraDims = naxis - minDim
        if extraDims < 0:
            print "Too few dimensions in data. "
            sys.exit(1)

        # Slice the data to strip the extra dims
        if extraDims == 0:
            xydata = data.copy()
        elif extraDims == 1:
            xydata = data[0].copy()
        elif extraDims == 2:
            xydata = data[0][0].copy()
        elif extraDims == 3:
            xydata = data[0][0][0].copy()
        else:
            print "Data array contains %s axes" % naxis 
            print "This script supports up to 5 axes only."
            sys.exit(1)
        del data
        
    # Strip unused dimensions from the header
    if not header is None:
        
        header = header.copy()
        naxis = header['NAXIS']
        
        stripHeadKeys = ['NAXIS','CRVAL','CRPIX','CDELT','CTYPE','CROTA',
                         'CD1_','CD2_', 'CUNIT']

        # Force a check on the maximum dimensions of relevant keywords
        if forceCheckDims>0:
            for key in stripHeadKeys:
                for i in range(forceCheckDims+1):
                    if key+str(i) in header:
                        if i > naxis:
                            naxis = i

        # Force a check on max dimensions of the PC keyword array
        if forceCheckDims>0:
            for i in range(1, forceCheckDims+1):
                for j in range(1, forceCheckDims+1):
                    if naxis < max([i,j]):
                        naxis = max([i,j])

        extraDims = naxis - minDim
        if extraDims < 0:
            print "Too few dimensions in data. "
            sys.exit(1)

        # Delete the keyword entries
        for i in range(minDim+1,naxis+1):
            for key in stripHeadKeys:
                if key+str(i) in header:
                    try:
                        del header[key+str(i)]
                    except Exception:
                        pass
                
        # Delete the PC array keyword entries
        for i in range(1,naxis+1):
            for j in range(1, naxis+1):
                key = "PC" + "%03d" % i + "%03d" % j
                if i>minDim or j>minDim:
                    if key in header:
                        try:
                            del header[key]
                        except Exception:
                            pass

        header['NAXIS'] = minDim

    # Return the relevant object(s)
    if not xydata is None and not header is None:
        return [xydata,header]
    elif not xydata is None and header is None:
        return xydata
    elif xydata is None and not header is None:
        return header
    else:
        print "Both header and data are 'Nonetype'."
        sys.exit(1)


#-----------------------------------------------------------------------------#
def fits_nans_to_blank(inFits, outFits=None, blankVal=-999.0, clobber=True,
                       inplace=False):
    """Convert NaNs in a FITS file to a blanking value"""
    
    if not os.path.exists(inFits):
        return 1

    inPath, inFile = os.path.split(inFits)
    inFitsRoot, inExt = os.path.splitext(inFile)

    # Set the output FITS file
    if inplace is False:
        if outFits is None:
            outFits = inFitsRoot + '.BK' + inExt
        if os.path.exists(outFits):
            if clobber == True:
                os.remove(outFits)
            else:
                return 1
    else:
        outFits = inFits

    # Open the files for reading and writing
    if inplace is False:
        hduLst = pf.open(inFits)
    else:
        hduLst = pf.open(inFits, mode='update')
    header = hduLst[0].header
    data = hduLst[0].data

    # Set the header blanking keyword
    oldBlankVal = None
    if 'BLANK' in header:
        oldBlankVal = header['BLANK']
    header['BLANK'] = blankVal
    if not 'BSCALE' in header:
        header['BSCALE'] = 1.0
    if not 'BZERO' in header:
        header['BZERO'], 0.0

    # Blank the relevant data
    hduLst[0].data = np.where(np.isnan(hduLst[0].data),
                              blankVal, hduLst[0].data)
    if not oldBlankVal is None:
        hduLst[0].data = np.where(np.equal(hduLst[0].data, oldBlankVal),
                                  blankVal, hduLst[0].data)

    # Write the changes to disk
    if inplace:
        hduLst.flush()
        return 0
    else:
        pf.writeto(outFits, data, header)
        return 0


#-----------------------------------------------------------------------------#
def fits_blank_to_nans(inFits, outFits=None, clobber=True,inplace=False):
    """Convert magic blanks in a FITS file to NaNs"""

    if not os.path.exists(inFits):
        return 1

    inPath,inFile = os.path.split(inFits)
    inFitsRoot,inExt = os.path.splitext(inFile)

    # Set the output FITS file
    if inplace is False:
        if outFits is None:
            outFits = inFitsRoot + '.BK' + inExt
        if os.path.exists(outFits):
            if clobber == True:
                os.remove(outFits)
            else:
                return 1
    else:
        outFits = inFits
        
    # Open the files for reading and writing
    if inplace is False:
        hduLst = pf.open(inFits)
    else:
        hduLst = pf.open(inFits, mode='update')
    header = hduLst[0].header
    data = hduLst[0].data

    # Check for pre-existing blank keywords
    oldBlankVal = None
    if 'BLANK' in header:
        oldBlankVal = header['BLANK']

    # Blank the data
    hduLst[0].data = np.where(np.equal(data, oldBlankVal),
                              np.nan, hduLst[0].data)

    # Update the header
    header['BLANK'] = -999.0
    if not 'BSCALE' in header:
        header['BSCALE'] = 1.0
    if not 'BZERO' in header:
        header['BZERO'] = 0.0

    # Write the changes to disk
    if inplace:
        hduLst.flush()
        return 0
    else:
        pf.writeto(outFits, data, header)
        return 0


#-----------------------------------------------------------------------------#
def fits_blank_zeros(inFitsFile, hdu=0):
    """Set values of exactly zero to NaNs in a FITS file"""
    
    head = pf.getheader(inFitsFile)
    data = pf.getdata(inFitsFile)
    data = np.where(data==0.0, np.nan, data)
    
    if os.path.exists(inFitsFile):
        os.remove(inFitsFile)
        
    pf.writeto(inFitsFile, data, head, output_verify='fix')
     

#-----------------------------------------------------------------------------#
def calc_beam_area(header):
    """Calculate the effective beam area in px given the Gaussian FWHMs"""

    # Read the beam and pixel dimensions
    try:
        bmaj = header['CLEANBMJ']
        bmin = header['CLEANBMN']
    except Exception:
        [bmaj, bmin, bpa] = get_beam_from_history(header)
        if bmaj is None:
            print "WARNING: Using hardcoded beam size (1.5'')."
            bmaj = 1.5 / 3600.0
            bmin = 1.5 / 3600.0
        
    pixscaleX = header['CDELT1']
    pixscaleX = header['CDELT2']

    # Convert FWHMs to sigma
    gfactor = 2.0 * m.sqrt(2.0 * m.log(2.0))  # FWHM=gfactor*sigma
    sigmaMaj = bmaj / gfactor
    sigmaMin = bmin / gfactor
    
    # Calculate the beam area in the input units and in pixels
    beamAreaSq = 2.0 * m.pi * sigmaMaj * sigmaMin
    pixAreaSq = abs(pixscaleX * pixscaleY)
    beamArea_pix = beamAreaSq /pixAreaSq
    
    return beamArea_pix



#-----------------------------------------------------------------------------#
def get_beam_from_history(header):
    """Search the FITS HISTORY fields for the AIPS clean beam string"""

    bmaj = None
    bmin = None
    bpa = None
    history = header.get_history()
    history.reverse()

    #'AIPS   CLEAN BMAJ=  4.3403E-04 BMIN=  3.1039E-04 BPA= -11.55'
    beamHistStr = 'AIPS\s+CLEAN\sBMAJ=\s+(\S+)\s+BMIN=\s+(\S+)\s+BPA=\s+(\S+)'
    bmHistPat = re.compile(beamHistStr)
    
    for i in range(len(history)):
        
        mch = bmHistPat.match(history[i])
        if mch:
            bmaj = float(mch.group(1))
            bmin = float(mch.group(2))
            bpa  = float(mch.group(3))
            break

    return bmaj, bmin, bpa


#-----------------------------------------------------------------------------#
def save_fits_1D_spec(outName, specArr, cardDict, clobber=True):
    """Save a FITS file containing a single simple 1D spectrum"""
    
    if len(specArr.shape) > 1:
        print "Err: spectrum to be saved is not a 1D array."
        sys.exit(1)
    try:
        hdu = pf.PrimaryHDU(specArr)
        h = hdu.header
        h.update('NAXIS1', len(specArr))
        for key, val in cardDict.iteritems():
            h.update(key, val)
        if os.path.exists(outName):
            if clobber:
                os.remove(outName)
            else:
                print "Err: file '%s' exists." % outName
                return(1)

        hdu.writeto(outName, output_verify='ignore')
        return 0
        
    except Exception:
        return 1
    
    
#-----------------------------------------------------------------------------#
def fits_make_lin_axis(head, axis=0):
    """Create an array with the values of a chosen linear FITS axis"""
    
    axis = int(axis)
    
    if head['NAXIS'] < axis + 1:
        return []
    
    i = str(int(axis) + 1)
    start = head['CRVAL' + i] + (1 - head['CRPIX' + i]) * head['CDELT' + i]
    stop = (head['CRVAL' + i] + (head['NAXIS' + i] + 1 - head['CRPIX' + i]) * 
            head['CDELT' + i])
    
    return np.arange(start, stop, head['CDELT' + i])


#-----------------------------------------------------------------------------#
def get_subfits(inFileName, x_deg, y_deg, radius_deg, zMin_w=None, zMax_w=None,
                do3D=False, wrapX180=False):
    """Cut out a sub-section of a FITS file (cube or map)"""

    # Enforce the sign convention whereby the coordinates range from
    # 180 to -180 degrees, rather than 0 to 360.This is necessary for files
    # which allow negative Galactic coordinates.
    if wrapX180 and  x_deg>180.0:
        x_deg -= 360.0

    # Some sanity checks
    if not os.path.exists(inFileName): 
        return None, None

    # Read the FITS file
    hduLst = pf.open(inFileName, 'readonly', memmap=True)       
    w = mkWCSDict(hduLst[0].header, wrapX180, 5)
    wcs2D = pw.WCS(w['header2D'])
    
    # Find the pixel at the center of the new sub-image
    [ [xCent_pix, yCent_pix] ] = wcs2D.wcs_sky2pix([(x_deg, y_deg)], 0)
    xCentInt_pix = round(xCent_pix)
    yCentInt_pix = round(yCent_pix)

    # Determine the X and Y bounds of the subarray within the main array
    radX_pix = m.ceil( radius_deg / abs(w['xdelt']) )
    radY_pix = m.ceil( radius_deg / abs(w['ydelt']) )
    xMin_pix = max(xCentInt_pix - radX_pix, 0)
    xMax_pix = min(xCentInt_pix + radX_pix, w['header2D']['NAXIS1']-1 )
    yMin_pix = max(yCentInt_pix - radY_pix, 0)
    yMax_pix = min(yCentInt_pix + radY_pix, w['header2D']['NAXIS2']-1 )
    zMin_chan = None
    zMax_chan = None
    if zMin_w and zMax_w and w['header3D']['NAXIS']>2:
        zMin_chan = max(m.ceil( world2chan(w, zMin_w) ), 0)
        zMax_chan = min(m.floor( world2chan(w, zMax_w) ) ,
                                w['header3D']['NAXIS3']-1)

    # Shift the refpix to reflect the new offset
    headSub = hduLst[0].header.copy()
    headSub['CRPIX1'] = hduLst[0].header['CRPIX1'] - xMin_pix
    headSub['CRPIX2'] = hduLst[0].header['CRPIX2'] - yMin_pix
    if zMin_chan and zMax_chan:
        headSub['CRPIX3'] = hduLst[0].header['CRPIX3'] - zMin_chan
      
    # Slice the data
    nAxis = len(hduLst[0].data.shape)
    if zMin_chan and zMax_chan:
        if nAxis == 2:
            dataSub = hduLst[0].data[yMin_pix:yMax_pix+1,
                                     xMin_pix:xMax_pix+1]
        elif nAxis == 3:
            dataSub = hduLst[0].data[zMin_chan:zMax_chan+1,
                                     yMin_pix:yMax_pix+1,
                                     xMin_pix:xMax_pix+1]
        elif nAxis == 4: 
            dataSub = hduLst[0].data[:,
                                     zMin_chan:zMax_chan+1,
                                     yMin_pix:yMax_pix+1,
                                     xMin_pix:xMax_pix+1]
        else:
            return None, None
    else:
        if nAxis == 2:
            dataSub = hduLst[0].data[yMin_pix:yMax_pix+1,
                                     xMin_pix:xMax_pix+1]
        elif nAxis == 3:
            dataSub = hduLst[0].data[:,
                                     yMin_pix:yMax_pix+1,
                                     xMin_pix:xMax_pix+1]
        elif nAxis == 4: 
            dataSub = hduLst[0].data[:,
                                     :,
                                     yMin_pix:yMax_pix+1,
                                     xMin_pix:xMax_pix+1]
        else:
            return None, None

    # Free the memory by closing the file
    hduLst.close()

    # Update the sub header
    try:
        headSub.update("DATAMIN", np.nanmin(dataSub.flatten()), '')
        headSub.update("DATAMAX", np.nanmax(dataSub.flatten()), '')
        headSub.update("NAXIS1", dataSub.shape[-1], '')
        headSub.update("NAXIS2", dataSub.shape[-2], '')
        if nAxis > 2:
            headSub.update("NAXIS3", dataSub.shape[-3], '')
    except Exception:
        pass
    
    return dataSub, headSub


#-----------------------------------------------------------------------------#
def get_cube_spec(inFitsFile, outFitsFile, x_deg, y_deg, radius_deg=None,
                  restfreq_Hz=None):
    """Extract a spectrum from a FITS cube"""

    # Open the files for reading
    cubeData      = pf.getdata(inFitsFile)
    cubeHead      = pf.getheader(inFitsFile)
    w             = get_wcs(cubeHead)
    x_pix, y_pix  = world2pix (w,x_deg,y_deg)
    x_pixInt      = int(round(x_pix))
    y_pixInt      = int(round(y_pix))
    radius_pix    = None
    if radius_deg:
        radius_pix = radius_deg / w['pixscale']

    # Extract the spectrum from a single pixel or a sum under
    # a circle of radius 'radius_deg'
    if not radius_deg:
        spectrum = cubeData[:, y_pixInt, x_pixInt]
    else:
        minXpx = int(m.floor(x_pix-radius_pix))
        maxXpx = int(m.ceil(x_pix+radius_pix))
        minYpx = int(m.floor(y_pix-radius_pix))
        maxYpx = int(m.ceil(y_pix+radius_pix))
        if minXpx<0:
            minXpx=0
        if minYpx<0:
            minYpx=0
        if maxXpx>(w['xnaxis']-1):
            maxXpx=(w['xnaxis']-1)
        if maxYpx>(w['ynaxis']-1):
            maxYpx=(w['ynaxis']-1)

        npix = 0
        spectrum = np.zeros_like(cubeData[:,y_pixInt,x_pixInt])
        for i in range(minXpx,maxXpx+1):
            for j in range(minYpx,maxYpx+1):
                sep = calc_linear_sep(x_pix,y_pix,i,j)
                if sep<=radius_pix:                
                    spectrum+=cubeData[:,j,i]
                    npix+=1
        print "done.\nSummed spectrum over %d pixels." % npix

    # Create an ASCII spectrum
    velChans = []
    for i in range(len(spectrum)):
        v = vox2world(w,i)
        velChans.append(v/1000.0)
    outRoot,extn =os.path.splitext(outFitsFile) 
    outDatFile = outRoot+'.dat'
    if os.path.exists(outDatFile): os.remove(outDatFile)
    print "Writing an ASCII spectrum to '%s'..." % outDatFile,
    DFILE = open(outDatFile,'w')
    for j in range(len(spectrum)):
        DFILE.write("%.4f " % (velChans[j]))
        DFILE.write("%.4f\n" % (spectrum[j]))
    DFILE.close()
    print "done."

    # Write the spectrum to a FITS file
    if os.path.exists(outFitsFile): os.remove(outFitsFile)
    print "Writing a FITS spectrum to '%s'..." % outFitsFile,
    specHDU = pf.PrimaryHDU(spectrum.reshape(1,1,spectrum.shape[0]))
    specHead = specHDU.header
    specHead.update("CTYPE1", cubeHead['CTYPE3'], '')
    specHead.update("CRPIX1", cubeHead['CRPIX3'], '')
    specHead.update("CRVAL1", cubeHead['CRVAL3'], '')
    specHead.update("CDELT1", cubeHead['CDELT3'], '')    
    specHead.update("CTYPE2", cubeHead['CTYPE1'], '')
    specHead.update("CRPIX2", cubeHead['CRPIX1'], '')
    specHead.update("CRVAL2", cubeHead['CRVAL1'], '')
    specHead.update("CDELT2", cubeHead['CDELT1'], '')    
    specHead.update("CTYPE3", cubeHead['CTYPE2'], '')
    specHead.update("CRPIX3", cubeHead['CRPIX2'], '')
    specHead.update("CRVAL3", cubeHead['CRVAL2'], '')
    specHead.update("CDELT3", cubeHead['CDELT2'], '')
    specHead.update("BUNIT", cubeHead['BUNIT'], '')
    specHead.update("EPOCH", cubeHead['EPOCH'], '')
    if not 'RESTFREQ' in header:
        specHead.update("RESTFREQ", restfreq_Hz, '')
    hduList = pf.HDUList([specHDU])
    hduList.writeto(outFitsFile, output_verify='fix')
    print "done."


#-----------------------------------------------------------------------------#
def biscuit_cutter(inFileName, outFileName, x_deg, y_deg, radius_deg,
                   vmin_w=None, vmax_w=None, do3D=False,
                   negXCoords=False):
    """Cut out a subimage from a FITS cube given a center and radius."""

    # Enforce the sign convention whereby the coordinates range from
    # 180 to -180 degrees, rather than 0 to 360.This is necessary for files
    # which allow negative Galactic coordinates.
    if negXCoords and  x_deg>180.0:
        x_deg -= 360.0
        
    # Some sanity checks
    if not os.path.exists(inFileName): 
        return 1
    
    # Read the header and data array from the input image
    header = pf.getheader(inFileName)
    data = pf.getdata(inFileName)

    # Strip unused dimensions from the array
    xydata = strip_data_dims(data, do3D)
    
    # Parse the header for the WCS information
    w = get_wcs(header)

    # Convert the requested centre coordinates to pixel coordinates
    [x_pix, y_pix] = world2pix(w, x_deg, y_deg)
    
    # Round the requested coordinates to the nearest 
    # pixel center and set as the reference pixel values,  
    cos_factor= m.cos(m.radians(y_deg)) 
    x_pix_int = int(x_pix+0.5)
    y_pix_int = int(y_pix+0.5)
    xrval_new, yrval_new = pix2world(w, x_pix_int, y_pix_int)
    
    # Determine the X-Y bounds of the sub-image in pixels.
    xbound=[]
    ybound=[]
    radius_x_px = int(radius_deg/w['xdelt'])
    radius_y_px = int(radius_deg/w['ydelt'])
    xbound.append(x_pix_int+radius_x_px)
    xbound.append(x_pix_int-radius_x_px)
    ybound.append(y_pix_int+radius_y_px)
    ybound.append(y_pix_int-radius_y_px)
    xmin=min(xbound)-1   # Convert 1-based indexing to 0-based indexing
    xmax=max(xbound)-1
    ymin=min(ybound)-1
    ymax=max(ybound)-1
    if xmin<0:
        xshift=xmin
        xmin=0
    else:
        xshift=0
    if ymin<0:
        yshift=ymin
        ymin=0
    else:
        yshift=0
    yrpix_new = (radius_y_px*(w['ydelt']/abs(w['ydelt'])))+1+yshift
    xrpix_new = (radius_x_px*(w['xdelt']/abs(w['xdelt'])))+1+xshift

    # Determine the Z bounds of the sub image in pixels
    if vmin_w and vmax_w:
        zmin = round(world2vox(w,vmin_w))
        zmax = round(world2vox(w,vmax_w))
        if zmin<1:
            zmin=1
        if zmax>w['znaxis']:
            zmax=w['znaxis']
        zrpix_new = round((zmax-zmin)/2.0)
        zrval_new = vox2world(w,zrpix_new+zmin-1)
        zmin=zmin-1 # Convert 1-based indexing to 0-based indexing
        zmax=zmax-1        
        
    # Slice the data into a subimage, according to the data array dimensions
    datasub=np.array([])
    if not do3D:
        datasub = xydata[ymin:ymax+1,xmin:xmax+1].copy()
    else:
        if vmin_w and vmax_w:
            datasub = xydata[zmin:zmax+1,ymin:ymax+1,xmin:xmax+1].copy()
        else:
            datasub = xydata[:,ymin:ymax+1,xmin:xmax+1].copy()

    # Update the header keywords
    headersub=header.copy()
    try:
        headersub.update("CRVAL"+str(w['x_indx']), xrval_new, '')
        headersub.update("CRPIX"+str(w['x_indx']), xrpix_new, '')
        headersub.update("CRVAL"+str(w['y_indx']), yrval_new, '')
        headersub.update("CRPIX"+str(w['y_indx']), yrpix_new, '')        
        headersub.update("DATAMIN", np.nanmin(datasub), '')
        headersub.update("DATAMAX", np.nanmax(datasub), '')
        if w['coord_type']=='GAL': headersub.update("EQUINOX", '1950.0  ', '')
        if vmin_w and vmax_w:
            headersub.update("CRVAL"+str(w['z_indx']), zrval_new, '')
            headersub.update("CRPIX"+str(w['z_indx']), zrpix_new, '')
    except Exception:
        del header; del data; del datasub
        return 1

    # Write the new image and header to the file.
    if os.path.exists(outFileName): os.remove(outFileName)
    try:
        pf.writeto(outFileName, datasub, headersub)
    except Exception:
        del header; del data; del datasub
        return 3

    # Clean up and return
    del header
    del data
    del datasub
    return 0


#-----------------------------------------------------------------------------#
def biscuit_cutter_swarp(inFiles, inDir, outFile, outDir, xCent_deg, yCent_deg,
                         xSide_deg, ySide_deg=None, outSystem='EQUATORIAL'):
    """Cutout a FITS image from a larger 2D FITS file using SWARP"""
    
    # Temporarally change to the temp dir
    curDir = os.getcwd()
    os.chdir(outDir)

    if ySide_deg is None:
        ySide_deg = xSide_deg

    header = pf.getheader(inDir + '/' + inFiles[0])
    w = mkWCSDict(header)
    xSidePix = abs(xSide_deg/w['xdelt'])
    ySidePix = abs(ySide_deg/w['xdelt'])
    del header
            
    # Setup the output paths
    fileRoot, fileExt = os.path.splitext(outFile)
    outWeightFits = fileRoot +"_wei" + fileExt
    #configFile = outDir + '/' + fileRoot + "_swarp.cfg"
    configFile = fileRoot + "_swarp.cfg"
    if os.path.exists(configFile):
        os.remove(configFile)
    if os.path.exists(outFile):
        os.remove(outFile)

    # Write a swarp config file
    CFILE = open(configFile,'w')
    CFILE.write( "IMAGEOUT_NAME          %s\n" % (outFile))
    CFILE.write( "WEIGHTOUT_NAME         %s\n" % (outWeightFits))
    CFILE.write( "CELESTIAL_TYPE         %s\n" % outSystem)
    CFILE.write( "PROJECTION_TYPE        SIN\n")
    CFILE.write( "CENTER_TYPE            MANUAL\n")
    CFILE.write( "CENTER                 %f %f\n" % (xCent_deg, yCent_deg))
    CFILE.write( "IMAGE_SIZE             %d,%d\n" % (xSidePix, ySidePix))
    CFILE.write( "SUBTRACT_BACK          N\n")
    CFILE.write( "BLANK_BADPIXELS        N\n")
    CFILE.write( "WRITE_XML              N\n")
    CFILE.write( "VERBOSE_TYPE           QUIET\n")
    CFILE.close()

    # Construct the SWARP command and run
    try:
        command = 'swarp -c %s ' % configFile
        for inFile in inFiles:
            command += ' %s/%s' % (inDir, inFile)
        commands.getstatusoutput(command)
    except Exception:            
        os.remove(configFile)
        os.chdir(curDir)
        return 1

    # Clean up
    os.remove(configFile)
    os.remove(outWeightFits)
    os.chdir(curDir)
    return 0


#-----------------------------------------------------------------------------#
def query_fits_map(header, data, x_deg, y_deg, doStrip=True, doApp=False,
                   r_deg=10.0/60.0):
    """Query the value in a 2D FITS map at a given coordinate in deg"""

    dataValue = None
    
    # Strip unused dimensions from the array
    if doStrip:
        data, header = strip_fits_dims(data, header, 2, 5)
    
    # Extract the data
    try:
        w = mkWCSDict(header)
        wcs = pw.WCS(w['header2D'])
        [[x_pix, y_pix]] =  wcs.wcs_sky2pix([(x_deg, y_deg)], 0)
        dataValue = data[int(round(y_pix)), int(round(x_pix))]
        r_pix = r_deg / w['pixscale']
        if doApp:
            xMax_pix = min(int(round(x_pix + r_pix)), w['xnaxis'])
            xMin_pix = max(int(round(x_pix - r_pix)), 0)
            yMax_pix = min(int(round(y_pix + r_pix)), w['ynaxis'])
            yMin_pix = max(int(round(y_pix - r_pix)), 0)        
            dataSub = data[yMin_pix:yMax_pix, xMin_pix:xMax_pix ]
            ms = calc_stats(dataSub)
        
    except Exception:
        print "Fail query FITS pixel."

    if doApp:
        return dataValue, ms
    else:        
        return dataValue
