#!/usr/bin/env python
#=============================================================================#
#                                                                             #
# NAME:     util_phot.py                                                      #
#                                                                             #
# PURPOSE:  Utility functions to perform photometry on FITS images.           #
#                                                                             #
# MODIFIED: 14-May-2015 by C. Purcell                                         #
#                                                                             #
# CONTENTS:                                                                   #
#                                                                             #
# polystr_to_arr    ... convert a polygon csv string to a numpy array 2XN     #
# polylst_parse     ... convert a list of strings into [(x1,y1),(x2,y2) ...]  #
# expand_polygon    ... expand a polygon aperture by a given buffer           #
# get_pix_in_poly   ... return the indices of the pixels inside a polygon     #
# point_in_poly     ... query if a point is inside the polygon                #
# point_in_poly_old ... old (loopy) version of point in polygon query         #
# poly_area         ... calculate the area of a polygon                       #
# polyarr2pix       ... convert array of vertices to pixel coordinates        #
# gen_polyell       ... generate a polyArray approximating an ellipse         #
# measure ellipse   ... measure image properties under an ellipse             #
# is_in_ellipse     ... test if a point is in an ellipse                      #
# measure_pix_props ... measure the properties of a 2D array                  #
# do_ap_phot        ... perform aperture photometry under an ellipse (OLD)    #
#                                                                             #
# TODO:                                                                       #
#                                                                             #
# * Convert get_pix_in_poly to vectorised routine                             #
# * Check and fix do_ap_phot                                                  #
#                                                                             #
#=============================================================================#
import os
import sys
import math as m
import numpy as np
try:
    import astropy.io.fits as pf
    import astropy.wcs.wcs as pw
except ImportError:
    import pyfits as pf
    import pywcs as pw
from shapely.geometry import Polygon as Polyshape
from matplotlib.path import Path

from util_misc import *
from util_fits import *


#-----------------------------------------------------------------------------#
def polystr_to_arr(polystr):
    """Convert the polygon strings to numpy arrays 2XN"""

    try: 
        polystrLst = polystr.split(',')
    except Exception:
        polystrLst = ['']
        
    if polystrLst==[''] :
        polyArr = np.array([])
    else:
        polyArr = np.array([float(x) for x in polystrLst]).reshape(-1,2)
        
        # Close the polygon if not already closed
        if not np.all(polyArr[0]== polyArr[-1]):
            polyArr = np.append(polyArr, [polyArr[0]], 0)
            
    return polyArr


#-----------------------------------------------------------------------------#
def polylst_parse(polyLst):
    """Convert a list of strings into a list of coordinate tuples"""

    polyLst = [float(x) for x in polyLst]
    polyLstTup = []
    for i in range(0, len(polyLst), 2):
        polyLstTup.append(tuple(polyLst[i:i+2]))

    return polyLstTup
    
    
#-----------------------------------------------------------------------------#
def expand_polygon(polyArr, offset):
    """Expand a polygon using the Shapely buffer function"""

    try:
        # Convert the polygon vertex array to a LineStr object
        polyShape = Polyshape(polyArr)
        
        # Offset the polygon by an ammount offset_deg
        polyShapeNew = polyShape.buffer(offset)
        
        # Convert the LineStr object back to a numpy array
        outPolyArr = np.array(polyShapeNew.exterior.xy).transpose()

        return outPolyArr
    
    except Exception:
        
        return np.array([])


#-----------------------------------------------------------------------------#
def get_pix_in_poly(poly):
    """Return the indices of the pixels inside a polygon
    TODO: Convert this to a vectored version using matplotlib.nxutils"""

    indices = []

    # Get the bounds of the poly vertices
    px = np.array(poly)[:,0]
    py =  np.array(poly)[:,1]
    px_min = int(m.floor(min(px)))
    px_max = int(m.ceil(max(px)))
    py_min = int(m.floor(min(py)))
    py_max = int(m.ceil(max(py)))
    
    for i in range(px_min, px_max+1):
        for j in range(py_min, py_max+1):
            if point_in_poly(i, j, poly):
                    indices.append((i, j))

    return indices


#-----------------------------------------------------------------------------#
def point_in_poly(i, j, poly):
    """Query if a point is inside the polygon"""
    
    return Path(poly).contains_point([i, j])


#-----------------------------------------------------------------------------#
def point_in_poly_old(px, py, poly):
    """Query if a point is inside the polygon (old, slow manual version)"""

    cn = 0

    i = -1
    j = 0
    while j < len(poly):
        qx, qy = poly[i]
        rx, ry = poly[j]

        if (px, py) == (qx, qy):
            return True

        if (    ((qy <= py) and (ry > py)) or \
                ((qy > py) and (ry <= py))    ):
            vt = (py - qy) / (ry - qy)
            if (px < qx + vt * (rx - qx)):
                cn += 1

        i = j
        j += 1

    return cn % 2 


#-----------------------------------------------------------------------------#
def poly_area(p):
    """Calculate the area of a polygon (from Stack Overflow)"""

    # Close the polygon
    if not p[0][0] == p[-1][0] and p[0][1] == p[-1][1]:
        p.append(p[0])

    try:
        p = p.tolist()
    except Exception:
        pass
    
    def segments(p):
        return zip(p, p[1:] + [p[0]])
    
    return 0.5 * abs(sum(x0*y1 - x1*y0
                         for ((x0, y0), (x1, y1)) in segments(p)))


#-----------------------------------------------------------------------------#
def polyarr2pix(polyArr, header):
    """Convert an array of WC verties to pixel-based coordinates"""
    
    w = mkWCSDict(header)
    wcs = pw.WCS(w['header2D'])
    
    return wcs.wcs_sky2pix(polyArr, 0)


#-----------------------------------------------------------------------------#
def gen_polyell(minAx, majAx, pa_deg=0.0, xCen=0.0, yCen=0.0, nPts=20):
    """Generate a polygon approximating an ellipse"""
    polyLst = []
    pa_rad = m.radians(pa_deg)
    dt = m.pi * 2.0 / nPts 
    for i in range(nPts+1):
        t = dt * i
        x = ( xCen
              + minAx * m.cos(t) * m.cos(pa_rad)
              - majAx * m.sin(t) * m.sin(pa_rad) )
        y = ( yCen
              + minAx * m.cos(t) * m.sin(pa_rad) 
              + majAx * m.sin(t) * m.cos(pa_rad) )
        polyLst.append( (x, y) )

    return np.array(polyLst)


#-----------------------------------------------------------------------------#
def measure_ellipse(header, data, x_deg, y_deg, maj_deg, min_deg, pa_deg):
    """Measure the image properties under an ellipse"""

    wcs = pw.WCS(header)
    w = mkWCSDict(header)

    # Convert inputs to pixel coordinates and calculate region bounds
    [[x_pix, y_pix]] =  wcs.wcs_sky2pix([(x_deg,y_deg)], 0)
    maj_pix = maj_deg / w['pixscale']
    min_pix = min_deg / w['pixscale']
    minX_pix = int(max(math.floor(x_pix - maj_pix), 0.0))
    maxX_pix = int(min(math.ceil(x_pix + maj_pix), header['NAXIS1']-1))
    minY_pix = int(max(math.floor(y_pix - maj_pix), 0.0))
    maxY_pix = int(min(math.ceil(y_pix + maj_pix), header['NAXIS2']-1))

    # Calculate a masked distance array and the indices of the pixels
    # inside a radius of maj_pix
    yi, xi = np.indices(data.shape)
    distArr = np.sqrt( np.power(xi - x_pix, 2) + np.power(yi - y_pix, 2) )
    distArr = np.where(distArr<=maj_pix, distArr, np.nan)
    indices = np.where(~np.isnan(distArr))

    # Measure the properties of the pixels
    ms = measure_pix_props (array, indices=indices, docent=True)


#-----------------------------------------------------------------------------#
def is_in_ellipse(x1, y1, x, y, a, b, paDeg):
    """Check if a point is in an ellipse"""

    x1 = float(x1)
    y1 = float(y1)
    x = float(x)
    y = float(y)
    a = float(a)
    b = float(b)
    paRad = math.radians(float(paDeg))

    dX = x - x1
    dY = y - y1
    
    inTest = ( math.pow((dX * math.cos(paRad) + dY * math.sin(paRad)) /b, 2.0)
             +math.pow((dX * math.sin(paRad) - dY * math.cos(paRad)) /a, 2.0))

    if inTest < 1.0:
        return True
    else:
        return False

            
#-----------------------------------------------------------------------------#
def measure_pix_props(dataArr, indices=None, doPositions=False,  clip=0.0,
                      nIter=10):
    """Calculate the statistics of selected pixels in an array.
    NOTE: 'indices' contains two columns ( arr(yIndices), arr(xIndices) )"""

    # If indices are provided mask off other pixels as NaNs
    if not indices is None:
        mask = np.ones_like(dataArr)
        mask[indices] = 0
        nanArr = np.where(np.equal(mask, 1.0), np.nan, dataArr)
    else:
        nanArr = dataArr
        indices = np.indices(nanArr.shape)

    # Perform the measurements
    ms = calc_clipped_stats(nanArr, clip, nIter)
    #ms['mean']
    #ms['median']
    #ms['stdev']
    #ms['madfm']
    #ms['max']
    #ms['min']
    #ms['centmax']
    #ms['npix'] 
    ms['summ'] = np.nansum(nanArr)
    ms['summ2'] = np.nansum(np.power(nanArr, 2.0))
    ms['centmax'].reverse()

    # Variables to contain positional measurements
    ms['cent'] = 0                # Geometric centroid
    ms['wcent'] = 0               # Weighted centroid
    ms['radWCentPix'] = 0.0       # Weighted radius
    
    ms['rSum'] = 0.0              # Radii sums, used in error calcs
    ms['rSum2'] = 0.0             # 
    ms['rAmpSum'] = 0.0           #
    
    ms['xSum'] = 0.0              # X & Y positional sums
    ms['xSum2'] = 0.0             #         
    ms['xAmpSum'] = 0.0           # 
    ms['ySum'] = 0.0              # 
    ms['ySum2'] = 0.0             # 
    ms['yAmpSum'] = 0.0           #
    
    ms['srcArr'] = []             #

    # If a 2D array, calculate the centroids and radius
    if doPositions:
        if not len(nanArr.shape)==2:
            print "WARNING: array is not 2D. Can't measure positions."
        else:

            # Calculate the weighted centroid
            weightSum = 0.0
            for j, i in np.array(indices).transpose():
                ms['xSum'] += i
                ms['xSum2'] += i**2.0
                ms['xAmpSum'] += ( i * dataArr[j,i] )
                ms['ySum'] += j
                ms['ySum2'] += j**2.0       
                ms['yAmpSum'] += ( j * dataArr[j,i] )
                ms['srcArr'].append( dataArr[j,i] )
                weightSum += dataArr[j,i]
            ms['wcent'] = ( ( ms['xAmpSum'] / weightSum),
                            ( ms['yAmpSum'] / weightSum) )
        
            # Calculate the geometric centroid
            ms['cent'] = (float(sum( indices[1]) ) / len(indices[1]),
                          float(sum( indices[0]) ) / len(indices[0]))
    
            # Calculate the brightness weighted radius (ref to wcent)
            weightSum = 0.0
            for j, i in np.array(indices).transpose():
                r = calc_linear_sep( ms['wcent'][0], ms['wcent'][1],
                                     float(i), float(j))
                ms['rSum'] += r
                ms['rSum2'] += r**2.0
                ms['rAmpSum'] += (dataArr[j,i] * r)
                weightSum += dataArr[j,i]
            ms['radWCentPix'] = ms['rAmpSum'] / weightSum
        ms['srcArr'] = np.array( ms['srcArr'] )
    ms['srcArr'] = np.array(ms['srcArr'])
    
    return ms


#-----------------------------------------------------------------------------#
def do_ap_phot(header, data, x_deg, y_deg, srcMaj_deg, srcMin_deg, srcPA_deg,
               innerRad_deg, outerRad_deg):
    """Perform aperture photometry
    TODO: Edit this for clarity and fix the errors"""

    wcs = pw.WCS(header)
    w = mkWCSDict(header)

    # Aperture sizes in pixels
    srcMaj_pix = srcMaj_deg/w['pixscale']
    srcMin_pix = srcMin_deg/w['pixscale']
    innerRad_pix = innerRad_deg/w['pixscale']
    outerRad_pix = outerRad_deg/w['pixscale']

    # Set the region bounds in *zero based coords*
    [[x_pix, y_pix]] =  wcs.wcs_sky2pix([(x_deg,y_deg)], 0)
    minX_pix = int(max(math.floor(x_pix - outerRad_pix), 0.0))
    maxX_pix = int(min(math.ceil(x_pix + outerRad_pix), header['NAXIS1']-1))
    minY_pix = int(max(math.floor(y_pix - outerRad_pix), 0.0))
    maxY_pix = int(min(math.ceil(y_pix + outerRad_pix), header['NAXIS2']-1))
    
    # Loop over the relevant pixels summing their values
    npixSrc    = 0
    npixSky    = 0
    skyDataLst = []
    srcDataLst = []
    srcIndxLst = []
    sigDataLst = []
    for i in range(minX_pix, maxX_pix+1):
        for j in range(minY_pix, maxY_pix+1):

            sep = calc_linear_sep(x_pix, y_pix, i, j)
            
            # Pick pixels from an elliptical source aperture
            a = srcMaj_pix
            b = srcMin_pix
            x = x_pix
            y = y_pix
            if is_in_ellipse(i, j, x, y, a, b, srcPA_deg):
                npixSrc += 1
                srcDataLst.append(data[j,i])
                srcIndxLst.append((float(i), float(j)))
                #dataTmp[j,i] = 2.0
            
            # Pick pixels from the sky annulus
            if sep >= innerRad_pix and sep <= outerRad_pix:
                if not np.isnan(data[j,i]):
                    npixSky += 1
                    skyDataLst.append(data[j,i])
                    #dataTmp[j,i] = 1.0

    # Measure the relevant source properties
    mSrc = calc_stats(srcDataLst)
    medianSrc  = mSrc['median']
    meanSrc    = mSrc['mean']
    madSrc     = mSrc['madfm']
    stdevSrc   = mSrc['stdev']
    maxSrc     = mSrc['max']
    minSrc     = mSrc['min']
    npixSrc    = mSrc['npix']
        
    # Measure the relevant sky properties
    mSky = calc_clipped_stats(skyDataLst, clip=3.0)
    medianSkyClp  = mSky['median']
    meanSkyClp    = mSky['mean']
    madSkyClp     = mSky['madfm']
    stdevSkyClp   = mSky['stdev']
    maxSkyClp     = mSky['max']
    minSkyClp     = mSky['min']
    npixSkyClp    = mSky['npix']
    
    # Calculate the integrated flux and error
    # F.Masci, IPAC: 'Flux-Uncertainty from Aperture Photometry'
    sumSrc = np.sum(srcDataLst)
    integFlux = (sumSrc - (medianSkyClp * npixSrc))
    t1 = (npixSrc * madSkyClp)**2.0
    t2 = ( (math.pi * npixSrc**2.0 * madSkyClp**2.0) /
           (2.0 * npixSkyClp) )
    dintegFlux = math.sqrt( t1 + t2)
    
    # Calculate the S/N
    maxSrcAperture = float(np.nanmax(srcDataLst))
    sigmaPix = float((maxSrcAperture - medianSkyClp)/madSkyClp)
    
    # Calculate the brightness weighted centroid
    xSum = 0.0
    ySum = 0.0
    wSum = 0.0
    for l in range(len(srcIndxLst)):
        if srcDataLst[l] > 0.0:
            i, j = srcIndxLst[l]
            xSum += (i*srcDataLst[l])
            ySum += (j*srcDataLst[l])
            wSum += srcDataLst[l]
    wCent_pix = ((xSum/wSum), (ySum/wSum))
    [[raWCent_deg, decWCent_deg]] = wcs.wcs_pix2sky([wCent_pix], 0)
            
    # Calculate the brightness weighted radius
    sepSum = 0.0
    intSum = 0.0
    for l in range(len(srcIndxLst)):
        if srcDataLst[l] >0.0:
            sep = calc_linear_sep(wCent_pix[0], wCent_pix[1],
                               srcIndxLst[l][0], srcIndxLst[l][1])
            sepSum += (srcDataLst[l] * sep)
            intSum += srcDataLst[l]
    radWCent_pix = sepSum / intSum
    radWCent_deg = radWCent_pix * w['pixscale']    
    nanFlag = 0
    if np.any(np.array(srcDataLst)==0.0):
        nanFlag = 1
    
    # Insert the measurements into the measurements dictionary
    measurements = {}
    measurements['SRC FLUX'] = integFlux
    measurements['SRC dFLUX'] = dintegFlux
    measurements['SRC RA_HR'] = x_deg/15.0
    measurements['SRC DEC_DEG'] = y_deg
    measurements['SRC RA_WCENT_HR'] = raWCent_deg/15.0
    measurements['SRC DEC_WCENT_DEG'] = decWCent_deg
    measurements['SRC RAD_WCENT_DEG'] = radWCent_deg
    measurements['SRC MAX'] = maxSrc
    measurements['SRC MIN'] = minSrc
    measurements['SRC MEDIAN'] = medianSrc
    measurements['SRC MEAN'] = meanSrc
    measurements['SRC MADFM'] = madSrc
    measurements['SRC STDEV'] = stdevSrc
    measurements['SRC NPTS'] = npixSrc
    measurements['SKY MAX'] = maxSkyClp
    measurements['SKY MIN'] = minSkyClp
    measurements['SKY MEDIAN'] = medianSkyClp
    measurements['SKY MEAN'] = meanSkyClp
    measurements['SKY MADFM'] = madSkyClp
    measurements['SKY STDEV'] = stdevSkyClp
    measurements['SKY NPTS'] = npixSkyClp
    measurements['SIGMA PIX'] = sigmaPix
    measurements['srcMaj_deg'] = srcMaj_deg
    measurements['srcMin_deg'] = srcMin_deg
    measurements['innerRad_deg'] = innerRad_deg
    measurements['outerRad_deg'] = outerRad_deg
    measurements['nanFlag'] = nanFlag

    return measurements

