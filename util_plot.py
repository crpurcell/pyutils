#!/usr/bin/env python
#=============================================================================#
#                                                                             #
# NAME:     util_plot.py                                                      #
#                                                                             #
# PURPOSE:  Common function for plotting data.                                #
#                                                                             #
# MODIFIED: 14-Oct-2015 by C. Purcell                                         #
#                                                                             #
# CONTENTS:                                                                   #
#                                                                             #
#  bin_xydata_errs ... bin X-Y data given a vector of bin edges               #
#  mk_hist_poly    ... create line-segments for an open histogram             #
# label_format_log ... format labels for a log plot                           #
#                                                                             #
#=============================================================================#
import copy
import os
import sys
import shutil
import math as m
import numpy as np


#-----------------------------------------------------------------------------#
def bin_xydata_errs(x, y, b, mode='mean'):
    """
    Bin X-Y data given a vector of bin edges. Return num and stdev in each 
    bin. Expects x and y vectors (same shape) and a vector of bin edges (N+1).
    """

    binGrid = []
    
    # Loop through lower edged of the first (nBin-1) bins
    for i in range(len(b)-2):

        e = []
        # loop through the input x array and assign y values to a bin
        for j in range(len(x)):
            if x[j] >= b[i] and x[j] < b[i+1]:
                e.append(y[j])
        binGrid.append(e)
                
    # Last bin
    i = len(b)-2
    e = []
    for j in range(len(x)):
        if x[j] >= b[i] and x[j] <= b[i+1]:
            e.append(y[j])
    binGrid.append(e)
    n = []
    e = []
    if mode ==  'sum':
        for i in range(len(b)-1):
            n.append(sum(binGrid[i]))
            e.append(0.0)
    elif mode ==  'mean':
        for i in range(len(b)-1):
            n.append(np.mean(binGrid[i]))
            e.append(np.std(binGrid[i]))
    elif mode == 'median':
        for i in range(len(b)-1):
            n.append(np.median(binGrid[i]))
            e.append(np.std(binGrid[i]))
    else:
        for i in range(len(b)-1):
            n.append(len(binGrid[i]))
            e.append(m.sqrt(n[-1]))

    return b, n, e


#-----------------------------------------------------------------------------#
def mk_hist_poly(bins, n, logScaleY=False, zeroPt=0.8, addZeroPt=True):
    """Create the line segments for the a polygon used to draw a histogram"""

    if logScaleY is True:
        for i in range(len(n)):
            if n[i] <= 0.0:
                n[i] = zeroPt 
    else:
        zeroPt = 0.0

    # Starting position
    polyCoordLst = []
    if addZeroPt:
        polyCoordLst.append([bins[0],zeroPt])

    # Form the line segments
    i = 0
    j = 0
    while i <= len(bins)-1:
        if j < len(n):
            polyCoordLst.append([bins[i],n[j]])
        if i == j:
            i += 1
        else:
            j += 1

    # Ground the polygon line and close
    if addZeroPt:
        polyCoordLst.append([bins[-1],zeroPt])
        polyCoordLst.append([bins[0],zeroPt])
    polyCoords = np.array(polyCoordLst)
                        
    return polyCoords


#-----------------------------------------------------------------------------#
def label_format_log(switchExp=3.0):
    """Return a function to format labels for log axes. Switches to power
    format for numbers greater than switchExp power of 10."""    

    def rfunc(num, pos=None):
        if num > 0.0:
            exponent = m.log10(num)
            if exponent>=0:
                if exponent<switchExp:
                    return "%.0f" % num    
                else:
                    return r'$10^{%i}$' % (m.log10(num))
            else:
                if exponent>-switchExp:
                    format_code = "%"+".%sf" % (str(abs(int(exponent))))
                    return format_code % num
                else:
                    return r'$10^{%i}$' % (m.log10(num))
    return rfunc
                
