#!/usr/bin/env python
#=============================================================================#
#                                                                             #
# NAME:     util_misc.py                                                      #
#                                                                             #
# PURPOSE:  Common function for inclusion in python scripts.                  #
#                                                                             #
# MODIFIED: 14-May-2015 by C. Purcell                                         #
#                                                                             #
# CONTENTS:                                                                   #
#                                                                             #
#  message           ... print a line of text delimited by string of ---      #
#  run_command       ... run a system command using os.system                 #
#  config_read       ... parse a key=value configuration file                 #
#  ascii_dat_read    ... read the columns from an ascii or csv file           #
#  lst_dir           ... list and filter the contents of the current dir      #
#  lst_to_float      ... convert the elements in a list to floats             #
#  dms2deg           ... convert sexigesimal to degrees                       #
#  deg2dms           ... convert degrees to sexigesimal                       #
#  calc_linear_sep   ... linear separation between two points                 #
#  get_separation    ... seperation on a curved surface                       #
#  check_sep         ... return True if sep < minSep_deg                      #
#  numeric_compare   ... with 'sort'  to sort numeric strings as numbers      #
#  num2str           ... format a number to a string                          #
#  nanmean           ... np.mean ignoring NaNs                                #
#  nanmedian         ... np.median ignoring NaNs                              #
#  nanstd            ... np.std ignoring NaNs                                 #
#  MAD               ... calculate the madfm                                  #
#  resample          ... resample an array by a factor                        #
#  calc_stats_old    ... calculate the statistics of an array (fails)         #
#  calc_stats        ... calculate the statistics of an array (uses masking)  #
#  calc_clipped_stats .. calculate the stats after sigma-clipping             #
#  polystr_to_arr    ... convert a string of coords to a polygon vertex array #
#  check_overlap     ... check for overlap of two squares on sky              #
#  sort_nicely       ... sort a list in human order                           #
#  select_into_namedlst . select from a DB into a named list                  #
#  rowlst_to_rowdict ... convert a named list into a named dictionary         #
#  env_var_app       ... append a string to an environment variable           #
#  env_var_pre       ... prepend a string to an environment variable          #
#  spawnDaemon       ... call an executable file as a detached procress       #
#  to_precision      ... return a formateed string representation of x        #
#                                                                             #
#=============================================================================#

# Import standard python modules
import os
import sys
import shutil
import copy
import re
import math as m
import numpy as np
import numpy.ma as ma
from scipy import median
from scipy import stats
from scipy.stats import norm
import pyslalib as sla


#-----------------------------------------------------------------------------#
# Feedback messages to the user delimited by two lines of '-' signs.          #
#-----------------------------------------------------------------------------#
def message (string):
    print "\n"+"-"*80
    print ">>> %s" % string
    print "-"*80


#-----------------------------------------------------------------------------#
# Execute system command (e.g., a miriad task)                                #
#-----------------------------------------------------------------------------#
def run_command(command, doPrint=True):

    # Clean up the spaces in the command
    spaces = re.compile('\s+')
    command = command.strip()
    command = spaces.sub(' ', command)

    # Print the command to screen and execute
    if doPrint:
        print "-"*80
        print ">", command
        print "-"*80
    os.system(command)


#-----------------------------------------------------------------------------#
# Read a configuration file and output a 'KEY=VALUE' dictionary               #
#-----------------------------------------------------------------------------#
def config_read(filename, delim='=', doValueSplit=True):
    
    configTable = {}
    CONFIGFILE = open(filename, "r")
    
    # Compile a few useful regular expressions
    spaces = re.compile('\s+')
    commaAndSpaces = re.compile(',\s+')
    commaOrSpace = re.compile('[\s|,]')
    brackets = re.compile('[\[|\]\(|\)|\{|\}]')
    comment = re.compile('#.*')
    quotes = re.compile('\'[^\']*\'')
    keyVal = re.compile('^.+' + delim + '.+')

    # Read in the input file, line by line
    for line in CONFIGFILE:

        valueLst=[]
        line = line.rstrip("\n\r")

        # Filter for comments and blank lines
        if not comment.match(line) and keyVal.match(line):

            # Weed out internal comments & split on 1st space
            line = comment.sub('',line)
            (keyword, value) = line.split(delim,1)

            # If the line contains a value
            keyword = keyword.strip()              # kill external whitespace
            keyword = spaces.sub('', keyword)      # kill internal whitespaces
            value = value.strip()                  # kill external whitespace
            value = spaces.sub(' ', value)         # shrink internal whitespace
            value = value.replace("'", '')         # kill quotes
            value = commaAndSpaces.sub(',', value) # kill ambiguous spaces

            # Split comma/space delimited value strings
            if doValueSplit:
                valueLst = commaOrSpace.split(value)
                if len(valueLst)<=1:
                    valueLst = valueLst[0]
                configTable[keyword] = valueLst
            else:
                configTable[keyword] = value

    return configTable


#-----------------------------------------------------------------------------#
# Read the columns in an ASCII or CSV file                                    #
#-----------------------------------------------------------------------------#
def ascii_dat_read(filename, delim=",", doNone=False, doFloatCols=[]):

    columnDict = {}
    configTable = {}

    DATFILE = open(filename, "r")
     
    # Compile a few useful regular expressions
    spaces = re.compile('\s+')
    comma_and_spaces = re.compile(',\s+')
    comma_or_space = re.compile('[\s|,]')
    brackets = re.compile('[\[|\]\(|\)|\{|\}]')
    comment = re.compile('#.*')
    quotes = re.compile('\'[^\']*\'')
    keyVal = re.compile('^.+=.+')
    words = re.compile('\S+')

    nValid = 0
    
    # Read in the input file, line by line
    for line in DATFILE:

        line = line.rstrip("\n\r")
        # Filter for comments and blank lines
        #if not comment.match(line) and words.match(line):
        if not comment.match(line):
            
            # Weed out internal comments 
            line = comment.sub('', line)
            line = line.strip()              # kill external whitespace
            line = spaces.sub(' ', line)     # shrink internal whitespace
            
            if keyVal.match(line):
                keyword, value = line.split('=',1)
                keyword = keyword.strip()          # kill external whitespace
                keyword = spaces.sub('', keyword)  # kill internal whitespaces
                if value:
                    configTable[keyword] = value
                
            else:

                # Split on commas
                line = line.split(delim)
                ncols = len(line)
                nValid += 1

                # Add the value to the column (padding if first instance)
                for i in range(len(line)):
                    if columnDict.has_key(i+1):
                        columnDict[i+1].append(line[i])
                    else:
                        columnDict[i+1] = ['']*(nValid-1) + [line[i]]
                        maxCols = ncols
                        
    # Pad the last columns if necessary
    for k in columnDict.keys():
        if len(columnDict[k]) < nValid:
            columnDict[k] += ['']*(nValid-len(columnDict[k]))

    # Convert to floats
    if len(doFloatCols)>0:
        for k in doFloatCols:
            for i in range(len(columnDict[k])):
                try:
                    columnDict[k][i] = float(columnDict[k][i])
                except Exception:
                    pass
        
    # Convert '' to None
    if doNone:
        for k in columnDict.keys():
            for i in range(len(columnDict[k])):
                if columnDict[k][i]=='':
                    columnDict[k][i] = None
            
    DATFILE.close()
    
    return columnDict, configTable


#-----------------------------------------------------------------------------#
# Convert the elements in a list to floats, catch errors                      #
#-----------------------------------------------------------------------------#
def lst_to_float(lst, blank='none'):
    for i in range(len(lst)):
        try:
            lst[i] = float(lst[i])
        except Exception:
            if blank=='none':
                lst[i] = None
            else:
                pass

    return lst
    
#-----------------------------------------------------------------------------#
# List the files in the current directory, filter according to a regexp       #
#-----------------------------------------------------------------------------#
def lst_dir(dir,regexp='^.*\..*$'):
    
    dirLst = os.listdir(dir)
    pattern = re.compile(regexp)
    lst = []
    for file_name in dirLst:
        if pattern.match(file_name):
            lst.append(file_name)
    lst.sort()

    return lst


#-----------------------------------------------------------------------------#
# Convert DD:MM:DD.D coordinates to decimal degrees.                          #
#-----------------------------------------------------------------------------#
def dms2deg(dms):
    
    try:
        delim = re.compile('[,| |:|h|d|m|s]')
        dmsLst = delim.split(dms)
        d = float(dmsLst[0])
        m = float(dmsLst[1])
        s = float(dmsLst[2])
        sign = 1
        if d != 0.0:
            sign = int( d/abs(d) )
        return sign * ( abs(d) + m / 60.0 + s / 3600.0 )
    except Exception:
        return None
    

#-----------------------------------------------------------------------------#
# Convert a float in degrees to 'dd mm ss' format                             #
#-----------------------------------------------------------------------------#
def deg2dms(deg, delim=':', doSign=False, nPlaces=2):

    try:
        angle = abs(deg)
        sign=1
        if angle!=0: sign = angle/deg
        
        # Calcuate the degrees, min and sec
        dd = int(angle)
        rmndr = 60.0*(angle - dd)
        mm = int(rmndr)
        ss = 60.0*(rmndr-mm)

        # If rounding up to 60, carry to the next term
        if float("%05.2f" % ss) >=60.0:
            mm+=1.0
            ss = ss - 60.0
        if float("%02d" % mm) >=60.0:
            dd+=1.0
            mm = mm -60.0
        if nPlaces> 0:
            formatCode = "%0" + "%s.%sf" % (str(2 + nPlaces + 1), str(nPlaces))
        else:
            formatCode = "%02.0f"
        if sign>0:
            if doSign:
                formatCode = "+%02d%s%02d%s" + formatCode
            else:
                formatCode = "%02d%s%02d%s" + formatCode
        else:
            formatCode = "-%02d%s%02d%s" + formatCode
        return formatCode % (dd, delim, mm, delim, ss)
        
    except Exception:
        return None

    
#-----------------------------------------------------------------------------#
# Calculate the linear separation between two points                          #
#-----------------------------------------------------------------------------#
def calc_linear_sep(x1, y1, x2, y2):
    
    x1 = float(x1)
    y1 = float(y1)
    x2 = float(x2)
    y2 = float(y2)
        
    return m.sqrt( m.pow( (x2 - x1), 2.0) + m.pow( (y2 -y1), 2.0) )


#-----------------------------------------------------------------------------#
# Calculate the seperation on a curved surface                                #
#-----------------------------------------------------------------------------#
def get_separation(x1_deg, y1_deg, x2_deg, y2_deg):

    sep_deg = None
    try:
        sep_deg = m.degrees(sla.slalib.sla_sep(m.radians(x1_deg),
                                               m.radians(y1_deg),
                                               m.radians(x2_deg),
                                               m.radians(y2_deg)))
    except Exception:
        pass
    
    return sep_deg

    # OLD INCORRECT FUNCTION. TODO: FIX
    #sep_deg = m.degrees( m.acos( m.cos( m.radians(90.0 - x1_deg) )
    #                             * m.cos( m.radians(90.0 - x2_deg) )
    #                             + m.sin( m.radians(90.0 - x1_deg) )
    #                             * m.sin( m.radians(90.0 - x2_deg) )
    #                             * m.cos( m.radians(y2_deg - y1_deg) )
    #                             ) )


#-----------------------------------------------------------------------------#
# Calculate the seperation between two points and check if same               #
#-----------------------------------------------------------------------------#
def check_sep(x1_deg, y1_deg, x2_deg, y2_deg, minSep_deg):

    sep_deg = get_separation(x1_deg, y1_deg, x2_deg, y2_deg)

    # Return True if the two points are closer than min_sep_deg
    if sep_deg <= minSep_deg:
        return True
    else:
        return False


#-----------------------------------------------------------------------------#
# Used with 'sort' to sort numeric strings as numbers                         #
#-----------------------------------------------------------------------------#
def numeric_compare(x, y):
    
    if float(x)>float(y):
        return 1
    elif float(x)==float(y):
        return 0
    else: # x<y
        return -1


#-----------------------------------------------------------------------------#
# Format a number to a string                                                 #
#-----------------------------------------------------------------------------#
def num2str(num, pre=2, post=2, dosign=False, do_zeros=False, doplus=True):
    
    num=float(num)
    if do_zeros:
        lead = '%0'
    else: lead = '%'
    
    if dosign:
        sign=1
        if not num==0:
            if abs(num)/num<0:
                sign=0
        if sign>0:
            if doplus and do_zeros:
                lead = '+' + lead
            else:
                lead = ' '+lead
            format_code = "%s%s.%sf" %(lead, str(pre+post+1), str(post))
        else:
            format_code = "%s%s.%sf" %(lead, str(pre+post+2), str(post))
    else:
        format_code = "%s%s.%sf" %(lead, str(pre+post+1), str(post))
        
    return format_code % num


#-----------------------------------------------------------------------------#
def nanmean(arr, **kwargs):
    """
    Returns mean ignoring NAN
    """
    
    return ma.mean( ma.masked_where(arr!=arr, arr), **kwargs )


#-----------------------------------------------------------------------------#
# Returns median ignoring NAN                                                 #
#-----------------------------------------------------------------------------#
def nanmedian(arr, **kwargs):
    """
    Returns median ignoring NAN
    """
    
    return ma.median( ma.masked_where(arr!=arr, arr), **kwargs )


#-----------------------------------------------------------------------------#
def nanstd(arr, **kwargs):
    """
    Returns stdev ignoring NAN
    """
    
    return ma.std( ma.masked_where(arr!=arr, arr), **kwargs )


#-----------------------------------------------------------------------------#
# Median Absolute Deviation From Median (MADFM) - Adam Ginsburg's code.       #
#-----------------------------------------------------------------------------#
def MAD(a, c=0.6745, axis=None):
    """
    Median Absolute Deviation along given axis of an array:
    median(abs(a - median(a))) / c
    c = 0.6745 is the constant to convert from MAD to std

    """
    
    a = ma.masked_where(a!=a, a)
    if a.ndim == 1:
        d = ma.median(a)
        m = ma.median(ma.fabs(a - d) / c)
    else:
        d = ma.median(a, axis=axis)
        # I don't want the array to change so I have to copy it?
        if axis > 0:
            aswp = ma.swapaxes(a,0,axis)
        else:
            aswp = a
        m = ma.median(ma.fabs(aswp - d) / c, axis=0)

    return m


#-----------------------------------------------------------------------------#
# Resample the data array by a factor. Copied from APLpy                      #
#-----------------------------------------------------------------------------#
def resample(array, factor):
    
    nx, ny = np.shape(array)
    nx_new = nx / factor
    ny_new = ny / factor
    array2 = np.empty((nx_new, ny))
    
    for i in range(nx_new-1):
        array2[i,:] = np.sum(array[i*factor:i*factor+1, :], axis=0)

    array3 = np.empty((nx_new, ny_new))
    for j in range(ny_new-1):
        array3[:, j] = np.sum(array2[:, j*factor:j*factor+1], axis=1)

    return array3


#-----------------------------------------------------------------------------#
# Calculate the statistics of an array                                        #
#-----------------------------------------------------------------------------#
def calc_stats_old(a, maskzero=False):
    
    statsDict = {}
    a = np.array(a)
    if maskzero:
        a = np.where( np.equal(a, 0.0), np.nan, a)

    # Check that array is not all NaNs
    statsDict['npix'] = int(np.sum(np.where(np.isnan(a),0.0,1.0)))
    if statsDict['npix']>=2:
        statsDict['stdev'] = float(stats.nanstd(a.flatten()))
        statsDict['mean'] = float(stats.nanmean(a.flatten()))
        statsDict['median'] = float(stats.nanmedian(a.flatten()))
        statsDict['max'] = float(np.nanmax(a))
        statsDict['min'] = float(np.nanmin(a))
        statsDict['centmax'] = list(np.unravel_index(np.nanargmax(a),
                                                     a.shape))
        statsDict['madfm'] = float(MAD(a.flatten()))
        statsDict['npix'] = int(np.sum(np.where(np.isnan(a),0.0,1.0)))
        statsDict['success'] = True
        
    else:
        statsDict['npix'] == 0
        statsDict['stdev']   = 0.0
        statsDict['mean']    = 0.0
        statsDict['median']  = 0.0
        statsDict['max']     = 0.0
        statsDict['min']     = 0.0
        statsDict['centmax'] = (0.0, 0.0)
        statsDict['madfm']   = 0.0
        statsDict['success'] = False
        
    return statsDict


#-----------------------------------------------------------------------------#
# Calculate the statistics of an array (masked version)                       #
#-----------------------------------------------------------------------------#
def calc_stats(a, maskzero=False):
    """
    Calculate the statistics of an array (masked version).
    """
    
    statsDict = {}
    a = np.array(a)

    # Mask off bad values and count valid pixels
    if maskzero:
        a = np.where( np.equal(a, 0.0), np.nan, a)
    am = ma.masked_invalid(a)
    statsDict['npix'] = np.sum(~am.mask)
    
    if statsDict['npix']>=2:
        statsDict['stdev'] = float(np.std(am))
        statsDict['mean'] = float(np.mean(am))
        statsDict['median'] = float(nanmedian(am))
        statsDict['max'] = float(np.max(am))
        statsDict['min'] = float(np.min(am))
        statsDict['centmax'] = list(np.unravel_index(np.argmax(am),
                                                     a.shape))
        statsDict['madfm'] = float(MAD(am.flatten()))
        statsDict['success'] = True
        
    else:
        statsDict['npix'] == 0
        statsDict['stdev']   = 0.0
        statsDict['mean']    = 0.0
        statsDict['median']  = 0.0
        statsDict['max']     = 0.0
        statsDict['min']     = 0.0
        statsDict['centmax'] = (0.0, 0.0)
        statsDict['madfm']   = 0.0
        statsDict['success'] = False
        
    return statsDict


#-----------------------------------------------------------------------------#
# Calculate the mean and stdev of an array given a sigma clip                 #
#-----------------------------------------------------------------------------#
def calc_clipped_stats(data, clip=3.0, nIter=10, maskzero=False):

    if maskzero:
        data = np.where( np.equal(data, 0.0), np.nan, data)
        
    ms = calc_stats(data)

    if clip>0 and nIter>0:
        convergeFlg = 0
        itCnt = 0
        while convergeFlg==0 and itCnt<nIter:
            meanOld, stdOld, madOld = ms['mean'], ms['stdev'], ms['madfm']
            minVal = ms['mean'] - (clip * ms['madfm'])
            maxVal = ms['mean'] + (clip * ms['madfm'])
            
            # Blank values outside the 3-sigma range
            dataMsk = np.where(np.greater(data, maxVal), np.nan, data)
            dataMsk = np.where(np.less(data, minVal), np.nan, dataMsk)
            
            # Measure the statistics
            ms = calc_stats(dataMsk)
            dataMsk = []
    
            if ms['mean'] == meanOld and ms['madfm'] == madOld:
                convergFlg = 1
            itCnt += 1

    return ms


#-----------------------------------------------------------------------------#
# Calculate the mean and stdev of an array given a sigma clip                 #
#-----------------------------------------------------------------------------#
def calc_clipped_stats_old(data, clip=3.0, nIter=10):
    
    data = np.array(data).flatten()
    
    mean = float(stats.nanmean(data))
    std = float(stats.nanstd(data))
    mad = float(MAD(data))

    if clip > 0.0:
        convergeFlg = 0
        itCnt = 0
        while convergeFlg==0 and itCnt<nIter:
            meanOld, stdOld, madOld = mean, std, mad
            minVal = mean - (clip * mad)
            maxVal = mean + (clip * mad)

            # Blank values outside the 3-sigma range
            dataMsk = np.where(np.greater(data, maxVal), np.nan, data)
            dataMsk = np.where(np.less(data, minVal), np.nan, dataMsk)

            # Measure the statistics
            mean = stats.nanmean(dataMsk)
            median = stats.nanmedian(dataMsk)
            std = stats.nanstd(dataMsk)
            mad = MAD(dataMsk)
            npix = np.sum(np.where(np.isnan(dataMsk),0.0,1.0))
            dataMsk = []
            
            if mean == meanOld and mad == madOld:
                convergFlg = 1
            itCnt += 1
            

    # Assemble the measurements into a dictionary
    m = {}
    m['mean'] = float(mean)
    m['median'] = float(median)
    m['stdev'] = float(std)
    m['madfm'] = float(mad)
    m['npix'] =int(npix)
    m['max'] = float(np.nanmax(data))
    m['min'] = float(np.nanmin(data))
    del data
    
    # If all nans
    if m['npix'] == 0:
        m['stdev'] = 0.0
        m['mean'] = 0.0
        m['median'] = 0.0
        m['max'] = 0.0
        m['min'] = 0.0
        m['centmax'] = (0.0,0.0)
        m['madfm'] = 0.0
        m['success'] = False
    else:
        m['success'] = True

    return m


#-----------------------------------------------------------------------------#
# Convert the polygon strings to numpy arrays 2XN                             #
#-----------------------------------------------------------------------------#
def polystr_to_arr(polystr):
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
# Convert the polygon strings to numpy arrays 2XN                             #
#-----------------------------------------------------------------------------#
def parse_poly_string(poly):
    return polystr_to_arr(poly)


#-----------------------------------------------------------------------------#
# Check for overlap of two squares in RA and Dec or l and b.                  #
#-----------------------------------------------------------------------------#
def check_overlap(x1, y1, width1, height1, x2, y2, width2, height2):

    xOverlap=False; yOverlap=False

    # Calculate the top and bottom Y
    y1_top=y1+height1/2.0        # Tile Y
    y1_bot=y1-height1/2.0
    y2_top=y2+height2/2.0        # Field Y
    y2_bot=y2-height2/2.0

    # Measure the positions of the top and bottom limits of the squares
    top2_lt_top1 = y2_top <= y1_top
    top2_gt_bot1 = y2_top >= y1_bot
    bot2_gt_bot1 = y2_bot >= y1_bot
    bot2_lt_top1 = y2_bot <= y1_top
    top1_lt_top2 = y1_top <= y2_top
    top1_gt_bot2 = y1_top >= y2_bot
    bot1_gt_bot2 = y1_bot >= y2_bot
    bot1_lt_top2 = y1_bot <= y2_top
    
    # Check for overlap of the line segments
    if top2_lt_top1 and top2_gt_bot1: yOverlap=True
    if bot2_gt_bot1 and bot2_lt_top1: yOverlap=True    
    if top1_lt_top2 and top1_gt_bot2: yOverlap=True
    if bot1_gt_bot2 and bot1_lt_top2: yOverlap=True

    # Test for X overlap, if a Y overlap is found
    if yOverlap:

    #
    #     |----|        
    # |---|-|  |  <--   Because X varies with Y, need to check 
    # |   |-|--|  <--   for overlap at the two middle Y boundaries
    # |     |
    # |-----|
    # 
        # Find the middle Y coordinates
        y_list = [y1_top,y1_bot,y2_top,y2_bot]
        y_list.sort()
        
        # 1) Top Y position:

        # Calculate the left and right X
        x1_left=x1+(width1/2.0)/m.cos(m.radians(y_list[2]))
        x1_right=x1-(width1/2.0)/m.cos(m.radians(y_list[2]))
        x2_left=x2+(width2/2.0)/m.cos(m.radians(y_list[2]))
        x2_right=x2-(width2/2.0)/m.cos(m.radians(y_list[2]))
        
        # Measure the positions 
        left2_lt_left1 = x2_left <= x1_left
        left2_gt_right1 = x2_left >= x1_right
        right2_gt_right1 = x2_right >= x1_right
        right2_lt_left1 = x2_right <= x1_left

        left1_lt_left2 = x1_left <= x2_left
        left1_gt_right2 = x1_left >= x2_right
        right1_gt_right2 = x1_right >= x2_right
        right1_lt_left2 = x1_right <= x2_left
        
        # Check for overlap of the line segments
        if left2_lt_left1 and left2_gt_right1: xOverlap=True
        if right2_gt_right1 and right2_lt_left1: xOverlap=True
        if left1_lt_left2 and left1_gt_right2: xOverlap=True
        if right1_gt_right2 and right1_lt_left2: xOverlap=True
        
        # 2) Bottom Y position:
        
        # Calculate the left and right X
        x1_left=x1+(width1/2.0)/m.cos(m.radians(y_list[1]))
        x1_right=x1-(width1/2.0)/m.cos(m.radians(y_list[1]))
        x2_left=x2+(width2/2.0)/m.cos(m.radians(y_list[1]))
        x2_right=x2-(width2/2.0)/m.cos(m.radians(y_list[1]))
        
        # Measure the positions 
        left2_lt_left1 = x2_left <= x1_left
        left2_gt_right1 = x2_left >= x1_right
        right2_gt_right1 = x2_right >= x1_right
        right2_lt_left1 = x2_right <= x1_left

        left1_lt_left2 = x1_left <= x2_left
        left1_gt_right2 = x1_left >= x2_right
        right1_gt_right2 = x1_right >= x2_right
        right1_lt_left2 = x1_right <= x2_left
        
        # Check for overlap of the line segments
        if left2_lt_left1 and left2_gt_right1: xOverlap=True
        if right2_gt_right1 and right2_lt_left1: xOverlap=True
        if left1_lt_left2 and left1_gt_right2: xOverlap=True
        if right1_gt_right2 and right1_lt_left2: xOverlap=True
        
    if (xOverlap and yOverlap):
        return True
    else:
        return False


#-----------------------------------------------------------------------------#
# Sort in human order                                                         #
#-----------------------------------------------------------------------------#
def sort_nicely( l ):
    
  convert = lambda text: int(text) if text.isdigit() else text 
  alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
  l.sort( key=alphanum_key ) 


#-----------------------------------------------------------------------------#
# Sort multiple lists                                                         #
#-----------------------------------------------------------------------------#
def sort_lists_by(lists, key_list=0, desc=False):
    return zip(*sorted(izip(*lists), reverse=desc,
                 key=lambda x: x[key_list]))


#-----------------------------------------------------------------------------#
# Store the results of a select statement into a list of dictionaries         #
#-----------------------------------------------------------------------------#
def select_into_namedlst(cursor, sql, args=[]):

    rowLst = []
    if args == []:
        cursor.execute(sql)
    else:
        cursor.execute(sql, tuple(args))
    columnNameLst = zip(*cursor.description)[0]
    rows = cursor.fetchall()
    for row in rows:
        e = {}
        for i in range(len(columnNameLst)):
            e[columnNameLst[i]] = row[i]
        rowLst.append(e)
    return rowLst


#-----------------------------------------------------------------------------#
# Convert a list of dictionaries into a dictionary of dictionaries            #
#-----------------------------------------------------------------------------#
def rowlst_to_rowdict(rowLst, keyName):
    rowDict = {}
    for e in rowLst:
        key = e.pop(keyName)
        rowDict[key] = e

    return rowDict


#-----------------------------------------------------------------------------#
# Append a string to an environment variable                                  #
#-----------------------------------------------------------------------------#
def env_var_app(key,value,delim=':'):
    if os.environ.has_key(key):
        os.environ[key]+=delim+str(value)
    else:
         os.environ[key]=str(value)


#-----------------------------------------------------------------------------#
# Prepend a string to an environment variable                                 #
#-----------------------------------------------------------------------------#
def env_var_pre(key,value,delim=':'):
    if os.environ.has_key(key):
        os.environ[key]=str(value)+delim+os.environ[key]
    else:
         os.environ[key]=str(value)

    
#-----------------------------------------------------------------------------#
# Call an executable file as a detached procress                              #
#-----------------------------------------------------------------------------#
def spawnDaemon(path_to_executable, *args):
    """Spawn a completely detached subprocess (i.e., a daemon).
    E.g.:
    spawnDaemon("../bin/producenotify.py", "producenotify.py", "xx")
    """

    # The standard I/O file descriptors are redirected to /dev/null.
    if (hasattr(os, "devnull")):
       REDIRECT_TO = os.devnull
    else:
       REDIRECT_TO = "/dev/null"
    
    # fork the first time (to make a non-session-leader child process)
    try:
        pid = os.fork()
    except OSError, e:
        raise RuntimeError("1st fork failed: %s [%d]" % (e.strerror, e.errno))
    if pid != 0:
        # parent (calling) process is all done
        return

    # detach from controlling terminal (to make child a session-leader)
    os.setsid()
    try:
        pid = os.fork()
    except OSError, e:
        raise RuntimeError("2nd fork failed: %s [%d]" % (e.strerror, e.errno))
        raise Exception, "%s [%d]" % (e.strerror, e.errno)
    if pid != 0:
        # child process is all done
        os._exit(0)
        
    # grandchild process now non-session-leader, detached from parent
    # grandchild process must now close all open files
    try:
        maxfd = os.sysconf("SC_OPEN_MAX")
    except (AttributeError, ValueError):
        maxfd = 1024

    for fd in range(0,maxfd):
        try:
           os.close(fd)
        except OSError: # ERROR, fd wasn't open to begin with (ignored)
           pass

    # redirect stdin, stdout and stderr to /dev/null
    os.open(REDIRECT_TO, os.O_RDWR) # standard input (0)
    os.dup2(0, 1)
    os.dup2(0, 2)

    # and finally let's execute the executable for the daemon!
    try:
       
       os.execv(path_to_executable, (path_to_executable,) +  tuple(args))
    except Exception, e:
       # oops, we're cut off from the world, let's just give up
       os._exit(255)


#-----------------------------------------------------------------------------#
def read_poly_file(fileName):

    polyLst = []
    
    # Compile a few useful regular expressions
    spaces = re.compile('\s+')
    commaAndSpaces = re.compile(',\s+')
    commaOrSpace = re.compile('[\s|,]')
    brackets = re.compile('[\[|\]\(|\)|\{|\}]')
    comment = re.compile('#.*')
    quotes = re.compile('\'[^\']*\'')
    keyVal = re.compile('^.+=.+')

    # Read in the input file, line by line
    FH = open(fileName, "r")
    for line in FH:
        
        line = line.rstrip("\n\r")

        
        # Filter for comments and blank lines
        if not comment.match(line) and keyVal.match(line):
        
            # Weed out internal comments & split
            line = comment.sub('',line)
            (keyword, value) = line.split('=',1)

            # Clean up syntax
            keyword = keyword.strip()              # kill external whitespace
            keyword = spaces.sub('', keyword)      # kill internal whitespaces
            value = value.strip()                  # kill external whitespace
            value = spaces.sub(' ', value)         # shrink internal whitespace
            value = value.replace("'", '')         # kill quotes
            value = commaAndSpaces.sub(',', value) # kill ambiguous spaces
            value = brackets.sub('', value)

            if keyword=='POLY':
                polyLst.append(value)
                
    return polyLst


#-----------------------------------------------------------------------------#
def to_precision(x, p):
    """
    returns a string representation of x formatted with a precision of p

    Based on the webkit javascript implementation taken from here:
    https://code.google.com/p/webkit-mirror/source/browse/JavaScriptCore/kjs/number_object.cpp

    Python code from:
    http://randlet.com/blog/python-significant-figures-format/
    """

    x = float(x)

    if x == 0.:
        return "0." + "0"*(p-1)

    out = []

    if x < 0:
        out.append("-")
        x = -x

    e = int(m.log10(x))
    tens = m.pow(10, e - p + 1)
    n = m.floor(x/tens)

    if n < m.pow(10, p - 1):
        e = e -1
        tens = m.pow(10, e - p+1)
        n = m.floor(x / tens)

    if abs((n + 1.) * tens - x) <= abs(n * tens -x):
        n = n + 1

    if n >= m.pow(10,p):
        n = n / 10.
        e = e + 1

    m1 = "%.*g" % (p, n)

    if e < -2 or e >= p:
        out.append(m1[0])
        if p > 1:
            out.append(".")
            out.extend(m1[1:p])
        out.append('e')
        if e > 0:
            out.append("+")
        out.append(str(e))
    elif e == (p -1):
        out.append(m1)
    elif e >= 0:
        out.append(m1[:e+1])
        if e+1 < len(m1):
            out.append(".")
            out.extend(m1[e+1:])
    else:
        out.append("0.")
        out.extend(["0"]*-(e+1))
        out.append(m1)

    return "".join(out)
