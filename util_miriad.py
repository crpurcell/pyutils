#!/usr/bin/python
#=============================================================================#
#                                                                             #
# NAME:     util_miriad.py                                                    #
#                                                                             #
# PURPOSE:  Common functions utilising MIIRAD commands.                       #
#                                                                             #
# MODIFIED: 14-May-2015 by C. Purcell                                         #
#                                                                             #
# CONTENTS:                                                                   #
#                                                                             #
# set_miriad_env     ... setup the MIRIAD environment                         #
# env_var_app        ... append a string to a system  environment variable    #
# env_var_pre        ... prepend a string to a system environment variable    #
# run_command        ... run a system command                                 #
# miriad_imcomb      ... use MIRIAD Imcomb to combine FITS images             #
# miriad_gal_j2000   ... use MIRIAD Regrid to convert GAL->J2000              #
# miriad_gal_j20001  ... use MIRIAD Regrid to convert GAL->J2000 (2012 ver)   #
# miriad_regrid_tin  ... use MIRIAD Regrid to match FITS file gridding        #
# miriad_regrid_tin1 ... as above, but 2012 version                           #
# fits_to_miriad     ... convert a FITS file to a MIRIAD dataset              #
# miriad_to_fits     ... convert a MIRIAD dataset to a FITS file              #
#                                                                             #
#=============================================================================#
import os
import sys
import re
import string
import shutil
import random
import commands


#-----------------------------------------------------------------------------#
def set_miriad_env(mirHome,arch='linux64'):
    """Set up the MIRIAD environment"""
    env_var_pre('MIRCAT',mirHome+'/cat')
    env_var_pre('MIRINC',mirHome+'/inc')
    env_var_pre('MIRPROG',mirHome+'/prog')
    env_var_pre('MIRSUBS',mirHome+'/subs')
    env_var_pre('MIRCAT',mirHome+'/cat')
    env_var_pre('MIRARCH',arch)
    env_var_pre('MIRBIN',mirHome+'/'+arch+'/bin')
    env_var_pre('MIRLIB',mirHome+'/'+arch+'/lib')
    env_var_pre('PATH',mirHome+'/'+arch+'/bin')


#-----------------------------------------------------------------------------#
def env_var_app(key,value,delim=':'):
    """Append values to (existing) environment variables"""
    
    if os.environ.has_key(key):
        os.environ[key]+=delim+str(value)
    else:
         os.environ[key]=str(value)


#-----------------------------------------------------------------------------#
def env_var_pre(key,value,delim=':'):
    """Prepend values to (existing) environment variables"""

    if os.environ.has_key(key):
        os.environ[key]=str(value)+delim+os.environ[key]
    else:
         os.environ[key]=str(value)


#-----------------------------------------------------------------------------#
def run_command(command, doPrint=True, verb=True):
    """Run a system command and print a message"""

    # Clean up the spaces in the command
    spaces = re.compile('\s+')
    command = command.strip()
    command = spaces.sub(' ', command)

    # Print the command to screen and execute
    if doPrint:
        print "-"*80
        print ">", command
        print "-"*80
    status,msg = commands.getstatusoutput(command)
    
    if verb:
        print msg

    return status, msg


#-----------------------------------------------------------------------------#
def miriad_imcomb(inFitsFiles, outFitsFile):
    """Use the MIRIAD imcomb task to combine multiple images"""
    
    outDir = ''
    inDirs=[]
    inFiles=[]
    stitchList=[]

    # Separate the paths and filenames
    outDir, outFile = os.path.split(outFitsFile)
    if outDir=='':
        outDir = '.'
    for i in range(len(inFitsFiles)):
        inDir, inFile = os.path.split(inFitsFiles[i])
        if inDir=='':
            inDir = '.'
        inDirs.append(inDir)
        inFiles.append(inFile)
        
    # Create a TMP directory to work in
    randStr = "IMCOMB"+"".join(random.sample(string.letters+string.digits, 4))
    tmpDir = outDir + '/' + randStr
    if os.path.exists(tmpDir):
        shutil.rmtree(tmpDir, True)
    
    # Change to the working directory
    startDir = os.getcwd()
    os.mkdir(tmpDir)
    os.chdir(tmpDir)

    # Load and regrid the files   
    try:
        for i in range(len(inFiles)):
                    
            # Make a local copy to operate on
            if os.path.exists(inFiles[i]):
                os.remove(inFiles[i])
            shutil.copy(startDir + '/' + inFitsFiles[i], inFiles[i])

            # Convert to MIRIAD format and regrid
            if os.path.exists(inFiles[i]+'.xy'): 
                shutil.rmtree(inFiles[i]+'.xy',True)
            command="fits in=%s op=xyin out=%s.xy" % (inFiles[i],inFiles[i])
            err, dummy = commands.getstatusoutput(command)
            if err:
                raise Exception()

            # Regrid the remaining files
            if i>0:
                if os.path.exists(inFiles[i]+'.rg'):
                    shutil.rmtree(inFiles[i]+'.rg',True)
                command="regrid in=%s.xy out=%s.rg " % (inFiles[i],inFiles[i])\
                         + "tin=%s axes=1,2 options=offset" % (stitchList[0])
                err,dummy = commands.getstatusoutput(command)
                if err:                    
                    raise Exception()
                stitchList.append(inFiles[i]+'.rg')                
            else:
                stitchList.append(inFiles[i]+'.xy')
    except Exception:
        if os.path.exists(tmpDir): shutil.rmtree(tmpDir,True)
        os.chdir(startDir)
        return 1

    # Stitch the files together
    try:
        if os.path.exists(outFile+'.xy'):
            shutil.rmtree(outFile+'.xy',True)
        stitchStr = (",".join(stitchList))
        command="imcomb in=%s out=%s.xy options=relax" % \
                 (stitchStr,outFile)
        err,dummy = commands.getstatusoutput(command)
        if err: raise Exception()
    except Exception:
        if os.path.exists(tmpDir): shutil.rmtree(tmpDir,True)
        os.chdir(startDir)
        return 2

    # Convert the output back to a FITS file
    try:
        command="fits in=%s.xy op=xyout out=%s" % (outFile,outFile)
        err,dummy = commands.getstatusoutput(command)
        if err:
            raise Exception()
        if os.path.exists(outFile):
            shutil.move(outFile, startDir + '/' + outDir + '/' + outFile)
    except Exception:
        if os.path.exists(tmpDir):
            shutil.rmtree(tmpDir,True)
        os.chdir(startDir)
        return 3
    
    # Remove the temporary directory
    if os.path.exists(tmpDir):
        shutil.rmtree(tmpDir,True)
    os.chdir(startDir)
    return 0
    

#-----------------------------------------------------------------------------#
def miriad_gal_j2000(infilename, outfilename, delMirFile=True):
    """Use the MIRIAD regrid task to convert an image from GAL to J2000"""
    
    # Strip the path from the in and outfile names
    inPath,inFile   = os.path.split(infilename)
    outPath,outFile = os.path.split(outfilename)

    # Make a local copy to operate on
    if not os.path.exists('TMP'): os.mkdir('TMP')
    if os.path.exists('TMP/'+inFile): os.remove('TMP/'+inFile)
    shutil.copy(infilename,'TMP/'+inFile)

    try:
        command="fits in=%s op=xyin out=%s.xy" % ('TMP/'+inFile,'TMP/'+inFile)
        dummy=os.system(command)
        command="regrid in=%s.xy options=galeqsw,equisw out=%s.xy" % \
            ('TMP/'+inFile,'TMP/'+outFile)
        dummy=os.system(command)
        command="fits in=%s.xy op=xyout out=%s" % ('TMP/'+outFile,
                                                   'TMP/'+outFile)
        dummy=os.system(command)
    except Exception:
        print "Regridding failed!"
        return 1

    if os.path.exists('TMP/'+outFile): shutil.move('TMP/'+outFile,outfilename)

    # Clean up
    if not delMirFile:
        if os.path.exists(outPath+'/'+outFile+'.xy'):
            shutil.rmtree(outPath+'/'+outFile+'.xy',True)        
        shutil.move('TMP/'+outFile+'.xy',outPath+'/'+outFile+'.xy')
    
    if os.path.exists('TMP'): shutil.rmtree('TMP',True)


#-----------------------------------------------------------------------------#
def miriad_gal_j20001(inFitsPathFile,outFitsPathFile):
    """Use the MIRIAD regrid task to convert an image from GAL to J2000"""
    
    # Separate the paths and filenames
    inDir,inFile = os.path.split(inFitsPathFile)    
    outDir,outFile = os.path.split(outFitsPathFile)
    if inDir=='': 
        inDir = '.'
    if outDir=='': 
        outDir = '.'
    if inFile==outFile:
        outRoot, outExt = os.path.splitext(outFile)
        outFileTmp = outRoot+'.TMP'+outExt
    else:
        outFileTmp = outFile        
    
    # Create a TMP directory to work in
    randStr = "REGRD"+"".join(random.sample(string.letters+string.digits, 4))
    tmpDir = outDir + '/' + randStr
    if os.path.exists(tmpDir): shutil.rmtree(tmpDir,True)
    
    # Change to the working directory
    startDir = os.getcwd()
    os.mkdir(tmpDir)
    os.chdir(tmpDir)

#    try:
    if True:   
        # Make a local copy to operate on
        if os.path.exists(inFile): os.remove(inFile)
        shutil.copy(inFitsPathFile,inFile)
        
        # Convert to miriad format
        if os.path.exists(inFile+'.xy'): shutil.rmtree(inFile+'.xy',True)
        command="fits in=%s op=xyin out=%s.xy" % (inFile,inFile)
        err,dummy = commands.getstatusoutput(command)
        if err: raise Exception()

        # Regrid to J2000 Equatorial coordinate frame
        if os.path.exists(inFile+'.rg'): shutil.rmtree(inFile+'.rg',True)
        command="regrid in=%s.xy out=%s.rg options=galeqsw,equisw"  % \
                 (inFile,inFile)
        err,dummy = commands.getstatusoutput(command)
        if err: raise Exception()

#    except Exception:
#        print "MEEP"
#        if os.path.exists(tmpDir): shutil.rmtree(tmpDir,True)
#        os.chdir(startDir)
#        return 1

    # Convert the output back to a FITS file
    try:
        command="fits in=%s.rg op=xyout out=%s" % (inFile,outFileTmp)
        err,dummy = commands.getstatusoutput(command)
        if err: raise Exception()
        if os.path.exists(outFileTmp):
            if os.path.exists(outDir + '/' + outFile):
                os.remove(outDir + '/' + outFile)
            shutil.move(outFileTmp,outDir + '/' + outFile)
        else:
            raise Exception()
    except Exception:
        if os.path.exists(tmpDir): shutil.rmtree(tmpDir,True)
        os.chdir(startDir)
        return dummy
    
    # Remove the temporary directory
    if os.path.exists(tmpDir): shutil.rmtree(tmpDir,True)
    os.chdir(startDir)
    return 0
        

#-----------------------------------------------------------------------------#
def miriad_regrid_tin(inFitsFile,outFitsFile,tinFitsFile):
    """Use MIRIAD to regrid a FITS file using a template"""

    # Strip the path from the in and outfile names
    inPath,inFile   = os.path.split(inFitsFile)
    tinPath,tinFile  = os.path.split(tinFitsFile)
    outPath,outFile = os.path.split(outFitsFile)

    # Make a local copies to operate on
    if not os.path.exists('TMP'): os.mkdir('TMP')

    retVal=0
    try:
        # Convert to miriad format
        if os.path.exists('TMP/'+inFile+'.m'):
            shutil.rmtree('TMP/'+inFile+'.m',True)
        command="fits in=%s op=xyin out=%s" % (inFitsFile,'TMP/'+inFile+'.m')
        dummy=os.system(command)
        if os.path.exists('TMP/'+tinFile+'.m'):
            shutil.rmtree('TMP/'+tinFile+'.m',True)
        command="fits in=%s op=xyin out=%s" % (tinFitsFile,'TMP/'+tinFile+'.m')
        dummy=os.system(command)

        # Regrid using the template
        if os.path.exists('TMP/'+outFile+'.m'):
            shutil.rmtree('TMP/'+outFile+'.m',True)
        command="regrid tin=%s in=%s out=%s" % \
                 ('TMP/'+tinFile+'.m','TMP/'+inFile+'.m','TMP/'+outFile+'.m')
        dummy=os.system(command)

        # Convert the output to FITS format
        if os.path.exists(outFitsFile): os.remove(outFitsFile)
        command="fits in=%s.m op=xyout out=%s" % ('TMP/'+outFile,outFitsFile)
        dummy=os.system(command)
        
    except Exception:
        print "Regridding failed!"
        retVal = 1

    # Clean up
    if os.path.exists('TMP'): shutil.rmtree('TMP',True)
    if retVal: raise Exception


#-----------------------------------------------------------------------------#
def miriad_regrid_tin1(inFitsFile,outFitsFile,tinFitsFile):
    """Use MIRIAD to regrid a FITS file using a template (2012 version)"""
    
    # Strip the path from the in and outfile names
    inPath,inFile   = os.path.split(inFitsFile)
    tinPath,tinFile  = os.path.split(tinFitsFile)
    outPath,outFile = os.path.split(outFitsFile)

    # Make a local copies to operate on
    if not os.path.exists('TMP'): 
        os.mkdir('TMP')

    retVal=0
    try:
        # Convert to miriad format
        if os.path.exists('TMP/'+inFile+'.m'):
            shutil.rmtree('TMP/'+inFile+'.m',True)
        command="fits in=%s op=xyin out=%s" % (inFitsFile,'TMP/'+inFile+'.m')
        dummy=os.system(command)
        if os.path.exists('TMP/'+tinFile+'.m'):
            shutil.rmtree('TMP/'+tinFile+'.m',True)
        command="fits in=%s op=xyin out=%s" % (tinFitsFile,'TMP/'+tinFile+'.m')
        dummy=os.system(command)

        # Regrid using the template
        if os.path.exists('TMP/'+outFile+'.m'):
            shutil.rmtree('TMP/'+outFile+'.m',True)
        command="regrid tin=%s in=%s out=%s " % \
                 ('TMP/'+tinFile+'.m','TMP/'+inFile+'.m','TMP/'+outFile+'.m')
        dummy=os.system(command)

        # Convert the output to FITS format
        if os.path.exists(outFitsFile): os.remove(outFitsFile)
        command="fits in=%s.m op=xyout out=%s" % ('TMP/'+outFile,outFitsFile)
        dummy=os.system(command)
        
    except Exception:
        print "Regridding failed!"
        retVal = 1

    # Clean up
    if os.path.exists('TMP'): shutil.rmtree('TMP',True)
    if retVal: 
        raise Exception
    

#-----------------------------------------------------------------------------#
def fits_to_miriad(inFitsFile, outMir=None, clobber=True, verb=True):
    """Convert a FITS file to a MIRIAD dataset"""
    
    if not os.path.exists(inFitsFile):
        return 1

    inPath,inFile = os.path.split(inFitsFile)
    inFitsRoot,inExt = os.path.splitext(inFile)
    
    if outMir is None:
        outMir = inFitsRoot + '.mir'
    if os.path.exists(outMir):
        if clobber == True:
            shutil.rmtree(outMir, True)
        else:
            return 1
    
    cmd = "fits op=xyin in=%s out=%s" % (inFitsFile,outMir)
    status,msg = run_command(cmd,verb=verb)

    return status


#-----------------------------------------------------------------------------#
def miriad_to_fits(inMir, outFitsFile=None, clobber=True, verb=True):
    """Convert MIRIAD dataset to a FITS file"""
    
    if not os.path.exists(inMir):
        return 1

    inPath,inMirData = os.path.split(inMir)
    inMirRoot,inExt = os.path.splitext(inMirData)
    if outFitsFile is None:
        outFitsFile = inMirRoot + '.fits'
    if os.path.exists(outFitsFile):
        if clobber == True:
            os.remove(outFitsFile)
        else:
            return 1
    
    cmd = "fits op=xyout in=%s out=%s" % (inMir,outFitsFile)
    status,msg = run_command(cmd,verb=verb)

    return status
