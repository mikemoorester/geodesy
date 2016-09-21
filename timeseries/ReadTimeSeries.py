#!/usr/bin/env python
from __future__ import division, print_function, absolute_import

import re
import numpy as np
import datetime as dt
import gpsTime as gt

#     0      1         2        3               4         5         6         7        8    9       10       11    12       13       14   15
# YYYYMMDD HHMNSC    DecYr     MJD              N         E         U        dN       +-    F       dE       +-    F        dU       +-   F
#  20100103 115900  2010.0068  55199.4993     -35.2     -32.9      11.8       0.8      1.6   0      -0.8      1.3   0      -3.5      4.7   0
def _parseTSFITepoch(line):
    data = line.split()
    epoch = {}
    epoch['YYYYMMDD'] = data[0]
    epoch['HHMMSS']   = data[1]
    epoch['DecYr']    = data[2]
    epoch['mjd']      = data[3]
    epoch['North']    = data[4]
    epoch['East']     = data[5]
    epoch['H']        = data[6]
    epoch['dN']       = data[7]
    epoch['sN']       = data[8]
    epoch['dE']       = data[10]
    epoch['sE']       = data[11]
    epoch['dU']       = data[13]
    epoch['sU']       = data[14]
    epoch['Flag']     = data[15]

    return epoch

def readTSFIT(tsfit_file):
    """
    Read the tsfit residual file:

    tsfit = readTSFIT(filename)

    """
    read_header = 0
    data = []
    epochs = []
    tsfit = {}
    tsfit['num_epochs'] = 0

    # Comment RGX 
    asterix_RGX = re.compile('^\#')
    #
    # Header regexs to get some meta data from the pbo time series
    #
    # 4-character ID: SEY1
    ssss_RGX = re.compile('4-character ID:')
    #
    # Station name  : Republic_of_Seyc
    name_RGX = re.compile('Station name  :')
    #
    # First Epoch   : 20001109 115900
    firstEpoch_RGX = re.compile('First Epoch   :')
    #
    # Last Epoch    : 20120902 115900
    lastEpoch_RGX = re.compile('Last Epoch    :')
    #
    # Release Data  :   130124 213047
    releaseData_RGX = re.compile('Release Data  :')
    #
    # XYZ Reference position :   3602870.54224  5238174.48216  -516275.39017 (Unknown)
    xyzReferencePos_RGX = re.compile('XYZ Reference position :')
    #
    # NEU Reference position :    -4.6737188338   55.4794047615  537.19792 (Unknown/WGS84)'
    neuReferencePos_RGX = re.compile('NEU Reference position :')
    #
    with open(tsfit_file,'r') as tsfitFH:
        for line in tsfitFH:
            # Haven;t finished reading the header
            if read_header == 0:
                if asterix_RGX.search(line):
                    #print('Found a comment:',line)
                    read_header = 1
                elif ssss_RGX.search(line):
                    #print('Found ssss',line)
                    tsfit['site'] = line[16:20]
                    #print('SITE:',pbo['site'],':SITE')
                elif name_RGX.search(line):
                    print('Found name')
                elif firstEpoch_RGX.search(line):
                    tsfit['firstEpoch'] = line[16:31]
                    print('Found first epoch',tsfit['firstEpoch'])
                elif lastEpoch_RGX.search(line):
                    tsfit['lastEpoch'] = line[16:31]
                    print('Found last epoch',tsfit['lastEpoch'])
                elif releaseData_RGX.search(line):
                    print('Found release data')
                elif xyzReferencePos_RGX.search(line):
                    print('Found xyz reference position')
                elif neuReferencePos_RGX.search(line):
                    print('Found neu reference position')
            else:
                # remove the new line character
                line = line.rstrip()
                #print('About to look at:',line)
                epoch = _parseTSFITepoch(line)
                epochs.append(epoch)
                tsfit['num_epochs'] = tsfit['num_epochs'] + 1
                #print("ANSWER:",YYYYMMDD, HHMMSS, X, Y, Z)
                #print("Epoch:",epoch)
                #print()

    tsfit['epochs'] = epochs
    return tsfit

def _parsePBOEpoch(line):
    #YYYYMMDD, HHMMSS, JJJJ, X, Y, Z, Sx, Sy, Sz, Rxy, Rxz, Ryz, Nlat, Elong, H, dN, dE, dU, sN, sE, sU, Rne, Rnu, Reu, Soln = line.split('\s+')
    #data = line.split('\s+')
    data = line.split()
    epoch = {}
    epoch['YYYYMMDD'] = data[0]
    epoch['HHMMSS']   = data[1]
    epoch['JJJJ']     = data[2]
    epoch['X']        = data[3]
    epoch['Y']        = data[4]
    epoch['Z']        = data[5]
    epoch['Sx']       = data[6]
    epoch['Sy']       = data[7]
    epoch['Sz']       = data[8]
    epoch['Rxy']      = data[9]
    epoch['Rxz']      = data[10]
    epoch['Ryz']      = data[11]
    epoch['Nlat']     = data[12]
    epoch['Elong']    = data[13]
    epoch['H']        = data[14]
    epoch['dN']       = data[15]
    epoch['dE']       = data[16]
    epoch['dU']       = data[17]
    epoch['sN']       = data[18]
    epoch['sE']       = data[19]
    epoch['sU']       = data[20]
    epoch['Rne']      = data[21]
    epoch['Rnu']      = data[22]
    epoch['Reu']      = data[23]
    epoch['Soln']     = data[24]

    return epoch

def readPBOTS(pbot_file):
    """
    Read the standard PBO time series data format:

    TsN, TsE, TsU, Stat = readPBOTS(filename)

    """
    read_header = 0
    data = []
    epochs = []
    pbo = {}
    pbo['num_epochs'] = 0

    # Comment RGX 
    asterix_RGX = re.compile('^\*')
    #
    # Header regexs to get some meta data from the pbo time series
    #
    # Format Version: 1.0.2
    pboVersion_RGX = re.compile('Format Version:')
    #
    # PBO Station Position Time Series. Reference Frame : Unknown
    referenceFrame_RGX = re.compile('PBO Station Position Time Series. Reference Frame :')
    #
    # 4-character ID: SEY1
    ssss_RGX = re.compile('4-character ID:')
    #
    # Station name  : Republic_of_Seyc
    name_RGX = re.compile('Station name  :')
    #
    # First Epoch   : 20001109 115900
    firstEpoch_RGX = re.compile('First Epoch   :')
    #
    # Last Epoch    : 20120902 115900
    lastEpoch_RGX = re.compile('Last Epoch    :')
    #
    # Release Data  :   130124 213047
    releaseData_RGX = re.compile('Release Data  :')
    #
    # XYZ Reference position :   3602870.54224  5238174.48216  -516275.39017 (Unknown)
    xyzReferencePos_RGX = re.compile('XYZ Reference position :')
    #
    # NEU Reference position :    -4.6737188338   55.4794047615  537.19792 (Unknown/WGS84)'
    neuReferencePos_RGX = re.compile('NEU Reference position :')
    #
    with open(pbot_file,'r') as pboFH:
        for line in pboFH:
            # Haven;t finished reading the header
            if read_header == 0:
                if asterix_RGX.search(line):
                    #print('Found a comment:',line)
                    read_header = 1
                elif pboVersion_RGX.search(line):
                    print('Found pbo Version')
                elif referenceFrame_RGX.search(line):
                    print('Found reference frame')
                elif ssss_RGX.search(line):
                    #print('Found ssss',line)
                    pbo['site'] = line[16:20]
                    #print('SITE:',pbo['site'],':SITE')
                elif name_RGX.search(line):
                    print('Found name')
                elif firstEpoch_RGX.search(line):
                    print('Found first epoch')
                elif lastEpoch_RGX.search(line):
                    print('Found last epoch')
                elif releaseData_RGX.search(line):
                    print('Found release data')
                elif xyzReferencePos_RGX.search(line):
                    print('Found xyz reference position')
                elif neuReferencePos_RGX.search(line):
                    print('Found neu reference position')
            else:
                # remove the new line character
                line = line.rstrip()
                #print('About to look at:',line)
                #YYYYMMDD, HHMMSS, JJJJ, X, Y, Z, Sx, Sy, Sz, Rxy, Rxz, Ryz, Nlat, Elong, H, dN, dE, dU, sN, sE, sU, Rne, Rnu, Reu, Soln = line.split('\s+')
                #data = line.split('\s+')
                #data = line.split()
                #YYYYMMDD = data[0]
                #HHMMSS = data[1]
                #X = data[2]
                #Y = data[3]
                #Z = data[4]
                epoch = _parsePBOEpoch(line)
                epochs.append(epoch)
                pbo['num_epochs'] = pbo['num_epochs'] + 1
                #print("ANSWER:",YYYYMMDD, HHMMSS, X, Y, Z)
                #print("Epoch:",epoch)
                #print()

    pbo['epochs'] = epochs
    return pbo

#% Conert jday back to DecYrs.  First convert MJD to MATLAB's serial date
#% number
#jday = jday + 678942;
#% Now covert.  Needs to be done one at a time due to leap years
#for i = 1:length(jday);
#    bdfull = datevecl(jday(i));
#    decyr = DateToYr(bdfull,6);
#    jday(i) = decyr;
#end
#% Now creat the DataN, DataE, DataU arrays
#DataN = [jday'; n'*1000; err_n'*1000];
#DataE = [jday'; e'*1000; err_e'*1000];
#DataU = [jday'; u'*1000; err_u'*1000];

def get_height(pbo):
    '''
    Given a PBO or TSFIT data structure return the height time series as a numpy array
    '''
    height = np.zeros(pbo['num_epochs'])
    ctr = 0
    for epoch in pbo['epochs']:
        height[ctr] = epoch['H']
        ctr = ctr + 1

    return height

def get_array(tso,parameter):
    '''
    Given a PBO or TSFIT data structure return the requested time series as a numpy array

    tso = timeseries data object, either tsfit or pbo

    parameter = parameter to return as an array
            'H' => hieght
            'dU' => mean removed heigth time series
            'mp_ts' = > convert the time series into matplotlib values

    '''
    npArray = np.zeros(tso['num_epochs'])
    ctr = 0
    for epoch in tso['epochs']:
        if parameter == 'mp_ts':
            year  = epoch['YYYYMMDD'][0:4]
            month = epoch['YYYYMMDD'][4:6] 
            dom   = epoch['YYYYMMDD'][6:8]
            #print("YYYY",epoch['YYYYMMDD'][0:4])
            npArray[ctr] = gt.ymdhms2mdt(year,month,dom,0,0,0.0)
        else:     
            npArray[ctr] = float( epoch[parameter] )
        ctr = ctr + 1


    # check for any zeros left in the mp_ts case, set to NaN 
    #if parameter == 'mp_ts' :
    #    #criterion = ( npArray[:] > -0.0001 ) & ( npArray[:] < 0.00001)
    #    criterion = npArray[:] < 0.00001
    #    ind = np.array(np.where(criterion))
    #    #print("number of NaN",np.size(ind))
    #    npArray[ind] = np.nan

    return npArray

def startEnd2dto(tso):
    '''

    Return the start and end epoch as a date time object

    '''

    year  = tso['firstEpoch'][0:4]
    month = tso['firstEpoch'][4:6] 
    dom   = tso['firstEpoch'][6:8]
    hh = 0 
    mm = 0
    sec = 0
    ms = 0
    sdto = dt.datetime(int(year),int(month),int(dom),int(hh),int(mm),int(sec),int(ms))

    year  = tso['lastEpoch'][0:4]
    month = tso['lastEpoch'][4:6] 
    dom   = tso['lastEpoch'][6:8]
    edto = dt.datetime(int(year),int(month),int(dom),int(hh),int(mm),int(sec),int(ms))
    return sdto,edto
