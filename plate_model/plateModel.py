#!/usr/bin/env python
from __future__ import print_function

import numpy as np
import sys
import json

from scipy.stats import t

#import matplotlib.pyplot as plt
#from mpl_toolkits.basemap import Basemap
#from matplotlib.patches import Polygon

import argparse
import pickle

#from shapely.geometry.polygon import LinearRing
#from shapely.geometry import Point
#from shapely.geometry import Polygon as Pg

import sinex
import geodetic as gt

#==============================================================================
 

#def plotPointsOnPlate(lon,lat,pmap):
#    """
#       plotPointsOnPlate(plate,lon,lat)
#
#       plot a polygon onto a basemap
#
#    """
#    print("lat:",lat)
#    x, y = pmap(np.array(lon), np.array(lat))
#    print("y",y)
#    pmap.plot(x,y,'bo',markersize=5)
#
#    return pmap

#def plotPolygon(lons,lats,pmap):
#    """
#       plotPolygon(lons,lats)
#
#       plot a polygon onto a basemap
#
#    """
#    # compute native map projection coordinates of lat/lon grid.
#    x, y = pmap(np.array(lons), np.array(lats))
#    xy = zip(x,y)
#    poly = Polygon(xy, facecolor='gray',alpha=0.1)
#    plt.gca().add_patch(poly)
#    return pmap
#    #plt.show()

#def plotShapeFile(pmap,shapefile):
#    """
#       plotPolygon(lons,lats)
#
#       plot a polygon onto a basemap
#
#    """
#    # compute native map projection coordinates of lat/lon grid.
#    pmap.readshapefile(shapefile,'anything',drawbounds=True, linewidth=2, color='b')
#
## an attempt at a decorator, neeed to check this 
#def checkdeg(func):
#    """
#    A simple bound check:
#
#    Check the degrees is less than 360. and greate than 0
#    """
#    def wrapper(*args):
#        if args[0] > 0.0 and args[0] < 360.0:
#            func(args)
#        else:
#            print("Degrees is out of bounds",args)
#
#    return wrapper

#@checkdeg
#==============================================================================
# def decdeg2dms(dd):
#     """
#     deg,mm,ss = decdeg2dms(dd)
# 
#     Convert decimal degree to degrees, minutes, seconds
# 
#      Input:
#              dd    decimal degrees (float)
# 
#      Output:
#              deg   degrees
#              mm    minutes
#              ss    seconds
#     """
#     mm,ss = divmod(dd*3600,60)
#     deg,mm = divmod(mm,60)
#    return deg,mm,ss
#==============================================================================

def setVelocityITRF2008AustraliaPlate(sites):
    """

    sites = setVelocity(sites)

    Set up a dictionary of values for sites used to determine the
    plate motion of the Australian plate from ITRF2008 derived velocities.

    Values were obtained from Alatmimi et al 2012 ITRF2008 PLATE Motion Model
    Journal of Geophysical research, vol 117, B07402, doi:10.1029/2011JB008930,2012
    Table A1.

     Input:  sites  a dictionary of sites, can be empty sites = {} or should have the same format as 
                    output dictionary

     Output: sites  a dictionary which will include the australian plate sites defined by ITRF2008
                    with the following structure:
                    sites['SSSS']                 SSSS is site id
                         ['SSSS']['Ve']           Ve velocity in east direction (mm/a)
                         ['SSSS']['Vn']           Vn velocity in north direction (mm/a)
                         ['SSSS']['sVe']          sVe sigma Ve (mm/a)
                         ['SSSS']['sVn']          sVn sigma Vn (mm/a)
                         ['SSSS']['long']         longitude in degrees
                         ['SSSS']['lat']          latitude in degrees


    """
    # Define the data structure if needed
    stations = ('YAR2','NNOR','KARR','DARW','CEDU','ALIC','ADE1','BUR1','TOW2','HOB2',
                'PARK','TIDB','STR1','SYDN','SUNM','KOUC','NOUM')

    for station in stations:
        if station not in sites: sites[station] = {}

    # set the values as published in the paper
    sites['YAR2']['Ve']  = 38.96
    sites['YAR2']['Vn']  = 57.79 
    sites['YAR2']['sVn'] = 0.06
    sites['YAR2']['sVe'] = 0.05
    sites['YAR2']['lon'] = 115.347
    sites['YAR2']['lat'] = -28.884

    sites['NNOR']['Ve']  = 38.66
    sites['NNOR']['Vn']  = 58.08 
    sites['NNOR']['sVn'] = 0.07
    sites['NNOR']['sVe'] = 0.08
    sites['NNOR']['lon'] = 116.193 
    sites['NNOR']['lat'] = -30.879

    sites['KARR']['Ve']  = 38.89
    sites['KARR']['Vn']  = 58.36  
    sites['KARR']['sVn'] = 0.06
    sites['KARR']['sVe'] = 0.05
    sites['KARR']['lon'] = 117.097 
    sites['KARR']['lat'] = -20.853

    sites['DARW']['Ve']  = 35.62 
    sites['DARW']['Vn']  = 59.38 
    sites['DARW']['sVn'] = 0.08
    sites['DARW']['sVe'] = 0.05
    sites['DARW']['lon'] = 131.133 
    sites['DARW']['lat'] = -12.761

    sites['CEDU']['Ve']  = 29.11 
    sites['CEDU']['Vn']  = 58.70 
    sites['CEDU']['sVn'] = 0.06
    sites['CEDU']['sVe'] = 0.06
    sites['CEDU']['lon'] = 133.810 
    sites['CEDU']['lat'] = -31.694

    sites['ALIC']['Ve']  = 32.01
    sites['ALIC']['Vn']  = 59.01 
    sites['ALIC']['sVn'] = 0.06
    sites['ALIC']['sVe'] = 0.05
    sites['ALIC']['lon'] = 133.886 
    sites['ALIC']['lat'] = -23.529

    sites['ADE1']['Ve']  = 24.76
    sites['ADE1']['Vn']  = 58.46 
    sites['ADE1']['sVn'] = 0.06
    sites['ADE1']['sVe'] = 0.06
    sites['ADE1']['lon'] = 138.647 
    sites['ADE1']['lat'] = -34.549

    sites['BUR1']['Ve']  = 15.54 
    sites['BUR1']['Vn']  = 57.29 
    sites['BUR1']['sVn'] = 0.18
    sites['BUR1']['sVe'] = 0.22
    sites['BUR1']['lon'] = 145.915 
    sites['BUR1']['lat'] = -40.860

    sites['TOW2']['Ve']  = 28.84 
    sites['TOW2']['Vn']  = 55.87 
    sites['TOW2']['sVn'] = 0.06
    sites['TOW2']['sVe'] = 0.05
    sites['TOW2']['lon'] = 147.056 
    sites['TOW2']['lat'] = -19.150

    sites['HOB2']['Ve']  = 14.12 
    sites['HOB2']['Vn']  = 55.80 
    sites['HOB2']['sVn'] = 0.06
    sites['HOB2']['sVe'] = 0.07
    sites['HOB2']['lon'] = 147.439 
    sites['HOB2']['lat'] = -42.613

    sites['PARK']['Ve']  = 19.40 
    sites['PARK']['Vn']  = 55.05 
    sites['PARK']['sVn'] = 0.20
    sites['PARK']['sVe'] = 0.21
    sites['PARK']['lon'] = 148.265 
    sites['PARK']['lat'] = -32.823

    sites['TIDB']['Ve']  = 18.23
    sites['TIDB']['Vn']  = 55.49 
    sites['TIDB']['sVn'] = 0.06
    sites['TIDB']['sVe'] = 0.07
    sites['TIDB']['lon'] = 148.980 
    sites['TIDB']['lat'] = -35.218

    sites['STR1']['Ve']  = 18.40 
    sites['STR1']['Vn']  = 55.52 
    sites['STR1']['sVn'] = 0.06
    sites['STR1']['sVe'] = 0.07
    sites['STR1']['lon'] = 149.010 
    sites['STR1']['lat'] = -35.134

    sites['SYDN']['Ve']  = 17.99 
    sites['SYDN']['Vn']  = 54.60 
    sites['SYDN']['sVn'] = 0.12
    sites['SYDN']['sVe'] = 0.14
    sites['SYDN']['lon'] = 151.150 
    sites['SYDN']['lat'] = -33.603

    sites['SUNM']['Ve']  = 21.57 
    sites['SUNM']['Vn']  = 53.65 
    sites['SUNM']['sVn'] = 0.20
    sites['SUNM']['sVe'] = 0.22
    sites['SUNM']['lon'] = 153.035 
    sites['SUNM']['lat'] = -27.328

    sites['KOUC']['Ve']  = 24.15
    sites['KOUC']['Vn']  = 47.81 
    sites['KOUC']['sVn'] = 0.08
    sites['KOUC']['sVe'] = 0.10
    sites['KOUC']['lon'] = 164.287 
    sites['KOUC']['lat'] = -20.432

    sites['NOUM']['Ve']  = 20.52
    sites['NOUM']['Vn']  = 46.12 
    sites['NOUM']['sVn'] = 0.07
    sites['NOUM']['sVe'] = 0.09
    sites['NOUM']['lon'] = 166.410 
    sites['NOUM']['lat'] = -22.135

    return sites

def nuvel1APlateModel():
    """
    NNR-NUVEL-1:(absolute plate motion, no-net rotation)

    Argus D.F., R.G. Gordon, No-net-rotation model of current plate velocities incorporating plate motion model NUVEL-1. Geophys. Res. Lett. (18) 2039-2042, 1991.

    NNR-NUVEL-1A:(absolute plate motion, no-net rotation)

    DeMets. C., R. G. Gordon, D. F. Argus, and S. Stein, Effect of recent revisions to the geomagnetic reversal time scale on estimate of current plate motions, Geophys. Res. Lett., vol. 21, no. 20, 2191-2194, 1994.

    """

    nuvel1A = {}

    plate_names = ('AMUR','ANTA','ARAB','AUST','CARB','EURA','INDI','NAZC','NOAM','NUBI',
                   'PCFC','SOAM','SOMA','SUND')

    for name in plate_names:
        nuvel1A[name] = {}

    #0 AUST
    #np.radians( ( plateModel['omega_X'] / 1000. ) / 3600.)
    # units are radians/Myr
    #nuvel1A['AUST']['omega_X']  = 0.009349 / 1000000.
    #nuvel1A['AUST']['omega_sX'] = 100.
    #nuvel1A['AUST']['omega_Y']  = 0.000284 / 1000000.
    #nuvel1A['AUST']['omega_sY'] = 100.
    #nuvel1A['AUST']['omega_Z']  = 0.016252 / 1000000.
    #nuvel1A['AUST']['omega_sZ'] = 100.
    #nuvel1A['AUST']['omega']    = 1.0744 / 1000000.
    #nuvel1A['AUST']['std_omega']= 100.

    # Obtained from satellite geodesy page 528 Seber
    # units are radians/Myr -> expecting units to be radians/year
    nuvel1A['AUST']['omega_X']  = 0.007839 / 1000000.
    nuvel1A['AUST']['omega_sX'] = 100.
    nuvel1A['AUST']['omega_Y']  = 0.005124 / 1000000.
    nuvel1A['AUST']['omega_sY'] = 100.
    nuvel1A['AUST']['omega_Z']  = 0.006282 / 1000000.
    nuvel1A['AUST']['omega_sZ'] = 100.
    nuvel1A['AUST']['omega']    = 0.00# / 1000000.
    nuvel1A['AUST']['std_omega']= 100.

    return nuvel1A

def plateModelITRF2008():
    """

    plateModel, Trate = plateModelITRF2008()

    Set up a dictionary of values for absolute plate rotation poles from 
    ITRF2008 derived velocities.

    Values were obtained from Alatmimi et al 2012 ITRF2008 PLATE Motion Model
    Journal of Geophysical research, vol 117, B07402, doi:10.1029/2011JB008930,2012
    Table 3.

     Input:  
     
     Output: 
             plateModels      a dictionary of plates, with the same format as:

                              plate[name][omega_X]        (mas/a)
                                         [omega_sX]       (mas/a)
                                         [omega_Y]        (mas/a)
                                         [omega_sY]       (mas/a)
                                         [omega_Z]        (mas/a)
                                         [omega_sZ]       (mas/a)
                                         [omega]          (deg/Ma)
                                         [std_omega]      (deg/Ma)

             Trate            a dictionary of the translation rate components:
     
                              Trate[Tx]     (mm/a)
                              Trate[sTx]    (mm/a)
                              Trate[Ty]     (mm/a)
                              Trate[sTy]    (mm/a)
                              Trate[Tz]     (mm/a)
                              Trate[sTz]    (mm/a)

    """

    Trate = {}
    Trate['Tx']  = 0.41
    Trate['sTx'] = 0.27
    Trate['Ty']  = 0.22
    Trate['sTy'] = 0.32
    Trate['Tz']  = 0.41
    Trate['sTz'] = 0.30

    plateModels = {}

    plate_names = ('AMUR','ANTA','ARAB','AUST','CARB','EURA','INDI','NAZC','NOAM','NUBI',
                   'PCFC','SOAM','SOMA','SUND')

    for name in plate_names:
        plateModels[name] = {}

    # AMUR
    plateModels['AMUR']['omega_X']  = -0.190
    plateModels['AMUR']['omega_sX'] = 0.040
    plateModels['AMUR']['omega_Y']  = -0.442
    plateModels['AMUR']['omega_sY'] = 0.051
    plateModels['AMUR']['omega_Z']  = 0.915
    plateModels['AMUR']['omega_sZ'] = 0.049
    plateModels['AMUR']['omega']    = 0.287
    plateModels['AMUR']['std_omega']= 0.008
    # ANTA
    plateModels['ANTA']['omega_X']  = -0.252
    plateModels['ANTA']['omega_sX'] = 0.008
    plateModels['ANTA']['omega_Y']  = -0.302
    plateModels['ANTA']['omega_sY'] = 0.006
    plateModels['ANTA']['omega_Z']  = 0.643
    plateModels['ANTA']['omega_sZ'] = 0.009
    plateModels['ANTA']['omega']    = 0.209
    plateModels['ANTA']['std_omega']= 0.003
    # ARAB
    plateModels['ARAB']['omega_X']  = 1.202
    plateModels['ARAB']['omega_sX'] = 0.082
    plateModels['ARAB']['omega_Y']  = -0.054
    plateModels['ARAB']['omega_sY'] = 0.100
    plateModels['ARAB']['omega_Z']  = 1.485
    plateModels['ARAB']['omega_sZ'] = 0.063
    plateModels['ARAB']['omega']    = 0.531
    plateModels['ARAB']['std_omega']= 0.027
    # AUST
    plateModels['AUST']['omega_X']  = 1.504
    plateModels['AUST']['omega_sX'] = 0.007
    plateModels['AUST']['omega_Y']  = 1.172
    plateModels['AUST']['omega_sY'] = 0.007
    plateModels['AUST']['omega_Z']  = 1.228
    plateModels['AUST']['omega_sZ'] = 0.007
    plateModels['AUST']['omega']    = 0.630
    plateModels['AUST']['std_omega']= 0.002
    # CARB
    plateModels['CARB']['omega_X']  = 0.049
    plateModels['CARB']['omega_sX'] = 0.201
    plateModels['CARB']['omega_Y']  = -1.088
    plateModels['CARB']['omega_sY'] = 0.417
    plateModels['CARB']['omega_Z']  = 0.664
    plateModels['CARB']['omega_sZ'] = 0.146
    plateModels['CARB']['omega']    = 0.354
    plateModels['CARB']['std_omega']= 0.122
    # EURA
    plateModels['EURA']['omega_X']  = -0.083
    plateModels['EURA']['omega_sX'] = 0.008
    plateModels['EURA']['omega_Y']  = -0.534
    plateModels['EURA']['omega_sY'] = 0.007
    plateModels['EURA']['omega_Z']  = 0.750
    plateModels['EURA']['omega_sZ'] = 0.008
    plateModels['EURA']['omega']    = 0.257
    plateModels['EURA']['std_omega']= 0.002
    # INDI
    plateModels['INDI']['omega_X']  = 1.232
    plateModels['INDI']['omega_sX'] = 0.031
    plateModels['INDI']['omega_Y']  = 0.303
    plateModels['INDI']['omega_sY'] = 0.128
    plateModels['INDI']['omega_Z']  = 1.540
    plateModels['INDI']['omega_sZ'] = 0.030
    plateModels['INDI']['omega']    = 0.554
    plateModels['INDI']['std_omega']= 0.017
    # NAZC
    plateModels['NAZC']['omega_X']  = -0.330
    plateModels['NAZC']['omega_sX'] = 0.011
    plateModels['NAZC']['omega_Y']  = -1.551
    plateModels['NAZC']['omega_sY'] = 0.029
    plateModels['NAZC']['omega_Z']  = 1.625
    plateModels['NAZC']['omega_sZ'] = 0.013
    plateModels['NAZC']['omega']    = 0.631
    plateModels['NAZC']['std_omega']= 0.005
    # NOAM
    plateModels['NOAM']['omega_X']  = 0.035
    plateModels['NOAM']['omega_sX'] = 0.008
    plateModels['NOAM']['omega_Y']  = -0.662
    plateModels['NOAM']['omega_sY'] = 0.009
    plateModels['NOAM']['omega_Z']  = -0.100
    plateModels['NOAM']['omega_sZ'] = 0.008
    plateModels['NOAM']['omega']    = 0.186
    plateModels['NOAM']['std_omega']= 0.002
    # NUBI
    plateModels['NUBI']['omega_X']  = 0.095
    plateModels['NUBI']['omega_sX'] = 0.009
    plateModels['NUBI']['omega_Y']  = -0.598
    plateModels['NUBI']['omega_sY'] = 0.007
    plateModels['NUBI']['omega_Z']  = 0.723
    plateModels['NUBI']['omega_sZ'] = 0.009
    plateModels['NUBI']['omega']    = 0.262
    plateModels['NUBI']['std_omega']= 0.003
    # PCFC
    plateModels['PCFC']['omega_X']  = -0.411
    plateModels['PCFC']['omega_sX'] = 0.007
    plateModels['PCFC']['omega_Y']  = 1.036
    plateModels['PCFC']['omega_sY'] = 0.007
    plateModels['PCFC']['omega_Z']  = -2.166
    plateModels['PCFC']['omega_sZ'] = 0.009
    plateModels['PCFC']['omega']    = 0.677
    plateModels['PCFC']['std_omega']= 0.002
    # SOAM
    plateModels['SOAM']['omega_X']  = -0.243
    plateModels['SOAM']['omega_sX'] = 0.009
    plateModels['SOAM']['omega_Y']  = -0.311
    plateModels['SOAM']['omega_sY'] = 0.010
    plateModels['SOAM']['omega_Z']  = -0.154
    plateModels['SOAM']['omega_sZ'] = 0.009
    plateModels['SOAM']['omega']    = 0.118
    plateModels['SOAM']['std_omega']= 0.002
    # SOMA
    plateModels['SOMA']['omega_X']  = -0.080
    plateModels['SOMA']['omega_sX'] = 0.028
    plateModels['SOMA']['omega_Y']  = -0.745
    plateModels['SOMA']['omega_sY'] = 0.030
    plateModels['SOMA']['omega_Z']  = 0.897
    plateModels['SOMA']['omega_sZ'] = 0.012
    plateModels['SOMA']['omega']    = 0.325
    plateModels['SOMA']['std_omega']= 0.007
    # SUND
    plateModels['SUND']['omega_X']  = 0.047
    plateModels['SUND']['omega_sX'] = 0.381
    plateModels['SUND']['omega_Y']  = -1.000
    plateModels['SUND']['omega_sY'] = 1.570
    plateModels['SUND']['omega_Z']  = 0.975
    plateModels['SUND']['omega_sZ'] = 0.045
    plateModels['SUND']['omega']    = 0.388
    plateModels['SUND']['std_omega']= 0.308

    return plateModels, Trate

def determineEulerPole(snx,Q,args): 
    """

    omega = determineEulerPole(sites)

    Determine the euler pole rates for a set of stations.

     Input: 
            sites   dictionary stations to calculate a plate motion model from

     Output:
            omega   vector of absolute plane rotation poles in cartesian form

                    omega[0] -> omega_x
                    omega[1] -> omega_y
                    omega[2] -> omega_z

            omega_V covariance matrix 

    """
    #============================================================================
    # A is the design matrix, partial derivatives with respect to the parameters 
    #       (number of estimates x number of observations)
    # n is the number of stations with horizontal velocities on the plate
    # Q is the a priori VCV matrix of observations
    # P inverse of the covariance matrix of observations
    # b is the misclose vector (l - e : Observed minus computed)
    # l is observation vector
    # e are the adjusted observations
    #
    # Solution of parameters: 
    #                         dX = (A^T P A)^-1 A^T P b
    # 
    #============================================================================
    m = np.shape(Q)[0]
    n = int(np.float(m) / 2.)
    #print("m,n",m,n)
    newQ = np.zeros((n,n))

    # Set up the matrices based on the number of stations
    A = np.zeros((n,3))
    # vector of observations
    l = np.zeros(n)
    # vector of estimates 
    e = np.zeros(n)
    # estimated (e) - observed (l)
    b = np.zeros(n)

    # work out the euler poles...
    omega = np.zeros(3)

    # check to see if we want to input a plate model as apriori values
    # will read in the JSON configuration file and apply the plate value e.g. AUST to the
    # omega values
    if args.apriori:
        nuvel1A = nuvel1APlateModel()
        omega[0] = nuvel1A['AUST']['omega_X']
        omega[1] = nuvel1A['AUST']['omega_Y']
        omega[2] = nuvel1A['AUST']['omega_Z']
        #tmp = readPlateModelFile(args)
        #print("apriori values to be set for:",args.apriori,tmp)
        #omega = np.radians( tmp / 3600. / 1000. )

    domega = np.ones(3)

    proc_sites = []
    oldIndices = []
    newIndices = []
    design_idx = 0
    for site in snx['stations']:
        #print("snx->",site) 
        # default behaviour is to only use one velocity per station
        ctr = 1
        if args.allVelocities:
            ctr = np.size(snx[site]['VELX'])
            
        for i in range(0,ctr):
        #for i in range(0,1): #np.size(snx[site]['VELX'])):
            #=========================================
            # Design matrix A should be:
            #                            |  0  Z -Y |
            #                            | -Z  0  X |
            #                            |  Y -X  0 |
            #  A(3*n,3) =                |----------|
            #                            |  .  .  . |
            #                            |  .  .  . |
            #                            |  .  .  . |
            #
            #=========================================
            #print("i",i,"site",site,snx[site]['VELX_index'][i])
            xidx = snx[site]['VELX_index'][i]
            yidx = snx[site]['VELY_index'][i]
            zidx = snx[site]['VELZ_index'][i]
            #print("xidx:",site, xidx,snx[site]['STAX'][i],snx[site]['STAY'][i],snx[site]['STAZ'][i] )

            A[design_idx,1]   =  snx[site]['STAZ'][i]
            A[design_idx,2]   = -1.*snx[site]['STAY'][i]

            A[design_idx+1,0] = -1.*snx[site]['STAZ'][i]
            A[design_idx+1,2] =  snx[site]['STAX'][i]

            A[design_idx+2,0] =  snx[site]['STAY'][i]
            A[design_idx+2,1] = -1.*snx[site]['STAX'][i]
            
            l[design_idx]   = snx[site]['VELX'][i] 
            l[design_idx+1] = snx[site]['VELY'][i]
            l[design_idx+2] = snx[site]['VELZ'][i]
            
            newQ[design_idx,design_idx]     = Q[xidx,xidx]
            newQ[design_idx+1,design_idx+1] = Q[yidx,yidx]
            newQ[design_idx+2,design_idx+2] = Q[zidx,zidx]
            
            # GET the cross-correlations of the matrix
            #print("Q[:,xidx]",xidx,Q[:,xidx])
            newIndices.append(design_idx)
            newIndices.append(design_idx+1)
            newIndices.append(design_idx+2)

            oldIndices.append(xidx)
            oldIndices.append(yidx)
            oldIndices.append(zidx)

            design_idx = design_idx + 3
            #print("Forming A matrix:",site,i)
    # Find the corr-correlations
    for i in range(0,np.size(newIndices)):
        for j in range(0,np.size(oldIndices)):
            newI = newIndices[i]
            newJ = newIndices[j]
            oldI = oldIndices[i]
            oldJ = oldIndices[j]

            newQ[newI,newJ] = Q[oldI,oldJ]
            newQ[newJ,newI] = Q[oldJ,oldI]
        
    del(Q)
    Q = newQ
    # What about iterating...
    # cofactor matrix (AtWA)
    P     = np.linalg.pinv(Q)
    
    AtP   = np.dot(A.T,P)
    N     = np.dot(AtP,A)

    #================================================
    # add tight constraint to the parameters so that the estimates do not change
    #================================================
    if args.constrain:
        C     = np.eye(np.shape(N)[0]) * 10E-32
        C_inv = np.linalg.pinv(C)
        N     = np.add(N,C_inv)

    N_inv = np.linalg.pinv(N)
    #
    # degrees of freedom:
    # n is the number of measuring stations
    # Caution: the third compnent of velocity vector depends on the first
    # two components, and therefore cannot be used as an independent measure
    # so this becomes 2*n not 3*n. is this correct??, we are not putting it in
    # in terms of radial components??
    #
    dof = (n/3 * 2) - 3

    sigma0 = np.sqrt( np.dot(np.dot(b.T,P),b) / dof)
    print("sigma0",sigma0,n,dof)

    # calculate the cofactor matrix of unknowns..  
    Qxx = N_inv * sigma0**2

    # uncertainties on the estimates..
    omega_sX = np.sqrt(Qxx[0,0])
    omega_sY = np.sqrt(Qxx[1,1])
    omega_sZ = np.sqrt(Qxx[2,2])
    #omega_sX = np.sqrt(Qxx[0,0]) * 3600. * 1000. 
    #omega_sY = np.sqrt(Qxx[1,1]) * 3600. * 1000.
    #omega_sZ = np.sqrt(Qxx[2,2]) * 3600. * 1000.
    #print("N_inv",omega_sX,omega_sY,omega_sZ)

    #================================================
    # check to see if we need to iterate the solution
    #================================================
    itr = 0
    MAXITRS = 5
    epsilon = np.linalg.norm(np.eye(3) * np.radians(0.001 / 3600. /1000.))
    #eps = np.radians(0.00001 / 3600. /1000.)
    eps = np.radians(0.001 / 3600. /1000.)
    #print("epsilon",epsilon,eps)

    while(np.linalg.norm(domega) > epsilon and itr < MAXITRS):
        print("iteration:",itr,"domega",np.linalg.norm(domega))

        # observed minus computed
        b = np.subtract(l,e) 

        # AtWb
        AtPb = np.dot(AtP,b)
        domega = np.dot(N_inv,AtPb)
        #print("domega",domega)

        # Apply the adjustment to the parameters
        omega = np.add(omega,domega)

        # now work out what the estimated velocity would be
        # with the plate model 
        tmpModels = {}
        tmpModels['AUST'] = {} 
        tmpModels['AUST']['omega_X'] = omega[0]
        tmpModels['AUST']['omega_Y'] = omega[1]
        tmpModels['AUST']['omega_Z'] = omega[2]
        tmpModels['AUST']['omega_sX'] = omega_sX 
        tmpModels['AUST']['omega_sY'] = omega_sY 
        tmpModels['AUST']['omega_sZ'] = omega_sZ 
        tmpSites = {}
        
        for site in snx['stations']:
            tmpSites[site] = {}
            tmpSites[site]['X'] = []
            tmpSites[site]['Y'] = []
            tmpSites[site]['Z'] = []

            for i in range(0,np.size(snx[site]['STAX'])):
                tmpSites[site]['X'].append(snx[site]['STAX'][i])
                tmpSites[site]['Y'].append(snx[site]['STAY'][i])
                tmpSites[site]['Z'].append(snx[site]['STAZ'][i])
                proc_sites.append(site+'_'+str(i)+'_X')
                proc_sites.append(site+'_'+str(i)+'_Y')
                proc_sites.append(site+'_'+str(i)+'_Z')
            
        #print("\ntmp plate model:",tmpPlateModel,"\n")
        #tmpSites = determineVelocities(tmpSites,tmpPlateModel)
        tmpSites = determineVelocities(tmpSites,tmpModels['AUST'])

        # TODO am I mainting the right order here??
        # is the e vector matching the same sites for the l vector???

        # set up the new vector of estimates based on the last set of
        # adjusted parameters..
        tctr = 0
        keylist = tmpSites.keys()
        keylist.sort()
        for tsite in keylist:
            #print("e tsite:",tsite)
            for i in range(0,np.size(tmpSites[tsite]['velX'])):
                e[tctr] = tmpSites[tsite]['velX'][i] 
                tctr = tctr + 1
                e[tctr] = tmpSites[tsite]['velY'][i]
                tctr = tctr + 1
                e[tctr] = tmpSites[tsite]['velZ'][i]
                tctr = tctr + 1

        sigma0 = np.sqrt( np.dot(np.dot(b.T,P),b) / dof)
        print("sigma0",sigma0,n,dof)

        # calculate the cofactor matrix of unknowns..  
        Qxx = N_inv * sigma0**2

        # uncertainties on the estimates..
        omega_sX = np.sqrt(Qxx[0,0])
        omega_sY = np.sqrt(Qxx[1,1])
        omega_sZ = np.sqrt(Qxx[2,2])

        # quality of corrected observations 
        # Ql = AQxA^T
        Ql = np.dot(np.dot(A,N_inv),A.T)
        P = np.linalg.pinv(Ql)
        #print("N_inv",omega_sX,omega_sY,omega_sZ)
        itr = itr +1
  
    b = np.subtract(l,e)
    # another method of displaying the results, will stick to cartesian representation 
    #
    # angular_velocity = np.degrees( np.sqrt(omega[0]**2 + omega[1]**2 + omega[2]**2) )*10**6
    # omega_lat = np.degrees(np.arctan( omega[2] / np.sqrt(omega[0]**2 + omega[1]**2)))
    # omega_lon = np.degrees(np.arctan( omega[1] / omega[0] ))
    # print("angular velocity, omega_lat, omega_lon:",angular_velocity,omega_lat,omega_lon)

    # convert omega from radians to mas:
    pomega = np.degrees(omega)*3600.*1000.
    pomega_sX = np.sqrt(Qxx[0,0]) * 3600. * 1000. 
    pomega_sY = np.sqrt(Qxx[1,1]) * 3600. * 1000.
    pomega_sZ = np.sqrt(Qxx[2,2]) * 3600. * 1000.

    print("")
    print("================================================")
    print("Plate Model Results")
    print("{0:.5f} +/-{1:.5f} {2:.5f} +/-{3:.5f} {4:.5f} +/-{5:.5f}".format(
          pomega[0], pomega_sX, pomega[1], pomega_sY, pomega[2], pomega_sZ))
    print("================================================")
    print("CATREF Results:")
    print("1.526 +/-0.015 1.185 +/-0.016 1.215 +/-0.016")
    print("================================================")
    print("ITRF2008 values:")
    print("1.504 +/-0.007 1.172 +/-0.007 1.228 +/-0.007")
    print("================================================")
    # calculate the cofactor matrix of residuals
    # Qvv = W^-1 - A N^-1 A^T
    Qvv = np.subtract(P,np.dot(np.dot(A,N_inv),A.T))

    # TODO move these test to show the results in V_n and V_e
    if args.baardsTest:
        btest = baardsTest(b,Qvv,P)
    if args.tauTest:
        ttest = tauTest(b,Qvv,sigma0,dof,proc_sites)

    # return omega_V in units of mas
    #Qxx = np.degrees(Qxx)  * 3600. * 1000.

    return omega, Qxx 

def tauTest(V,Qvv,sigma0,dof,stns):
    """
    tau test

    tauTest(V,Qvv,sigma0)

    Input:
            V       vector of residuals
            Qvv     cofactor matrix of residuals
            P       apripro sigma of observations


    Output:

    """

    v_standard = np.zeros(np.size(V))
    v_test     = np.zeros(np.size(V))

    # work out the standard residuals
#   print("=========================================")
#   #print("residual,normalized_residual,test_residual,pval")
#   #print("size of residuals:",np.size(V))
#   for i in range(0,np.size(V)):
#       tau_calc = V[i] / ( sigma0 * np.sqrt(Qvv[i,i]) )
#       cval = t.isf(0.05/2.,dof)
#       tau_crit = cval * np.sqrt(dof) / np.sqrt(dof -1 + cval**2 ) 
#       #v_standard[i] = V[i] / np.sqrt(Qvv[i,i])
#       #v_test[i] = v_standard[i] / sigma0
#       if tau_crit > tau_calc and tau_calc > (tau_crit*-1.) :
#           print("{0} {1:6.3f} ".format(stns[i],V[i]))
#       else:
#           print("{0} {1:6.3f} ***".format(stns[i],V[i]))

#   print("=========================================")
    return v_test
    
def baardsTest(V,Qvv,P):
    """
    baards data snooping algorithm

    baardsTest(V,Qvv)

    Input:
            V       vector of residuals
            Qvv     cofactor matrix of residuals
            P       apripro sigma of observations

    Output:

    """

    v_standard = np.zeros(np.size(V))
    v_test     = np.zeros(np.size(V))
    # work out the standard residuals
    for i in range(0,np.size(V)):
        #print(i,Qvv[i,i])
        v_standard[i] = V[i] / np.sqrt(Qvv[i,i])
        v_test[i] = v_standard[i] / np.sqrt(P[i,i])

    return v_test

def determineVelocities(sites,plateModel,Trate = [0.0,0.0,0.0]):
    """

    sites = determineVelocities(sites,plateModel)
    sites = determineVelocities(sites,plateModel,Trate)

    Work out the velocity of the station based upon the plate model, supplied

     Input:
            sites
            plateModel
            Trate        vector of translation, default is [0. 0. 0.]
                         in mm.

    """
    #print("In determine velocities",sites)
    # Trate is set in mm, need to
    # convert Trate to m, if greater than 0
    if np.linalg.norm(Trate) > 0.00001:
        if Trate[0] > 0.00001 : Trate[0] = Trate[0]/1000.
        if Trate[1] > 0.00001 : Trate[1] = Trate[1]/1000.
        if Trate[2] > 0.00001 : Trate[2] = Trate[2]/1000.

    for site in sites:
        #print("Calculating velocity for site:",site)
        # make sure an XYZ is defined, if not work it out from the lat, lon, h 
        if ('X' not in sites[site]) or ('Y' not in sites[site]) or ('Z' not in sites[site]):
            if 'lat' not in sites[site] or 'lon' not in sites[site]:
                print(site,"does not have any coordinates, unable to compute a velocity")
                continue

            # TODO check the validity of this approach
            # this is using an ellipse but euler is assuming a sphere...
            # work out the XYZ of the station..
            if 'h' not in sites[site]:
                sites[site]['h'] = 0.0

            #XYZ = llh2xyz(ell2sph(sites[site]['lat']),sites[site]['lon'],sites[site]['h'])
            XYZ = gt.llh2xyz(sites[site]['lat'],sites[site]['lon'],sites[site]['h'])
            #print("SITES:",site,sites[site])
            sites[site]['X'] = XYZ[0]
            sites[site]['Y'] = XYZ[1]
            sites[site]['Z'] = XYZ[2]
            #print(sites[site])
        
        # Convert from mas to radians 
        #omega_X  = np.radians( ( plateModel['omega_X'] / 1000. ) / 3600.) 
        #omega_Y  = np.radians( ( plateModel['omega_Y'] / 1000. ) / 3600.) 
        #omega_Z  = np.radians( ( plateModel['omega_Z'] / 1000. ) / 3600.) 
        #omega_sX = np.radians( ( plateModel['omega_sX'] / 1000. ) / 3600.) 
        #omega_sY = np.radians( ( plateModel['omega_sY'] / 1000. ) / 3600.) 
        #omega_sZ = np.radians( ( plateModel['omega_sZ'] / 1000. ) / 3600.) 

        omega_X  = plateModel['omega_X'] 
        omega_Y  = plateModel['omega_Y'] 
        omega_Z  = plateModel['omega_Z'] 
        omega_sX = plateModel['omega_sX'] 
        omega_sY = plateModel['omega_sY'] 
        omega_sZ = plateModel['omega_sZ'] 

        sites[site]['velX'] = [] 
        sites[site]['velY'] = [] 
        sites[site]['velZ'] = [] 
        sites[site]['vel'] = [] 
        sites[site]['velsX'] = [] 
        sites[site]['velsY'] = [] 
        sites[site]['velsZ'] = [] 

        #print("Translation rate:",Trate)
        for i in range(0,np.size(sites[site]['X'])):
            # Now work out the velocity
            sites[site]['velX'].append( ( sites[site]['Z'][i] * omega_Y 
                                        - sites[site]['Y'][i] * omega_Z )
                                        + Trate[0] )

            sites[site]['velY'].append( ( sites[site]['X'][i] * omega_Z -  
                                             sites[site]['Z'][i] * omega_X )
                                           + Trate[1] )
        
        
            sites[site]['velZ'].append( ( sites[site]['Y'][i] * omega_X -  
                                             sites[site]['X'][i] * omega_Y )
                                           + Trate[2] )

            sites[site]['vel'].append( np.sqrt( sites[site]['velX'][i]**2 
                                           + sites[site]['velY'][i]**2 
                                           + sites[site]['velZ'][i]**2 ) )

            # sigmas of velocities:
            sites[site]['velsX'].append( ( sites[site]['Z'][i] * omega_sY -
                                        sites[site]['Y'][i] * omega_sZ ) )# + sTrate[0]
            sites[site]['velsY'].append( ( sites[site]['X'][i] * omega_sZ -
                                        sites[site]['Z'][i] * omega_sX ) )# + sTrate[1]
            sites[site]['velsZ'].append( ( sites[site]['Y'][i] * omega_sX -
                                        sites[site]['X'][i] * omega_sY ) )# + sTrate[2]

    return sites

def readConfig(args):
    """

    args = readConfig(args)

    read in the configuration file in JSON format

    """
    with open(args.config) as json_data_file:
        config = json.load(json_data_file)
    #print(config)

    # check to see if we want to exclude stations, or process on a limited set
    if 'station_options' in config:
        if 'exclude' in config['station_options']:
            args.exclude=config['station_options']['exclude']
            print("args.exclude",args.exclude)
        
        if 'include' in config['station_options']:
            args.include=config['station_options']['include']

    return args 

def readPlateModelFile(args):
    """
    omega = readPlateModelFile(args)
    """
    omega = np.zeros(3)
    # plate model should be in JSON format
    print("Trying to open:",args.model)
    with open(args.model) as json_data_file:
        config = json.load(json_data_file)

    if 'plateModels' in config:
        #print(config['plateModels'],args.apriori)#if 'AUST' in config['plateModels']:
        if args.apriori[0] == 'AUST' and 'AUST' in config['plateModels']:
            omega[0]=config['plateModels']['AUST']['omega_X']
            omega[1]=config['plateModels']['AUST']['omega_Y']
            omega[2]=config['plateModels']['AUST']['omega_Z']
            #omega[3]=config['plateModels']['AUST']['omega_sX']
            #omega[4]=config['plateModels']['AUST']['omega_sY']
            #omega[5]=config['plateModels']['AUST']['omega_sZ']
            #omega[6]=config['plateModels']['AUST']['omega']
            #omega[7]=config['plateModels']['AUST']['std_omega']
        
    return omega



#==============================================================================
if __name__ == "__main__":
    #=======
    parser = argparse.ArgumentParser(description='Calculate the plate model from a GPS velocity SINEX file',
                                     formatter_class=argparse.RawTextHelpFormatter, 
                                     epilog='''\


    Example:

    1)  To calculate the plate euler pole for the Australian plate, from the ITRF2008 GPS velocity solution, compute
        the tau statistics, and plot the stations used:

        > ./plateModel.py -f ./data/ITRF2008-TRF.SNX -p data/Australian.plt --tau --plot  
            
    2)  If you plan to include exclude stations from the solution, it will be faster and more efficient to convert 
        SINEX file into a binary file:
        
        > ./plateModel.py -f ./data/ITRF2008-TRF.SNX --save itrf2008

        Then calculate the euler models from the binary file with:
        
        > ./plateModel.py -b itrf2008 -p data/Australian.plt --tau --plot

        Then exclude a station from the list, and recompute the euler vector:
        
        > ./plateModel.py -b itrf2008 -p data/Australian.plt --tau --plot --exclude COCO TONG

    3)  It can be useful to import a plate model and then check to see if the GPS velocities fit to that model.
        To this we need to tightly constrain the plate model parameters.
        
        > ./plateModel.py -b itrf2008 -p data/Australian.plt --model ./config/plateModelsITRF2008.cfg --apriori AUST --constrain --tau

    4)  Compare with ITRF 2008 results:

        > ./plateModel.py -b itrf2008 -p ./data/Australian.plt --config ./config/itrf2008_australian_plate.cfg


            ''')

    parser.add_argument('--plot',dest='plot',default=False,action='store_true',help="Plot the plate boundaries")
    parser.add_argument('--gmt',dest='gmt',default=True,action='store_true',help="Output the velocity residuals into a GMT plot file for psvelo")

    #parser.add_argument('--save',dest='save',default=False,action='store_true',help="Store the parsed sinex file as a pickle data structure, and the covariance matrix as a numpy binary file")
    parser.add_argument('-f',dest='sinexFile',default='',help="Sinex File to read to obtain estimates")
    parser.add_argument('-p',dest='plateFiles',default='',nargs='+',help="Plate boundary files based on Nu-vel")
    parser.add_argument('--save',dest='save',default='',help="Store the parsed sinex file as a pickle data structure, and the covariance matrix as a numpy binary file")
    parser.add_argument('-b',dest='binaryData',default='',help="Import a previously saved binary files")

    parser.add_argument('--exclude',dest='exclude',nargs='+',help="Stations to exclude eg MAC1 or MAC1_0")
    parser.add_argument('--include',dest='include',nargs='+',help="Only include the following list of stations/velocities eg MAC1 or MAC1_0")
    parser.add_argument('--AV','--allVelocities',dest='allVelocities',default=False,action='store_true',help="Include multiple velocity solutions from each station")
    parser.add_argument('--config',dest='config',default='',help="Configuration file")

#   parser.add_argument('--model',dest='model',default='',help="Configuration file for plate models")
    parser.add_argument('--apriori',dest='apriori',nargs='+',help="plates to add apriori to setup as apriori values eg. Australia")
    parser.add_argument('--constrain',dest='constrain',default=False,action='store_true',help="Constrain the plate model estimates")

    parser.add_argument('--tau',dest='tauTest',default=False,action='store_true',help="Carry out a Tau test")
    parser.add_argument('--baard',dest='baardsTest',default=False,action='store_true',help="Carry out a baards test")

    args = parser.parse_args()
    #======
    lon = []
    lat = []
    codes = []

    # read in the configuration file in JSON format
    if args.config:
        args = readConfig(args)

# Need to install shapely module for this to work
#   if args.plateFiles:
#       for pf in args.plateFiles:
#           #print("Found the following:",pf)
#           plate = np.genfromtxt(pf,skip_header=1)

#           # define a polygon in shapely
#           platePolygon = Pg(plate)

    if args.sinexFile:
        snx,cova = sinex.readSINEX(args.sinexFile)
        #print("COVA:",np.shape(cova))

        if args.save:
            with open(args.save+'.pkl','wb') as pklID:
                pickle.dump(snx,pklID,2)
            # now save the covariance matrix as a numpy compressed file    
            np.savez_compressed(args.save,cova=cova)
            sys.exit(0)

    if args.binaryData:
        npzfile = np.load(args.binaryData+'.npz')
        cova = npzfile['cova']

        # Now read the pickle file
        with open(args.binaryData+'.pkl','rb') as pklID:
            snx = pickle.load(pklID)
      
    # check to see which stations are within the plate boundary
    if 'stations' in snx:
        for station in snx['stations']:
            codes.append(station)

        # Need shapely module for this to work (comment out the two lines above to restore logic)
#       for station in snx['stations']:
#           a = Point(snx[station]['lon'],snx[station]['lat'])
#           if a.within(platePolygon):
#               #print("Station is within the plate:",station,snx[station]['lon'],snx[station]['lat'])
#               #lon.append(snx[station]['lon'])
#               #lat.append(snx[station]['lat'])
#               codes.append(station)
#           else:
#               if not args.exclude :
#                   args.exclude = []
#               args.exclude.append(station)
#   else :
#       for station in snx:
#           a = Point(snx[station]['lon'],snx[station]['lat'])
#           if a.within(platePolygon):
#               #print("Station is within the plate:",station,snx[station]['lon'],snx[station]['lat'])
#               #lon.append(snx[station]['lon'])
#               #lat.append(snx[station]['lat'])
#               codes.append(station)
#           else:
#               if not args.exclude :
#                   args.exclude = []
#               args.exclude.append(station)
    #==========================================================================
    # formup the matrices based on the plate(s) being solved for
    # indices to keep
    # only interested in the velocity estimates
    #==========================================================================
    if args.exclude:
        snx, cova = sinex.removeSolution(args.exclude,snx,cova)

    omega,omega_V = determineEulerPole(snx,cova,args)

    #==========================================================================
    # work out the velocity residuals
    #  - define the new estimated plate
    #  - work out the velocities based on these estimates
    #==========================================================================

    sites = {}
    for site in snx['stations']:
        sites[site] = {}
        sites[site]['X'] = []
        sites[site]['Y'] = []
        sites[site]['Z'] = []

        for i in range(0,np.size(snx[site]['STAX'])):
            sites[site]['X'].append(snx[site]['STAX'][i])
            sites[site]['Y'].append(snx[site]['STAY'][i])
            sites[site]['Z'].append(snx[site]['STAZ'][i])

    #==========================================================================
    print("")
    print("Velocity residuals, compared to Nuvel1aNNR:")
    print("======================")
    print("Site     Ve      Vn")
    print("        (mm)    (mm)")
    print("======================")
    #==========================================================================
    nuvel1A = nuvel1APlateModel()
    nvel = determineVelocities(sites,nuvel1A['AUST'])

    dns = []
    des = []
    dvec = []
    stations = []
    for site in sites:
        for i in range(0,np.size(snx[site]['VELX'])):
            vxyz = np.zeros(3)
            vxyz[0] = snx[site]['VELX'][i]
            vxyz[1] = snx[site]['VELY'][i]
            vxyz[2] = snx[site]['VELZ'][i]

            lat_p,lon_p,h = gt.xyz2llh(snx[site]['STAX'][i],snx[site]['STAY'][i],snx[site]['STAZ'][i])
            neu = gt.geo2topo(lat_p, lon_p, vxyz )

            # Velocities as determined from a pre-defined plate model
            lat_e,lon_e,h = gt.xyz2llh(nvel[site]['X'][i], nvel[site]['Y'][i], nvel[site]['Z'][i])
            vxyz_est = np.zeros(3)
            vxyz_est[0] = nvel[site]['velX'][i]
            vxyz_est[1] = nvel[site]['velY'][i]
            vxyz_est[2] = nvel[site]['velZ'][i]
            neu_est = gt.geo2topo(lat_e,lon_e, vxyz_est )
            lat.append(lat_e)
            lon.append(lon_e)

            dn = (neu[0] - neu_est[0]) * 1000.
            de = (neu[1] - neu_est[1]) * 1000.
            key = site+"_"+str(i)
            dns.append(dn)
            des.append(de)
            dvec.append(np.sqrt(dn**2 + de**2))
            stations.append(key)

    ind = np.argsort(dvec)
    for i in ind:
        print("{0} {1:>7.3f} {2:>7.3f} ".format(stations[i],des[i],dns[i]))

    print("======================")

    if args.gmt:
        print("\n****----> GMT plot file -> velocity_residuals_nuvel.gmt\n")
        with open("velocity_residuals_nuvel.gmt","w") as vres:
            for i in ind:
                print("{0:>7.3f} {1:>7.3f} {2:>7.3f} {3:>7.3f} {4}".format(lon[i],lat[i],des[i],dns[i],stations[i]),file=vres)

    #==========================================================================
    print("")
    print("===========================================")
    print("Velocity residuals SINEX - Velocities")
    print("derived from estimated plate model")
    print("===========================================")
    print("Site     Ve      Vn")
    print("        (mm)    (mm)")
    print("===========================================")
    #==========================================================================
    plateModels = {}
    plateModels['AUST'] = {} 
    plateModels['AUST']['omega_X'] = omega[0]
    plateModels['AUST']['omega_Y'] = omega[1]
    plateModels['AUST']['omega_Z'] = omega[2]
    plateModels['AUST']['omega_sX'] = omega_V[0,0] 
    plateModels['AUST']['omega_sY'] = omega_V[1,1]
    plateModels['AUST']['omega_sZ'] = omega_V[2,2]

    #print("estimated    nuvel")
    #print("omega_X {0:>15.15f} {1:>15.15f}".format(plateModels['AUST']['omega_X'], nuvel1A['AUST']['omega_X']))
    #print("omega_Y {0:>15.15f} {1:>15.15f}".format(plateModels['AUST']['omega_Y'], nuvel1A['AUST']['omega_Y']))
    #print("omega_Z {0:>15.15f} {1:>15.15f}".format(plateModels['AUST']['omega_Z'], nuvel1A['AUST']['omega_Z']))

    sites_est = determineVelocities(sites,plateModels['AUST'])

    dns = []
    des = []
    dvec = []
    stations = []

    for site in sites:
        for i in range(0,np.size(snx[site]['VELX'])):
            vxyz = np.zeros(3)
            vxyz[0] = snx[site]['VELX'][i]
            vxyz[1] = snx[site]['VELY'][i]
            vxyz[2] = snx[site]['VELZ'][i]

            lat_p,lon_p,h = gt.xyz2llh(snx[site]['STAX'][i],snx[site]['STAY'][i],snx[site]['STAZ'][i])
            neu = gt.geo2topo(lat_p, lon_p, vxyz )

            # Velocities as determined from a pre-defined plate model
            vxyz_est = np.zeros(3)
            vxyz_est[0] = sites_est[site]['velX'][i]
            vxyz_est[1] = sites_est[site]['velY'][i]
            vxyz_est[2] = sites_est[site]['velZ'][i]

            lat_e,lon_e,h = gt.xyz2llh(sites_est[site]['X'][i],sites_est[site]['Y'][i],sites_est[site]['Z'][i])
            neu_est = gt.geo2topo(lat_e,lon_e, vxyz_est )
            lat.append(lat_e)
            lon.append(lon_e)

            dn = (neu[0] - neu_est[0]) * 1000.
            de = (neu[1] - neu_est[1]) * 1000.

            key = site+"_"+str(i)
            dns.append(dn)
            des.append(de)
            dvec.append(np.sqrt(dn**2 + de**2))
            stations.append(key)

    ind = np.argsort(dvec)
    for i in ind:
        print("{0} {1:>7.3f} {2:>7.3f} ".format(stations[i],des[i],dns[i]))

    print("===========================================")
    if args.gmt:
        print("\n****----> GMT plot file -> velocity_residuals_plate.gmt\n")
        with open("velocity_residuals_plate.gmt","w") as vres:
            for i in ind:
                print("{0:>9.5f} {1:>9.5f} {2:>9.5f} {3:>9.5f} {4}".format(lon[i],lat[i],des[i],dns[i],stations[i]),file=vres)
    print("======================")

    #==========================================================================
    if args.plot:
        import matplotlib.pyplot as plt
        from mpl_toolkits.basemap import Basemap
        from matplotlib.patches import Polygon
        # global map 
        #pmap = Basemap(projection='robin',area_thresh = 1000.0,lat_0=-45,lon_0=10,resolution='l')
        # australian map
        # lat_ts latitude of true scale
        # lat_0,lon_0  central point

        #pmap = Basemap(width=13000000,height=9000000,projection='laea',
        #        lat_ts=-30, lat_0=-30, lon_0=125,
        #        resolution='l')
        #fig = plt.figure()
        #ax1 = fig.add_subplot(111)
        pmap = Basemap(width=13000000,height=9000000,projection='lcc',
                lat_1=15, lat_2=-60,
                lat_0=-30, lon_0=125,
                resolution='l')
        #       resolution='l',ax=ax1)

        # draw coastlines, country boundaries, fill continents.
        pmap.drawcoastlines(linewidth=0.25)
        pmap.drawcountries(linewidth=0.25)
        #pmap.fillcontinents(color='gray',lake_color='aqua')
        pmap.shadedrelief()
        # draw the edge of the map projection region (the projection limb)
        #pmap.drawmapboundary(fill_color='aqua')

        # draw lat/lon grid lines every 30 degrees.
        pmap.drawmeridians(np.arange(0,360,30))
        pmap.drawparallels(np.arange(-90,90,30))
        x, y = pmap(lon, lat)
        pmap.plot(x,y,'bo',markersize=5,picker=True)
        plotPolygon(plate[:,0],plate[:,1],pmap)

#       if 1:
#           def onpick3(event):
#               ind = event.ind
#               print 'onpick3 scatter:', ind, np.take(x, ind), np.take(y, ind)

#           x, y = pmap(lon, lat)
#           col = pmap.plot(x,y,'bo',markersize=5,picker=True)
#           plotPolygon(plate[:,0],plate[:,1],pmap)
#           fig.canvas.mpl_connect('pick_event', onpick3)

        plt.show()

