#!/usr/bin/env python
from __future__ import division, print_function, absolute_import

import numpy as np

def refell(type):
    '''
    Useage:  [a,b,e2,finv]=refell(type)
    Input:   type - reference ellipsoid type (char)
                    CLK66 = Clarke 1866
                    GRS67 = Geodetic Reference System 1967
                    GRS80 = Geodetic Reference System 1980
                    WGS72 = World Geodetic System 1972
                    WGS84 = World Geodetic System 1984
    Output:  a    - major semi-axis (m)
             b    - minor semi-axis (m)
             e2   - eccentricity squared
             finv - inverse of flattening
    '''
    type.upper()

    if (type=='CLK66' or type=='NAD27'):
        a=6378206.4 
        finv=294.9786982 
    elif type=='GRS67':
        a=6378160.0 
        finv=298.247167427 
    elif (type=='GRS80' or type=='NAD83'):
        a=6378137.0 
        finv=298.257222101 
    elif (type=='WGS72'):
        a=6378135.0 
        finv=298.26 
    elif (type=='WGS84'):
        a=6378137.0 
        finv=298.257223563 

    f = 1./finv 
    b = a*(1.-f) 
    e2 = 1.-(1.-f)**2 
    #print("e2",e2)
    return a,b,e2,finv

def llh2xyz(lat,lon,h,datum='GRS80') :
    '''
[x,y,z] = llh2xyz(lat,lon,h,a,finv)
#
#    [x,y,z] = llh2xyz(lat,lon,h) defaults to GRS80
#
#    Converts ellipsoidal coordinates to cartesian.
#
#     Input:   lat  - vector of ellipsoidal latitudes (degrees)
#              lon  - vector of ellipsoidal E longitudes (degrees)
#              h    - vector of ellipsoidal heights (m)
#              a    - semi-axis (m) (default GRS80)
#              finv - inverse flattening of ellipsoid (default GRS80)
#
#     Output:  [x,y,z] cartesian coordinates (m)
#
#     Ref: Geodesy: Wolfgang Torge pg. 51-52
#          eq 3.26 and 3.34
    '''
    
    a,b,e2,finv = refell(datum)
    v=a/np.sqrt(1.-e2*np.sin(np.radians(lat))*np.sin(np.radians(lat)))
    x=(v+h)*np.cos(np.radians(lat))*np.cos(np.radians(lon))
    y=(v+h)*np.cos(np.radians(lat))*np.sin(np.radians(lon))
    z=(v*(1.-e2)+h)*np.sin(np.radians(lat));

    return ([x,y,z])


def xyz2llh(X,Y,Z,datum='GRS80'):
    """
    calculate geodetic coordinates
    latitude, longitude, height given Cartesian
    coordinates X,Y,Z, and reference ellipsoid
    values semi-major axis (a) and the inverse
    of flattening (finv).
    """
    a,b,e2,finv = refell(datum)
    elat = 1.e-12
    eht  = 1.e-5

    # Initial values
    p = np.sqrt(X*X+Y*Y)
    lat = np.arctan2(Z,p*(1.-e2))

    h = 0.
    dh = 1.
    dlat = 1.
    maxit = 8

    for i in range(0,maxit):
        lat0 = lat
        h0   = h

        v = a/np.sqrt(1-e2*np.sin(lat)*np.sin(lat))
        h = p/np.cos(lat)-v

        lat = np.arctan2(Z, p*(1.-e2*v/(v+h)))
        dlat = np.abs(lat-lat0)
        dh = np.abs(h-h0)
        if dlat > elat and dh > eht:
            #print("Completed iteration",i)
            break

    return(np.degrees(lat),np.degrees(np.arctan2(Y,X)), h)
    
#def topocent(X,dx) :
#    """
#    %TOPOCENT  Transformation of vector dx into topocentric coordinate
#    %          system with origin at X.
#    %          Both parameters are 3 by 1 vectors.
#    %          Output: D    vector length in units like the input
#    %                  Az   azimuth from north positive clockwise, degrees
#    %                  El   elevation angle, degrees
#
#    """
#    dtr = np.pi/180.
#    phi,lambda,h = togeod(6378137,298.257223563,X(1),X(2),X(3))
#    cl = cos(lambda*dtr); sl = sin(lambda*dtr)
#    cb = cos(phi*dtr); sb = sin(phi*dtr)
#    F = [-sl -sb*cl cb*cl;
#      cl -sb*sl cb*sl;
#       0    cb   sb];
#
#    local_vector = np.transpose(F)*dx
#
#    E = local_vector(1)
#    N = local_vector(2)
#    U = local_vector(3)
#
#    hor_dis = np.sqrt(E**2 + N**2)
#    if hor_dis < 1.e-20:
#        Az = 0
#        El = 90
#    else:
#        Az = atan2(E,N)/dtr;
#        El = atan2(U,hor_dis)/dtr;
#    
#    if Az < 0:
#        Az = Az+360
#
#    D = np.sqrt(dx(1)**2+dx(2)**2+dx(3)**2)
#
##    return(Az,El,D)

def decdeg2dms(dd):
    """
    deg,mm,ss = decdeg2dms(dd)

    Convert decimal degree to degrees, minutes, seconds

     Input:
             dd    decimal degrees (float)

     Output:
             deg   degrees
             mm    minutes
             ss    seconds
    """
    mm,ss = divmod(dd*3600,60)
    deg,mm = divmod(mm,60)
    return deg,mm,ss
  
def geo2topo(lat,lon,xyz):
    """
    
    Vxyz = geo2topo(lat,long,Vxyz)

    Convert a geocentric (XYZ to neu) 
   
     Input:
      lat      latitude (degrees)
      long     longitude (degrees)
      xyz      numpy array where: 
                     xyz[0] = X velocity
                     xyz[1] = Y velocity
                     xyz[2] = Z velocity

     Output:
      neu      numpy array with shape (3,)
       
    """
    # get the rotation matrix for converting XYZ to NEU 
    r = rotationMatrix(lat,lon)
    neu = np.dot(r,xyz)
    #print("lat,lon",lat,lon)
    #print("r",r)
    #print("xyz:",xyz)
    #print("neu:",neu)
    #print("Now invert the process:")
    #XYZ = np.dot(r.T,neu)
    #print("STARTED with xyz:",xyz)
    #print("ENDED   with xyz:",XYZ)
    return neu

def topo2geo(lat,lon,neu):
    """
    
    xyz = topo2geo(lat,long,neu)

    Convert a geocentric (XYZ to neu) 
   
     Input:
      lat      latitude (degrees)
      long     longitude (degrees)
      neu      numpy array where: 
                     neu[0] = North
                     neu[1] = East
                     neu[2] = Up

     Output:
      xyz      numpy array with shape (3,)
       
    """
    # get the rotation matrix for converting XYZ to NEU 
    r = rotationMatrix(lat,lon)
    xyz = np.dot(r.T,neu)
    return xyz 

def rotationMatrix(lat,lon):
    """

    R = rotationMatrix(lat,lon)

    Set up a rotation matrix for converting XYZ velocities into
    neu velocities

    R = | -sin(lat)cos(lon) -sin(lat)sin(lon) cos(lat) |
        | -sin(lon)         cos(lon)          0        |
        | cos(lat)cos(lon)  cos(lat)sin(lon)  sin(lat) |

     Input:
      lat     latitude in degrees
      lon     longitude in degrees

     Output:
      R       R(3,3) a rotation matrix to take a velocity in x,y,z to neu


    """
    # set up a rotation matrix for converting XYZ to neu:
    R = np.zeros((3,3))

    lat = np.radians(lat)
    lon = np.radians(lon)

    R[0,0] = -np.sin(lat) * np.cos(lon) 
    R[0,1] = -np.sin(lat) * np.sin(lon)
    R[0,2] = np.cos(lat)

    R[1,0] = -np.sin(lon) 
    R[1,1] = np.cos(lon)
    R[1,2] = 0.

    R[2,0] = np.cos(lat) * np.cos(lon)
    R[2,1] = np.cos(lat) * np.sin(lon)
    R[2,2] = np.sin(lat)

    return R
 
def lat_ell2sph(elat,datum='GRS80'):
    """
    slat = ell2sph(elat)
 
    Convert ellipsoidal latitude to spherical latitude.
 
     Input:
         elat   geodetic latitude (degrees)
         datum  default = GRS80 
 
     Output:
         slat   spherical latitude (degrees)

    """
    a,b,e2,finv = refell(datum)
    slat =  np.degrees( np.arctan( (1. - e2) * np.tan(np.radians(elat)) ) )

    return slat   
    
if __name__ == "__main__":
    print("Not much to do here")
