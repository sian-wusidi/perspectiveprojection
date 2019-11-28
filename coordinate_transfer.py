 # -*- coding: utf-8 -*-
'''
Created in August, 2019
@author: Sidi Wu
'''

import math
import numpy as np

#a = 6378.137
#b = 6356752.314245
#b = a
#f = (a - b) / a
#f =1/298.257223563
#e_sq = f * (2-f)
#e_sq= 6.69437999014e-3

a        =  6378.137;
#f         =  1.0/298.257223563;
f=0;
b         =  a * ( 1.0 - f );
f        =  1-b/a;
eccsq    =  1 - b*b/(a*a);
ecc      =  math.sqrt(eccsq);
EARTH_A  =  a;
EARTH_B  =  b;
EARTH_F  =  f;
EARTH_Ecc=  ecc;
EARTH_Esq=  eccsq;
def geodetic_to_ecef(lat, lon, h):
    # (lat, lon) in WSG-84 degrees
    # h in meters
    lamb = math.radians(lat)
    phi = math.radians(lon)
    s = math.sin(lamb)
    N = a / math.sqrt(1 - e_sq * s * s)

    sin_lambda = math.sin(lamb)
    cos_lambda = math.cos(lamb)
    sin_phi = math.sin(phi)
    cos_phi = math.cos(phi)

    x = (h + N) * cos_lambda * cos_phi
    y = (h + N) * cos_lambda * sin_phi
    z = (h + (1 - e_sq) * N) * sin_lambda

    return x, y, z

def ecef_to_enu(x, y, z, lat0, lon0, h0):
    lamb = math.radians(lat0)
    phi = math.radians(lon0)
    s = math.sin(lamb)
    N = a / math.sqrt(1 - e_sq * s * s)

    sin_lambda = math.sin(lamb)
    cos_lambda = math.cos(lamb)
    sin_phi = math.sin(phi)
    cos_phi = math.cos(phi)

    x0 = (h0 + N) * cos_lambda * cos_phi
    y0 = (h0 + N) * cos_lambda * sin_phi
    z0 = (h0 + (1 - e_sq) * N) * sin_lambda
    
    xd = x - x0
    yd = y - y0
    zd = z - z0

    xEast = -sin_phi * xd + cos_phi * yd
    yNorth = -cos_phi * sin_lambda * xd - sin_lambda * sin_phi * yd + cos_lambda * zd
    zUp = cos_lambda * cos_phi * xd + cos_lambda * sin_phi * yd + sin_lambda * zd

    return xEast, yNorth, zUp

def ecef_to_enu(x, y, z, x0,y0,z0):
    lat0,lon0,h0=ECEF_to_geodetic(x0,y0,z0)
    lamb = math.radians(lat0)
    phi = math.radians(lon0)
    #s = math.sin(lamb)
    #N = a / math.sqrt(1 - e_sq * s * s)

    sin_lambda = math.sin(lamb)
    cos_lambda = math.cos(lamb)
    sin_phi = math.sin(phi)
    cos_phi = math.cos(phi)

    #x0 = (h0 + N) * cos_lambda * cos_phi
    #y0 = (h0 + N) * cos_lambda * sin_phi
    #z0 = (h0 + (1 - e_sq) * N) * sin_lambda
    
    xd = x - x0
    yd = y - y0
    zd = z - z0
    #R=np.array([[-sin_lambda,cos_lambda,0],[-sin_phi*cos_lambda,-sin_phi*sin_lambda,cos_phi],[cos_phi*cos_lambda,cos_phi*sin_lambda,sin_phi]])
    xEast = -sin_phi * xd + cos_phi * yd
    yNorth = -cos_phi * sin_lambda * xd - sin_lambda * sin_phi * yd + cos_lambda * zd
    zUp = cos_lambda * cos_phi * xd + cos_lambda * sin_phi * yd + sin_lambda * zd

    return xEast, yNorth, zUp
    #return R

def enu_to_ecef(xE,yN,zU,x0,y0,z0):
    lat0,lon0,h0=ECEF_to_geodetic(x0,y0,z0)
    lamb = math.radians(lon0)
    phi = math.radians(lat0)
    #s = math.sin(lamb)
    #N = a / math.sqrt(1 - e_sq * s * s)

    sin_lambda = math.sin(lamb)
    cos_lambda = math.cos(lamb)
    sin_phi = math.sin(phi)
    cos_phi = math.cos(phi)
    R=np.array([[-sin_lambda,-sin_phi*cos_lambda,cos_phi*cos_lambda],[cos_lambda,-sin_phi*sin_lambda,cos_phi*sin_lambda],[0,cos_phi,sin_phi]])
    X=np.dot(R,np.array([[xE],[yN],[zU]]))
    return X


def geodetic_to_enu(lat, lon, h, lat_ref, lon_ref, h_ref):
    x, y, z = geodetic_to_ecef(lat, lon, h)
    #s = math.sin(lamb)
    #N = a / math.sqrt(1 - e_sq * s * s)

    
    
    return ecef_to_enu(x, y, z, lat_ref, lon_ref, h_ref)


def ECEF_to_geodetic(x,y,z):

    rp     = math.sqrt ( x*x + y*y + z*z );
    dtr =  math.pi/180.0;
    flatgc = math.asin ( z / rp )/dtr;
    testval= np.abs(x) + np.abs(y);
    flon=0.0;
    if ( testval < 1.0e-10):
        flon = 0.0 
    else:
        flon = math.atan2 (y,x)/dtr
    if (flon < 0.0 ):
        flon = flon + 360.0
    p      =  math.sqrt( x*x + y*y );
    llhvec=np.zeros(3)
    if ( p < 1.0e-10 ): 
          flat = 90.0
          if ( z < 0.0 ):
              flat = -90.0
          altkm = rp - rearth(flat);
          llhvec[0]  = flat;
          llhvec[1]  = flon;
          llhvec[2]  = altkm;

          return  llhvec
    
    rnow  =  rearth(flatgc);
    altkm =  rp - rnow;
    flat  =  gc2gd (flatgc,altkm);
          
    rrnrm =  radcur(flat);
    rn    =  rrnrm[1];
    esq = EARTH_Esq;
    for kount in range(5):
        slat  =  math.sin(dtr*flat);
        tangd =  ( z + rn*esq*slat ) / p;
        flatn =  math.atan(tangd)/dtr;
    
        dlat  =  flatn - flat;
        flat  =  flatn;
        clat  =  math.cos( dtr*flat );
    
        rrnrm =  radcur(flat);
        rn    =  rrnrm[1];
    
        altkm =  (p/clat) - rn;
    
        if (np.abs(dlat) < 1.0e-12): break

      
     
    llhvec[0]  = flat;
    llhvec[1]  = flon;
    llhvec[2]  = altkm;

    return  llhvec 
      
def rearth (lat):
    rrnrm =  radcur ( lat );
    r     =  rrnrm[0];

    return  r 

def radcur(lat):
    #compute the radius at the geodetic latitude lat (in degrees)
    #output r,rn,rm in km
    rrnrm=np.zeros(3)
    dtr =  math.pi/180.0;

    asq   = a*a;
    bsq   = b*b;
    eccsq  =  1 - bsq/asq;
    ecc = math.sqrt(eccsq);
    clat  =  math.cos(dtr*lat);
    slat  =  math.sin(dtr*lat);

    dsq   =  1.0 - eccsq * slat * slat;
    d     =  math.sqrt(dsq);

    rn    =  a/d;
    rm    =  rn * (1.0 - eccsq ) / dsq;

    rho   =  rn * clat;
    z     =  (1.0 - eccsq ) * rn * slat;
    rsq   =  rho*rho + z*z;
    r     =  math.sqrt( rsq );

    rrnrm[0]  =  r;
    rrnrm[1]  =  rn;
    rrnrm[2]  =  rm;
    
    return rrnrm

def gc2gd(flatgc, altkm):
#    geocentric latitude to geodetic latitude
#
#     Input:
#               flatgc    geocentric latitude deg.
#               altkm     altitide in km
#     ouput:
#               flatgd    geodetic latitude in deg
    dtr =  math.pi/180.0;
    rtd   = 1/dtr;
    rrnrm = np.zeros(3)
    ecc   =  EARTH_Ecc;
    esq   =  ecc*ecc;
    altnow  =  altkm;

    rrnrm   =  radcur (flatgc);
    rn      =  rrnrm[1];
     
    ratio   = 1 - esq*rn/(rn+altnow);

    tlat    = math.tan(dtr*flatgc) / ratio;
    flatgd  = rtd * math.atan(tlat);

#       now use this approximation for gd-lat to get rn etc.

    rrnrm   =  radcur ( flatgd );
    rn      =  rrnrm[1];

    ratio   =  1  - esq*rn/(rn+altnow)
    tlat    =  math.tan(dtr*flatgc)/ratio;
    flatgd  =  rtd * math.atan(tlat);

    return  flatgd