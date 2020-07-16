# DEMSI_XY2GEODETIC - Convert X and Y coordinates on stereographic grid to lat-long
#
# [lat,lon,SLAT]=demsi_xy2geodetic(X,Y,SGN)
#
# This function converts x and y coordinates on a polar stereographic grid
# to latitudes and longitudes, and is essentially a JPL algorithm rewritten
# to work in Python.
#
# Input:
# X   - polar stereographic X coordinate (m)
# Y   - polar stereographic Y coordinate (m)
# SGN - Hemisphere specifier: 1=Northern Hemisphere, -1=Southern Hemisphere
#
# Where (X,Y)=(0,0) is at the pole.
#
# Output:
# lat   - Latitude (radians)
# lon   - Longitude (radians)
# SLAT  - Latitude of true distance
#
# -------------------------------------------------------------------------
# This converts from Polar Stereographic (X,Y) coordinates
# to geodetic latitude and longitude for the polar regions. The equations
# are from Snyder, J. P., 1982,  Map Projections Used by the U.S.
# Geological Survey, Geological Survey Bulletin 1532, U.S. Government
# Printing Office.  See JPL Technical Memorandum 3349-8401 for details.
# CODE CONVERTED FROM A JPL SCRIPT BY C. S. Morris and  V. J. Troisi
# -------------------------------------------------------------------------
#
# See also: demsi_geodetic2xy
#
# Code by Andrew Roberts and Travis Davis December 6 2017
# Department of Oceanography, Naval Postgraduate School

import numpy as np
import scipy

def demsi_xy2geodetic(X, Y, SGN):

    # set defaults
    E2 = 0.006693883
    E = np.sqrt(E2)     # eccentricity of Hughes ellipsoid
    RE = 6378273.0      # mean earth radius in m
    SLAT = 70.          # latitude of true distance
    PI = np.pi          # PI
    SL  = SLAT*PI/180.
    RHO = np.sqrt(np.power(X,2.)+np.power(Y,2.))
    CM  = np.divide(np.cos(SL),np.sqrt(1.-E2*np.power(np.sin(SL),2.)))
    T   = np.divide(np.tan(PI/4.-SL/2.),np.power(np.divide(1.-E*np.sin(SL), \
          1.+E*np.sin(SL)),E/2.))

    # special case of true latitude being at the pole
    if np.abs((SLAT-90.))<1.E-5:
        T = np.dot(RHO,np.sqrt(np.dot(np.power(1.+E,1.+E), \
            np.power(1.-E,1-E))))/(2*RE)
    else:
        T = np.divide(np.dot(RHO,T),RE*CM)

    CHI = PI/2.-2.*np.arctan(T)
    ALAT  = CHI+(E2/2.+5.0/24.*np.power(E2,2.)+np.power(E2,3.)/12.)* \
            np.sin(2.*CHI)+(7.0/48.*np.power(E2,2.)+ \
            29.0/240.*np.power(E2,3.))*np.sin(4.*CHI)+ \
            7.0/120.*np.power(E2,3.)*np.sin(6.*CHI)
    ALAT  = SGN*ALAT
    ALONG = np.arctan2(SGN*X,-SGN*Y)
    ALONG = SGN*ALONG

    # special case of being incredibly close to the pole
    ALAT  = np.where(RHO < 100.0, PI*SGN/2, ALAT)
    ALONG = np.where(RHO < 100.0, 0.0, ALONG)

    #-------------------------------------------------------------------------*
    lon = ALONG
    lat = ALAT

    return (lat, lon)
