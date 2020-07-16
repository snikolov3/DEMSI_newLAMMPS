# DEMSI_ATEST_STEREOGRAPHIC - Test demsi_xy2geodetic and demsi_geodetic2xy
#
# This script runs a quick test on the two DEMSI stereographic projection
# routines, demsi_xy2geodetic and demsi_geodetic2xy.
#
# Code by Travis Davis and Andrew Roberts December 6 2017
# Department of Oceanography, Naval Postgraduate School

import numpy as np
import scipy

from demsi_geodetic2xy import demsi_geodetic2xy
from demsi_xy2geodetic import demsi_xy2geodetic

# Test Case #1
lat1 = 90
lon1 = 20
SGN1 = 1

# Test Case #2
lat2 = 70.
lon2 = 45.
SGN2 = 1

# Test Case #3
lat3 = -60.
lon3 = 45.
SGN3 = -1

# Test Evaluation
(X1,Y1) = demsi_geodetic2xy(lat1,lon1,SGN1)
(X2,Y2) = demsi_geodetic2xy(lat2,lon2,SGN2)
(X3,Y3) = demsi_geodetic2xy(lat3,lon3,SGN3)

(lat1r,lon1r)    = demsi_xy2geodetic(X1,Y1,SGN1)
(lat2r,lon2r)    = demsi_xy2geodetic(X2,Y2,SGN2)
(lat3r,lon3r)    = demsi_xy2geodetic(X3,Y3,SGN3)

print("-------------")
print(X1,Y1)
print(X2,Y2)
print(X3,Y3)
print("-------------")
print(lat1r,lat1,lon1r,lon1r)
print(lat2r,lat2,lon2r,lon2r)
print(lat3r,lat3,lon3r,lon3r)
print("-------------")
