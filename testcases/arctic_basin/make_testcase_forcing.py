import numpy as np
from netCDF4 import Dataset
import netCDF4
import sys, math, glob, os

sys.path.append("../../utils/stereoprojection/")
from demsi_xy2geodetic import demsi_xy2geodetic
from demsi_geodetic2xy import demsi_geodetic2xy

from scipy.interpolate import RegularGridInterpolator

import matplotlib.pyplot as plt

#-------------------------------------------------------------------------------

def convert_CORE_latslons_to_vertices(latsIn, lonsIn):

    nLatsIn = len(latsIn)
    nLonsIn = len(lonsIn)

    latsOut = np.zeros(nLatsIn+1)
    lonsOut = np.zeros(nLonsIn+1)

    for iLat in range(1,nLatsIn):
        latsOut[iLat] = 0.5 * (latsIn[iLat-1] + latsIn[iLat])

    latsOut[0]  = -90.0
    latsOut[-1] =  90.0

    lonsOut[0:nLonsIn] = lonsIn[0:nLonsIn]
    lonsOut[-1] = 360.0

    degreesToRadian = math.pi / 180.0

    latsOut[:] = latsOut[:] * degreesToRadian
    lonsOut[:] = lonsOut[:] * degreesToRadian

    for iLon in range(0,nLonsIn+1):
        if (lonsOut[iLon] < 0.0):
            lonsOut[iLon] = lonsOut[iLon] + 2.0 * math.pi

    return latsOut, lonsOut

#-------------------------------------------------------------------------------

def convert_CORE_data_to_vertices(arrayIn):

    nLatsIn = arrayIn.shape[0]
    nLonsIn = arrayIn.shape[1]

    arrayOut = np.zeros((nLatsIn+1,nLonsIn+1))

    for iLat in range(1,nLatsIn):
        arrayOut[iLat,:-1] = 0.5 * np.add(arrayIn[iLat-1,:],arrayIn[iLat,:])

    arrayOut[0,:] = np.mean(arrayIn[1,:])
    arrayOut[-1,:] = np.mean(arrayIn[-2,:])

    arrayOut[:,-1] = arrayOut[:,0]

    return arrayOut

#-------------------------------------------------------------------------------

def interpolate_CORE_scalar_field(filenameExpressionCORE, varname, varOut):

    filenameCORE = glob.glob(filenameExpressionCORE)[0]
    fileCORE = Dataset(filenameCORE,"r")

    nTimes = len(fileCORE.dimensions["TIME"])

    prefix = "%s:" %(varname)

    for iTime in range(0,nTimes):

        arrayCOREIn = fileCORE.variables[varname][iTime,:,:]

        arrayCORE = convert_CORE_data_to_vertices(arrayCOREIn)

        interp = RegularGridInterpolator((latCORE, lonCORE), arrayCORE)

        varOut[iTime,:,:] = interp(latlonDEMSI).reshape((nx,ny))

    fileCORE.close()

#-------------------------------------------------------------------------------

def interpolate_CORE_vector_field(filenameCOREU, filenameCOREV, varnameU, varnameV, varOutU, varOutV):

    fileCOREU = Dataset(filenameCOREU,"r")
    fileCOREV = Dataset(filenameCOREV,"r")

    nTimes = len(fileCOREU.dimensions["TIME"])

    prefix = "%s, %s:" %(varnameU, varnameV)

    for iTime in range(0,nTimes):

        arrayCOREUIn = fileCOREU.variables[varnameU][iTime,:,:]
        arrayCOREVIn = fileCOREV.variables[varnameV][iTime,:,:]

        arrayCOREURot = np.add(np.multiply(arrayCOREUIn,alpha[:,:]),np.multiply(arrayCOREVIn,gamma[:,:]))
        arrayCOREVRot = np.add(np.multiply(arrayCOREUIn,beta [:,:]),np.multiply(arrayCOREVIn,delta[:,:]))

        arrayCOREU = convert_CORE_data_to_vertices(arrayCOREURot)
        arrayCOREV = convert_CORE_data_to_vertices(arrayCOREVRot)

        interpU = RegularGridInterpolator((latCORE, lonCORE), arrayCOREU)
        interpV = RegularGridInterpolator((latCORE, lonCORE), arrayCOREV)

        varOutU[iTime,:,:] = interpU(latlonDEMSI).reshape((nx,ny))
        varOutV[iTime,:,:] = interpV(latlonDEMSI).reshape((nx,ny))

    fileCOREU.close()
    fileCOREV.close()

#-------------------------------------------------------------------------------

def increment_time(month,day,hour,minute,second):

    days_in_month = [31,28,31,30,31,30,31,31,30,31,30,31]

    hour = hour + 6
    if (hour >= 24):
        hour = hour - 24
        day = day + 1
        if (day > days_in_month[month-1]):
            day = 1
            month = month + 1

    return month,day,hour,minute,second

#-------------------------------------------------------------------------------

# data input location
dataDir = os.environ['DEMSI_DATA_DIR']

# read in grid
print("Read in grid")
filenameGrid = "grid.nc"

fileGrid = Dataset(filenameGrid,"r")

nx = len(fileGrid.dimensions["nx"])
ny = len(fileGrid.dimensions["ny"])

resolution = fileGrid.variables["resolution"][:]
dx = resolution[0]
dy = resolution[1]

latDEMSI = fileGrid.variables["latitude"][:]
lonDEMSI = fileGrid.variables["longitude"][:]
lonDEMSI = np.where(lonDEMSI < 0.0, lonDEMSI+2.0*math.pi, lonDEMSI)

fileGrid.close()

Lx = dx * float(nx-1)
Ly = dy * float(ny-1)




# forcing
yearStart = 1999
yearEnd = 2001


# flattened DEMSI positions
latDEMSIFlatten = latDEMSI.flatten()
lonDEMSIFlatten = lonDEMSI.flatten()
latlonDEMSI = np.zeros((nx*ny,2))
for ij in range(0,nx*ny):
    latlonDEMSI[ij,0] = latDEMSIFlatten[ij]
    latlonDEMSI[ij,1] = lonDEMSIFlatten[ij]




# core lat lon
print("Calculate interpolation weights")
filenameCORE = dataDir + "/CORE-II/u_10.2000.nc"
fileCORE = Dataset(filenameCORE,"r")

nLatsCORE = len(fileCORE.dimensions["LAT"])
nLonsCORE = len(fileCORE.dimensions["LON"])

latCOREIn = fileCORE.variables["LAT"][:]
lonCOREIn = fileCORE.variables["LON"][:]
latCORE, lonCORE = convert_CORE_latslons_to_vertices(latCOREIn, lonCOREIn)

nTimes = len(fileCORE.dimensions["TIME"])

fileCORE.close()


# DEMSI basis vectors in CORE space
print("DEMSI basis in CORE space")
alpha = np.zeros((nLatsCORE,nLonsCORE))
beta  = np.zeros((nLatsCORE,nLonsCORE))
gamma = np.zeros((nLatsCORE,nLonsCORE))
delta = np.zeros((nLatsCORE,nLonsCORE))

theta = -25.0
thetaRadians = np.radians(theta)
c, s = np.cos(thetaRadians), np.sin(thetaRadians)

for iLat in range(0,nLatsCORE):
    for iLon in range(0,nLonsCORE):

        x0,    y0    = demsi_geodetic2xy(latCORE[iLat],          lonCORE[iLon],          1)
        uhatx, uhaty = demsi_geodetic2xy(latCORE[iLat],          lonCORE[iLon] + 1.0e-8, 1)
        vhatx, vhaty = demsi_geodetic2xy(latCORE[iLat] + 1.0e-8, lonCORE[iLon],          1)

        x0r =  c * x0 + s * y0
        y0r = -s * x0 + c * y0

        uhatxr =  c * uhatx + s * uhaty
        uhatyr = -s * uhatx + c * uhaty

        vhatxr =  c * vhatx + s * vhaty
        vhatyr = -s * vhatx + c * vhaty

        vmag = math.sqrt(math.pow(uhatxr - x0r,2) + math.pow(uhatyr - y0r,2))
        uhatxr = (uhatxr - x0r) / vmag
        uhatyr = (uhatyr - y0r) / vmag

        vmag = math.sqrt(math.pow(vhatxr - x0r,2) + math.pow(vhatyr - y0r,2))
        vhatxr = (vhatxr - x0r) / vmag
        vhatyr = (vhatyr - y0r) / vmag

        alpha[iLat,iLon] = uhatxr
        beta [iLat,iLon] = uhatyr
        gamma[iLat,iLon] = vhatxr
        delta[iLat,iLon] = vhatyr




print("Perform interpolation")
for year in range(yearStart, yearEnd+1):

    print("Year:", year)

    # file out
    filenameSixHourlyOut = "forcing_CORE_six_hourly.%4.4i.nc" %(year)
    fileOut = Dataset(filenameSixHourlyOut,"w",format="NETCDF3_CLASSIC")

    fileOut.climatology = 0;

    strLen = 64

    fileOut.createDimension("strLen", strLen)
    fileOut.createDimension("Time", nTimes)
    fileOut.createDimension("nx", nx)
    fileOut.createDimension("ny", ny)

    # times
    varOutTime = fileOut.createVariable("Time","c",dimensions=["Time","strLen"])
    month = 1
    day = 1
    hour = 3
    minute = 0
    second = 0
    for iTime in range(0,nTimes):
        timeStr = "%4.4i-%2.2i-%2.2i_%2.2i:%2.2i:%2.2i" %(year,month,day,hour,minute,second)
        varOutTime[iTime,0:19] = netCDF4.stringtochar(np.array(timeStr, 'S19'))
        month,day,hour,minute,second = increment_time(month,day,hour,minute,second)

    # define fields
    varOutU = fileOut.createVariable("xAtmWind","d",dimensions=["Time","nx","ny"])
    varOutV = fileOut.createVariable("yAtmWind","d",dimensions=["Time","nx","ny"])
    varOutT = fileOut.createVariable("airTemperature","d",dimensions=["Time","nx","ny"])
    varOutQ = fileOut.createVariable("specificHumidity","d",dimensions=["Time","nx","ny"])

    # uAirVelocity, vAirVelocity
    filenameCOREU = dataDir + "/CORE-II/u_10.%4.4i.nc" %(year)
    filenameCOREV = dataDir + "/CORE-II/v_10.%4.4i.nc" %(year)
    interpolate_CORE_vector_field(filenameCOREU, filenameCOREV, "U_10_MOD", "V_10_MOD", varOutU, varOutV)
    #varOutU[:] = 100.0
    #varOutV[:] = 0.0

    # airTemperature, specificHumidity
    filenameCORET = dataDir + "/CORE-II/t_10.%4.4i.nc" %(year)
    filenameCOREQ = dataDir + "/CORE-II/q_10.%4.4i.nc" %(year)
    interpolate_CORE_scalar_field(filenameCORET, "T_10_MOD", varOutT)
    interpolate_CORE_scalar_field(filenameCOREQ, "Q_10_MOD", varOutQ)

    # close output file
    fileOut.close()




# atmosphere monthly file out
filenameMonthlyOut = "forcing_CORE_monthly_clim_atmos.0000.nc"
fileOut = Dataset(filenameMonthlyOut,"w",format="NETCDF3_CLASSIC")

fileOut.climatology = 1;

strLen = 64

nTimes = 12

fileOut.createDimension("strLen", strLen)
fileOut.createDimension("Time", nTimes)
fileOut.createDimension("nx", nx)
fileOut.createDimension("ny", ny)

# times
varOut = fileOut.createVariable("Time","c",dimensions=["Time","strLen"])
year = 0
day = 15
hour = 0
minute = 0
second = 0
for iTime in range(0,nTimes):
    month = iTime + 1
    timeStr = "%4.4i-%2.2i-%2.2i_%2.2i:%2.2i:%2.2i" %(year,month,day,hour,minute,second)
    varOutTime[iTime,0:19] = netCDF4.stringtochar(np.array(timeStr, 'S19'))

# define fields
varOutC = fileOut.createVariable("cloudFraction","d",dimensions=["Time","nx","ny"])
varOutP = fileOut.createVariable("precipitationRate","d",dimensions=["Time","nx","ny"])

# cloudiness
filenameCOREC = dataDir + "/CORE-II/cldf.omip.nc"
interpolate_CORE_scalar_field(filenameCOREC, "cldf", varOutC)

# precipitation
filenameCOREP = dataDir + "/CORE-II/prec.nmyr.nc"
interpolate_CORE_scalar_field(filenameCOREP, "prec", varOutP)

# close output file
fileOut.close()




# ocean monthly file out
filenameMonthlyOut = "forcing_CORE_monthly_clim_ocean.0000.nc"
fileOut = Dataset(filenameMonthlyOut,"w",format="NETCDF3_CLASSIC")

fileOut.climatology = 1;

strLen = 64

nTimes = 12

fileOut.createDimension("strLen", strLen)
fileOut.createDimension("Time", nTimes)
fileOut.createDimension("nx", nx)
fileOut.createDimension("ny", ny)

# times
varOut = fileOut.createVariable("Time","c",dimensions=["Time","strLen"])
year = 0
day = 15
hour = 0
minute = 0
second = 0
for iTime in range(0,nTimes):
    month = iTime + 1
    timeStr = "%4.4i-%2.2i-%2.2i_%2.2i:%2.2i:%2.2i" %(year,month,day,hour,minute,second)
    varOutTime[iTime,0:19] = netCDF4.stringtochar(np.array(timeStr, 'S19'))

# define fields
varOut = fileOut.createVariable("seaSurfaceSalinity","d",dimensions=["Time","nx","ny"])
varOut[:] = 34.0

varOut = fileOut.createVariable("seaSurfaceTemperature","d",dimensions=["Time","nx","ny"])
varOut[:] = -1.836

varOut = fileOut.createVariable("xOcnCurrents","d",dimensions=["Time","nx","ny"])
varOut[:] = 0.0

varOut = fileOut.createVariable("yOcnCurrents","d",dimensions=["Time","nx","ny"])
varOut[:] = 0.0

varOut = fileOut.createVariable("seaSurfaceTiltU","d",dimensions=["Time","nx","ny"])
varOut[:] = 0.0

varOut = fileOut.createVariable("seaSurfaceTiltV","d",dimensions=["Time","nx","ny"])
varOut[:] = 0.0

varOut = fileOut.createVariable("oceanMixedLayerDepth","d",dimensions=["Time","nx","ny"])
varOut[:] = 50.0

varOut = fileOut.createVariable("oceanHeatFluxConvergence","d",dimensions=["Time","nx","ny"])
varOut[:] = 0.0

# close output file
fileOut.close()
