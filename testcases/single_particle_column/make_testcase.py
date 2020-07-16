import sys
from netCDF4 import Dataset
import netCDF4
import numpy as np
import os
import string
import math
sys.path.append("../../utils/testcases/")
from testcase_utils import write_out_grid, write_particle_file, add_init_connectivity

#-------------------------------------------------------------------------------

def convert_time(timeIn):

    days_in_month = [31,28,31,30,31,30,31,31,30,31,30,31]

    date = timeIn.split("_")[0]
    time = timeIn.split("_")[1]

    year  = int(date.split("-")[0])
    month = int(date.split("-")[1])
    day   = int(date.split("-")[2])

    hour   = int(time.split(":")[0])
    minute = int(time.split(":")[1])
    second = int(time.split(":")[2])

    if (hour == 24):
        hour = 0
        day = day + 1
        if (day > days_in_month[month-1]):
            day = 1
            month = month + 1
            if (month > 12):
                month = 1
                year = year + 1

    timeOut = "%4.4i-%2.2i-%2.2i_%2.2i:%2.2i:%2.2i" %(year,month,day,hour,minute,second)
    return timeOut

#-------------------------------------------------------------------------------

# MPAS-Seaice one dimensional test case

# data input location
dataDir = os.environ['DEMSI_DATA_DIR']

# grid
latTestCase = 71.35
lonTestCase = -156.5

degreesToRadians = math.pi / 180.0

nx=2
ny=2

lx=1000
ly=1000

dx, dy = lx/(nx-1), ly/(ny-1)

lat = np.zeros((nx,ny))
lon = np.zeros((nx,ny))

lat[:] = latTestCase * degreesToRadians
lon[:] = lonTestCase * degreesToRadians

write_out_grid("grid.nc", nx, ny, dx, dy, lat, lon)

# particles
nParticles = 1

nTypes = 1

x = np.zeros(nParticles)
y = np.zeros(nParticles)
r = np.zeros(nParticles)
d = np.zeros(nParticles, dtype="i")
t = np.ones(nParticles, dtype="i")
f = np.ones(nParticles)
h = np.ones(nParticles)

radius = 500.0

x[0] = radius
y[0] = radius

r[0] = radius

d[0] = 1


maxVertices = 4

nv = np.zeros((nParticles), dtype="i")
nv[:] = 4
xv = np.zeros((nParticles,maxVertices), dtype="i")
yv = np.zeros((nParticles,maxVertices), dtype="i")

xv[0][:] = [0.0, 2.0*radius, 2.0*radius, 0.0       ]
yv[0][:] = [0.0, 0.0       , 2.0*radius, 2.0*radius]

filenameOut = "particles_in.nc"
write_particle_file(filenameOut, nParticles, d, nTypes, t, x, y, r, f, h, maxVertices, nv, xv, yv)

add_init_connectivity(filenameOut)


# forcing
years = [1989,1990,1991]

nTimesTotal = 0
for year in years:

    # file in
    filenameSixHourlyIn = dataDir + "/mpas-seaice/sc_71.35_-156.5/LYq_six_hourly.%4.4i.nc" %(year)
    fileIn = Dataset(filenameSixHourlyIn,"r")

    nTimes = len(fileIn.dimensions["Time"])
    nTimesTotal = nTimesTotal + nTimes

    xtime = fileIn.variables["xtime"][:]

    airTemperature = fileIn.variables["airTemperature"][:]
    airSpecificHumidity = fileIn.variables["airSpecificHumidity"][:]
    uAirVelocity = fileIn.variables["uAirVelocity"][:]
    vAirVelocity = fileIn.variables["vAirVelocity"][:]

    if (year == years[0]):
        xtimeTotal               = xtime.copy()
        airTemperatureTotal      = airTemperature.copy()
        airSpecificHumidityTotal = airSpecificHumidity.copy()
        uAirVelocityTotal        = uAirVelocity.copy()
        vAirVelocityTotal        = vAirVelocity.copy()
    else:
        xtimeTotal               = np.concatenate((xtimeTotal,xtime))
        airTemperatureTotal      = np.concatenate((airTemperatureTotal,airTemperature))
        airSpecificHumidityTotal = np.concatenate((airSpecificHumidityTotal,airSpecificHumidity))
        uAirVelocityTotal        = np.concatenate((uAirVelocityTotal,uAirVelocity))
        vAirVelocityTotal        = np.concatenate((vAirVelocityTotal,vAirVelocity))

    fileIn.close()

# file out
filenameSixHourlyOut = "forcing_CORE_six_hourly.nc"
fileOut = Dataset(filenameSixHourlyOut,"w",format="NETCDF3_CLASSIC")

fileOut.climatology = 0;

strLen = 64
nx = 2
ny = 2

fileOut.createDimension("strLen", strLen)
fileOut.createDimension("Time", nTimesTotal)
fileOut.createDimension("nx", nx)
fileOut.createDimension("ny", ny)

varOut = fileOut.createVariable("Time","c",dimensions=["Time","strLen"])
for iTime in range(0,nTimesTotal):
    timeStr = convert_time(''.join(map(bytes.decode, xtimeTotal[iTime,0:19])))
    varOut[iTime,0:19] = netCDF4.stringtochar(np.array(timeStr, 'S19'))

varOut = fileOut.createVariable("airTemperature","d",dimensions=["Time","nx","ny"])
for iTime in range(0,nTimesTotal):
    for ix in range(0,nx):
        for iy in range(0,ny):
            varOut[iTime,ix,iy] = airTemperatureTotal[iTime,0]

varOut = fileOut.createVariable("specificHumidity","d",dimensions=["Time","nx","ny"])
for iTime in range(0,nTimesTotal):
    for ix in range(0,nx):
        for iy in range(0,ny):
            varOut[iTime,ix,iy] = airSpecificHumidityTotal[iTime,0]

varOut = fileOut.createVariable("xAtmWind","d",dimensions=["Time","nx","ny"])
for iTime in range(0,nTimesTotal):
    for ix in range(0,nx):
        for iy in range(0,ny):
            varOut[iTime,ix,iy] = uAirVelocityTotal[iTime,0]

varOut = fileOut.createVariable("yAtmWind","d",dimensions=["Time","nx","ny"])
for iTime in range(0,nTimesTotal):
    for ix in range(0,nx):
        for iy in range(0,ny):
            varOut[iTime,ix,iy] = vAirVelocityTotal[iTime,0]

fileOut.close()


# atmos monthly climatology
filenameMonthlyClimIn = dataDir + "/mpas-seaice/sc_71.35_-156.5/LYq_monthly.nc"
fileIn = Dataset(filenameMonthlyClimIn,"r")

nTimes = len(fileIn.dimensions["Time"])

xtime = fileIn.variables["xtime"][:]

cloudFraction = fileIn.variables["cloudFraction"][:]
rainfallRate  = fileIn.variables["rainfallRate"][:]

fileIn.close()

filenameMonthlyClimOut = "forcing_CORE_monthly_clim_atmos.0000.nc"
fileOut = Dataset(filenameMonthlyClimOut,"w",format="NETCDF3_CLASSIC")

fileOut.climatology = 1;

strLen = 64
nx = 2
ny = 2

fileOut.createDimension("strLen", strLen)
fileOut.createDimension("Time", nTimes)
fileOut.createDimension("nx", nx)
fileOut.createDimension("ny", ny)

varOut = fileOut.createVariable("Time","c",dimensions=["Time","strLen"])
for iTime in range(0,nTimes):
    timeStr = ''.join(map(bytes.decode, xtime[iTime,0:19]))
    varOut[iTime,0:19] = netCDF4.stringtochar(np.array(timeStr, 'S19'))

varOut = fileOut.createVariable("cloudFraction","d",dimensions=["Time","nx","ny"])
for iTime in range(0,nTimes):
    for ix in range(0,nx):
        for iy in range(0,ny):
            varOut[iTime,ix,iy] = cloudFraction[iTime,0]

varOut = fileOut.createVariable("precipitationRate","d",dimensions=["Time","nx","ny"])
for iTime in range(0,nTimes):
    for ix in range(0,nx):
        for iy in range(0,ny):
            varOut[iTime,ix,iy] = rainfallRate[iTime,0]

fileOut.close()

# ocean monthly climatology
filenameMonthlyClimIn = dataDir + "/mpas-seaice/sc_71.35_-156.5/oceanmixed_ice_depth_sc_71.35_-156.5.nc"
fileIn = Dataset(filenameMonthlyClimIn,"r")

nTimes = len(fileIn.dimensions["Time"])

xtime = fileIn.variables["xtime"][:]

seaSurfaceSalinity       = fileIn.variables["seaSurfaceSalinity"][:]
seaSurfaceTemperature    = fileIn.variables["seaSurfaceTemperature"][:]
uOceanVelocity           = fileIn.variables["uOceanVelocity"][:]
vOceanVelocity           = fileIn.variables["vOceanVelocity"][:]
seaSurfaceTiltU          = fileIn.variables["seaSurfaceTiltU"][:]
seaSurfaceTiltV          = fileIn.variables["seaSurfaceTiltV"][:]
oceanMixedLayerDepth     = fileIn.variables["oceanMixedLayerDepth"][:]
oceanHeatFluxConvergence = fileIn.variables["oceanHeatFluxConvergence"][:]

fileIn.close()

filenameMonthlyClimOut = "forcing_CORE_monthly_clim_ocean.0000.nc"
fileOut = Dataset(filenameMonthlyClimOut,"w",format="NETCDF3_CLASSIC")

fileOut.climatology = 1;

strLen = 64
nx = 2
ny = 2

fileOut.createDimension("strLen", strLen)
fileOut.createDimension("Time", nTimes)
fileOut.createDimension("nx", nx)
fileOut.createDimension("ny", ny)

varOut = fileOut.createVariable("Time","c",dimensions=["Time","strLen"])
for iTime in range(0,nTimes):
    timeStr = ''.join(map(bytes.decode, xtime[iTime,0:19]))
    varOut[iTime,0:19] = netCDF4.stringtochar(np.array(timeStr, 'S19'))

varOut = fileOut.createVariable("seaSurfaceSalinity","d",dimensions=["Time","nx","ny"])
for iTime in range(0,nTimes):
    for ix in range(0,nx):
        for iy in range(0,ny):
            varOut[iTime,ix,iy] = seaSurfaceSalinity[iTime,0]

varOut = fileOut.createVariable("seaSurfaceTemperature","d",dimensions=["Time","nx","ny"])
for iTime in range(0,nTimes):
    for ix in range(0,nx):
        for iy in range(0,ny):
            varOut[iTime,ix,iy] = seaSurfaceTemperature[iTime,0]

varOut = fileOut.createVariable("xOcnCurrents","d",dimensions=["Time","nx","ny"])
for iTime in range(0,nTimes):
    for ix in range(0,nx):
        for iy in range(0,ny):
            varOut[iTime,ix,iy] = uOceanVelocity[iTime,0]

varOut = fileOut.createVariable("yOcnCurrents","d",dimensions=["Time","nx","ny"])
for iTime in range(0,nTimes):
    for ix in range(0,nx):
        for iy in range(0,ny):
            varOut[iTime,ix,iy] = vOceanVelocity[iTime,0]

varOut = fileOut.createVariable("seaSurfaceTiltU","d",dimensions=["Time","nx","ny"])
for iTime in range(0,nTimes):
    for ix in range(0,nx):
        for iy in range(0,ny):
            varOut[iTime,ix,iy] = seaSurfaceTiltU[iTime,0]

varOut = fileOut.createVariable("seaSurfaceTiltV","d",dimensions=["Time","nx","ny"])
for iTime in range(0,nTimes):
    for ix in range(0,nx):
        for iy in range(0,ny):
            varOut[iTime,ix,iy] = seaSurfaceTiltV[iTime,0]

varOut = fileOut.createVariable("oceanMixedLayerDepth","d",dimensions=["Time","nx","ny"])
for iTime in range(0,nTimes):
    for ix in range(0,nx):
        for iy in range(0,ny):
            varOut[iTime,ix,iy] = oceanMixedLayerDepth[iTime,0]

varOut = fileOut.createVariable("oceanHeatFluxConvergence","d",dimensions=["Time","nx","ny"])
for iTime in range(0,nTimes):
    for ix in range(0,nx):
        for iy in range(0,ny):
            varOut[iTime,ix,iy] = oceanHeatFluxConvergence[iTime,0]

fileOut.close()
