from netCDF4 import Dataset
import netCDF4
from scipy.interpolate import griddata
import math
import numpy as np

#-------------------------------------------------------------------------------

def interpolate_MPAS_variable(nTimes, nx, ny, xDEMSI, yDEMSI, xMPAS, yMPAS, arrayMPAS):

    arrayDEMSI = np.zeros((nTimes, nx, ny))

    for iTime in range(0,nTimes):

        print("    Time: ", iTime, " of ", nTimes)

        arrayDEMSI[iTime, :, :] = griddata((xMPAS, yMPAS), arrayMPAS[iTime,:], (xDEMSI.flatten(), yDEMSI.flatten()), method='linear').reshape(xDEMSI.shape)

    return arrayDEMSI

#-------------------------------------------------------------------------------

# MPAS mesh
print("MPAS mesh data...")
filenameMesh = "seaice_QU_120km.nc"

fileMesh = Dataset(filenameMesh, "r")

nCells = len(fileMesh.dimensions["nCells"])

latCell = fileMesh.variables["latCell"][:]
lonCell = fileMesh.variables["lonCell"][:]

fileMesh.close()

xMPAS = []
yMPAS = []
for iCell in range(0,nCells):
    if (latCell[iCell] > 0.0):
        xMPAS.append(math.cos(lonCell[iCell]) * math.cos(latCell[iCell]))
        yMPAS.append(math.sin(lonCell[iCell]) * math.cos(latCell[iCell]))
xMPAS = np.array(xMPAS)
yMPAS = np.array(yMPAS)

nCellsReduced = len(xMPAS)




# MPAS ocean data
print("MPAS ocean data...")
filenameIn = "oceanmixed_ice_depth_QU120km.nc"

fileIn = Dataset(filenameIn, "r")

nTimes = len(fileIn.dimensions["Time"])

seaSurfaceSalinityIn       = fileIn.variables["seaSurfaceSalinity"][:]
seaSurfaceTemperatureIn    = fileIn.variables["seaSurfaceTemperature"][:]
uOceanVelocityIn           = fileIn.variables["uOceanVelocity"][:]
vOceanVelocityIn           = fileIn.variables["vOceanVelocity"][:]
seaSurfaceTiltUIn          = fileIn.variables["seaSurfaceTiltU"][:]
seaSurfaceTiltVIn          = fileIn.variables["seaSurfaceTiltV"][:]
oceanMixedLayerDepthIn     = fileIn.variables["oceanMixedLayerDepth"][:]
oceanHeatFluxConvergenceIn = fileIn.variables["oceanHeatFluxConvergence"][:]

fileIn.close()

seaSurfaceSalinity       = np.zeros((nTimes,nCellsReduced))
seaSurfaceTemperature    = np.zeros((nTimes,nCellsReduced))
uOceanVelocity           = np.zeros((nTimes,nCellsReduced))
vOceanVelocity           = np.zeros((nTimes,nCellsReduced))
seaSurfaceTiltU          = np.zeros((nTimes,nCellsReduced))
seaSurfaceTiltV          = np.zeros((nTimes,nCellsReduced))
oceanMixedLayerDepth     = np.zeros((nTimes,nCellsReduced))
oceanHeatFluxConvergence = np.zeros((nTimes,nCellsReduced))

for iTime in range(0,nTimes):
    iCellReduced = 0
    for iCell in range(0,nCells):
        if (latCell[iCell] > 0.0):

            seaSurfaceSalinity      [iTime,iCellReduced] = seaSurfaceSalinityIn      [iTime,iCell]
            seaSurfaceTemperature   [iTime,iCellReduced] = seaSurfaceTemperatureIn   [iTime,iCell]
            uOceanVelocity          [iTime,iCellReduced] = uOceanVelocityIn          [iTime,iCell]
            vOceanVelocity          [iTime,iCellReduced] = vOceanVelocityIn          [iTime,iCell]
            seaSurfaceTiltU         [iTime,iCellReduced] = seaSurfaceTiltUIn         [iTime,iCell]
            seaSurfaceTiltV         [iTime,iCellReduced] = seaSurfaceTiltVIn         [iTime,iCell]
            oceanMixedLayerDepth    [iTime,iCellReduced] = oceanMixedLayerDepthIn    [iTime,iCell]
            oceanHeatFluxConvergence[iTime,iCellReduced] = oceanHeatFluxConvergenceIn[iTime,iCell]

            iCellReduced += 1




# DEMSI grid file
print("DEMSI grid data...")
filenameGrid = "grid.nc"

fileGrid = Dataset(filenameGrid, "r")

nxDEMSI = len(fileGrid.dimensions["nx"])
nyDEMSI = len(fileGrid.dimensions["ny"])

latitudeDEMSI  = fileGrid.variables["latitude"][:]
longitudeDEMSI = fileGrid.variables["longitude"][:]

fileGrid.close()

xDEMSI = np.zeros((nxDEMSI, nyDEMSI))
yDEMSI = np.zeros((nxDEMSI, nyDEMSI))
for i in range(0, nxDEMSI):
    for j in range(0, nyDEMSI):
        xDEMSI[i, j] = math.cos(longitudeDEMSI[i, j]) * math.cos(latitudeDEMSI[i, j])
        yDEMSI[i, j] = math.sin(longitudeDEMSI[i, j]) * math.cos(latitudeDEMSI[i, j])




# interpolation
print("Interpolation...")
print("  seaSurfaceSalinity")       ; seaSurfaceSalinityDEMSI       = interpolate_MPAS_variable(nTimes, nxDEMSI, nyDEMSI, xDEMSI, yDEMSI, xMPAS, yMPAS, seaSurfaceSalinity)
print("  seaSurfaceTemperature")    ; seaSurfaceTemperatureDEMSI    = interpolate_MPAS_variable(nTimes, nxDEMSI, nyDEMSI, xDEMSI, yDEMSI, xMPAS, yMPAS, seaSurfaceTemperature)
print("  uOceanVelocity")           ; uOceanVelocityDEMSI           = interpolate_MPAS_variable(nTimes, nxDEMSI, nyDEMSI, xDEMSI, yDEMSI, xMPAS, yMPAS, uOceanVelocity)
print("  vOceanVelocity")           ; vOceanVelocityDEMSI           = interpolate_MPAS_variable(nTimes, nxDEMSI, nyDEMSI, xDEMSI, yDEMSI, xMPAS, yMPAS, vOceanVelocity)
print("  seaSurfaceTiltU")          ; seaSurfaceTiltUDEMSI          = interpolate_MPAS_variable(nTimes, nxDEMSI, nyDEMSI, xDEMSI, yDEMSI, xMPAS, yMPAS, seaSurfaceTiltU)
print("  seaSurfaceTiltV")          ; seaSurfaceTiltVDEMSI          = interpolate_MPAS_variable(nTimes, nxDEMSI, nyDEMSI, xDEMSI, yDEMSI, xMPAS, yMPAS, seaSurfaceTiltV)
print("  oceanMixedLayerDepth")     ; oceanMixedLayerDepthDEMSI     = interpolate_MPAS_variable(nTimes, nxDEMSI, nyDEMSI, xDEMSI, yDEMSI, xMPAS, yMPAS, oceanMixedLayerDepth)
print("  oceanHeatFluxConvergence") ; oceanHeatFluxConvergenceDEMSI = interpolate_MPAS_variable(nTimes, nxDEMSI, nyDEMSI, xDEMSI, yDEMSI, xMPAS, yMPAS, oceanHeatFluxConvergence)




# create output file
print("Output file...")
filenameOut = "forcing_CORE_monthly_clim_ocean.0000.nc"

fileOut = Dataset(filenameOut,"w",format="NETCDF3_CLASSIC")

fileOut.climatology = 1;

fileOut.createDimension("strLen", 64)
fileOut.createDimension("Time", nTimes)
fileOut.createDimension("nx", nxDEMSI)
fileOut.createDimension("ny", nyDEMSI)

# times
varOutTime = fileOut.createVariable("Time","c",dimensions=["Time","strLen"])
year = 0
day = 15
hour = 0
minute = 0
second = 0
for iTime in range(0,nTimes):
    month = iTime + 1
    timeStr = "%4.4i-%2.2i-%2.2i_%2.2i:%2.2i:%2.2i" %(year,month,day,hour,minute,second)
    varOutTime[iTime,0:19] = netCDF4.stringtochar(np.array(timeStr, 'S19'))

var = fileOut.createVariable("seaSurfaceSalinity", "d", dimensions=["Time","nx","ny"])
var[:] = seaSurfaceSalinityDEMSI[:]

var = fileOut.createVariable("seaSurfaceTemperature", "d", dimensions=["Time","nx","ny"])
var[:] = seaSurfaceTemperatureDEMSI[:]

var = fileOut.createVariable("xOcnCurrents", "d", dimensions=["Time","nx","ny"])
var[:] = uOceanVelocityDEMSI[:]

var = fileOut.createVariable("yOcnCurrents", "d", dimensions=["Time","nx","ny"])
var[:] = vOceanVelocityDEMSI[:]

var = fileOut.createVariable("seaSurfaceTiltU", "d", dimensions=["Time","nx","ny"])
var[:] = seaSurfaceTiltUDEMSI[:]

var = fileOut.createVariable("seaSurfaceTiltV", "d", dimensions=["Time","nx","ny"])
var[:] = seaSurfaceTiltVDEMSI[:]

var = fileOut.createVariable("oceanMixedLayerDepth", "d", dimensions=["Time","nx","ny"])
var[:] = oceanMixedLayerDepthDEMSI[:]

var = fileOut.createVariable("oceanHeatFluxConvergence", "d", dimensions=["Time","nx","ny"])
var[:] = oceanHeatFluxConvergenceDEMSI[:]

fileOut.close()
