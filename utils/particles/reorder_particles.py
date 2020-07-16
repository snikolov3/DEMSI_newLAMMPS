from netCDF4 import Dataset
import numpy as np
import argparse
import sys

#-------------------------------------------------------------------------------

def reorder_variable(fileIn, fileOut, globalID, varname):

    if (len(fileIn.variables[varname].dimensions) == 1):
        reorder_variable_1D(fileIn, fileOut, globalID, varname)
    elif (len(fileIn.variables[varname].dimensions) == 2):
        reorder_variable_2D(fileIn, fileOut, globalID, varname)

#-------------------------------------------------------------------------------

def reorder_variable_1D(fileIn, fileOut, globalID, varname):

    varIn = fileIn.variables[varname][:]

    varSorted = [var for _,var in sorted(zip(globalID,varIn))]

    varOut = fileOut.createVariable(varname,"d",dimensions=["nParticles"])
    varOut[:] = varSorted[:]

#-------------------------------------------------------------------------------

def reorder_variable_2D(fileIn, fileOut, globalID, varname):

    varIn = fileIn.variables[varname][:]

    nParticles = len(fileIn.dimensions["nParticles"])
    dimName = fileIn.variables[varname].dimensions[1]
    n = len(fileIn.dimensions[dimName])

    varSorted  = np.zeros((nParticles,n))

    for i in range(0,n):
        varSorted[:,i] = [var for _,var in sorted(zip(globalID,varIn[:,i]))]

    varOut = fileOut.createVariable(varname,"d",dimensions=["nParticles","nCategories"])
    varOut[:] = varSorted[:]

#-------------------------------------------------------------------------------

def reorder_particle_files(filenameIn, filenameOut):

    # open file in
    fileIn = Dataset(filenameIn,"r")

    globalID = fileIn.variables["globalID"][:]
    globalIDSorted = sorted(globalID)

    strLen = len(fileIn.dimensions["strLen"])
    nParticles = len(fileIn.dimensions["nParticles"])
    nCategories = len(fileIn.dimensions["nCategories"])

    Time = fileIn.variables["Time"][:]

    # open file out
    fileOut = Dataset(filenameOut,"w",format="NETCDF3_CLASSIC")

    fileOut.createDimension("strLen",strLen)
    fileOut.createDimension("nParticles",nParticles)
    fileOut.createDimension("nCategories",nCategories)

    TimeOut = fileOut.createVariable("Time","c",dimensions=["strLen"])
    TimeOut[:] = Time[:]

    globalIDOut = fileOut.createVariable("globalID","i",dimensions=["nParticles"])
    globalIDOut[:] = globalIDSorted[:]

    # iterate over variables
    for varname in fileIn.variables:
        if (varname != "globalID" and
            "nParticles" in fileIn.variables[varname].dimensions):
            reorder_variable(fileIn, fileOut, globalID, varname)

    # close files
    fileIn.close()
    fileOut.close()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Reorder a particle file by globalID.')

    parser.add_argument('-i', dest='filenameIn', help='Input file to reorder.')
    parser.add_argument('-o', dest='filenameOut', help='Output reordered file.')

    args = parser.parse_args()

    if (args.filenameIn.strip() == args.filenameOut.strip()):
        print("Output file must be different to input file.")
        sys.exit()

    reorder_particle_files(args.filenameIn, args.filenameOut)
