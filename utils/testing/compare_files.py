#!/usr/bin/env python

from __future__ import print_function
from netCDF4 import Dataset
import argparse
import numpy as np

#-------------------------------------------------------------------------------

def check_attributes(variable1, variable2, globalAtt=False):

    if (globalAtt):
        variableName = "global attributes"
    else:
        variableName = variable1.name;

    # global attributes:
    globalAttributes1 = set(variable1.ncattrs())
    globalAttributes2 = set(variable2.ncattrs())

    attsIn1not2 = globalAttributes1 - globalAttributes2
    if (len(attsIn1not2) > 0):
        print("Attributes in file1 but not file2 for {}: ".format(variableName))
        for attname in sorted(attsIn1not2):
            print("  -- ",attname)

    attsIn2not1 = globalAttributes2 - globalAttributes1
    if (len(attsIn2not1) > 0):
        print("Attributes in file2 but not file1 for {}: ".format(variableName))
        for attname in sorted(attsIn2not1):
            print("  -- ",attname)

    # check global attributes
    globalAttributesInBoth = globalAttributes1 & globalAttributes2
    globalAttributeMaxLen=0
    for globalAttribute in globalAttributesInBoth:
        globalAttributeMaxLen = max(globalAttributeMaxLen,len(globalAttribute))
    for globalAttribute in sorted(globalAttributesInBoth):
        if (getattr(file1, globalAttribute) != getattr(file1, globalAttribute)):
            print("Global Attribute {:{}} differs for {}: file1: {}, file2: {}".format(globalAttribute, globalAttributeMaxLen, variableName, getattr(file1, globalAttribute), getattr(file1, globalAttribute)))

#-------------------------------------------------------------------------------

def compare_files(filename1, filename2):

    file1 = Dataset(filename1,"r")
    file2 = Dataset(filename2,"r")

    # check dimensions present
    dimensions1 = set(file1.dimensions.keys())
    dimensions2 = set(file2.dimensions.keys())

    dimsIn1not2 = dimensions1 - dimensions2
    if (len(dimsIn1not2) > 0):
        print("Dimensions in file1 but not file2: ")
        for dimname in sorted(dimsIn1not2):
            print("  -- ",dimname)

    dimsIn2not1 = dimensions2 - dimensions1
    if (len(dimsIn2not1) > 0):
        print("Dimensions in file2 but not file1: ")
        for dimname in sorted(dimsIn2not1):
            print("  -- ",dimname)

    # check dimensions
    dimensionsInBoth = dimensions1 & dimensions2
    dimensionMaxLen=0
    for dimension in dimensionsInBoth:
        dimensionMaxLen = max(dimensionMaxLen,len(dimension))
    for dimension in dimensionsInBoth:
        len1 = len(file1.dimensions[dimension])
        len2 = len(file2.dimensions[dimension])
        if (len1 != len2):
            print("Dimension size differs for {:{}}: file1: {:5}, file2: {:5}".format(dimension, dimensionMaxLen, len1, len2))

    # check variables present
    varnames1 = set(file1.variables.keys())
    varnames2 = set(file2.variables.keys())

    varsIn1not2 = varnames1 - varnames2
    if (len(varsIn1not2) > 0):
        print("Variables in file1 but not file2: ")
        for varname in sorted(varsIn1not2):
            print("  -- ",varname)

    varsIn2not1 = varnames2 - varnames1
    if (len(varsIn2not1) > 0):
        print("Variables in file2 but not file1: ")
        for varname in sorted(varsIn2not1):
            print("  -- ",varname)

    # check variables
    varnamesInBoth = varnames1 & varnames2
    varnameMaxLen=0
    for varname in varnamesInBoth:
        varnameMaxLen = max(varnameMaxLen,len(varname))
    globalIDIndices1 = np.argsort(file1.variables["globalID"][:])
    globalIDIndices2 = np.argsort(file2.variables["globalID"][:])
    for varname in sorted(varnamesInBoth):
        try:
            diff = file1.variables[varname][globalIDIndices1] - file2.variables[varname][globalIDIndices2]

            mindiff  = np.amin(diff)
            maxdiff  = np.amax(diff)
            meandiff = np.mean(diff)

            if (mindiff != 0.0 or maxdiff != 0.0 or meandiff != 0.0):
                print("Variable {:{}} differs: min: {:15e}, max: {:15e}, mean: {:15e}".format(varname, varnameMaxLen, mindiff, maxdiff, meandiff))

        except:
            print("Could not compare arrays for variable: ", varname)

    # variable attributes
    for varname in sorted(varnamesInBoth):
        check_attributes(file1.variables[varname], file2.variables[varname])

    # global attributes
    check_attributes(file1, file2, True)

    file1.close()
    file2.close()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    # command line arguments
    parser = argparse.ArgumentParser(description='Check differences between two netcdf files')

    parser.add_argument('--f1', dest='filename1',
                        help='Name of the first file to check.')
    parser.add_argument('--f2', dest='filename2',
                        help='Name of the second file to check.')

    args = parser.parse_args()

    compare_files(args.filename1, args.filename2)
