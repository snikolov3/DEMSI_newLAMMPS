import sys, os
from netCDF4 import Dataset
from reorder_particles import reorder_particle_files
import argparse
import numpy as np

parser = argparse.ArgumentParser(description='Reorder a particle file by globalID.')

parser.add_argument('--f1', dest='filenameIn1', help='First particle file to compare.')
parser.add_argument('--f2', dest='filenameIn2', help='Second particle file to compare.')
parser.add_argument('-v', dest='verbose', action='store_true', help='Output cell differences.')

args = parser.parse_args()

if (os.path.isfile("tmp_compare_1.nc")):
    os.remove("tmp_compare_1.nc")
if (os.path.isfile("tmp_compare_2.nc")):
    os.remove("tmp_compare_2.nc")

reorder_particle_files(args.filenameIn1, "tmp_compare_1.nc")
reorder_particle_files(args.filenameIn2, "tmp_compare_2.nc")

file1 = Dataset("tmp_compare_1.nc", "r")
file2 = Dataset("tmp_compare_2.nc", "r")

nParticles1 = len(file1.dimensions["nParticles"])
nParticles2 = len(file2.dimensions["nParticles"])

if (nParticles1 != nParticles2):
    print("nParticles1 != nParticles2")
    sys.exit()

varnamesDifferent = []
minDifferences = []
maxDifferences = []

for varname in file1.variables:

    var1 = file1.variables[varname][:]
    var2 = file2.variables[varname][:]

    if (not np.array_equal(var1, var2)):
        varnamesDifferent.append(varname)

        diff = var1 - var2
        minDifferences.append(np.amin(diff))
        maxDifferences.append(np.amax(diff))

        if (args.verbose and var1.ndim == 1):

            for i in range(0,var1.size):
                if (var1[i] != var2[i]):
                    print(i, varname, var1[i], var2[i], var1[i]-var2[i])

if (len(varnamesDifferent) > 0):
    print(args.filenameIn1, " different: ", varnamesDifferent, minDifferences, maxDifferences)
else:
    print(args.filenameIn1, " same")

file1.close()
file2.close()
