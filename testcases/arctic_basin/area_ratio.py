import glob
from netCDF4 import Dataset
import numpy as np

filenames = sorted(glob.glob("./output/particles_arctic_basin.*"))

for filename in filenames:

    filein = Dataset(filename,"a")

    nParticles = len(filein.dimensions["nParticles"])

    effectiveElementArea = filein.variables["effectiveElementArea"][:]
    areaInit = filein.variables["areaInit"][:]

    areaRatio = np.divide(effectiveElementArea,areaInit)

    var = filein.createVariable("areaRatio","d",dimensions=["nParticles"])
    var[:] = areaRatio[:]

    filein.close()
