from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import glob

#times = [
#    "0001-01-02_00:00:00",
#    "0001-01-06_00:05:00",
#    "0001-01-16_00:05:00"]

#labels = [
#    "1 day",
#    "5 days",
#    "16 days"]

times = []
linetypes = []
daysInMonth = [31,28]
#daysInMonth = [2,28]
for month in range(1,2):
    for day in range(0,daysInMonth[month-1]):
        for hour in range(0,24,24):
            times.append("0001-%2.2i-%2.2i_%2.2i:00:00" %(month,day+1,hour))
            linetypes.append("solid")

#times = [
#    "0001-01-02_00:00:00",
#    "0001-01-03_00:00:00"]

#labels = [
#    "1 day",
#    "5 days"]

#linetypes = [
#    "solid",
#    "solid"]


for time, linetype in zip(times, linetypes):

    filename = "./output/particles_out.%s.nc" %(time)

#filenames = sorted(glob.glob("./output/particles_out.*.nc"))

#for filename in filenames:

    fileIn = Dataset(filename,"r")

    print(filename)

    nParticles  = len(fileIn.dimensions["nParticles"])
    nCategories = len(fileIn.dimensions["nCategories"])

    typeIn = fileIn.variables["type"][:]
    xIn = fileIn.variables["x"][:,0]
    xIn[:] = xIn[:] / 1000.0
    iceAreaCategoryIn   = fileIn.variables["iceAreaCategory"][:]
    iceVolumeCategoryIn = fileIn.variables["iceVolumeCategory"][:]

    x = []
    iceAreaCategory = []
    iceVolumeCategory = []
    for iParticle in range(0,nParticles):
        if (typeIn[iParticle] == 1):
            x.append(xIn[iParticle])
            iceAreaCategory.append(iceAreaCategoryIn[iParticle,:])
            iceVolumeCategory.append(iceVolumeCategoryIn[iParticle,:])
    x = np.array(x)
    iceAreaCategory = np.array(iceAreaCategory)
    iceVolumeCategory = np.array(iceVolumeCategory)

    nParticles -= 1

    xx = [x[0]-5.0]
    h = [0.0]

    for iParticle in range(0,nParticles):
        iceVolume = 0.0
        iceArea = 0.0
        for iCategory in range(0,nCategories):
            iceVolume += iceVolumeCategory[iParticle,iCategory]
            iceArea   += iceAreaCategory  [iParticle,iCategory]
        h.append(iceVolume)
        xx.append(x[iParticle])

    plt.plot(xx, h, c="black", ls=linetype, lw=0.3)

    fileIn.close()

axis = plt.gca()
axis.set_xlim(0.0,1000.0)
axis.set_ylim(0.0,5.0)
axis.set_xlabel("Position (km)")
axis.set_ylabel("Ice thickness (m)")
#plt.gca().legend(labels)

plt.tight_layout()
plt.savefig("thicknesses.eps")
