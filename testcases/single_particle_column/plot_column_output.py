from netCDF4 import Dataset
import glob
import matplotlib.pyplot as plt

filenames = sorted(glob.glob("./output/*"))

iceAreaCategory = []
iceVolumeCategory = []
snowVolumeCategory = []
surfaceTemperature = []

for filename in filenames:

    fileIn = Dataset(filename,"r")

    iceAreaCategory.append(fileIn.variables["iceAreaCategory"][0,0])
    iceVolumeCategory.append(fileIn.variables["iceVolumeCategory"][0,0])
    snowVolumeCategory.append(fileIn.variables["snowVolumeCategory"][0,0])
    surfaceTemperature.append(fileIn.variables["surfaceTemperature"][0,0])

    fileIn.close()

fig, ax1 = plt.subplots()

ax1.plot(surfaceTemperature,color="green")
ax1.set_ylabel("Temperature (C)")
ax1.set_xlabel("Time step")
ax1.set_xlim(0,4700)
ax1.set_ylim(None,0)
ax1.set_title("DEMSI")

ax2 = ax1.twinx()

ax2.plot(iceVolumeCategory,color="red")
ax2.plot(snowVolumeCategory,color="blue")
ax2.set_ylabel("Thickness (m)")
ax2.set_ylim(0,None)

plt.savefig("single_particle_column_output.png")
