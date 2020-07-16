from netCDF4 import Dataset
import glob
import matplotlib.pyplot as plt

filenames = sorted(glob.glob("./output/*"))

x = []
ux = []

for filename in filenames:

    filein = Dataset(filename,"r")

    x.append(filein.variables["x"][0])
    ux.append(filein.variables["ux"][0])

    filein.close()

fig, axis = plt.subplots()
axis.plot(x)
axis.set_xlabel("Time")
axis.set_ylabel("x position (m)")
plt.savefig("position.png")

fig, axis = plt.subplots()
axis.plot(ux)
axis.set_xlabel("Time")
axis.set_ylabel("x velocity (m/s)")
plt.savefig("velocity.png")
