from netCDF4 import Dataset
import glob
import matplotlib.pyplot as plt
import numpy as np

# all filenames
filenames = sorted(glob.glob("./output/*"))

# initial positions of elements
fileIn = Dataset(filenames[0],"r")
x0 = fileIn.variables["x"][:]
y0 = fileIn.variables["y"][:]
r0 = fileIn.variables["radius"][:]
fileIn.close()

# final positions
fileIn = Dataset(filenames[-1],"r")
xf = fileIn.variables["x"][:]
yf = fileIn.variables["y"][:]
fileIn.close()

# position change since start of simulation
for i in range(0,100):

    xTime = []
    for filename in filenames:

        fileIn = Dataset(filename,"r")

        x = fileIn.variables["x"][:]
        xTime.append(x[i]-x0[i])

        fileIn.close()

    plt.plot(xTime)

plt.show()

# y position
changed = False
for filename in filenames:

    fileIn = Dataset(filename,"r")

    y = fileIn.variables["y"][:]

    if (not np.array_equal(y0,y)):
        print("y arrays changed: ", filename, np.amin(y-y0), np.amax(y-y0))
        changed = True

    fileIn.close()

if (changed):
    print("FAIL: y values changed")
else:
    print("PASS: y values unchanged")

# plot final position
plt.plot(xf-x0)
plt.show()


# max overlap
overlaps = []
for filename in filenames:

    fileIn = Dataset(filename,"r")

    x = fileIn.variables["x"][:]

    overlap = x[1:] - x[:-1] - r0[1:] - r0[:-1]

    overlaps.append(np.amax(overlap))

plt.plot(overlaps)
plt.show()
