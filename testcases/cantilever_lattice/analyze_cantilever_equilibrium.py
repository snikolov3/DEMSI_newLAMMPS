from netCDF4 import Dataset
import glob
from math import sqrt, pow
import matplotlib.pyplot as plt
import numpy as np

#-------------------------------------------------------------------------------

def cantilever_deflection(x, L, h, E, P):

    # https://en.wikipedia.org/wiki/Euler-Bernoulli_beam_theory

    I = pow(h,3) / 12.0

    w = (P * pow(x,2) * (3.0 * L - x)) / (6.0 * E * I)

    return w

#-------------------------------------------------------------------------------

# loop over cantilever output files
filenameExpression = "./output/*"

filenames = sorted(glob.glob(filenameExpression))

velocities = []

iFile = -1
for filename in filenames:
    iFile = iFile + 1

    fileIn = Dataset(filename,"r")

    nParticles = len(fileIn.dimensions["nParticles"])

    x = fileIn.variables["x"][:]
    y = fileIn.variables["y"][:]
    ux = fileIn.variables["ux"][:]
    uy = fileIn.variables["uy"][:]

    fileIn.close()

    velocity = np.sqrt(np.power(ux,2) + np.power(uy,2))
    velocities.append(np.mean(velocity))

plt.semilogy(velocities)
plt.title("Mean velocity of cantilever elements over time.")
plt.xlabel("Time (output index)")
plt.ylabel("Mean velocity (m/s)")
plt.show()


# last file
fileIn = Dataset(filenames[-1],"r")

nParticles = len(fileIn.dimensions["nParticles"])

x = fileIn.variables["x"][:]
y = fileIn.variables["y"][:]
ux = fileIn.variables["ux"][:]
uy = fileIn.variables["uy"][:]
radius = fileIn.variables["radius"][:]

fileIn.close()

plt.scatter(x[:], y[:])

xMin =  1e30
yMin = -1e30
carryOn = True
xBottom = []
yBottom = []
while(True):
    validIdx = np.where(np.logical_and(x < xMin,y > yMin))[0]
    if (len(validIdx) > 0):
        iMin = validIdx[y[validIdx].argmin()]
        xMin = x[iMin]
        yMin = y[iMin]
        xBottom.append(x[iMin])
        yBottom.append(y[iMin])
    else:
        break

plt.scatter(xBottom, yBottom, c="red")

xMin = -1e30
yMin =  1e30
carryOn = True
xTop = []
yTop = []
while(True):
    validIdx = np.where(np.logical_and(x > xMin,y < yMin))[0]
    if (len(validIdx) > 0):
        iMin = validIdx[y[validIdx].argmax()]
        xMin = x[iMin]
        yMin = y[iMin]
        xTop.append(x[iMin])
        yTop.append(y[iMin])
    else:
        break

plt.scatter(xTop, yTop, c="green")

plt.title("Cantilever element positions for final output.")
plt.xlabel("x position (m)")
plt.ylabel("y position (m)")
plt.show()

plt.plot(np.subtract(xBottom,np.amin(xBottom)),np.subtract(yBottom,np.amax(yBottom)),color="green",marker='x')
plt.plot(np.subtract(xTop,np.amin(xTop)),np.subtract(yTop,np.amax(yTop)),color="blue",marker='+')

# analytical: https://en.wikipedia.org/wiki/Euler-Bernoulli_beam_theory

# first file
fileIn = Dataset(filenames[0],"r")

nParticles = len(fileIn.dimensions["nParticles"])

xFirst = fileIn.variables["x"][:]
yFirst = fileIn.variables["y"][:]
radiusFirst = fileIn.variables["radius"][:]

fileIn.close()

r = np.mean(radiusFirst)
L = np.amax(xFirst) - np.amin(xFirst)
h = np.amax(yFirst) - np.amin(yFirst) + 2.0 * r

E = 1e9
P = -1e9 * 3.0

n = 100
x = []
y = []
for i in range(0,n):
    xx = (float(i) / float(n-1)) * L
    x.append(xx)
    w = cantilever_deflection(xx, L, h, E, P)
    y.append(w)

plt.plot(x, y, color="red")

plt.title("Cantilever deflection compared to analytical formula")
plt.xlabel("x position (m)")
plt.ylabel("Vertical beam deflection")
plt.legend(["Bottom cantilever elements","Top cantilever elements","Analytical formula"])
plt.show()
