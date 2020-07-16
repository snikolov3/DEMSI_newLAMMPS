from netCDF4 import Dataset
import glob
import matplotlib.pyplot as plt
import numpy as np


nBondsArray = []
totalBondLengths = []

filenames = sorted(glob.glob("./output/*"))

for filename in filenames:

    filein = Dataset(filename,"r")

    try:
        nBonds = len(filein.dimensions["nBonds"])
        bondCrackFraction = filein.variables["bondCrackFraction"][:,:]
        bondLength = filein.variables["bondLength"][:]
        bond_indices = np.where(bondCrackFraction[:,1] > bondCrackFraction[:,0])[0]
        nBondsArray.append(len(bond_indices))
        totalBondLengths.append(np.sum(np.multiply(bondCrackFraction[bond_indices,1]-bondCrackFraction[bond_indices,0], bondLength[bond_indices])))
    except:
        nBondsArray.append(0)
        totalBondLengths.append(0.0)

    filein.close()

# plot
fig, ax1 = plt.subplots()
ax1.plot(nBondsArray, '-x', color="red")
ax1.set_xlabel("Time")
ax1.set_ylabel("nBonds")
ax1.tick_params('y', colors="red")
ax2 = ax1.twinx()
ax2.plot(totalBondLengths, color="blue")
ax2.set_ylabel("Total bond length")
ax2.tick_params('y', colors="blue")
plt.savefig("bond_timeseries.png")
