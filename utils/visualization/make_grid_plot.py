import argparse
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np

#-------------------------------------------------------------------------------

def plot_gridded_data(filenameIn, filenameOut, varname):

    fileIn = Dataset(filenameIn,"r")

    array = fileIn.variables[varname][:]
    if (array.ndim == 3):
        array = array[0,:,:]

    fileIn.close()

    fig, axis = plt.subplots()

    cax = axis.imshow(np.transpose(array),origin="lower")
    fig.colorbar(cax)

    plt.savefig(filenameOut,dpi=150)
    plt.cla()
    plt.close(fig)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Make a plot of DEMSI gridded data.')

    parser.add_argument('-i', dest='filenameIn', help='File to plot')
    parser.add_argument('-o', dest='filenameOut', help='Output image')
    parser.add_argument('-v', dest='varname', help='Variable to plot')

    args = parser.parse_args()

    plot_gridded_data(args.filenameIn, args.filenameOut, args.varname)
