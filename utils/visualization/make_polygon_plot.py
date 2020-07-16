from netCDF4 import Dataset
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import argparse
import numpy as np
import matplotlib.cm as cm

#-------------------------------------------------------------------------------

def plot_polygons(filenameIn, filenameOut, xplotmin=None, xplotmax=None, yplotmin=None, yplotmax=None, gridFilename=None, varname=None):

    # read in data
    fileIn = Dataset(filenameIn,"r")

    nParticles = len(fileIn.dimensions["nParticles"])
    maxVertices = len(fileIn.dimensions["maxVertices"])

    nVerticesOnCell = fileIn.variables["nVerticesOnCell"][:]
    xVertex = fileIn.variables["xVertex"][:]
    yVertex = fileIn.variables["yVertex"][:]

    # plot
    fig, axis = plt.subplots()

    if (varname is not None):

        try:
            v = fileIn.variables[varname][:]
        except:
            print("Unknown varname: ", varname)
            raise

        cmap = cm.viridis

        patches = []
        colors = []

        for iParticle in range(0,nParticles):

            vertices = []
            for iVertex in range(0,nVerticesOnCell[iParticle]):
                vertices.append((xVertex[iParticle,iVertex],yVertex[iParticle,iVertex]))

            patches.append(Polygon(vertices,True))
            colors.append(v[iParticle])

        pc = PatchCollection(patches, cmap=cmap)
        pc.set_array(np.array(colors))

        axis.add_collection(pc)

        fig.colorbar(pc)

    else:

        patches = []

        for iParticle in range(0,nParticles):

            vertices = []
            for iVertex in range(0,nVerticesOnCell[iParticle]):
                vertices.append((xVertex[iParticle,iVertex],yVertex[iParticle,iVertex]))

            patches.append(Polygon(vertices,True))

        pc = PatchCollection(patches, facecolors="white", edgecolors="black", linewidths=0.5)

        axis.add_collection(pc)


    # plot domain
    if (gridFilename is not None):
        gridFile = Dataset(gridFilename,"r")
        nx = len(gridFile.dimensions["nx"])
        ny = len(gridFile.dimensions["ny"])
        resolution = gridFile.variables["resolution"][:]
        gridFile.close()
        xplotminUse = 0.0
        xplotmaxUse = (nx-1) * resolution[0]
        yplotminUse = 0.0
        yplotmaxUse = (ny-1) * resolution[1]
    else:
        if (xplotmin is None or \
            xplotmax is None or \
            yplotmin is None or \
            yplotmax is None):
            print("Must specify plot domain")
            sys.exit(-1)
        else:
            xplotminUse = xplotmin
            xplotmaxUse = xplotmax
            yplotminUse = yplotmin
            yplotmaxUse = yplotmax

    axis.set_xlim((xplotminUse, xplotmaxUse))
    axis.set_ylim((yplotminUse, yplotmaxUse))
    axis.set_aspect('equal')

    # create plot
    plt.savefig(filenameOut,dpi=300)
    plt.cla()
    plt.close(fig)

    fileIn.close()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Make a plot of DEMSI particle polygons.')

    parser.add_argument('-i', dest='filenameIn', help='File to plot')
    parser.add_argument('-o', dest='filenameOut', help='Output image')
    parser.add_argument('--x0', dest='xplotmin', type=float, default=None, help='Minimum plot x')
    parser.add_argument('--x1', dest='xplotmax', type=float, default=None, help='Maximum plot x')
    parser.add_argument('--y0', dest='yplotmin', type=float, default=None, help='Minimum plot y')
    parser.add_argument('--y1', dest='yplotmax', type=float, default=None, help='Maximum plot y')
    parser.add_argument('-g', dest='gridFilename', default=None, help='Grid file to get domain from')
    parser.add_argument('-v', dest='varname', help='Variable to plot')

    args = parser.parse_args()

    plot_polygons(args.filenameIn, args.filenameOut, args.xplotmin, args.xplotmax, args.yplotmin, args.yplotmax, args.gridFilename, args.varname)
