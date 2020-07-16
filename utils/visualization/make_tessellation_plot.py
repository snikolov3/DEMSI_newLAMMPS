from netCDF4 import Dataset
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import argparse
import numpy as np
import matplotlib.cm as cm

#-------------------------------------------------------------------------------

def plot_tessellation(filemeshIn, filenameIn, filenameOut, varname, xplotmin=None, xplotmax=None, yplotmin=None, yplotmax=None, gridFilename=None, cmin=None, cmax=None):

    # read in mesh
    fileIn = Dataset(filemeshIn,"r")

    nParticles = len(fileIn.dimensions["nParticles"])
    maxVertices = len(fileIn.dimensions["maxVertices"])

    t = fileIn.variables["type"][:]
    nVerticesOnCell = fileIn.variables["nVerticesOnCell"][:]
    xVertex = fileIn.variables["xVertex"][:]
    yVertex = fileIn.variables["yVertex"][:]

    fileIn.close()

    # read in data
    fileIn = Dataset(filenameIn,"r")

    globalIndexInit = fileIn.variables["globalIndexInit"][:]
    v = fileIn.variables[varname][0,:]

    fileIn.close()

    # plot
    fig, axis = plt.subplots()

    cmap = cm.viridis

    patches = []
    colors = []

    for iParticle in range(0,nParticles):
        iParticleMesh = globalIndexInit[iParticle]

        if (t[iParticleMesh] != 2):

            vertices = []
            for iVertex in range(0,nVerticesOnCell[iParticleMesh]):
                vertices.append((xVertex[iParticleMesh,iVertex],yVertex[iParticleMesh,iVertex]))

            patches.append(Polygon(vertices,True))
            colors.append(v[iParticle])

    pc = PatchCollection(patches, cmap=cmap)
    pc.set_array(np.array(colors))
    pc.set_clim(cmin, cmax)

    axis.add_collection(pc)

    fig.colorbar(pc)

    # coastline
    patches = []

    for iParticle in range(0,nParticles):
        iParticleMesh = globalIndexInit[iParticle]

        if (t[iParticleMesh] == 2):

            vertices = []
            for iVertex in range(0,nVerticesOnCell[iParticleMesh]):
                vertices.append((xVertex[iParticleMesh,iVertex],yVertex[iParticleMesh,iVertex]))

            patches.append(Polygon(vertices,True,fc="grey"))

    pc = PatchCollection(patches, match_original=True)

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

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Make a plot of DEMSI initial tessellation.')

    parser.add_argument('-m', dest='filemeshIn', help='Mesh file to plot')
    parser.add_argument('-i', dest='filenameIn', help='Data file to plot')
    parser.add_argument('-o', dest='filenameOut', help='Output image')
    parser.add_argument('-v', dest='varname', help='Variable to plot')
    parser.add_argument('--x0', dest='xplotmin', type=float, default=None, help='Minimum plot x')
    parser.add_argument('--x1', dest='xplotmax', type=float, default=None, help='Maximum plot x')
    parser.add_argument('--y0', dest='yplotmin', type=float, default=None, help='Minimum plot y')
    parser.add_argument('--y1', dest='yplotmax', type=float, default=None, help='Maximum plot y')
    parser.add_argument('-g', dest='gridFilename', default=None, help='Grid file to get domain from')
    parser.add_argument('--c0', dest='cmin', type=float, default=None, help='Minimum plot value')
    parser.add_argument('--c1', dest='cmax', type=float, default=None, help='Maximum plot value')

    args = parser.parse_args()

    plot_tessellation(args.filemeshIn, args.filenameIn, args.filenameOut, args.varname, args.xplotmin, args.xplotmax, args.yplotmin, args.yplotmax, args.gridFilename, args.cmin, args.cmax)
