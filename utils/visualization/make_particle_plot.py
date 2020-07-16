from netCDF4 import Dataset
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
from matplotlib.collections import LineCollection
import glob, os, sys
import argparse
import numpy as np
import matplotlib.cm as cm

#-------------------------------------------------------------------------------

def plot_particles_figure(fig, axis, args):

    # read in data
    fileIn = Dataset(args["filenameIn"],"r")

    nParticles = len(fileIn.dimensions["nParticles"])

    x = fileIn.variables["x"][:,0]
    y = fileIn.variables["x"][:,1]
    r = fileIn.variables["radius"][:]
    try:
        t = fileIn.variables["type"][:]
    except:
        pass

    try:

        nBonds = len(fileIn.dimensions["nBonds"])

        b1p1 = fileIn.variables["bondEndPoint1Particle1"][:]
        b2p1 = fileIn.variables["bondEndPoint2Particle1"][:]
        b1p2 = fileIn.variables["bondEndPoint1Particle2"][:]
        b2p2 = fileIn.variables["bondEndPoint2Particle2"][:]

        bondCrackFraction = fileIn.variables["bondCrackFraction"][:]

    except:
        nBonds = 0

    # plot

    if (args["varname"] is not None):

        try:
            v = fileIn.variables[args["varname"]][:]
        except:
            print("Unknown varname: ", args["varname"])
            raise

        cmap = cm.viridis

        patches = []
        colors = []

        for iParticle in range(0, nParticles):
            if (t[iParticle] == 1):
                patches.append(Circle((x[iParticle], y[iParticle]), r[iParticle]))
                colors.append(v[iParticle])

        pc = PatchCollection(patches, cmap=cmap)
        pc.set_array(np.array(colors))
        pc.set_clim([args["cmin"], args["cmax"]])

        axis.add_collection(pc)
        fig.colorbar(pc)

        # coastline
        if (args["coastline"]):
            patches = []
            for iParticle in range(0, nParticles):
                if (t[iParticle] == 2):
                    patches.append(Circle((x[iParticle], y[iParticle]), r[iParticle], fc=args["coastcolor"]))

            pc = PatchCollection(patches, match_original=True)

            axis.add_collection(pc)

    elif (args["plotType"]):

        cmap = cm.viridis

        patches = []

        for iParticle in range(0, nParticles):
            patches.append(Circle((x[iParticle], y[iParticle]), r[iParticle]))

        pc = PatchCollection(patches, cmap=cmap)
        pc.set_array(t)

        axis.add_collection(pc)

        # legend
        if (args["useLegend"]):
            norm = matplotlib.colors.Normalize(vmin=np.amin(t), vmax=np.amax(t))

            patchList = []
            for uniqueType in set(t):
                legendKey = "Type %i" %(uniqueType)
                patchList.append(mpatches.Patch(color=cmap(norm(uniqueType)), label=legendKey))
            plt.legend(handles=patchList)

    else:

        patches = []

        for iParticle in range(0, nParticles):
            patches.append(Circle((x[iParticle], y[iParticle]), r[iParticle], fill=False, edgecolor='k', linewidth=0.3))

        pc = PatchCollection(patches, match_original=True)

        axis.add_collection(pc)

    # bonds
    if (args["plotBonds"]):

        for iBond in range(0,nBonds):

            if (bondCrackFraction[iBond,0] != bondCrackFraction[iBond,1]):

                axis.plot([b1p1[iBond,0],b2p1[iBond,0]],[b1p1[iBond,1],b2p1[iBond,1]],color="darkcyan")
                axis.plot([b1p2[iBond,0],b2p2[iBond,0]],[b1p2[iBond,1],b2p2[iBond,1]],color="darkcyan")

                xEndPoint1 = 0.5 * (b1p1[iBond,0] + b1p2[iBond,0])
                yEndPoint1 = 0.5 * (b1p1[iBond,1] + b1p2[iBond,1])

                xEndPoint2 = 0.5 * (b2p1[iBond,0] + b2p2[iBond,0])
                yEndPoint2 = 0.5 * (b2p1[iBond,1] + b2p2[iBond,1])

                axis.arrow(xEndPoint1,yEndPoint1,xEndPoint2-xEndPoint1,yEndPoint2-yEndPoint1, head_width=250, head_length=500.0, length_includes_head=True)

                xCrackPoint1 = xEndPoint1 + bondCrackFraction[iBond,0] * (xEndPoint2 - xEndPoint1)
                yCrackPoint1 = yEndPoint1 + bondCrackFraction[iBond,0] * (yEndPoint2 - yEndPoint1)

                xCrackPoint2 = xEndPoint2 - (1.0-bondCrackFraction[iBond,1]) * (xEndPoint2 - xEndPoint1)
                yCrackPoint2 = yEndPoint2 - (1.0-bondCrackFraction[iBond,1]) * (yEndPoint2 - yEndPoint1)

                axis.plot([xEndPoint1,xCrackPoint1],[yEndPoint1,yCrackPoint1],color="red",linewidth=4)
                axis.plot([xEndPoint2,xCrackPoint2],[yEndPoint2,yCrackPoint2],color="red",linewidth=4)

    # plot domain
    if (args["gridFilename"] is not None and \
        args["xplotmin"] is None or \
        args["xplotmax"] is None or \
        args["yplotmin"] is None or \
        args["yplotmax"] is None):
        gridFile = Dataset(args["gridFilename"],"r")
        nx = len(gridFile.dimensions["nx"])
        ny = len(gridFile.dimensions["ny"])
        resolution = gridFile.variables["resolution"][:]
        gridFile.close()
        xplotminUse = 0.0
        xplotmaxUse = (nx-1) * resolution[0]
        yplotminUse = 0.0
        yplotmaxUse = (ny-1) * resolution[1]
    elif (args["gridFilename"] is None and \
          args["xplotmin"] is not None or \
          args["xplotmax"] is not None or \
          args["yplotmin"] is not None or \
          args["yplotmax"] is not None):
        xplotminUse = args["xplotmin"]
        xplotmaxUse = args["xplotmax"]
        yplotminUse = args["yplotmin"]
        yplotmaxUse = args["yplotmax"]
    else:
        print("Error: Either '-g' must be specified or all of '--x0', '--x1', '--y0', and '--y1'")
        sys.exit(-1)

    axis.set_xlim((xplotminUse, xplotmaxUse))
    axis.set_ylim((yplotminUse, yplotmaxUse))
    axis.set_aspect('equal')

    # remove tick marks and labels
    if (args["removeTicks"]):
        axis.ticklabel_format(style='plain')
        axis.tick_params( \
            axis='x', \
            which='both', \
            bottom=False, \
            top=False, \
            labelbottom=False)
        axis.tick_params( \
            axis='y', \
            which='both', \
            left=False, \
            right=False, \
            labelleft=False)

    fileIn.close()

    return fig, axis

#-------------------------------------------------------------------------------

def plot_particles(args):

    fig, axis = plt.subplots()

    fig, axis = plot_particles_figure(fig, axis, args)

    # create plot
    plt.savefig(args["filenameOut"],dpi=300)
    plt.cla()
    plt.close(fig)

#-------------------------------------------------------------------------------

def add_plot_particle_args(parser):

    parser.add_argument('--x0', dest='xplotmin', type=float, default=None, help='Minimum plot x')
    parser.add_argument('--x1', dest='xplotmax', type=float, default=None, help='Maximum plot x')
    parser.add_argument('--y0', dest='yplotmin', type=float, default=None, help='Minimum plot y')
    parser.add_argument('--y1', dest='yplotmax', type=float, default=None, help='Maximum plot y')
    parser.add_argument('-g', dest='gridFilename', default=None, help='Grid file to get domain from')
    parser.add_argument('-t', dest='plotType', action='store_true', help='Plot the particle type')
    parser.add_argument('-l', dest='useLegend', action='store_true', help='Include a legend of plot types')
    parser.add_argument('--removeticks', dest='removeTicks', action='store_true', help='Remove tick marks and labels from plot')
    parser.add_argument('-b', dest='plotBonds', action='store_true', help='Plot the particle bonds')
    parser.add_argument('-v', dest='varname', default=None, help='Variable to plot')
    parser.add_argument('--c0', dest='cmin', type=float, default=None, help='Minimum plot value')
    parser.add_argument('--c1', dest='cmax', type=float, default=None, help='Maximum plot value')
    parser.add_argument('-c','--coastline', dest='coastline', action='store_true', help='Plot the coastline with an input variable')
    parser.add_argument('--coastcolor', dest='coastcolor', default="grey", help='color to plot coastline elements')

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Make a plot of DEMSI particles.')

    parser.add_argument('-i', dest='filenameIn', required=True, help='File to plot')
    parser.add_argument('-o', dest='filenameOut', required=True, help='Output image')

    add_plot_particle_args(parser)

    args = parser.parse_args()

    plot_particles(vars(args))
