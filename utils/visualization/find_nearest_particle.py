import matplotlib.pyplot as plt
import argparse
from make_particle_plot import plot_particles_figure, add_plot_particle_args
import numpy as np
from netCDF4 import Dataset

#-------------------------------------------------------------------------------

def print_nearest_particle_index(px, py, filenameIn):

    # find nearest particle
    fileIn = Dataset(filenameIn,"r")

    nParticles = len(fileIn.dimensions["nParticles"])

    x = fileIn.variables["x"][:]
    y = fileIn.variables["y"][:]
    globalID = fileIn.variables["globalID"][:]

    fileIn.close()

    distance = np.sqrt(np.add(np.power(x-px, 2), np.power(y-py, 2)))
    iParticleNearest = np.argmin(distance)

    print("Nearest element: ", iParticleNearest, ", globalID: ", globalID[iParticleNearest], ", distance: ", np.amin(distance), ", position: ", x[iParticleNearest], y[iParticleNearest])

#-------------------------------------------------------------------------------

def find_nearest_particle(args):

    fig, axis = plt.subplots()

    fig, axis = plot_particles_figure(fig, axis, args)

    axis.plot([args["px"]], [args["py"]], marker="+", color="red")

    # create plot
    plt.savefig(args["filenameOut"],dpi=300)
    plt.cla()
    plt.close(fig)

    # print the index of the nearest particle to the point
    print_nearest_particle_index(args["px"], args["py"], args["filenameIn"])

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Make a plot of DEMSI particles and display a marker at specified point. Also return index of nearest particle to that marker.')

    parser.add_argument('--px', dest='px', type=float, required=True, help='x position of marker point.')
    parser.add_argument('--py', dest='py', type=float, required=True, help='y position of marker point.')
    parser.add_argument('-i', dest='filenameIn', required=True, help='File to plot')
    parser.add_argument('-o', dest='filenameOut', required=True, help='Output image')

    add_plot_particle_args(parser)

    args = parser.parse_args()

    find_nearest_particle(vars(args))
