import matplotlib.pyplot as plt
import argparse
from make_particle_plot import plot_particles_figure, add_plot_particle_args
import numpy as np
from netCDF4 import Dataset

#-------------------------------------------------------------------------------

def add_marker_to_particle(fig, axis, pi, filenameIn):

    # find nearest particle
    fileIn = Dataset(filenameIn,"r")

    x = fileIn.variables["x"][:]
    y = fileIn.variables["y"][:]

    fileIn.close()

    axis.plot([x[pi]], [y[pi]], marker="+", color="red")

    return fig, axis

#-------------------------------------------------------------------------------

def locate_particle(args):

    fig, axis = plt.subplots()

    fig, axis = plot_particles_figure(fig, axis, args)

    fig, axis = add_marker_to_particle(fig, axis, args["pi"], args["filenameIn"])

    # create plot
    plt.savefig(args["filenameOut"],dpi=300)
    plt.cla()
    plt.close(fig)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Make a plot of DEMSI particles and display a marker at a specified element.')

    parser.add_argument('--pi', dest='pi', type=int, required=True, help='index of particle to add marker to')
    parser.add_argument('-i', dest='filenameIn', required=True, help='File to plot')
    parser.add_argument('-o', dest='filenameOut', required=True, help='Output image')

    add_plot_particle_args(parser)

    args = parser.parse_args()

    locate_particle(vars(args))
