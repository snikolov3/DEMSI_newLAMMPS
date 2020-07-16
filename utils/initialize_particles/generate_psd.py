from  __future__ import print_function, division, absolute_import
# Python script to generate a histogram corresponding to the desired
# particle size distribution.

import numpy as np
import argparse
import sys

#-------------------------------------------------------------------------------

def write_histogram(binedges, binvalues, filename):
    """
    Write a histogram in the format expected by the particle placement code

    Parameters
    ----------
    binedges : (N+1) X 1 numpy array
      These are the edges of the bins of the histogram that define
      the particle size distribution. Units are in meters, and values
      correspond to particle diameter
    binvalues : N X 1 numpy array
      Values of frequency for each bin. These do not need to be normalized,
      as the particle placement code will do that.
    filename : string
      Name of file to write histogram to
    """
    if len(binedges) != len(binvalues)+1:
        print("Size of bin edges must be one greater than size of bin values\n")
        return

    with open(filename, "w") as f:
        binedges.tofile(f, sep="\n")
        f.write("\nValues:\n")
        binvalues.tofile(f, sep="\n")

#-------------------------------------------------------------------------------

def generate_psd(type, mean=None, stddev=None, min_diameter=None, minval=None, maxval=None, diameter=None):

    if (type == "gaussian"):

        # Gaussian distribution
        if (mean is None):
            print("Must specify mean diameter (-m) for gaussian option.")
            sys.exit()
        if (stddev is None):
            print("Must specify standard deviation (-s) for gaussian option.")
            sys.exit()
        if (min_diameter is None):
            print("Must specify minimum diameter (-n) for gaussian option.")
            sys.exit()

        values = np.random.normal(mean, stddev, 100000)
        if np.any(values < min_diameter):
            print("Warning: some radius values were less than min. diameter, setting these to the minimum diameter\n")
            values[values < min_diameter] = min_diameter

        binvalues, binedges = np.histogram(values, 100)

        outname = "gaussian_mean%i_stdev%i.txt" %(mean, stddev)
        write_histogram(binedges, binvalues, outname)

    elif (type == "uniform"):

        # uniform distribution
        if (minval is None):
            print("Must specify minimum diameter (-d1) for uniform option.")
            sys.exit()
        if (maxval is None):
            print("Must specify maximum diameter (-d2) for uniform option.")
            sys.exit()

        values = np.random.uniform(low=minval,high=maxval,size=100)

        binvalues, binedges = np.histogram(values, 100)

        outname = "uniform_min%i_max%i.txt" %(minval, maxval)
        write_histogram(binedges, binvalues, outname)

    elif (type == "monodisperse"):

        # monodisperse distribution
        if (diameter is None):
            print("Must specify element diameter (-d) for monodisperse option.")
            sys.exit()

        outname = "mono_%i.txt" %(diameter)
        eps = np.finfo(np.float).eps*10
        with open(outname, "w") as f:
            f.write(str(diameter-eps)+"\n"+str(diameter+eps)+"\nValues:\n1")

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Make a particle size distribution')

    parser.add_argument('-t', '--type', dest='type', help='Type of distribution', choices=['gaussian','uniform','monodisperse'])

    # gaussian options
    parser.add_argument('-m', '--mean', dest='mean', help='[gaussian] Mean diameter (m).', type=float)
    parser.add_argument('-s', '--stddev', dest='stddev', help='[gaussian] Standard deviation diameter (m).', type=float)
    parser.add_argument('-n', '--min', dest='min_diameter', help='[gaussian] Minimum diameter (m).', type=float)

    # uniform options
    parser.add_argument('-d1', '--minval', dest='minval', help='[uniform] Minimum diameter (m).', type=float)
    parser.add_argument('-d2', '--maxval', dest='maxval', help='[uniform] Maximum diameter (m).', type=float)

    # monodisperse options
    parser.add_argument('-d', '--diameter', dest='diameter', help='[monodisperse] Element diameter (m).', type=float)

    args = parser.parse_args()

    generate_psd(args.type, args.mean, args.stddev, args.min_diameter, args.minval, args.maxval, args.diameter)
