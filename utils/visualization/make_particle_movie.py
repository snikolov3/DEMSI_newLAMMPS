from make_particle_plot import plot_particles, add_plot_particle_args
import glob, os, sys
import argparse

parser = argparse.ArgumentParser(description='Make a movie of DEMSI particles.')

parser.add_argument('-f', dest='filenameTemplate', required=True, help='Files to create movie of')
parser.add_argument('-o', dest='movieName', required=True, help='Output movie name')

add_plot_particle_args(parser)

args = parser.parse_args()

# check tmp directory exists
if (not os.path.isdir("./tmp")):
    print("tmp directory not found")
    sys.exit()

# test if tmp empty
if (not os.listdir("./tmp") == []):
    print("tmp directory not empty")
    sys.exit()

filenamesIn = sorted(glob.glob(args.filenameTemplate))

iPlot = 0
for filenameIn in filenamesIn:
    iPlot = iPlot + 1

    print("Making frame %i of %i..." %(iPlot,len(filenamesIn)))

    args.filenameIn = filenameIn
    args.filenameOut = "./tmp/tmp_%5.5i.png" %(iPlot)

    plot_particles(vars(args))

cmd = "ffmpeg -i tmp/tmp_%%05d.png -vcodec h264 %s" %(args.movieName)
os.system(cmd)
