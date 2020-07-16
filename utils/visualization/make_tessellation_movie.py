from make_tessellation_plot import plot_tessellation
import glob, os, sys
import argparse

parser = argparse.ArgumentParser(description='Make a movie of DEMSI initial tessellation.')

parser.add_argument('-m', dest='filemeshIn', help='Mesh file to plot')
parser.add_argument('-f', dest='filenameTemplate', required=True, help='Files to create movie of')
parser.add_argument('-o', dest='movieName', required=True, help='Output movie name')
parser.add_argument('-v', dest='varname', help='Variable to plot')
parser.add_argument('--x0', dest='xplotmin', type=float, default=None, help='Minimum plot x')
parser.add_argument('--x1', dest='xplotmax', type=float, default=None, help='Maximum plot x')
parser.add_argument('--y0', dest='yplotmin', type=float, default=None, help='Minimum plot y')
parser.add_argument('--y1', dest='yplotmax', type=float, default=None, help='Maximum plot y')
parser.add_argument('-g', dest='gridFilename', default=None, help='Grid file to get domain from')
parser.add_argument('--c0', dest='cmin', type=float, default=None, help='Minimum plot value')
parser.add_argument('--c1', dest='cmax', type=float, default=None, help='Maximum plot value')

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

    filenameOut = "./tmp/tmp_%5.5i.png" %(iPlot)

    plot_tessellation(args.filemeshIn, filenameIn, filenameOut, args.varname, args.xplotmin, args.xplotmax, args.yplotmin, args.yplotmax, args.gridFilename, args.cmin, args.cmax)

cmd = "ffmpeg -i tmp/tmp_%%05d.png -vcodec h264 %s" %(args.movieName)
os.system(cmd)
