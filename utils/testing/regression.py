#!/usr/bin/env python

from __future__ import print_function
import sys
import subprocess
import argparse
import os
from modify_xml import modify_config_file
from compare_files import compare_files

# command line arguments
parser = argparse.ArgumentParser(description='Perform a smoke test')

parser.add_argument('-n', '--testName', dest='testName', 
                    help='Name of the test.')
parser.add_argument('-b', '--baseline', dest='baselineDir', 
                    help='Path to the DEMSI baseline directory.')
parser.add_argument('-t', '--testdir', dest='testDir', 
                    help='Path to the DEMSI test directory.')
parser.add_argument('-c', '--cmpfile', dest='cmpFile', 
                    help='Output filename to compare.')
parser.add_argument('-x', '--configFile', dest='demsiConfigFile',
                    help='Path to the DEMSI config file to test.')
parser.add_argument('-p', '--nProcs', dest='nProcs', type=int ,
                    help='Number of processors to test with.')
parser.add_argument('--nThreads', dest='nThreads', type=int , default=4,
                    help='Number of threads to test with.')
parser.add_argument('-d', '--simulationDuration', dest='simulationDuration', default=None,
                    help='Path to the DEMSI config file to test.')
parser.add_argument('--configChanges', dest='configChanges', default=None, nargs='*', metavar=('[ConfigName=ConfigValue]'),
                    help='Additonal modifications to default config file.')

args = parser.parse_args()

# create test output directory
outputDirectory = "./output_%s_testdir/" %(args.testName)
if not os.path.exists(outputDirectory):
    os.makedirs(outputDirectory)

# filename for modified xml config file
demsiConfigFileTest = os.path.splitext(os.path.basename(args.demsiConfigFile))[0] + "_" + args.testName + "_test.xml"

# changes to config file
configNames = []
configValues = []

# output particle directory to change in xml config file
configNames.append("particleWriteOutputDirectory")
configValues.append(outputDirectory)

# output ocean directory to change in xml config file
configNames.append("oceanWriteOutputDirectory")
configValues.append(outputDirectory)

# simulation duration to change in xml config file
if (args.simulationDuration is not None):
    configNames.append("simulationStopType")
    configValues.append("simulationDuration")
    configNames.append("simulationDuration")
    configValues.append(args.simulationDuration)

# additional changes to config
if (args.configChanges is not None):
    for configChange in args.configChanges:
        configNames.append(configChange.split("=")[0])
        configValues.append(configChange.split("=")[1])

modify_config_file(args.demsiConfigFile, demsiConfigFileTest, configNames, configValues)

# run smoke test
returnCode = subprocess.call(["mpirun", "-np", str(args.nProcs), args.testDir+"/demsi", demsiConfigFileTest,"-kokkos","threads","%d"%args.nThreads,"numa","1"])

if (returnCode != 0):
    exit(returnCode)

# comparison files
if (args.baselineDir != ""):

    relTestDir = os.path.relpath(os.getcwd(), args.testDir)

    baseFilename = args.baselineDir + "/" + relTestDir + "/" + outputDirectory + args.cmpFile
    testFilename = args.testDir     + "/" + relTestDir + "/" + outputDirectory + args.cmpFile

    # compare with baseline run
    returnCode = subprocess.call(["cmp", baseFilename, testFilename])

    if (returnCode == 0):
        # files identical
        print("Comparison files identical")
        exit(0)

    elif (returnCode == 1):
        # files differ
        print("Comparison files differ")
        compare_files(baseFilename,testFilename)
        exit(1)

    else:
        # other error
        print ("cmp: error code: ", returnCode)
        exit(1)

else:

    print ("No baseline defined")
    exit(1)
