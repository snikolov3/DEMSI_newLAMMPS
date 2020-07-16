#!/usr/bin/env python

import sys
import subprocess
import argparse
import os
from modify_xml import modify_config_file

# command line arguments
parser = argparse.ArgumentParser(description='Perform a check test')

parser.add_argument('-n', '--testName', dest='testName',
                    help='Name of the test.')
parser.add_argument('-e', '--executable', dest='demsiExecutable',
                    help='Path to the DEMSI executable to test.')
parser.add_argument('-x', '--configFile', dest='demsiConfigFile',
                    help='Path to the DEMSI config file to test.')
parser.add_argument('-p', '--nProcs', dest='nProcs', type=int ,
                    help='Number of processors to test with.')
parser.add_argument('-c', '--checkscript', dest='checkscript',
                    help='Name of script to check test results.')

args = parser.parse_args()

# create test output directory
outputDirectory = "./output_%s_testdir/" %(args.testName)
if not os.path.exists(outputDirectory):
    os.makedirs(outputDirectory)

# filename for modified xml config file
demsiConfigFileTest = os.path.splitext(os.path.basename(args.demsiConfigFile))[0] + "_" + args.testName + "_test.xml"

# output directory to change in xml config file
configNames = ["particleWriteOutputDirectory"]
configValues = [outputDirectory]

# modify the xml config file
modify_config_file(args.demsiConfigFile, demsiConfigFileTest, configNames, configValues)

# run check test
returnCode = subprocess.call(["mpirun", "-np", str(args.nProcs), args.demsiExecutable, demsiConfigFileTest])

if (returnCode == 1):
    # DEMSI process failed
    print("DEMSI process failed.")
    exit(returnCode)
else:
    # call check script
    returnCode = subprocess.call(["python", args.checkscript, outputDirectory])
    exit(returnCode)
