#!/usr/bin/env python

from __future__ import print_function
import xml.etree.ElementTree as ET
import argparse
import os.path

#-------------------------------------------------------------------------------

def modify_config(root, configName, configNewValue):

    for config in root.iter("Config"):

        if (config.attrib["name"] == configName):

            config.attrib["value"] = configNewValue

    return root

#-------------------------------------------------------------------------------

def modify_config_file(filenameIn, filenameOut, configNames, configNewValues):

    if (not os.path.isfile(filenameIn)):
        print("modify_config_file: Input file does not exist")
        exit(1)

    if (len(configNames) != len(configNewValues)):
        print("modify_config_file: configNames and configNewValues not same length")
        exit(1)

    tree = ET.parse(filenameIn)
    root = tree.getroot()

    for configName, configNewValue in zip(configNames, configNewValues):

        root = modify_config(root, configName, configNewValue)

    tree.write(filenameOut)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Modify a DEMSI XML input file')
    parser.add_argument('-i', '--input', dest='filenameIn', required=True,
                        help='Path to the input XML file.')
    parser.add_argument('-o', '--output', dest='filenameOut', required=True,
                        help='Path to the output XML file.')
    parser.add_argument('-n', '--names', dest='configNames', required=True, nargs='+',
                        help='Config names to modify.')
    parser.add_argument('-v', '--values', dest='configNewValues', required=True, nargs='+',
                        help='Values to modify configs with.')

    args = parser.parse_args()

    modify_config_file(args.filenameIn, args.filenameOut, args.configNames, args.configNewValues)
