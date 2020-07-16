import xml.etree.ElementTree as ET
import sys
import argparse
from six.moves.urllib.parse import urlparse
import os, os.path

#-------------------------------------------------------------------------------

def download_file(url, destination, filenameIn, filenameOut, proxy):

    import subprocess

    print(url, destination, filenameIn, filenameOut, proxy)

    if (proxy == "none"):
        process = subprocess.Popen(["wget", "-O", "%s%s" %(destination,filenameOut), "%s%s" %(url,filenameIn)], stdout=subprocess.PIPE)
    else:
        process = subprocess.Popen(["wget", "-O", "%s%s" %(destination,filenameOut), "-e", "use_proxy=yes", "-e", proxy, "%s%s" %(url,filenameIn)], stdout=subprocess.PIPE)

    while process.poll() is None:
        line = process.stdout.readline() # This blocks until it receives a newline.
        print(line)
    print(process.stdout.read())

#-------------------------------------------------------------------------------

# This script downloads data needed for various DEMSI test cases

# Arctic test case
print("Downloading data from arctic test case...")

parser = argparse.ArgumentParser(description='Download data needed for DEMSI test cases.')

parser.add_argument('-d', '--download', dest='downloadFile', help='Download configuration file', default="download.xml")
parser.add_argument('-p', '--proxytype', dest='proxytype', help='Proxy type', choices=['none','lanl'], default="none")

args = parser.parse_args()

# proxies
proxies = {"none": {"http":  "none",
                    "https": "none"},
           "lanl": {"http":  "http_proxy=http://proxyout.lanl.gov:8080",
                    "https": "https_proxy=http://proxyout.lanl.gov:8080"}}

# iterate over downloads
tree = ET.parse(args.downloadFile)
downloads = tree.getroot()

for download in downloads:

    url = download.find('url').text
    destination = download.find('destination').text

    if not os.path.exists(destination):
        os.makedirs(destination)

    filenames = download.find('filenames')
    for filename in filenames:

        inputFilename  = filename.find('input').text
        outputFilename = filename.find('output').text

        proxy = proxies[args.proxytype][urlparse(url).scheme]
        download_file(url, destination, inputFilename, outputFilename, proxy)

        if (os.path.splitext(outputFilename)[-1] == ".gz"):
            os.system("gunzip %s%s" %(destination,outputFilename))

# environment variable to set
print("Set the DEMSI data environment variable:")
print("export DEMSI_DATA_DIR=%s" %(os.getcwd()))
