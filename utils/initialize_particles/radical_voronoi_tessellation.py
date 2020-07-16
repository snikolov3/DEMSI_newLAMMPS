from netCDF4 import Dataset
import numpy as np
import argparse
import os

#-------------------------------------------------------------------------------

# function to get unique values
def unique(list1):

    # insert the list to the set
    list_set = set(list1)
    # convert the set to the list
    unique_list = (list(list_set))

    return

#-------------------------------------------------------------------------------

def clockwise(x, y):

    cx = np.mean(x)
    cy = np.mean(y)
    a = np.arctan2(y - cy, x - cx)

    order = a.ravel().argsort()

    x = x[order]
    y = y[order]

    return np.vstack([x,y])

#-------------------------------------------------------------------------------

def radical_voronoi_tessellation(filenameIn, filenameOut, voroExecutable, xmin, xmax, ymin, ymax):

    # particle in file
    fileDEMSI = Dataset(filenameIn,"r")

    nParticles = len(fileDEMSI.dimensions["nParticles"])

    ids = fileDEMSI.variables["globalID"][:]
    x = fileDEMSI.variables["x"][:,0]
    y = fileDEMSI.variables["x"][:,1]
    r = fileDEMSI.variables["radius"][:]

    fileDEMSI.close()

    maxRadius = np.amax(r)

    zmin = -maxRadius
    zmax =  maxRadius


    # create voro input file
    fileVORO = open("voro.tmp","w")

    for iParticle in range(0,nParticles):

        fileVORO.write("%i %e %e %e %e\n" %(ids[iParticle],x[iParticle],y[iParticle],0.0,r[iParticle]))

    fileVORO.close()


    # run voro
    cmd = "%s -r -o -g -c \"%%i %%w %%P\" %f %f %f %f %f %f voro.tmp" %(voroExecutable, xmin, xmax, ymin, ymax, zmin, zmax)
    os.system(cmd)


    # read in voro data
    filein = open("voro.tmp.vol","r")
    lines = filein.readlines()
    filein.close()


    nCells = len(lines)

    nVerticesArray = []
    for line in lines:

        lineSplit = line.split()
        nVerticesArray.append(int(lineSplit[1]))



    maxVertices = int(np.amax(nVerticesArray) / 2)


    xVertexOut = np.zeros((nCells,maxVertices))
    yVertexOut = np.zeros((nCells,maxVertices))

    cellIDs = []
    nVerticesArray = []
    iCell = -1
    for line in lines:
        iCell = iCell + 1

        lineSplit = line.split()
        cellIDs.append(int(lineSplit[0]))

        vertices = []
        for i in range(2,len(lineSplit)):

            vertexStr = lineSplit[i].strip(")").strip("(").split(",")
            vertices.append([float(vertexStr[0]),float(vertexStr[1])])

        vertices = np.array(vertices)
        vertices = np.unique(vertices,axis=0)

        nVertices = vertices.shape[0]

        nVerticesArray.append(nVertices)

        xVertex = vertices[:,0]
        yVertex = vertices[:,1]

        vertices = clockwise(xVertex, yVertex)

        for iVertex in range(0,nVertices):
            xVertexOut[iCell,iVertex] = vertices[0][iVertex]
            yVertexOut[iCell,iVertex] = vertices[1][iVertex]



    # write out DEMSI compatable netcdf file with tesselation
    fileOut = Dataset(filenameOut,"w",format="NETCDF3_CLASSIC")

    fileOut.createDimension("nParticles",nCells)
    fileOut.createDimension("maxVertices",maxVertices)

    var = fileOut.createVariable("globalID","i",dimensions=["nParticles"])
    var[:] = cellIDs[:]
    var.units = "-"

    var = fileOut.createVariable("nVerticesOnCell","i",dimensions=["nParticles"])
    var[:] = nVerticesArray[:]
    var.units = "-"

    var = fileOut.createVariable("xVertex","d",dimensions=["nParticles","maxVertices"])
    var[:] = xVertexOut[:]
    var.units = "m"

    var = fileOut.createVariable("yVertex","d",dimensions=["nParticles","maxVertices"])
    var[:] = yVertexOut[:]
    var.units = "m"

    fileOut.close()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Perform the radical tessellation of a DEMSI particle distribution.')

    parser.add_argument('-i', dest='filenameIn', help='Input DEMSI particle file to tessellate')
    parser.add_argument('-o', dest='filenameOut', help='Output file of tessellation')
    parser.add_argument('-v', dest='voroExecutable', help='Path to the voro++ executable')
    parser.add_argument('--x0', dest='xmin', type=float, default=None, help='Minimum x to tessellate')
    parser.add_argument('--x1', dest='xmax', type=float, default=None, help='Maximum x to tessellate')
    parser.add_argument('--y0', dest='ymin', type=float, default=None, help='Minimum y to tessellate')
    parser.add_argument('--y1', dest='ymax', type=float, default=None, help='Maximum y to tessellate')

    args = parser.parse_args()

    run_voro(args.filenameIn, args.gridFilename, args.filenameOut, args.voroExecutable, args.xmin, args.xmax, args.ymin, args.ymax)
