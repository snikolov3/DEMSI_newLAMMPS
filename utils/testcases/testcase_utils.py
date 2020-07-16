from netCDF4 import Dataset
import netCDF4
import numpy as np
import math

#-------------------------------------------------------------------------------

def write_out_grid(filenameOut, nx, ny, dx, dy, lat, lon):

    fileOut = Dataset(filenameOut,"w",format="NETCDF3_CLASSIC")

    fileOut.createDimension("TWO",2)
    fileOut.createDimension("nx",nx)
    fileOut.createDimension("ny",ny)

    resolutionVar = fileOut.createVariable("resolution","d",dimensions=["TWO"])
    resolutionVar[:] = [dx,dy]

    dataVar = fileOut.createVariable("latitude","d",dimensions=["nx","ny"])
    dataVar[:] = lat[:]

    dataVar = fileOut.createVariable("longitude","d",dimensions=["nx","ny"])
    dataVar[:] = lon[:]

    fileOut.close()

#-------------------------------------------------------------------------------

def write_out_fixed_forcing(filenameOut, nx, ny, ugx, ugy):

    fileOut = Dataset(filenameOut,"w",format="NETCDF3_CLASSIC")

    fileOut.climatology = 1

    fileOut.createDimension("strLen", 64)
    fileOut.createDimension("Time", 1)
    fileOut.createDimension("nx",nx)
    fileOut.createDimension("ny",ny)

    var = fileOut.createVariable("Time","c",dimensions=["Time","strLen"])

    year = 0
    month = 1
    day = 1
    hour = 0
    minute = 0
    second = 0

    timeStr = "%4.4i-%2.2i-%2.2i_%2.2i:%2.2i:%2.2i" %(year,month,day,hour,minute,second)
    var[0,0:19] = netCDF4.stringtochar(np.array(timeStr, 'S19'))

    varx = fileOut.createVariable("xAtmWind","d",dimensions=["Time","nx","ny"])
    vary = fileOut.createVariable("yAtmWind","d",dimensions=["Time","nx","ny"])

    varx[0,:,:] = ugx[:,:]
    vary[0,:,:] = ugy[:,:]

    fileOut.close()

#-------------------------------------------------------------------------------

def write_out_time_varying_forcing(filenameOut, nx, ny, ugx, ugy):

    # time varying forcing
    fileOut = Dataset(filenameOut,"w",format="NETCDF3_CLASSIC")

    fileOut.climatology = 0

    fileOut.createDimension("strLen", 64)
    fileOut.createDimension("Time", 10)
    fileOut.createDimension("nx",nx)
    fileOut.createDimension("ny",ny)

    var = fileOut.createVariable("Time","c",dimensions=["Time","strLen"])

    year = 1
    month = 1
    hour = 0
    minute = 0
    second = 0

    for day in range(0,10):
        timeStr = "%4.4i-%2.2i-%2.2i_%2.2i:%2.2i:%2.2i" %(year,month,day+1,hour,minute,second)
        var[day,0:19] = netCDF4.stringtochar(np.array(timeStr, 'S19'))

    varx = fileOut.createVariable("xAtmWind","d",dimensions=["Time","nx","ny"])
    vary = fileOut.createVariable("yAtmWind","d",dimensions=["Time","nx","ny"])
    for day in range(0,10):
        varx[day,:,:] = ugx[:,:] * (float(2-day) / 2.0)
        vary[day,:,:] = ugy[:,:] * (float(2-day) / 2.0)

    fileOut.close()

#-------------------------------------------------------------------------------

def add_init_polygons(filename, maxCellSize, maxVertices, nVerticesOnCell, xVertex, yVertex):

    # add polygon to output file
    fileParticlesOut = Dataset(filename,"a")

    fileParticlesOut.maxCellSize = maxCellSize

    fileParticlesOut.createDimension("maxVertices",maxVertices)

    var = fileParticlesOut.createVariable("nVerticesOnCell","i",dimensions=["nParticles"])
    var[:] = nVerticesOnCell[:]
    var.units = "-"

    var = fileParticlesOut.createVariable("xVertex","d",dimensions=["nParticles","maxVertices"])
    var[:] = xVertex[:]
    var.units = "m"

    var = fileParticlesOut.createVariable("yVertex","d",dimensions=["nParticles","maxVertices"])
    var[:] = yVertex[:]
    var.units = "m"

    fileParticlesOut.close()

#-------------------------------------------------------------------------------

def add_init_connectivity(filenamein):

    debug = False

    filein = Dataset(filenamein,"a")

    nParticles  = len(filein.dimensions["nParticles"])
    maxVertices = len(filein.dimensions["maxVertices"])

    x = filein.variables["x"][:,0]
    y = filein.variables["x"][:,1]

    nVerticesOnCell = filein.variables["nVerticesOnCell"][:]
    xVertex         = filein.variables["xVertex"][:]
    yVertex         = filein.variables["yVertex"][:]

    if (debug): print("nParticles: ", nParticles)

    # find all unique edges
    edgesNonUnique = []
    edgesCell = []
    xEdgesCell = {}
    yEdgesCell = {}

    for iParticle in range(0,nParticles):
        edgesCellCell = []
        for iVertexOnCell in range(0,nVerticesOnCell[iParticle]):

            iVertexOnCell1 = iVertexOnCell
            iVertexOnCell2 = (iVertexOnCell + 1) % nVerticesOnCell[iParticle]

            xEdge = 0.5 * (xVertex[iParticle,iVertexOnCell1] + xVertex[iParticle,iVertexOnCell2])
            yEdge = 0.5 * (yVertex[iParticle,iVertexOnCell1] + yVertex[iParticle,iVertexOnCell2])

            if (xVertex[iParticle,iVertexOnCell1] > xVertex[iParticle,iVertexOnCell2]):

                edge = (((xVertex[iParticle,iVertexOnCell1],yVertex[iParticle,iVertexOnCell1]),
                         (xVertex[iParticle,iVertexOnCell2],yVertex[iParticle,iVertexOnCell2])))

            elif (xVertex[iParticle,iVertexOnCell1] < xVertex[iParticle,iVertexOnCell2]):

                edge = (((xVertex[iParticle,iVertexOnCell2],yVertex[iParticle,iVertexOnCell2]),
                         (xVertex[iParticle,iVertexOnCell1],yVertex[iParticle,iVertexOnCell1])))

            else:

                if (yVertex[iParticle,iVertexOnCell1] > yVertex[iParticle,iVertexOnCell2]):

                    edge = (((xVertex[iParticle,iVertexOnCell1],yVertex[iParticle,iVertexOnCell1]),
                             (xVertex[iParticle,iVertexOnCell2],yVertex[iParticle,iVertexOnCell2])))

                elif (yVertex[iParticle,iVertexOnCell1] < yVertex[iParticle,iVertexOnCell2]):

                    edge = (((xVertex[iParticle,iVertexOnCell2],yVertex[iParticle,iVertexOnCell2]),
                             (xVertex[iParticle,iVertexOnCell1],yVertex[iParticle,iVertexOnCell1])))

                else:
                    print("vertices are the same")
                    sys.exit()

            edgesCellCell.append(edge)
            edgesNonUnique.append(edge)

            if (edge not in xEdgesCell):
                xEdgesCell[edge] = xEdge
                yEdgesCell[edge] = yEdge


        edgesCell.append(edgesCellCell)

    edges = set(edgesNonUnique)
    nEdges = len(edges)

    if (debug): print("nEdges: ", nEdges)

    edges = list(edges)

    edgesDict = {}

    iEdge = 0
    for edge in edges:

        edgesDict[edge] = iEdge
        iEdge += 1

    edgesOnCell = np.zeros((nParticles,maxVertices),dtype="i")
    edgesOnCell[:,:] = -1


    for iParticle in range(0,nParticles):
        for iVertexOnCell in range(0,nVerticesOnCell[iParticle]):

            edgesOnCell[iParticle,iVertexOnCell] = edgesDict[edgesCell[iParticle][iVertexOnCell]]


    if (debug): print(edgesOnCell)

    cellsOnEdgesDict = {}

    xEdges = np.zeros(nEdges)
    yEdges = np.zeros(nEdges)


    for iParticle in range(0,nParticles):
        for iVertexOnCell in range(0,nVerticesOnCell[iParticle]):

            iEdge = edgesOnCell[iParticle,iVertexOnCell]

            if (iEdge in cellsOnEdgesDict):

                cellsOnEdgesDict[iEdge].append(iParticle)

            else:

                cellsOnEdgesDict[iEdge] = [iParticle]

            xEdges[iEdge] = xEdgesCell[edgesCell[iParticle][iVertexOnCell]]
            yEdges[iEdge] = yEdgesCell[edgesCell[iParticle][iVertexOnCell]]


    if (debug): print(xEdges)
    if (debug): print(yEdges)
    if (debug): print(cellsOnEdgesDict)

    cellsOnEdge = np.zeros((nEdges,2),dtype="i")
    cellsOnEdge[:] = -1

    for iEdge in cellsOnEdgesDict:

        if (len(cellsOnEdgesDict[iEdge]) == 2):
            cellsOnEdge[iEdge,:] = cellsOnEdgesDict[iEdge][:]
        elif (len(cellsOnEdgesDict[iEdge]) == 1):
            cellsOnEdge[iEdge,0] = cellsOnEdgesDict[iEdge][0]
        else:
            print("too many cells on edge")
            sys.exit()

    if (debug): print(cellsOnEdge)

    cellsOnCell = np.zeros((nParticles,maxVertices),dtype="i")
    for iParticle in range(0,nParticles):
        for iVertexOnCell in range(0,nVerticesOnCell[iParticle]):
            iEdge = edgesOnCell[iParticle,iVertexOnCell]
            iCell1 = cellsOnEdge[iEdge,0]
            iCell2 = cellsOnEdge[iEdge,1]
            if (iParticle == iCell1):
                cellsOnCell[iParticle,iVertexOnCell] = iCell2
            else:
                cellsOnCell[iParticle,iVertexOnCell] = iCell1



    # gnuplot output
    if (debug):

        fileout = open("polygons.txt","w")

        for iParticle in range(0,nParticles):
            for iVertexOnCell in range(0,nVerticesOnCell[iParticle]):

                fileout.write("%e %e\n" %(xVertex[iParticle,iVertexOnCell],yVertex[iParticle,iVertexOnCell]))

            fileout.write("%e %e\n" %(xVertex[iParticle,0],yVertex[iParticle,0]))

            fileout.write("\n")

        fileout.close()

        fileout = open("edges.txt","w")

        for iEdge in range(0,nEdges):
            fileout.write("%e %e\n" %(xEdges[iEdge],yEdges[iEdge]))

        fileout.close()

        fileout = open("edges_single.txt","w")

        for iEdge in range(0,nEdges):
            if (cellsOnEdge[iEdge,0] == -1 or cellsOnEdge[iEdge,1] == -1):
                fileout.write("%e %e\n" %(xEdges[iEdge],yEdges[iEdge]))

        fileout.close()

        fileout = open("edges_duplicated.txt","w")

        for iEdge in range(0,nEdges):
            if (cellsOnEdge[iEdge,0] == cellsOnEdge[iEdge,1]):
                fileout.write("%e %e\n" %(xEdges[iEdge],yEdges[iEdge]))

        fileout.close()

        fileout = open("cell_to_cell.txt","w")

        for iEdge in range(0,nEdges):
            if (cellsOnEdge[iEdge,0] != -1 and cellsOnEdge[iEdge,1] != -1):
                fileout.write("%e %e\n" %(x[cellsOnEdge[iEdge,0]],y[cellsOnEdge[iEdge,0]]))
                fileout.write("%e %e\n" %(x[cellsOnEdge[iEdge,1]],y[cellsOnEdge[iEdge,1]]))
            fileout.write("\n")

        fileout.close()

        fileout = open("cells_on_cell.txt","w")

        for iParticle in range(0,nParticles):
            for iVertexOnCell in range(0,nVerticesOnCell[iParticle]):

                iCell1 = iParticle
                iCell2 = cellsOnCell[iParticle,iVertexOnCell]

                if (iCell1 != -1 and iCell2 != -1):

                    fileout.write("%e %e\n" %(x[iCell1],y[iCell1]))
                    fileout.write("%e %e\n" %(x[iCell2],y[iCell2]))
                    fileout.write("\n")

        fileout.close()

    # output connectivity
    filein.createDimension("nEdges",nEdges)
    try:
        filein.createDimension("TWO",2)
    except:
        # probably already exists
        pass

    var = filein.createVariable("xEdge","d",dimensions=["nEdges"])
    var[:] = xEdges[:]
    var.units = "m"

    var = filein.createVariable("yEdge","d",dimensions=["nEdges"])
    var[:] = yEdges[:]
    var.units = "m"

    var = filein.createVariable("cellsOnEdge","i",dimensions=["nEdges","TWO"])
    var[:] = cellsOnEdge[:]
    var.units = "-"

    var = filein.createVariable("edgesOnCell","i",dimensions=["nParticles","maxVertices"])
    var[:] = edgesOnCell[:]
    var.units = "-"

    var = filein.createVariable("cellsOnCell","i",dimensions=["nParticles","maxVertices"])
    var[:] = cellsOnCell[:]
    var.units = "-"

    filein.close()

#-------------------------------------------------------------------------------

def plot_particles_debug(particles, array):

    import matplotlib.pyplot as plt
    import matplotlib.patches as mp
    import matplotlib.collections as mc
    import matplotlib.cm as cm

    fig, axis = plt.subplots()

    patches = []
    colors = []

    for iParticle in range(0, particles.num_particles):
        patches.append(mp.Circle((particles.x[iParticle], particles.y[iParticle]), particles.radius[iParticle]))
        colors.append(array[iParticle])

    pc = mc.PatchCollection(patches, cmap=cm.viridis)
    pc.set_array(np.array(colors))

    axis.add_collection(pc)

    fig.colorbar(pc)

    axis.set_xlim((np.amin(np.subtract(particles.x,particles.radius)), np.amax(np.add(particles.x,particles.radius))))
    axis.set_ylim((np.amin(np.subtract(particles.y,particles.radius)), np.amax(np.add(particles.y,particles.radius))))
    axis.set_aspect('equal')

    plt.show()

#-------------------------------------------------------------------------------

def write_particle_file(filenameOut, nParticles, d, nTypes, t, x, y, r, f=None, h=None, maxVertices=None, nv=None, xv=None, yv=None):

    fileOut = Dataset(filenameOut,"w",format="NETCDF3_CLASSIC")

    fileOut.maxRadius = np.amax(r)
    fileOut.nTypes = nTypes

    fileOut.createDimension("nParticles",nParticles)
    fileOut.createDimension("TWO",2)

    var = fileOut.createVariable("globalID","i",dimensions=["nParticles"])
    var.units = "-"
    var[:] = d[:]

    var = fileOut.createVariable("type","i",dimensions=["nParticles"])
    var.units = "-"
    var[:] = t[:]

    var = fileOut.createVariable("x","d",dimensions=["nParticles","TWO"])
    var.units = "m"
    var[:,0] = x[:]
    var[:,1] = y[:]

    var = fileOut.createVariable("radius","d",dimensions=["nParticles"])
    var.units = "m"
    var[:] = r[:]

    if (f is not None):
        var = fileOut.createVariable("iceFraction","d",dimensions=["nParticles"])
        var.units = "-"
        var[:] = f[:]

    if (h is not None):
        var = fileOut.createVariable("iceThickness","d",dimensions=["nParticles"])
        var.units = "m"
        var[:] = h[:]

    if (maxVertices is not None):
        fileOut.createDimension("maxVertices",maxVertices)

    if (nv is not None):
        var = fileOut.createVariable("nVerticesOnCell","i",dimensions=["nParticles"])
        var.units = "-"
        var[:] = nv[:]

    if (xv is not None):
        var = fileOut.createVariable("xVertex","d",dimensions=["nParticles","maxVertices"])
        var.units = "m"
        var[:] = xv[:]

    if (yv is not None):
        var = fileOut.createVariable("yVertex","d",dimensions=["nParticles","maxVertices"])
        var.units = "m"
        var[:] = yv[:]

    if (nv is not None and xv is not None and yv is not None):
        maxCellSize = 0.0
        for iParticle in range(0,nParticles):
            for iVertex in range(1,nv[iParticle]):
                distance = math.sqrt(math.pow(xv[iParticle,iVertex]-xv[iParticle,0],2) +
                                     math.pow(yv[iParticle,iVertex]-yv[iParticle,0],2))
                maxCellSize = max(maxCellSize,distance)

        fileOut.maxCellSize = maxCellSize

    fileOut.close()

#-------------------------------------------------------------------------------
