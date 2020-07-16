import random
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Polygon
from matplotlib.collections import PatchCollection
import os, sys, math
from numpy.random import RandomState
from scipy.optimize import minimize
from netCDF4 import Dataset
from particle_distribution_quality import number_of_bonds, fraction_covered_area
import time
import argparse

timeVoronoi = 0.0
timeInscribed = 0.0
timeInPolygon = 0.0
timeOptimize = 0.0

#-------------------------------------------------------------------------------

def get_cell_coords(pt, a):
    """Get the coordinates of the cell that pt = (x,y) falls in."""

    return int(pt[0] // a), int(pt[1] // a)

#-------------------------------------------------------------------------------

def get_neighbours(cells, nx, ny, coords):
    """Return the indexes of points in cells neighbouring cell at coords.

    For the cell at coords = (x,y), return the indexes of points in the cells
    with neighbouring coordinates illustrated below: ie those cells that could
    contain points closer than r.

                                     ooo
                                    ooooo
                                    ooXoo
                                    ooooo
                                     ooo

    """

    dxdy = [(-1,-2),(0,-2),(1,-2),(-2,-1),(-1,-1),(0,-1),(1,-1),(2,-1),
            (-2,0),(-1,0),(1,0),(2,0),(-2,1),(-1,1),(0,1),(1,1),(2,1),
            (-1,2),(0,2),(1,2),(0,0)]
    neighbours = []
    for dx, dy in dxdy:
        neighbour_coords = coords[0] + dx, coords[1] + dy
        if not (0 <= neighbour_coords[0] < nx and
                0 <= neighbour_coords[1] < ny):
            # We're off the grid: no neighbours here.
            continue
        neighbour_cell = cells[neighbour_coords]
        if neighbour_cell is not None:
            # This cell is occupied: store this index of the contained point.
            neighbours.append(neighbour_cell)
    return neighbours

#-------------------------------------------------------------------------------

def point_valid(r, samples, cells, nx, ny, pt, a):
    """Is pt a valid point to emit as a sample?

    It must be no closer than r from any other point: check the cells in its
    immediate neighbourhood.

    """

    cell_coords = get_cell_coords(pt, a)
    for idx in get_neighbours(cells, nx, ny, cell_coords):
        nearby_pt = samples[idx]
        # Squared distance between or candidate point, pt, and this nearby_pt.
        distance2 = (nearby_pt[0]-pt[0])**2 + (nearby_pt[1]-pt[1])**2
        if distance2 < r**2:
            # The points are too close, so pt is not a candidate.
            return False
    # All points tested: if we're here, pt is valid
    return True

#-------------------------------------------------------------------------------

def get_point(samples, cells, nx, ny, a, r, width, height, k, refpt):
    """Try to find a candidate point relative to refpt to emit in the sample.

    We draw up to k points from the annulus of inner radius r, outer radius 2r
    around the reference point, refpt. If none of them are suitable (because
    they're too close to existing points in the sample), return False.
    Otherwise, return the pt.

    """
    i = 0
    while i < k:
        rho, theta = np.random.uniform(r, 2*r), np.random.uniform(0, 2*np.pi)
        pt = refpt[0] + rho*np.cos(theta), refpt[1] + rho*np.sin(theta)
        if not (0 <= pt[0] < width and 0 <= pt[1] < height):
            # This point falls outside the domain, so try again.
            continue
        if point_valid(r, samples, cells, nx, ny, pt, a):
            return pt
        i += 1
    # We failed to find a suitable point in the vicinity of refpt.
    return False

#-------------------------------------------------------------------------------

def poisson_disc_sampling(maxRadius, width, height):

    # https://scipython.com/blog/poisson-disc-sampling-in-python/
    # Choose up to k points around each reference point as candidates for a new
    # sample point
    k = 30

    # Minimum distance between samples
    #r = 1.7
    r = maxRadius

    #width, height = 60, 45

    # Cell side length
    a = r/np.sqrt(2)
    # Number of cells in the x- and y-directions of the grid
    nx, ny = int(width / a) + 1, int(height / a) + 1

    # A list of coordinates in the grid of cells
    coords_list = [(ix, iy) for ix in range(nx) for iy in range(ny)]
    # Initilalize the dictionary of cells: each key is a cell's coordinates, the
    # corresponding value is the index of that cell's point's coordinates in the
    # samples list (or None if the cell is empty).
    cells = {coords: None for coords in coords_list}

    # Pick a random point to start with.
    pt = (np.random.uniform(0, width), np.random.uniform(0, height))
    samples = [pt]
    # Our first sample is indexed at 0 in the samples list...
    cells[get_cell_coords(pt,a)] = 0
    # ... and it is active, in the sense that we're going to look for more points
    # in its neighbourhood.
    active = [0]

    nsamples = 1
    # As long as there are points in the active list, keep trying to find samples.
    while active:
        # choose a random "reference" point from the active list.
        idx = np.random.choice(active)
        refpt = samples[idx]
        # Try to pick a new point relative to the reference point.
        pt = get_point(samples, cells, nx, ny, a, r, width, height, k, refpt)
        if pt:
            # Point pt is valid: add it to the samples list and mark it as active
            samples.append(pt)
            nsamples += 1
            active.append(len(samples)-1)
            cells[get_cell_coords(pt,a)] = len(samples) - 1
        else:
            # We had to give up looking for valid points near refpt, so remove it
            # from the list of "active" points.
            active.remove(idx)

    return samples

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

def radical_voronoi_tessellation(nParticles, ids, x, y, r, voroExecutable, xmin, xmax, ymin, ymax, zmin, zmax):

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


    return maxVertices, cellIDs, nVerticesArray, xVertexOut, yVertexOut

#-------------------------------------------------------------------------------

def plot_voronoi_tessellation(filenameOut, nParticles, maxVertices, nVerticesArray, xVertex, yVertex, xmin, xmax, ymin, ymax):

    fig, axis = plt.subplots()


    patches = []
    for iParticle in range(0,nParticles):

        polygon = []
        for iVertex in range(0,nVerticesArray[iParticle]):
            polygon.append((xVertex[iParticle,iVertex],yVertex[iParticle,iVertex]))

        patches.append(Polygon(polygon,fill=False,edgecolor="red"))

    pc = PatchCollection(patches, match_original=True)
    axis.add_collection(pc)



    # plot options
    axis.set_xlim((xmin,xmax))
    axis.set_ylim((ymin,ymax))

    axis.set_aspect('equal')

    plt.savefig(filenameOut)

#-------------------------------------------------------------------------------

def plot_voronoi_tessellation_inscribed(filenameOut, nParticles, maxVertices, nVerticesArray, xVertex, yVertex, xmin, xmax, ymin, ymax, nParticlesCircle, xnew, ynew, rinscribed):

    fig, axis = plt.subplots()


    patches = []
    for iParticle in range(0,nParticles):

        polygon = []
        for iVertex in range(0,nVerticesArray[iParticle]):
            polygon.append((xVertex[iParticle,iVertex],yVertex[iParticle,iVertex]))

        patches.append(Polygon(polygon,fill=False,edgecolor="red"))

    pc = PatchCollection(patches, match_original=True)
    axis.add_collection(pc)

    patches = []
    for iParticle in range(0,nParticlesCircle):

        patches.append(Circle((xnew[iParticle], ynew[iParticle]), rinscribed[iParticle], fill=False, edgecolor='k', linewidth=0.3))

    pc = PatchCollection(patches, match_original=True)
    axis.add_collection(pc)

    # plot options
    axis.set_xlim((xmin,xmax))
    axis.set_ylim((ymin,ymax))

    axis.set_aspect('equal')

    plt.savefig(filenameOut)

#-------------------------------------------------------------------------------

def polygon_centroid(nVertices, xVertex, yVertex):

    area = 0.0

    for iVertex in range(0,nVertices):

        iVertexP1 = iVertex + 1
        if (iVertexP1 >= nVertices): iVertexP1 = 0

        area = area + xVertex[iVertex] * yVertex[iVertexP1] - xVertex[iVertexP1] * yVertex[iVertex]

    area = 0.5 * area

    xCentroid = 0.0
    yCentroid = 0.0

    for iVertex in range(0,nVertices):

        iVertexP1 = iVertex + 1
        if (iVertexP1 >= nVertices): iVertexP1 = 0

        xCentroid = xCentroid + (xVertex[iVertex] + xVertex[iVertexP1]) * (xVertex[iVertex] * yVertex[iVertexP1] - xVertex[iVertexP1] * yVertex[iVertex])
        yCentroid = yCentroid + (yVertex[iVertex] + yVertex[iVertexP1]) * (xVertex[iVertex] * yVertex[iVertexP1] - xVertex[iVertexP1] * yVertex[iVertex])

    xCentroid = xCentroid / (6.0 * area)
    yCentroid = yCentroid / (6.0 * area)

    return (xCentroid, yCentroid)

#-------------------------------------------------------------------------------

def is_left(x0, y0, x1, y1, x2, y2):

    return ( (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0) )

#-------------------------------------------------------------------------------

def in_polygon(x, y, nVertices, xVertex, yVertex):

    # http://geomalgorithms.com/a03-_inclusion.html

    windingNumber = 0

    for i0 in range(0, nVertices):

        i1 = i0 + 1
        if (i1 == nVertices): i1 = 0

        if (yVertex[i0] <= y):
            if (yVertex[i1] > y and is_left(xVertex[i0], yVertex[i0], xVertex[i1], yVertex[i1], x, y) > 0.0):
                windingNumber += 1

        else:
            if (yVertex[i1] <= y and is_left(xVertex[i0], yVertex[i0], xVertex[i1], yVertex[i1], x, y) < 0.0):
                windingNumber -= 1

    return (windingNumber != 0)

#-------------------------------------------------------------------------------

def distance_node_to_edge_2(node, xVertex1, yVertex1, xVertex2, yVertex2, edgeLength):

    distance = math.fabs((yVertex2-yVertex1)*node[0] - (xVertex2-xVertex1)*node[1] + xVertex2*yVertex1 - yVertex2*xVertex1)

    distance = distance / edgeLength

    return distance

#-------------------------------------------------------------------------------

def get_inscribed_circle_2(node, nVertices, xVertex, yVertex, edgeLengths):

    inscribed_radius = 1e30

    for iVertex1 in range(0,nVertices):

        iVertex2 = (iVertex1 + 1) % nVertices

        distance_to_node = distance_node_to_edge(node, xVertex[iVertex1], yVertex[iVertex1], xVertex[iVertex2], yVertex[iVertex2], edgeLengths[iVertex1])

        if (distance_to_node < inscribed_radius):
            inscribed_radius = distance_to_node

    return inscribed_radius

#-------------------------------------------------------------------------------

def largest_inscribed_circle_2(nVertices, xVertex, yVertex):

    rndState = RandomState(1234)

    min_x = np.amin(xVertex)
    max_x = np.amax(xVertex)
    min_y = np.amin(yVertex)
    max_y = np.amax(yVertex)

    maximum_inscribed_radius = 0.0
    best_node = None

    accuracy = 1e30
    minimum_accuracy = 1.0
    maximum_consecutive_misses = 5

    # precalc edge lengths
    edgeLengths = []
    for iVertex1 in range(0,nVertices):

        iVertex2 = (iVertex1 + 1) % nVertices

        edgeLengths.append(math.sqrt(math.pow((yVertex[iVertex2]-yVertex[iVertex1]),2) + math.pow((xVertex[iVertex2]-xVertex[iVertex1]),2)))

    edgeLengths = np.array(edgeLengths)

    it = 0
    # iterate random sampling region
    while (accuracy > minimum_accuracy):
        it += 1

        # begin loop through nodes
        consecutive_misses = 0

        while (consecutive_misses < maximum_consecutive_misses):

            # select coordinates at random within bounds
            x = rndState.uniform(min_x, max_x)
            y = rndState.uniform(min_y, max_y)
            node = (x, y)

            if (in_polygon(x, y, nVertices, xVertex, yVertex)):

                inscribed_radius = get_inscribed_circle(node, nVertices, xVertex, yVertex, edgeLengths)

                if (inscribed_radius > maximum_inscribed_radius):

                    best_node = node
                    maximum_inscribed_radius = inscribed_radius

                    consecutive_misses = 0

                else:

                    consecutive_misses += 1

        accuracy = min(max_x - min_x, max_y - min_y)

        min_x = best_node[0] - (max_x - min_x) / (math.sqrt(2.0) * 2.0)
        max_x = best_node[0] + (max_x - min_x) / (math.sqrt(2.0) * 2.0)
        min_y = best_node[1] - (max_y - min_y) / (math.sqrt(2.0) * 2.0)
        max_y = best_node[1] + (max_y - min_y) / (math.sqrt(2.0) * 2.0)

    return best_node, maximum_inscribed_radius

#-------------------------------------------------------------------------------

def distance_node_to_edge(node, xVertex1, yVertex1, xVertex2, yVertex2, edgeLength):

    distance = math.fabs((yVertex2-yVertex1)*node[0] - (xVertex2-xVertex1)*node[1] + xVertex2*yVertex1 - yVertex2*xVertex1)

    distance = distance / edgeLength

    return distance

#-------------------------------------------------------------------------------

def get_inscribed_circle(node, nVertices, xVertex, yVertex, edgeLengths):

    if (in_polygon(node[0], node[1], nVertices, xVertex, yVertex)):

        inscribed_radius = 1e30

        for iVertex1 in range(0,nVertices):

            iVertex2 = (iVertex1 + 1) % nVertices

            distance_to_node = distance_node_to_edge(node, xVertex[iVertex1], yVertex[iVertex1], xVertex[iVertex2], yVertex[iVertex2], edgeLengths[iVertex1])

            if (distance_to_node < inscribed_radius):
                inscribed_radius = distance_to_node

    else:

        inscribed_radius = 0.0

    return -inscribed_radius

#-------------------------------------------------------------------------------

def largest_inscribed_circle(nVertices, xVertex, yVertex):

    global timeOptimize

    # precalc edge lengths
    edgeLengths = []
    for iVertex1 in range(0,nVertices):
        iVertex2 = (iVertex1 + 1) % nVertices
        edgeLengths.append(math.sqrt(math.pow((yVertex[iVertex2]-yVertex[iVertex1]),2) + math.pow((xVertex[iVertex2]-xVertex[iVertex1]),2)))
    edgeLengths = np.array(edgeLengths)

    # guess initial as centroid
    node0 = np.array(polygon_centroid(nVertices, xVertex, yVertex))

    startTimeOptimize = time.time()
    res = minimize(get_inscribed_circle, node0, args=(nVertices, xVertex, yVertex, edgeLengths), method='nelder-mead', options={'xtol': 1, 'disp': False})
    endTimeOptimize = time.time()
    timeOptimize += endTimeOptimize - startTimeOptimize

    return res.x, -res.fun

#-------------------------------------------------------------------------------

def pack_discs(filenameOut, voroExecutable, xmin, xmax, ymin, ymax, rmin, rmax, nIterationsCentroidal, nIterationsInscribed):

    global timeVoronoi
    global timeInscribed
    global timeInPolygon

    samples = poisson_disc_sampling(rmax, xmax-2.0*rmax, ymax-2.0*rmax)

    nParticles = len(samples)

    ids = []
    x = []
    y = []
    r = []

    iParticle = 1
    for particle in samples:

        ids.append(iParticle)
        x.append(particle[0]+rmax)
        y.append(particle[1]+rmax)
        r.append(random.uniform(rmin*0.7,rmax*0.7))

        iParticle += 1


    ids = np.array(ids)
    x = np.array(x)
    y = np.array(y)
    r = np.array(r)

    maxRadius = np.amax(r)


    print("nParticles: ", nParticles)



    nParticlesVoronoi = 9 * nParticles

    outfile = open("quality.txt","w",buffering=1)


    startTimeTotal = time.time()

    # lloyds iteration
    nIterations = nIterationsCentroidal + nIterationsInscribed
    for it in range(0,nIterations):

        # replicate domain 9x9 for periodicity
        xVoronoi = []
        yVoronoi = []
        rVoronoi = []
        idsVoronoi = []
        id = 1
        for i in range(-1,2):
            for j in range(-1,2):
                for iParticle in range(0,nParticles):
                    xVoronoi.append(x[iParticle]+float(i)*xmax)
                    yVoronoi.append(y[iParticle]+float(j)*ymax)
                    rVoronoi.append(r[iParticle])
                    idsVoronoi.append(id)
                    id += 1
        xVoronoi = np.array(xVoronoi)
        yVoronoi = np.array(yVoronoi)
        rVoronoi = np.array(rVoronoi)
        idsVoronoi = np.array(idsVoronoi)



        # perform radical Voronoi tessellation
        startTimeVoronoi = time.time()
        maxVertices, cellIDs, nVerticesArray, xVertex, yVertex = radical_voronoi_tessellation(nParticlesVoronoi, idsVoronoi, xVoronoi, yVoronoi, rVoronoi, voroExecutable, -xmax, 2*xmax, -ymax, 2*ymax, -maxRadius, maxRadius)
        endTimeVoronoi = time.time()
        timeVoronoi = endTimeVoronoi - startTimeVoronoi

        # find large inscribed circle per polygon
        x = []
        y = []
        r = []
        rinscribed = []
        rinscribedVoronoi = []

        iParticleVoronoi = 0
        for i in range(-1,2):
            for j in range(-1,2):
                for iParticle in range(0,nParticles):

                    if (it < nIterationsCentroidal):
                        node = polygon_centroid(nVerticesArray[iParticleVoronoi], xVertex[iParticleVoronoi,:], yVertex[iParticleVoronoi,:])
                        radius = 0.0
                    else:
                        startTimeInscribed = time.time()
                        node, radius = largest_inscribed_circle(nVerticesArray[iParticleVoronoi], xVertex[iParticleVoronoi,:], yVertex[iParticleVoronoi,:])
                        endTimeInscribed = time.time()
                        timeInscribed += endTimeInscribed-startTimeInscribed

                    rinscribedVoronoi.append(radius)

                    startTimeInPolygon = time.time()
                    if (not in_polygon(node[0], node[1], nVerticesArray[iParticleVoronoi], xVertex[iParticleVoronoi,:], yVertex[iParticleVoronoi,:])):
                        for iVertex in range(0,nVerticesArray[iParticleVoronoi]):
                            print(xVertex[iParticleVoronoi,iVertex], yVertex[iParticleVoronoi,iVertex])
                        print()
                        print(polygon_centroid(nVerticesArray[iParticleVoronoi], xVertex[iParticleVoronoi,:], yVertex[iParticleVoronoi,:]))
                        print("not in polygon: ", i, j, iParticle, node, radius)
                        sys.exit()
                    endTimeInPolygon = time.time()
                    timeInPolygon += endTimeInPolygon-startTimeInPolygon


                    if (node[0] > xmin and node[0] <= xmax and
                        node[1] > ymin and node[1] <= ymax):

                        x.append(node[0])
                        y.append(node[1])
                        r.append(rVoronoi[iParticleVoronoi])
                        rinscribed.append(radius)

                    iParticleVoronoi += 1

        x = np.array(x)
        y = np.array(y)
        r = np.array(r)
        rinscribed = np.array(rinscribed)
        rinscribedVoronoi = np.array(rinscribedVoronoi)

        endTimeInscribed = time.time()

        nbond = number_of_bonds(xVoronoi, yVoronoi, rinscribedVoronoi, xmin, xmax, ymin, ymax, -100.0)
        fraction, overlap = fraction_covered_area(xVoronoi, yVoronoi, rinscribedVoronoi, xmin, xmax, ymin, ymax)

        print("Iteration: ", it+1, " of: ", nIterations, ", nParticles: ", nParticles, ", new particles: ", x.size)

        if (it >= nIterationsCentroidal):
            outfile.write("%f %f %f\n" %(nbond, fraction, overlap))


    endTimeTotal = time.time()
    timeTotal = endTimeTotal - startTimeTotal

    print("Times: ", timeTotal, timeVoronoi, timeInscribed, timeInPolygon, timeOptimize)

    outfile.close()

    # output
    fileOut = Dataset(filenameOut,"w")

    fileOut.periodic = "YES"
    fileOut.xRange = xmax
    fileOut.yRange = ymax

    fileOut.createDimension("nParticles",nParticles)

    var = fileOut.createVariable("x","d",dimensions=["nParticles"])
    var[:] = x[:]

    var = fileOut.createVariable("y","d",dimensions=["nParticles"])
    var[:] = y[:]

    var = fileOut.createVariable("radius","d",dimensions=["nParticles"])
    var[:] = rinscribed[:]

    fileOut.close()

    plot_voronoi_tessellation_inscribed("particles_lloyds.png", nParticlesVoronoi, maxVertices, nVerticesArray, xVertex, yVertex, xmin, xmax, ymin, ymax, nParticles, x, y, rinscribed)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Create a periodic particle distribution with uniform radius distribution with a variant of Lloyds algorithm.')

    parser.add_argument('-o', dest='filenameOut', required=True, help='Output filename.')
    parser.add_argument('-v', dest='voroExecutable', required=True, help='Location of the voro++ executable.')
    parser.add_argument('--x0', dest='xmin', required=True, type=float, help='x minimum of output domain.')
    parser.add_argument('--x1', dest='xmax', required=True, type=float, help='x maximum of output domain.')
    parser.add_argument('--y0', dest='ymin', required=True, type=float, help='y minimum of output domain.')
    parser.add_argument('--y1', dest='ymax', required=True, type=float, help='y maximum of output domain.')
    parser.add_argument('--r0', dest='rmin', required=True, type=float, help='radius minimum of particle distribution.')
    parser.add_argument('--r1', dest='rmax', required=True, type=float, help='radius maximum of particle distribution.')
    parser.add_argument('--nitc', dest='nIterationsCentroidal', default=10, type=int, help='Number of centroidal iterations.')
    parser.add_argument('--niti', dest='nIterationsInscribed',  default=100, type=int, help='Number of inscribed iterations.')

    args = parser.parse_args()

    pack_discs(args.filenameOut, args.voroExecutable, args.xmin, args.xmax, args.ymin, args.ymax, args.rmin, args.rmax, args.nIterationsCentroidal, args.nIterationsInscribed)
