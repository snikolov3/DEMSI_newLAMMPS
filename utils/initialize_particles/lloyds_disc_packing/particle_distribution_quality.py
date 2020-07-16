from netCDF4 import Dataset
import shapely.geometry as shapelyGeometry
import shapely.ops as shapelyOps
import numpy as np
import scipy.spatial
import matplotlib.pyplot as plt
import matplotlib.patches as mplPatches
import matplotlib.collections as mplCollections



# quality metrics

# 1. Fraction of covered area (maximize)
# 2. Area of overlapping elements (minimize)
# 3. Mean number of bonds per element (maximize)
# 4. Mean overlap (minimize)

#-------------------------------------------------------------------------------

def fraction_covered_area(xAll, yAll, rAll, xmin, xmax, ymin, ymax):

    # get elements from subregion
    maxRadius = np.amax(rAll)

    desiredElements = np.where((xAll > xmin - maxRadius) &
                               (xAll < xmax + maxRadius) &
                               (yAll > ymin - maxRadius) &
                               (yAll < ymax + maxRadius),
                               True, False)

    xs = np.extract(desiredElements, xAll)
    ys = np.extract(desiredElements, yAll)
    rs = np.extract(desiredElements, rAll)

    # fraction of area covered
    particles = []
    for x, y, r in zip(xs, ys, rs):
        particles.append(shapelyGeometry.Point((x,y)).buffer(r))

    particlesUnion = shapelyOps.cascaded_union(particles)

    testRegion = shapelyGeometry.Polygon([[xmin, ymin], [xmax, ymin], [xmax, ymax], [xmin, ymax]])

    areaOfParticles = testRegion.intersection(particlesUnion)

    testRegionArea = testRegion.area
    particlesArea = areaOfParticles.area

    # area of all particles
    totalAreaOfParticles = 0.0
    for particle in particles:
        totalAreaOfParticles += testRegion.intersection(particle).area
    areaOfOverlap = totalAreaOfParticles - areaOfParticles.area

    return particlesArea / testRegionArea, areaOfOverlap / testRegionArea

#-------------------------------------------------------------------------------

def plot_region(region):

    fig, ax = plt.subplots(figsize=(8, 8))

    for polygon in region:
        mpl_poly = mplPatches.Polygon(np.array(polygon.exterior), facecolor="g", lw=0, alpha=0.4)
        ax.add_patch(mpl_poly)

    ax.relim()
    ax.autoscale()

    plt.show()

#-------------------------------------------------------------------------------

def number_of_bonds(xAll, yAll, rAll, xmin, xmax, ymin, ymax, minOverlap = 0.0):

    maxRadius = np.amax(rAll)

    desiredElements = np.where((xAll > xmin - 4.0*maxRadius) &
                               (xAll < xmax + 4.0*maxRadius) &
                               (yAll > ymin - 4.0*maxRadius) &
                               (yAll < ymax + 4.0*maxRadius),
                               True, False)

    xs = np.extract(desiredElements, xAll)
    ys = np.extract(desiredElements, yAll)
    rs = np.extract(desiredElements, rAll)

    maxRadius = np.amax(rs)
    points = np.column_stack((xs, ys))

    tree = scipy.spatial.KDTree(points)
    distance_matrix = tree.sparse_distance_matrix(tree, 2.0*maxRadius)
    bonds = [k[0] for k,v in distance_matrix.items()
                 if rs[k[0]] + rs[k[1]] - v >= minOverlap and
                 ((xs[k[0]] > xmin and xs[k[0]] < xmax and ys[k[0]] > ymin and ys[k[0]] < ymax) or
                  (xs[k[1]] > xmin and xs[k[1]] < xmax and ys[k[1]] > ymin and ys[k[1]] < ymax))]

    desiredElements = np.where((xAll > xmin) &
                               (xAll < xmax) &
                               (yAll > ymin) &
                               (yAll < ymax),
                               True, False)

    return float(len(bonds))/float(np.sum(desiredElements))

#-------------------------------------------------------------------------------

def plot_elements(nParticles, x, y, r, xmin, xmax, ymin, ymax):

    fig, axis = plt.subplots()

    patches = []
    for iParticle in range(0,nParticles):

        patches.append(mplPatches.Circle((x[iParticle], y[iParticle]), r[iParticle], fill=False, edgecolor='k', linewidth=0.3))

    pc = mplCollections.PatchCollection(patches, match_original=True)
    axis.add_collection(pc)

    # plot options
    axis.set_xlim((xmin,xmax))
    axis.set_ylim((ymin,ymax))

    axis.set_aspect('equal')

    plt.show()

#-------------------------------------------------------------------------------

def particle_distribution_quality(filenamein, xmin, xmax, ymin, ymax, minOverlap = 0.0):

    filein = Dataset(filenamein,"r")

    nParticles = len(filein.dimensions["nParticles"])

    xAll = filein.variables["x"][:]
    yAll = filein.variables["y"][:]
    rAll = filein.variables["radius"][:]

    filein.close()

    nbond = number_of_bonds(xAll, yAll, rAll, xmin, xmax, ymin, ymax, minOverlap)

    fraction, overlap = fraction_covered_area(xAll, yAll, rAll, xmin, xmax, ymin, ymax)

    return nbond, fraction, overlap

#-------------------------------------------------------------------------------
