from  __future__ import print_function, division, absolute_import

import numpy as np
import scipy.spatial
import sys
import re
from netCDF4 import Dataset

class Particle:
    pass

class Bond:
    pass

def main():
    if len(sys.argv) != 3:
        print("Usage: dump_process.py <name of dump file> <time step>")
        return

    file_name=sys.argv[1]
    requested_step=int(sys.argv[2])

    found_timestep = False
    particles = []

    #Loop over dump file
    # Look for requested time step
    # Get list of particles at that time step
    with open(file_name, "r") as input_file:
        for line in input_file:
            if "ITEM: TIMESTEP" in line:
                time_step = int(input_file.readline())
                if time_step == requested_step:
                    found_timestep = True
                    if "ITEM: NUMBER OF ATOMS" not in input_file.readline():
                        print("Missing expected 'ITEM: NUMBER OF ATOMS' line in dump file\n")
                        return -1
                    num_atoms = int(input_file.readline())
                    if "ITEM: BOX BOUNDS" not in input_file.readline():
                        print("Missing expected 'ITEM: BOX BOUNDS' line in dump file\n")
                        return -1
                    xlo, xhi = map(float, input_file.readline().split())
                    ylo, yhi = map(float, input_file.readline().split())
                    zlo, zhi = map(float, input_file.readline().split())
                    if "ITEM: ATOMS" not in input_file.readline():
                        print("Missing expected 'ITEM: ATOMS' line in dump file\n")
                        return -1
                    for i in range(0, num_atoms):
                        words = input_file.readline().split()
                        p = Particle()
                        p.id = int(words[0])
                        p.type = int(words[1])
                        p.mass = float(words[2])
                        p.diam = 2*float(words[3])
                        p.x, p.y, p.z = map(float, words[4:7])
                        particles.append(p)

        if not found_timestep:
            print("Could not find requested time step "+str(requested_step)+" in dump file")
            return -1

    print("Read a total of "+str(len(particles))+" particles\n")
    #Output file base name
    basename = sys.argv[1][0:sys.argv[1].find(".dump")]
    outname = basename+"_step"+str(requested_step)

    #KD-tree search for all neighbor distances
    points = np.array([[p.x, p.y] for p in particles])
    radii = np.array([0.5*p.diam for p in particles])
    maxrad = np.max(radii)

    tree = scipy.spatial.KDTree(points)
    distance_matrix = tree.sparse_distance_matrix(tree, 2*maxrad)

    #Form bonds between any contacting particles
    bonds = {(particles[k[0]].id, particles[k[1]].id):v for k,v in distance_matrix.items() if particles[k[0]].diam + particles[k[1]].diam > 2*v and k[0] < k[1]}

    print("Adding a total of "+str(len(bonds))+" bonds\n")

    # creating lammps file
    with open(outname+".data","w") as output_file:
        output_file.write("LAMMPS data file with bonds added. Created from dump file "+str(file_name)+", step "+str(requested_step)+"\n")
        output_file.write(str(num_atoms)+" atoms\n")
        output_file.write(str(len(bonds))+" bonds\n\n")
        output_file.write("1 atom types\n")
        output_file.write("1 bond types\n\n")
        output_file.write(str(xlo)+" "+str(xhi)+" xlo xhi\n")
        output_file.write(str(ylo)+" "+str(yhi)+" ylo yhi\n")
        output_file.write(str(zlo)+" "+str(zhi)+" zlo zhi\n\n")
        output_file.write("Atoms # sphere\n\n")
        for p in particles:
            density = p.mass/(4/3*np.pi*(0.5*p.diam)**3)
            output_file.write(str(p.id)+" "+str(p.type)+" "+str(p.diam)+" "+str(density)+" "+str(p.x)+" "+str(p.y)+" "+str(p.z)+"\n")
        output_file.write("\nBonds\n\n")
        for i,b in enumerate(bonds.keys()):
            output_file.write(str(i+1)+" 1 "+str(b[0])+" "+str(b[1])+"\n")

    # create netcdf file
    fileOut = Dataset(outname+".nc","w",format="NETCDF3_CLASSIC")

    nParticles = len(particles)
    nCategories = 1
    nBonds = len(bonds)
    fileOut.createDimension("nParticles",nParticles)
    fileOut.createDimension("nBonds",nBonds)
    fileOut.createDimension("TWO",2)

    globalID = np.zeros(nParticles,dtype="i")
    x = np.zeros(nParticles)
    y = np.zeros(nParticles)
    radius = np.zeros(nParticles)
    for iParticle in range(0,nParticles):
        globalID[iParticle] = particles[iParticle].id
        x[iParticle] = particles[iParticle].x
        y[iParticle] = particles[iParticle].y
        radius[iParticle] = particles[iParticle].diam / 2.0
    var = fileOut.createVariable("globalID","i",dimensions=["nParticles"])
    var[:] = globalID[:]
    var.units = "-"
    var = fileOut.createVariable("x","d",dimensions=["nParticles"])
    var[:] = x[:]
    var.units = "m"
    var = fileOut.createVariable("y","d",dimensions=["nParticles"])
    var[:] = y[:]
    var.units = "m"
    var = fileOut.createVariable("radius","d",dimensions=["nParticles"])
    var[:] = radius[:]
    var.units = "m"

    bondsNetcdf = np.zeros((nBonds,2),dtype="i")
    for i,b in enumerate(bonds.keys()):
        bondsNetcdf[i,0] = b[0]
        bondsNetcdf[i,1] = b[1]

    var = fileOut.createVariable("bonds","i",dimensions=["nBonds","TWO"])
    var[:] = bondsNetcdf[:]
    var.units = "-"

    fileOut.close()


if __name__=="__main__":
    main()
