from __future__ import print_function, division, absolute_import

# Script to modify a lammps data file based on a grayscale image,
# as follows:
#
# 1. Any LAMMPS particles with centers in pixels that have values 
#  between 1 and 254 (inclusive) are removed.
# 2. Any LAMMPS particles that overlap any pixels with value 0
#  are flagged as 'immobile'. This is simply done by changing
#  their atom type to a value of 2, which is then used in the 
#  lammps script to keep them immobilized.
#
# To align the image to the lammps data file, the top left corner of
# the image is assumed to be at xlo, ylo. The input pixel size is then used to 
# scale the rest of the image.
#
# Usage: set_immobile_particles_image.py <LAMMPS input data file name> 
# <image file name> <size of pixel in x> <size of pixel in y>


import numpy as np
import sys
import scipy.ndimage

class LammpsData:
    pass
class Particle:
    pass
class Bond:
    pass

def read_lmp_data(filename):    
    lmp = LammpsData()
    with open(filename,"r") as f:
        lmp.comment = f.readline()
        words = f.readline().split()
        if words[1] == "atoms":
            lmp.natoms = int(words[0])
        else:
            print("Missing expected \'# atoms\' line in data file!\n")
            return

        words = f.readline().split()
        if len(words) > 1 and words[1] == "bonds":
            lmp.nbonds = int(words[0])
            next(f)
        else:
            lmp.nbonds = 0
        
        words = f.readline().split()
        if words[2] == "types":            
            lmp.n_atom_types = int(words[0])
        else:
            print("Missing expected '\# atom types\' line in data file!\n")
            return

        words = f.readline().split()
        if len(words) > 0:
            if words[1] == "bond" and words[2] == "types":
                lmp.n_bond_types = int(words[0])
                next(f)

        #Probably should add more error checks here..
        lmp.xlo, lmp.xhi = map(float, f.readline().split()[0:2])
        lmp.ylo, lmp.yhi = map(float, f.readline().split()[0:2])
        lmp.zlo, lmp.zhi = map(float, f.readline().split()[0:2])
        for _ in range(3):
            next(f)
        lmp.particles = []
        for _ in range(lmp.natoms):
            words = f.readline().split()
            p = Particle()
            p.id = int(words[0])
            p.type = int(words[1])
            p.diam = float(words[2])
            p.rad = 0.5*p.diam
            p.dens = float(words[3])
            p.x, p.y, p.z = map(float, words[4:7])
            p.remove_flag = False
            lmp.particles.append(p)
        if lmp.nbonds > 0:
            lmp.bonds = []
            for _ in range(0,3):
                next(f)
            for _ in range(lmp.nbonds):
                b = Bond()
                words = f.readline().split()
                b.id = int(words[0])
                b.type = int(words[1])
                b.id1 = int(words[2])
                b.id2 = int(words[3])
                b.remove_flag = False
                lmp.bonds.append(b)
    return lmp

def write_data(filename, lmp, add_comment=""):
    with open(filename, "w") as f:
        f.write(lmp.comment.rstrip()+add_comment.rstrip()+"\n")
        f.write(str(lmp.natoms)+" atoms\n")
        if lmp.nbonds > 0:
            f.write(str(lmp.nbonds)+" bonds\n")
        f.write("\n")
        f.write(str(lmp.n_atom_types)+" atom types\n")
        if lmp.nbonds > 0:
            f.write("1 bond types\n")
        f.write("\n")
        f.write(str(lmp.xlo)+" "+str(lmp.xhi)+" xlo xhi\n")
        f.write(str(lmp.ylo)+" "+str(lmp.yhi)+" ylo yhi\n")
        f.write(str(lmp.zlo)+" "+str(lmp.zhi)+" zlo zhi\n")
        f.write("\nAtoms\n\n")
        for p in lmp.particles:
            f.write(str(p.id)+" "+str(p.type)+" "+str(p.diam)+" "+str(p.dens)+" "+
                    str(p.x)+" "+str(p.y)+" "+str(p.z)+"\n")
        if lmp.nbonds > 0:
            f.write("\nBonds\n\n")
            for i,b in enumerate(lmp.bonds):
                f.write(str(i+1)+" 1 "+str(b.id1)+" "+str(b.id2)+"\n")


if __name__=="__main__":
    #Check arguments
    if (len(sys.argv) != 5):
        print("Usage: python set_immobile_particles_image.py <LAMMPS input data file name> <image file name> <size of pixel in x> <size of pixel in y>")
        sys.exit()

    #Read LAMMPS data file
    try:
        lmp = read_lmp_data(sys.argv[1])
    except:
        print("Error reading lammps data file "+sys.argv[1])
        sys.exit()

    #Read image
    try:
        img = scipy.ndimage.imread(sys.argv[2])
    except:
        print("Could not read image file "+sys.argv[2])
        sys.exit()

    #Parse pixel size
    try:
        dx, dy = map(float, sys.argv[3:5])
    except:
        print("Could not read pixel size")
        sys.exit()

    #Remove any particles outside the image
    for p in lmp.particles:
        xc = int((p.x - lmp.xlo)//dx)
        yc = int((p.y - lmp.ylo)//dy)
        if xc < 0 or yc < 0 or xc >= img.shape[0] or yc >= img.shape[1]:
            p.remove_flag = True

    #Lists of pixel coords corresponding to particle centers
    x_pixels = [int((p.x-lmp.xlo)//dx) for p in lmp.particles if not p.remove_flag]
    y_pixels = [int((p.y-lmp.ylo)//dy) for p in lmp.particles if not p.remove_flag]

    values_at_centers = img[x_pixels, y_pixels]

    #Indices of particles to delete
    remove_indices = np.where((values_at_centers > 0) & (values_at_centers < 255))[0].tolist()

    #Indices of particles to make immobile
    type2_indices = np.where(values_at_centers == 0)[0].tolist()

    for i in type2_indices:
        lmp.particles[i].type = 2

    for i in remove_indices:
        lmp.particles[i].remove_flag = True

    remove_ids = [p.id for p in lmp.particles if p.remove_flag]

    #Also remove bonds that involve any particles that are being removed
    if lmp.nbonds > 0:
        for b in lmp.bonds:
            if (b.id1 in remove_ids) or (b.id2 in remove_ids):
                b.remove_flag = True

    #Update lmp object and write to file
    lmp.particles = [p for p in lmp.particles if not p.remove_flag]
    lmp.natoms = len(lmp.particles)

    if lmp.nbonds > 0:
        lmp.bonds = [b for b in lmp.bonds if not b.remove_flag]
        lmp.nbonds = len(lmp.bonds)

    lmp.n_atom_types = len(np.unique([p.type for p in lmp.particles]))

    basename = sys.argv[1][0:sys.argv[1].find(".data")]
    outname = basename+"_imgcrop.data"
    write_data(outname, lmp)


