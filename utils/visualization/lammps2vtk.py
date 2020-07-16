from  __future__ import print_function, division, absolute_import
import numpy as np
import mmap
import re
import sys

class Particle:
    pass

class Frame:
    pass

class Trajectory:
    pass

def write_frame_vtk(particles, basename, istep):
    f = open(basename+"_"+str(istep)+".vtk", "w")
    f.write("# vtk DataFile Version 2.0\n")
    f.write("Time step "+str(istep)+"\n")
    f.write("ASCII\n")
    f.write("DATASET UNSTRUCTURED_GRID\n")
    f.write("POINTS "+str(len(particles))+" float\n")
    for p in particles:
        f.write(str(p.x)+" "+str(p.y)+" 0\n")
    f.write("POINT_DATA "+str(len(particles))+"\n")
    f.write("SCALARS Radius float 1\n")
    f.write("LOOKUP_TABLE default\n")
    for p in particles:
        f.write(str(p.rad)+" ")
    f.close()

def main():
    if len(sys.argv) != 3:
        print("Usage: dump_process.py <name of dump file> <frequency of skips>")
        return

    fname=sys.argv[1]
    freq=int(sys.argv[2])

    try:
        f=open(sys.argv[1],"r")
    except IOError:
        print("Unable to open file {}".format(sys.argv[1]))
        return
    
    #Output file base name
    basename=sys.argv[1][0:sys.argv[1].find(".dump")]
    
    #Could also define these as user inputs
    startframe = 0
    stopframe = -1

    #Mode: PDB or VTK (currently only VTK implemented)
    mode="VTK"

    #Dump column indices corresponding to x, y, diameter 
    # (can easily add more per-particle values, see p.x,p.y,p.rad line below)
    xcol = 4
    ycol = 5
    rcol = 3

    #Quick and dirty: assumes entire dump file fits into memory. 
    #Eventually do this via buffering scheme
    fmap = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)

    #Get list of positions for starts of all steps
    steps = []
    while(1):
        i=fmap.find(b"ITEM: TIMESTEP")
        if i < 0:
            break
        steps.append(i)
        fmap.seek(i+1,0)

    if stopframe < 0:
        lastframe = len(steps)+stopframe

    for istep,i in enumerate(range(startframe,lastframe)):
        print("Processing step ",istep+1," of ",lastframe-startframe)
        if istep % freq != 0:
            continue
        fmap.seek(steps[i], 0)
        for line in iter(fmap.readline, ""):
            if "ITEM: NUMBER OF ATOMS" in line.decode("utf8"):
                natoms = int(fmap.readline().decode("utf8"))
            if "ITEM: BOX BOUNDS" in line.decode("utf8"):
                xlo, xhi = map(float,fmap.readline().decode("utf8").split(" "))
                ylo, yhi = map(float,fmap.readline().decode("utf8").split(" "))
                zlo, zhi = map(float,fmap.readline().decode("utf8").split(" "))
            if "ITEM: ATOMS" in line.decode("utf8"):
                particles = []
                for j in range(0, natoms):
                    atomline = fmap.readline().decode("utf8")
                    if not atomline:
                        break
                    p=Particle()
                    p.x,p.y,p.rad = np.array(atomline.split(" "))[[xcol,ycol,rcol]].astype(float)
                    particles.append(p)
                write_frame_vtk(particles, basename, istep)
                break
    
   
    fmap.close()

if __name__=="__main__":
    main()
