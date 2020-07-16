#!/usr/bin/env python

import os

testcases = [
    "compression",
    "compression_shear",
    "tension",
    "tension_shear",
    "shear",
    "rotate"]

for testcase in testcases:

    print testcase
    cmd = "rm output/* ; ../../demsi config_bonded_%s.xml ; rm tmp/* ;  ../../utils/visualization/make_particle_movie.py -g grid.nc -f \"./output/particles_out.*\" -b -o two_particles_%s.mp4" %(testcase,testcase)
    print cmd
    os.system(cmd)
