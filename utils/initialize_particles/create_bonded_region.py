#Several ways to import pydemsi, most likely to be useful are:
#1. At the top of the script:
#import sys
#sys.path.append("your-path-to-where-pydemsi.py is")
#2. Add the path where pydemsi.py is to your PYTHONPATH environment variable
#3. Just keep pydemsi.py in the same directory as your importing script


import pydemsi
myregion = pydemsi.Particles.from_lattice(50,20,5000)
myregion.create_bonds()
#Optionally, plot particles with bonds as red lines (note: could do a lot better on bond visualization, e.g. color by damage, thickness, etc, plot from particle edges rather than centers)
#myregion.plot(draw_bonds = True) 

#To manipulate various particle attributes, numpy's 'where'
#function is handy.
#E.g. to set types for all particles with x < 20000 to type 2:
#  myregion.type[np.where(myregion.x < 20000)[0]] = 2
#
#To create a 'hole' of size 20000 in the middle of the region:
#  xc = np.mean(myregion.x)
#  yc = np.mean(myregion.y)
#  remove_indices = np.where(np.sqrt((myregion.x -xc)**2 + (myregion.y - yc)**2) < 20000)[0]
#  myregion.remove_particles(remove_indices)

myregion.to_particle_netcdf("bonded_rectangle_region.nc")



