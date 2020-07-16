import pydemsi
# Script to modify a lammps data file based on a grayscale image,
# as follows:
#
# 1. Any particles with centers in pixels that have values 
#  0 are removed.
# 2. Any particles with centers in pixels with values 1
#  are flagged as 'immobile'. This is simply done by changing
#  their atom type to a value of 2, which is then used in the 
#  lammps script to keep them immobilized.
#
# To align the image to the lammps data file, the bottom left corner of
# the image is assumed to be at xlo, ylo. The input pixel size is then used to 
# scale the rest of the image.
#
# Usage: set_immobile_particles_image.py <LAMMPS input data file name> 
# <image file name> <size of pixel in x> <size of pixel in y>


if __name__=="__main__":
    if (len(sys.argv) != 5):
        print("Usage: python set_immobile_particles_image.py <LAMMPS input data file name> <image file name> <size of pixel in x> <size of pixel in y>")
        sys.exit()

    basename = sys.argv[1][0:sys.argv[1].find(".data")]
    particles = pydemsi.Particles(sys.argv[1])
    particles.types_from_image(sys.argv[2], types=[-1,2], pixelvalues=[0,1], dx = float(sys.argv[3]), dy = float(sys.argv[4]))
    particles.to_lammps_data(basename+"_image_processed.data")
    particles.to_particle_netcdf(basename+"_image_processed.nc")
    
                             

                          


