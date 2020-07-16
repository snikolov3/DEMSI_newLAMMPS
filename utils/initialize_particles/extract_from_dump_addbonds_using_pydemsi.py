import pydemsi
import sys

def main():
    if len(sys.argv) != 3:
        print("Required arguments are <name of dump file> <time step>")
        return

    dumpname = sys.argv[1]
    dumpstep = int(sys.argv[2])
    myparticles = pydemsi.Particles.from_lammps_dump(dumpname, dumpstep)
    myparticles.create_bonds()

    basename = dumpname[0:dumpname.find(".dump")]
    myparticles.to_lammps_data(basename+".data")
    myparticles.to_particle_netcdf(basename+".nc")


if __name__=="__main__":
    main()
