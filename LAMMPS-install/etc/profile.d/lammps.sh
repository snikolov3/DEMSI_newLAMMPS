# set environment for LAMMPS and msi2lmp executables
# to find potential and force field files
LAMMPS_POTENTIALS=${LAMMPS_POTENTIALS-/ascldap/users/snikolo/DEMSI_cpu3/DEMSI/LAMMPS-install/share/lammps/potentials}
MSI2LMP_LIBRARY=${MSI2LMP_LIBRARY-/ascldap/users/snikolo/DEMSI_cpu3/DEMSI/LAMMPS-install/share/lammps/frc_files}
export LAMMPS_POTENTIALS MSI2LMP_LIBRARY
