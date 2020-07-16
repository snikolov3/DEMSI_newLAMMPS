/*
  DEMSI particle io structure

 */

#ifndef DEMSI_PARTICLES_IO_H
#define DEMSI_PARTICLES_IO_H

#include "demsi_particles.h"
#include "demsi_contacts.h"
#include "demsi_logging.h"
#include "demsi_time.h"
#include "demsi_configs.h"
#include "demsi_column.h"
#include "demsi_column_variables.h"

#include <Kokkos_Core.hpp>

/*! \file   demsi_particles_io.h
    \brief  Header file for the DEMSI::ParticlesRead and DEMSI::ParticlesWrite classes
 */

namespace DEMSI {

//------------------------------------------------------------------------------
/*! \class ParticlesColumnRead
    \brief Class to read in column variables.

    This class allows column variable to be read in from a netcdf file.
 */
//------------------------------------------------------------------------------
class ParticlesColumnRead {
public:

  /*! \brief Constructor for the ParticlesColumnRead class.
      \param partitionIn Pointer to the partition object
      \param logIn Pointer to the log object
      \param particlesIn Pointer to the particles object
      \param filenameIn Filename of the netcdf file to read column variables from
   */
  ParticlesColumnRead(DEMSI::Partition* partitionIn, DEMSI::Log* logIn, DEMSI::Particles* particlesIn, const std::string filenameIn);

  /*! \brief Destructor for the ParticlesColumnRead class. */
  ~ParticlesColumnRead();

  /*! \brief Determine if the column variable exists in the netcdf file
      \param varname Name of variable to test existance of
      \return If true column variable exists in the netcdf file
  */
  bool has_variable(const std::string varname);

  /*! \brief Read a column variable from the netcdf file
      \param varname Name of variable to read in
      \param columnVariable Column variable object to read data to
  */
  void get_variable(const std::string varname, DEMSI::ColumnVariable<double>* columnVariable);

private:

  /*! Pointer to the partition object with mpi info. */
  DEMSI::Partition* partition;

  /*! Pointer to the logging object. */
  DEMSI::Log* log;

  /*! Pointer to the Particles object. */
  DEMSI::Particles* particles;

  /*! Filename of input netcdf file */
  std::string filename;

  /*! Netcdf file identifier */
  int ncID;

  /*! Number of particles in input file including type 0 */
  int nParticlesFile;

  /*! Total number of particles across all processors */
  int nParticlesTotal;

  /*! Array of number of particles per processor */
  int *nParticlePerProc;

  /*! Array of particle number dispalcements for MPI calls */
  int *particleDispls;

  /*! Map from global IDs on the file to on processors */
  std::map<int, int> globalIDMap;

};

//------------------------------------------------------------------------------
/*! \class ParticlesWrite
    \brief Writes particle data in a DEMSI:Particles object to a netcdf
    file.

    MPI is used to gather particle data onto the zero processor before writing
    this data to a netcdf file. One file is created per time, since nParticles
    is not fixed. Filename are determined from a filename template and the
    current simulation time.
 */
//------------------------------------------------------------------------------
class ParticlesWrite {
public:

  /*! \brief Constructor for the ParticlesWrite class.
      \param partitionIn Pointer to the partition object with mpi info.
      \param logIn Pointer to the logging object.
      \param contactsIn Pointer to the DEMSI::Contacts object containing bond info.
      \param particlesIn Pointer to the DEMSI::Particles object to write from.
      \param columnIn Pointer to the DEMSI::Column object.
      \param streamNameIn Stream name.
      \param filenameTemplateIn Filename template for particle output file.
      \param clockIn Simulation clock to use to create write alarm.
      \param alarmStartTime Start time for particle write alarm.
      \param alarmInterval Time interval for particle write alarm.
      \param clobberIn If true clobber exisiting files.
      \param outputFieldNamesIn Vector of output column field names.
   */
  ParticlesWrite(DEMSI::Partition* partitionIn, DEMSI::Log* logIn, DEMSI::Contacts* contactsIn, DEMSI::Particles* particlesIn, DEMSI::Column* columnIn, const std::string streamNameIn, const std::string filenameTemplateIn, DEMSI::Clock* clockIn, DEMSI::Time alarmStartTime, DEMSI::TimeInterval alarmInterval, const bool clobberIn, std::vector<std::string> outputFieldNamesIn);

  /*! \brief Destructor for the ParticlesWrite class. */
  ~ParticlesWrite() = default;

  /*! \brief Write the particles data to a netcdf file at the current time. */
  void write(std::string postfix = "");

  /*! \brief Stream name. */
  std::string name(void);

private:

  /*! Pointer to the partition object with mpi info. */
  DEMSI::Partition* partition;

  /*! Pointer to the logging object. */
  DEMSI::Log* log;

  /*! Pointer to the DEMSI::Contacts object. */
  DEMSI::Contacts* contacts;

  /*! Pointer to the DEMSI::Particles object to be written to file. */
  DEMSI::Particles* particles;

  /*! Pointer to the DEMSI::Column object. */
  DEMSI::Column* column;

  /*! Stream name. */
  std::string streamName;

  /*! Filename template for output particle netcdf files. */
  std::string filenameTemplate;

  /*! Alarm to signal writing particle data to a netcdf file. */
  DEMSI::Alarm outputAlarm;

  /*! Pointer to the simulation clock for writing output times. */
  DEMSI::Clock* clock;

  /*! Flag if write is active. */
  bool active;

  /*! If true clobber existing files */
  bool clobber;

  /*! Vector of column output field names for this stream */
  std::vector<std::string> outputFieldNames;

  /*! \brief Write a two dimensional lammps particle variable to a netcdf file.
      \param fieldOut LAMMPS vector particle variable to write.
      \param varname Netcdf output file variable name.
      \param ncID Netcdf file ID for the output file.
      \param nParticlesTotal Total number of particles across all processors.
      \param nParticlesPerProc Array of number of particles per processor.
      \param particleDispls Displacement positions for each processor.
   */
  template <typename kokkos_view_type_2d>
  void write_particle_variable_2d(kokkos_view_type_2d fieldOut, std::string varname, const int ncID, const int nParticlesTotal, const int* nParticlesPerProc, const int* particleDispls);

  /*! \brief Write a one dimensional lammps particle variable to a netcdf file.
      \param fieldOut LAMMPS scalar particle variable to write.
      \param varname Netcdf output file variable name.
      \param ncID Netcdf file ID for the output file.
      \param nParticlesTotal Total number of particles across all processors.
      \param nParticlesPerProc Array of number of particles per processor.
      \param particleDispls Displacement positions for each processor.
   */
  void write_particle_variable(kokkos_view_type_1d_float fieldOut, std::string varname, const int ncID, const int nParticlesTotal, const int* nParticlesPerProc, const int* particleDispls);

  /*! \brief Write a one dimensional lammps globalID particle variable to a netcdf file.
      \param fieldOut LAMMPS globalID particle variable to write.
      \param varname Netcdf output file variable name.
      \param ncID Netcdf file ID for the output file.
      \param nParticlesTotal Total number of particles across all processors.
      \param nParticlesPerProc Array of number of particles per processor.
      \param particleDispls Displacement positions for each processor.
   */
  void write_particle_variable(kokkos_view_type_tagint_host fieldOut, std::string varname, const int ncID, const int nParticlesTotal, const int* nParticlesPerProc, const int* particleDispls);

  /*! \brief Write a double column variable to a netcdf file.
      \param columnVariable Double column variable to write out.
      \param ncID Netcdf file ID for the output file.
      \param nParticlesTotal Total number of particles across all processors.
      \param nParticlesPerProc Array of number of particles per processor.
      \param particleDispls Displacement positions for each processor.
   */
  void write_particle_variable(DEMSI::ColumnVariable<double>* columnVariable, const int ncID, const int nParticlesTotal, const int* nParticlesPerProc, const int* particleDispls);

  /*! \brief Write an int column variable to a netcdf file.
      \param columnVariable Int column variable to write out.
      \param ncID Netcdf file ID for the output file.
      \param nParticlesTotal Total number of particles across all processors.
      \param nParticlesPerProc Array of number of particles per processor.
      \param particleDispls Displacement positions for each processor.
   */
  void write_particle_variable(DEMSI::ColumnVariable<int>* columnVariable, const int ncID, const int nParticlesTotal, const int* nParticlesPerProc, const int* particleDispls);

  /*! \brief Write a vector of size 2 bond global ID arrays to a netcdf file.
      \param fieldOut Vector of bond global ID arrays.
      \param varname Netcdf output file bond variable name.
      \param ncID Netcdf file ID for the output file.
      \param nBondsTotal Total number of bonds across all processors.
      \param nBondsPerProc Array of number of bonds per processor.
      \param bondDispls Displacement positions for each processor.
  */
  void write_bond_variable(Kokkos::View<LAMMPS_NS::tagint*[2]>fieldOut, std::string varname, const int ncID, const int nBondsTotal, const int* nBondsPerProc, const int* bondDispls);

  /*! \brief Write a vector bond doubles to a netcdf file.
      \param fieldOut Vector of bond double variables.
      \param varname Netcdf output file bond variable name.
      \param ncID Netcdf file ID for the output file.
      \param nBondsTotal Total number of bonds across all processors.
      \param nBondsPerProc Array of number of bonds per processor.
      \param bondDispls Displacement positions for each processor.
  */
  void write_bond_variable(Kokkos::View<double*> fieldOut, std::string varname, const int ncID, const int nBondsTotal, const int* nBondsPerProc, const int* bondDispls);

  /*! \brief Write a vector of size 2 bond double arrays to a netcdf file.
      \param fieldOut Vector of bond double arrays.
      \param varname Netcdf output file bond variable name.
      \param ncID Netcdf file ID for the output file.
      \param nBondsTotal Total number of bonds across all processors.
      \param nBondsPerProc Array of number of bonds per processor.
      \param bondDispls Displacement positions for each processor.
  */
  void write_bond_variable(Kokkos::View<double*[2]> fieldOut, std::string varname, const int ncID, const int nBondsTotal, const int* nBondsPerProc, const int* bondDispls);

};

//------------------------------------------------------------------------------
/*! \class ParticlesOutputStreams
    \brief Class containing all the particle output streams and code to set up
    those streams from the configs file.
 */
//------------------------------------------------------------------------------
class ParticlesOutputStreams {
public:

  /*! \brief Constructor for the ParticlesOutputStreams class.
      \param partition Pointer to the partition object with mpi info.
      \param log Pointer to the logging object.
      \param configs Pointer to the configs object.
      \param contacts Pointer to the DEMSI::Contacts object containing bond info.
      \param particles Pointer to the DEMSI::Particles object to write from.
      \param column Pointer to the DEMSI::Column object.
      \param clock Simulation clock to use to create write alarm.
      \param alarmStartTime Start time for particle write alarm.
   */
  ParticlesOutputStreams(DEMSI::Partition* partition, DEMSI::Log* logIn, DEMSI::Configs* configs, DEMSI::Contacts* contacts, DEMSI::Particles* particles, DEMSI::Column* column, DEMSI::Clock* clock, DEMSI::Time alarmStartTime);

  /*! \brief Destructor for the ParticlesOutputStreams class. */
  ~ParticlesOutputStreams() = default;

  /*! \brief Write out all particle output streams */
  void write(void);

  /*! Pointer to the logging object. */
  DEMSI::Log* log;

private:

  /*! vector of pointers to the particle write streams */
  std::vector<DEMSI::ParticlesWrite*> streams;

};

} // namespace DEMSI

#endif
