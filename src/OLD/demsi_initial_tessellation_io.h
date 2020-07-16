#ifndef INITIAL_TESSELLATION_IO_H_
#define INITIAL_TESSELLATION_IO_H_

/*! \file   demsi_initial_tessellation_io.h
    \brief  Header file for the initial tessellation io classes.
*/

#include "demsi_partition.h"
#include "demsi_logging.h"
#include "demsi_initial_tessellation.h"
#include "demsi_time.h"

namespace DEMSI {

//------------------------------------------------------------------------------
/*! \class TessellationWrite
    \brief Writes ocean fields to a netcdf file.

    MPI is used to gather ocean data onto the zero processor before writing
    this data to a netcdf file. Filename are determined from a filename
    template and the current simulation time.
 */
//------------------------------------------------------------------------------
class TessellationWrite {
public:

  /*! \brief Constructor for the TessellationWrite class.
      \param partitionIn Pointer to the partition object with mpi info.
      \param logIn Pointer to the logging object.
      \param tessellationIn Pointer to the tessellation object.
      \param filenameTemplateIn Filename template for particle output file.
      \param clockIn Simulation clock to use to create write alarm.
      \param alarmStartTime Start time for particle write alarm.
      \param writeInterval Time interval string for particle write alarm.
      \param clobberIn If true clobber exisiting files.
      \param outputFieldNamesIn Vector of output column field names.
   */
  TessellationWrite(DEMSI::Partition* partitionIn, DEMSI::Log* logIn, DEMSI::Tessellation* tessellation, const std::string streamName, const std::string filenameTemplateIn, DEMSI::Clock* clockIn, DEMSI::Time alarmStartTime, const std::string writeInterval, const bool clobberIn, std::vector<std::string> outputFieldNamesIn);

  /*! \brief Destructor for the TessellationWrite class. */
  ~TessellationWrite() = default;

  /*! \brief Write the ocean data to a netcdf file at the current time. */
  void write(const bool writeNow);

  /*! \brief Name of stream. */
  std::string name(void);

private:

  /*! Pointer to the partition object with mpi info. */
  DEMSI::Partition* partition;

  /*! Pointer to the logging object. */
  DEMSI::Log* log;

  /*! Pointer to the tessellation object. */
  DEMSI::Tessellation* tessellation;

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

  /*! \brief Write a single ocean field to a netcdf file.
      \param field Pointer to the field to write.
      \param nc_id Netcdf file ID to write to.
      \param iTime Time index to write to.
      \param nParticlesInitAllProcs Total number of ocean cells
      \param nParticlesInitOnProcs Number of ocean cells in each processor
      \param displacements displacement vector
  */
  void write_field(std::pair<double*,std::string> field, const int nc_id, const size_t iTime, int nParticlesInitAllProcs, int* nParticlesInitOnProcs, int* displacements);

  /*! \brief Write the globalIndexInit field to a netcdf file.
      \param nc_id Netcdf file ID to write to.
      \param nParticlesInitAllProcs Total number of ocean cells
      \param nParticlesInitOnProcs Number of ocean cells in each processor
      \param displacements displacement vector
  */
  void write_globalIndexInit_field(const int nc_id, int nParticlesInitAllProcs, int* nParticlesInitOnProcs, int* displacements);

};

//------------------------------------------------------------------------------
/*! \class TessellationOutputStreams
    \brief Class containing all the ocean output streams and code to set up
    those streams from the configs file.
 */
//------------------------------------------------------------------------------
class TessellationOutputStreams {
public:

  /*! \brief Constructor for the TessellationOutputStreams class.
      \param partition Pointer to the partition object with mpi info.
      \param log Pointer to the logging object.
      \param configs Pointer to the configs object.
      \param tessellation Pointer to the tessellation object.
      \param clock Simulation clock to use to create write alarm.
      \param alarmStartTime Start time for particle write alarm.
   */
  TessellationOutputStreams(DEMSI::Partition* partition, DEMSI::Log* log, DEMSI::Configs* configs, DEMSI::Tessellation* tessellation, DEMSI::Clock* clock, DEMSI::Time alarmStartTime);

  /*! \brief Destructor for the TessellationOutputStreams class. */
  ~TessellationOutputStreams() = default;

  /*! \brief Write out all ocean output streams */
  void write(void);

  /*! \brief Write out particular stream */
  void write(const std::string streamName, const bool writeNow);

private:

  /*! vector of pointers to the ocean write streams */
  std::vector<DEMSI::TessellationWrite*> streams;

};

} // namespace DEMSI

#endif
