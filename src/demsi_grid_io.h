#ifndef GRID_IO_H_
#define GRID_IO_H_

/*! \file   demsi_grid_io.h
    \brief  Header file for the DEMSI::GridWrite class.
*/

#include "demsi_logging.h"
#include "demsi_partition.h"
#include "demsi_grid.h"
#include "demsi_time.h"
#include "demsi_configs.h"

#include <Kokkos_Core.hpp>
#include <vector>


namespace DEMSI {

//------------------------------------------------------------------------------
/*! \class GridWrite
    \brief Writes gridded fields to a netcdf file.

    MPI is used to gather gridded data onto the zero processor before writing
    this data to a netcdf file. Filename are determined from a filename
    template and the current simulation time.
 */
//------------------------------------------------------------------------------
class GridWrite {
public:

  /*! \brief Constructor for the GridWrite class.
      \param partitionIn Pointer to the partition object with mpi info.
      \param logIn Pointer to the logging object.
      \param gridIn Pointer to the DEMSI::Grid object.
      \param filenameTemplateIn Filename template for particle output file.
      \param clockIn Simulation clock to use to create write alarm.
      \param alarmStartTime Start time for particle write alarm.
      \param alarmInterval Time interval for particle write alarm.
      \param clobberIn If true clobber exisiting files.
      \param outputFieldNamesIn Vector of output column field names.
   */
  GridWrite(DEMSI::Partition* partitionIn, DEMSI::Log* logIn, DEMSI::Grid* gridIn, const std::string filenameTemplateIn, DEMSI::Clock* clockIn, DEMSI::Time alarmStartTime, DEMSI::TimeInterval alarmInterval, const bool clobberIn, std::vector<std::string> outputFieldNamesIn);

  /*! \brief Destructor for the GridWrite class. */
  ~GridWrite() = default;

  /*! \brief Write the gridded data to a netcdf file at the current time. */
  void write(void);

private:

  /*! Pointer to the partition object with mpi info. */
  DEMSI::Partition* partition;

  /*! Pointer to the logging object. */
  DEMSI::Log* log;

  /*! Pointer to the DEMSI::Grid object. */
  DEMSI::Grid* grid;

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

  /*! \brief Write a single Eulerian field to a netcdf file.
      \param field Pointer to the field to write.
      \param nc_id Netcdf file ID to write to.
      \param iTime Time index to write to.
  */
  void write_field(Kokkos::View<double**>* field, const int nc_id, const size_t iTime);

};

//------------------------------------------------------------------------------
/*! \class GridOutputStreams
    \brief Class containing all the grid output streams and code to set up
    those streams from the configs file.
 */
//------------------------------------------------------------------------------
class GridOutputStreams {
public:

  /*! \brief Constructor for the GridOutputStreams class.
      \param partition Pointer to the partition object with mpi info.
      \param log Pointer to the logging object.
      \param configs Pointer to the configs object.
      \param grid Pointer to the DEMSI::Grid object.
      \param clock Simulation clock to use to create write alarm.
      \param alarmStartTime Start time for particle write alarm.
   */
  GridOutputStreams(DEMSI::Partition* partition, DEMSI::Log* log, DEMSI::Configs* configs, DEMSI::Grid* grid, DEMSI::Clock* clock, DEMSI::Time alarmStartTime);

  /*! \brief Destructor for the GridOutputStreams class. */
  ~GridOutputStreams() = default;

  /*! \brief Write out all grid output streams */
  void write(void);

private:

  /*! vector of pointers to the grid write streams */
  std::vector<DEMSI::GridWrite*> streams;

};

} // namespace DEMSI

#endif
