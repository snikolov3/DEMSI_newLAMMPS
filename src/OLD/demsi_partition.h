/*
  DEMSI partition class

 */

#ifndef DEMSI_PARTITION_H
#define DEMSI_PARTITION_H

#include <mpi.h>

/*! \file   demsi_partition.h
    \brief  Header file for the DEMSI::Partition class
*/

namespace DEMSI {

//------------------------------------------------------------------------------
/*! \class Partition
    \brief Class describing the MPI information on this partition.
*/
//------------------------------------------------------------------------------
class Partition {
public:

  /*! \brief Partition class constructor.
      \param narg Pointer to the main nargs argument.
      \param arg Pointer to the main arg argument.
  */
  Partition(int *narg, char ***arg);

  /*! \brief Partition class default destructor */
  ~Partition() = default;

  /*! \brief Returns the process ID of this process in this partition.
      \return Processor ID of this process.
  */
  int this_proc(void) const;

  /*! \brief Returns the total number of processes in this partition.
      \return Total number of processes in this partition.
  */
  int nprocs(void) const;

  /*! \brief Returns the process ID of the master process.
      \return Process ID of the master process.
  */
  int master_proc(void) const;

  /*! \brief Returns whether this process is the master process.
      \return Whether this process is the master process.
  */
  bool on_master_proc(void) const;

  /*! \brief Returns the MPI communicator for this partition.
      \return MPI communicator for this partition.
  */
  MPI_Comm comm(void) const;

private:

  /*! MPI communicator for this partition. */
  MPI_Comm communicator;

  /*! Process ID for this process. */
  int thisProc;

  /*! Process ID of the master process. */
  int masterProc;

  /*! Number of processes on this partition. */
  int nProcs;

};

} // namespace DEMSI

#endif
