#include "demsi_partition.h"

#include <mpi.h>

namespace DEMSI {

// Partition class constructor.
Partition::Partition(int *narg, char ***arg) {

  MPI_Init(narg, arg);

  communicator = MPI_COMM_WORLD;

  MPI_Comm_rank(communicator, &thisProc);
  MPI_Comm_size(communicator, &nProcs);

  masterProc = 0;

}

// Returns the processor ID of this processor.
int Partition::this_proc(void) const {
  return thisProc;
}

// Returns the total number of processes in this partition.
int Partition::nprocs(void) const {
  return nProcs;
}

// Returns the process ID of the master process.
int Partition::master_proc(void) const {
  return masterProc;
}

// Returns whether this process is the master process.
bool Partition::on_master_proc(void) const {
  return (thisProc == masterProc);
}

// Returns the MPI communicator for this partition.
MPI_Comm Partition::comm(void) const {
  return communicator;
}

} // namespace DEMSI
