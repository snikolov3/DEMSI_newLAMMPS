#ifndef DEMSI_KOKKOSPARSER_H_
#define DEMSI_KOKKOSPARSER_H_

#include <mpi.h>
#include <Kokkos_Core.hpp>
#include "demsi_logging.h"

namespace DEMSI {

/*! \class KokkosParser
    \brief Class handling Kokkos command line arguments and returning parameters.
*/
class KokkosParser {

public:

  // Pointer to logging object
  DEMSI::Log* log;

  int num_threads;
  int numa;
  int device;
  int ngpu;

  KokkosParser(int, char **, MPI_Comm, Log* logIn);
  ~KokkosParser() {};

  int getNumberOfThreads() const { return num_threads; }
  int getNuma() const { return numa; }
  int getDeviceID() const { return device; }
  int getNumberOfGPUs() const { return ngpu; }

  Kokkos::InitArguments getKokkosInitArguments() const;

  void logStatus() const;

private:
  KokkosParser() {};                         // prohibit using the default constructor
  KokkosParser(const KokkosParser &) {};     // prohibit using the copy constructor

};

} // namespace DEMSI

#endif
