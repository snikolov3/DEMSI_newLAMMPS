#include "demsi_kokkosparser.h"
#include "demsi_logging.h"
#include <mpi.h>
#include <Kokkos_Core.hpp>

using namespace DEMSI;

KokkosParser::KokkosParser(int narg, char **arg, MPI_Comm communicator, Log* logIn) {
  // This function reimplements LAMMPS parsing from LAMMPS/src/lammps.cpp 
  // and LAMMPS/src/KOKKOS/kokkos.cpp 
  /* ----------------------------------------------------------------------
     LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
     http://lammps.sandia.gov, Sandia National Laboratories
     Steve Plimpton, sjplimp@sandia.gov
  
     Copyright (2003) Sandia Corporation.  Under the terms of Contract
     DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
     certain rights in this software.  This software is distributed under
     the GNU General Public License.
  
     See the README file in the top-level LAMMPS directory.
  ------------------------------------------------------------------------- */

  this->log = logIn;

  // process any command-line args that invoke Kokkos settings
  
#ifdef KOKKOS_HAVE_CUDA
  // with Cuda, we assume 1 gpu requested, even without a command line argument
  // this means the device with no argument being sent would be 0
  // if a 'g' or 'gpus' flag is given, it will overwrite this default
  ngpu = 1;
#else
  ngpu = 0;
#endif
  device = 0;

  // default number of threads if not specified
  num_threads = 1;
  numa = 1;

  // find kokkos arguments, beginning with second argument given (skips first one, the executable)
  int iarg = 1;
  int kkfirst = 0, kklast = 0;
  bool kokkos_args_found = false;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"-kokkos") == 0 ||
               strcmp(arg[iarg],"-k") == 0) {
      log->check(iarg+2 <= narg, "Invalid command-line argument");//, element_location_str(elementLocationIn));
      iarg += 1;
      // delimit any extra args for the Kokkos instantiation
      kkfirst = iarg;
      while (iarg < narg && arg[iarg][0] != '-') iarg++;
      kklast = iarg;
      kokkos_args_found = true;
      break;
    } else {
      iarg += 1;
    }
  }

  if (kokkos_args_found) { 

    // reset narg and arg just for kokkos arguments
    narg = kklast-kkfirst;
    arg = &arg[kkfirst];


    iarg = 0;
    while (iarg < narg) {
      if (strcmp(arg[iarg],"d") == 0 || strcmp(arg[iarg],"device") == 0) {
        log->check(iarg+2 <= narg, "Invalid Kokkos command-line args");
        device = atoi(arg[iarg+1]);
        iarg += 2;

      } else if (strcmp(arg[iarg],"g") == 0 ||
                 strcmp(arg[iarg],"gpus") == 0) {
#ifndef KOKKOS_HAVE_CUDA
        log->check(false, "GPUs are requested but Kokkos has not been compiled for CUDA");
#endif
        log->check(iarg+2 <= narg, "Invalid Kokkos command-line args");
        ngpu = atoi(arg[iarg+1]);
  
        int skip_gpu = 9999;
        if (iarg+2 < narg && isdigit(arg[iarg+2][0])) {
          skip_gpu = atoi(arg[iarg+2]);
          iarg++;
        }
        iarg += 2;
  
        char *str;
        if ((str = getenv("SLURM_LOCALID"))) {
          int local_rank = atoi(str);
          device = local_rank % ngpu;
          if (device >= skip_gpu) device++;
        }
        if ((str = getenv("MV2_COMM_WORLD_LOCAL_RANK"))) {
          int local_rank = atoi(str);
          device = local_rank % ngpu;
          if (device >= skip_gpu) device++;
        }
        if ((str = getenv("OMPI_COMM_WORLD_LOCAL_RANK"))) {
          int local_rank = atoi(str);
          device = local_rank % ngpu;
          if (device >= skip_gpu) device++;
        }
  
      } else if (strcmp(arg[iarg],"t") == 0 ||
                 strcmp(arg[iarg],"threads") == 0) {
        num_threads = atoi(arg[iarg+1]);
        iarg += 2;
  
      } else if (strcmp(arg[iarg],"n") == 0 ||
                 strcmp(arg[iarg],"numa") == 0) {
        numa = atoi(arg[iarg+1]);
        iarg += 2;
  
      } else log->check(false, "Invalid Kokkos command-line args");
    }
  } // else use defaults since no arguments were given

// This check catches if the user specifically sets gpus to 0
// otherwise, the default is set to 1
#ifdef KOKKOS_HAVE_CUDA
  log->check(ngpu > 0, "Kokkos has been compiled for CUDA but no GPUs are requested"); 
#endif
}

Kokkos::InitArguments KokkosParser::getKokkosInitArguments() const {
  Kokkos::InitArguments args;
  args.num_threads = num_threads;
  args.num_numa = numa;
  args.device_id = device;
  return args;
}

void KokkosParser::logStatus() const {
#ifdef KOKKOS_HAVE_CUDA
  (*log)() << "KOKKOS mode is enabled on GPU with nthreads: " << getNumberOfThreads() << " numa: " << getNuma() << " device_id: " << getDeviceID() << std::endl;
#else
  (*log)() << "KOKKOS mode is enabled on CPU with nthreads: " << getNumberOfThreads() << " numa: " << getNuma() << std::endl;
#endif
}

