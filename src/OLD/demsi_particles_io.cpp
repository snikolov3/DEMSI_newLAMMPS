#include "demsi_particles_io.h"

#include "demsi_particles.h"
#include "demsi_file_utils.h"
#include "demsi_time.h"
#include "demsi_contacts.h"
#include "demsi_configs.h"
#include "demsi_column.h"
#include "demsi_column_variables.h"

#include <netcdf.h>
#include <Kokkos_Core.hpp>
#include <string>
#include <iostream>
#include <fstream>
#include <mpi.h>
#include <vector>
#include <set>
#include <array>

namespace DEMSI {

//------------------------------------------------------------------------------
// class ParticlesColumnRead
//------------------------------------------------------------------------------

// constructor for the particles column read object
ParticlesColumnRead::ParticlesColumnRead(DEMSI::Partition* partitionIn, DEMSI::Log* logIn, DEMSI::Particles* particlesIn, const std::string filenameIn) {

  partition = partitionIn;
  log = logIn;
  particles = particlesIn;

  (*log)(DEMSI::LOG::DEBUG) << "   ...ParticlesColumnRead constructor" << std::endl;

  filename = filenameIn;

  // error handling
  int err;
  char mpiErrBuffer[MPI_MAX_ERROR_STRING];
  int mpiErrLen;

  // nParticles per proc
  if (partition->on_master_proc()) nParticlePerProc = new int[partition->nprocs()];
  int nParticlesSend = *(particles->nParticles);
  err = MPI_Gather(&nParticlesSend, 1, MPI_INT, &nParticlePerProc[0], 1, MPI_INT, partition->master_proc(), partition->comm());
  MPI_Error_string(err, mpiErrBuffer, &mpiErrLen);
  log->check(err == MPI_SUCCESS, "Problem gathering nParticles: ", (std::string) mpiErrBuffer);

  if (partition->on_master_proc()) {

    // total particles
    nParticlesTotal = 0;
    for (int iProc = 0 ; iProc<partition->nprocs() ; iProc++) {
      nParticlesTotal = nParticlesTotal + nParticlePerProc[iProc];
    }

    // particle displacements
    particleDispls = new int[partition->nprocs()];
    particleDispls[0] = 0;
    for (int iProc = 1 ; iProc<partition->nprocs() ; iProc++) {
      particleDispls[iProc] = particleDispls[iProc-1] + nParticlePerProc[iProc-1];
    }
  }

  // particle IDs on this processor
  int particleIDsLocal[std::max(*(particles->nParticles),1)];
  for (int iParticle=0 ; iParticle<*(particles->nParticles) ; iParticle++) {
    particleIDsLocal[iParticle] = (int) particles->globalID[iParticle];
  }

  // particleIDs on all processors
  int* particleIDs;
  if (partition->on_master_proc()) particleIDs = new int[nParticlesTotal];
  err = MPI_Gatherv(&particleIDsLocal[0], *(particles->nParticles), MPI_INT, &particleIDs[0], &nParticlePerProc[0], &particleDispls[0], MPI_INT, partition->master_proc(), partition->comm());
  MPI_Error_string(err, mpiErrBuffer, &mpiErrLen);
  log->check(err == MPI_SUCCESS, "Problem gathering particle ids: ", (std::string) mpiErrBuffer);

  // open the netcdf file
  if (partition->on_master_proc()) {

    err = nc_open(filename.c_str(), 0, &ncID);
    log->check(err == NC_NOERR, "Problem opening file: ", nc_strerror(err), ", for file: ", filename);

    // read total number of particles on the file
    int dimidNParticles;
    err = nc_inq_dimid(ncID, "nParticles", &dimidNParticles);
    log->check(err == NC_NOERR, "Problem getting dimid for nParticles: ", nc_strerror(err), ", for file: ", filename);

    size_t nParticlesFileIn;
    err = nc_inq_dimlen(ncID, dimidNParticles, &nParticlesFileIn);
    log->check(err == NC_NOERR, "Problem getting nParticles: ", nc_strerror(err), ", for file: ", filename);
    nParticlesFile = (int) nParticlesFileIn;

    // read type
    int varidType;
    err = nc_inq_varid(ncID, "type", &varidType);
    log->check(err == NC_NOERR, "Problem getting varid for type: ", nc_strerror(err), ", for file: ", filename);

    int typeIn[nParticlesFile];
    err = nc_get_var_int(ncID, varidType, &typeIn[0]);
    log->check(err == NC_NOERR, "Problem reading data for type: ", nc_strerror(err), ", for file: ", filename);

    // read globalids
    int varidGlobalID;
    err = nc_inq_varid(ncID, "globalID", &varidGlobalID);
    log->check(err == NC_NOERR, "Problem getting varid for globalID: ", nc_strerror(err), ", for file: ", filename);

    int globalID[nParticlesFile];
    err = nc_get_var_int(ncID, varidGlobalID, &globalID[0]);
    log->check(err == NC_NOERR, "Problem reading data for globalID: ", nc_strerror(err), ", for file: ", filename);

    // calculate particle index mapping between file and lammps
    for (int iParticleLocal = 0 ; iParticleLocal < nParticlesTotal ; iParticleLocal++) {
      for (int iParticleFile = 0 ; iParticleFile < nParticlesFile ; iParticleFile++) {
        if (particleIDs[iParticleLocal] == globalID[iParticleFile]) {
          globalIDMap[iParticleLocal] = iParticleFile;
        }
      } // iParticleFile
    } // iParticleLocal

    delete [] particleIDs;

  } // master proc

}

// destructor for the particles column read object
ParticlesColumnRead::~ParticlesColumnRead() {

  (*log)(DEMSI::LOG::DEBUG) << "   ...ParticlesColumnRead destructor" << std::endl;

  if (partition->on_master_proc()) {

    // close the netcdf file
    int err;
    err = nc_close(ncID);
    log->check(err == NC_NOERR, "Problem closing file: ", nc_strerror(err), ", for file: ", filename);

    // deallocate supporting info
    delete [] nParticlePerProc;
    delete [] particleDispls;

  } // master proc

}

// determine if netcdf input file has a particular variable
bool ParticlesColumnRead::has_variable(const std::string varname) {

  (*log)(DEMSI::LOG::DEBUG) << "   ...ParticlesColumnRead has_variable" << std::endl;

  int hasVariable;
  int err;
  char mpiErrBuffer[MPI_MAX_ERROR_STRING];
  int mpiErrLen;

  if (partition->on_master_proc()) {

    int varid;
    err = nc_inq_varid(ncID, varname.c_str(), &varid);
    if (err == NC_ENOTVAR) {
      hasVariable = 0;
    } else if (err == NC_NOERR) {
      hasVariable = 1;
    } else {
      log->abort("Problem finding variable ", varname, ": ", nc_strerror(err), ", for file: ", filename);
    }

  } // master proc

  err = MPI_Bcast(&hasVariable, 1, MPI_INT, partition->master_proc(), partition->comm());
  MPI_Error_string(err, mpiErrBuffer, &mpiErrLen);
  log->check(err == MPI_SUCCESS, "Problem broadcasting hasVariable: ", (std::string) mpiErrBuffer);
  return (bool) hasVariable;

}

// read in a particular column variable from the input netcdf file
void ParticlesColumnRead::get_variable(const std::string varname, DEMSI::ColumnVariable<double>* columnVariable) {

  (*log)(DEMSI::LOG::DEBUG) << "   ...ParticlesColumnRead get_variable" << std::endl;

  // error handling
  int err;
  char mpiErrBuffer[MPI_MAX_ERROR_STRING];
  int mpiErrLen;

  int displacements[partition->nprocs()];
  double* arraySend;

  // check that the input file and variable are compatible
  if (partition->on_master_proc()) {

    // get variable id
    int varid;
    err = nc_inq_varid(ncID, varname.c_str(), &varid);
    log->check(err == NC_NOERR, "Problem get variable id for ", varname, ": ", nc_strerror(err), ", for file: ", filename);

    // get number of dimensions
    int nDims;
    err = nc_inq_varndims(ncID, varid, &nDims);
    log->check(err == NC_NOERR, "Problem getting nDims for variable ", varname, ": ", nc_strerror(err), ", for file: ", filename);

    // get dimension ids
    int dimids[nDims];
    err = nc_inq_vardimid(ncID, varid, &dimids[0]);
    log->check(err == NC_NOERR, "Problem getting dimids for variable ", varname, ": ", nc_strerror(err), ", for file: ", filename);

    // check first dimension is "nParticles"
    char dimName[NC_MAX_NAME+1];
    err = nc_inq_dimname(ncID, dimids[0], &dimName[0]);
    log->check(err == NC_NOERR, "Problem getting first dimname for variable ", varname, ": ", nc_strerror(err), ", for file: ", filename);

    // check other dimension sizes
    for (int iDim = 1 ; iDim < nDims ; iDim++) {
      size_t dimLen;
      err = nc_inq_dimlen(ncID, dimids[iDim], &dimLen);
      log->check(err == NC_NOERR, "Problem getting dimsize for variable ", varname, ": ", nc_strerror(err), ", for file: ", filename);
    }

    // read whole array
    double* arrayRead = new double[nParticlesFile * columnVariable->size_per_particle()];

    err = nc_get_var_double(ncID, varid, &arrayRead[0]);
    log->check(err == NC_NOERR, "Problem reading variable: ", nc_strerror(err), ", for variable: ", varname);

    arraySend = new double[nParticlesTotal * columnVariable->size_per_particle()];
    for (int iParticleSend = 0 ; iParticleSend < nParticlesTotal ; iParticleSend++) {
      for (int j = 0 ; j < columnVariable->size_per_particle() ; j++) {
        int iParticleRead = globalIDMap.find(iParticleSend)->second;
        int iSend = iParticleSend*columnVariable->size_per_particle() + j;
        int iRead = iParticleRead*columnVariable->size_per_particle() + j;
        arraySend[iSend] = arrayRead[iRead];
      } // j
    } // iParticle

    delete [] arrayRead;

    // displacements
    for (int iProc = 0 ; iProc<partition->nprocs() ; iProc++) {
      displacements[iProc] = particleDispls[iProc] * columnVariable->size_per_particle();
    }

  } // master proc

  // scatter
  double* arrayRecv = new double[*(particles->nParticles) * columnVariable->size_per_particle()];
  err = MPI_Scatterv(&arraySend[0], &nParticlePerProc[0], &displacements[0], MPI_DOUBLE, &arrayRecv[0], *(particles->nParticles), MPI_DOUBLE, partition->master_proc(), partition->comm());
  MPI_Error_string(err, mpiErrBuffer, &mpiErrLen);
  log->check(err == MPI_SUCCESS, "Problem with MPI scatter: ", (std::string) mpiErrBuffer, ", for variable: ", varname, " and file: ", filename);

  if (partition->on_master_proc()) {
    delete [] arraySend;
  }

  if (*(particles->nParticles) > 0) {
    double* columnVariablePtr = columnVariable->get(0);
    for (int iParticle = 0 ; iParticle < *(particles->nParticles) ; iParticle++) {
      for (int j = 0 ; j < columnVariable->size_per_particle() ; j++) {
	*columnVariablePtr = arrayRecv[iParticle*columnVariable->size_per_particle() + j];
	columnVariablePtr++;
      } // j
    } // iParticle
  }

  delete [] arrayRecv;

}

//------------------------------------------------------------------------------
// class ParticlesWrite
//------------------------------------------------------------------------------

// Constructor for the ParticlesWrite class.
  ParticlesWrite::ParticlesWrite(DEMSI::Partition* partitionIn, DEMSI::Log* logIn, DEMSI::Contacts* contactsIn, DEMSI::Particles* particlesIn, DEMSI::Column* columnIn, const std::string streamNameIn, const std::string filenameTemplateIn, DEMSI::Clock* clockIn, DEMSI::Time alarmStartTime, DEMSI::TimeInterval alarmInterval, const bool clobberIn, std::vector<std::string> outputFieldNamesIn) {

  partition = partitionIn;
  log = logIn;
  contacts = contactsIn;
  particles = particlesIn;
  column = columnIn;
  streamName = streamNameIn;
  filenameTemplate = filenameTemplateIn;
  clobber = clobberIn;
  outputFieldNames = outputFieldNamesIn;

  // check if filename is none or variant
  std::string filenameTemplateUpper = filenameTemplate;
  std::transform(filenameTemplateUpper.begin(), filenameTemplateUpper.end(), filenameTemplateUpper.begin(), toupper);
  if (filenameTemplateUpper == "NONE") {
    active = false;
  } else {
    active = true;
    outputAlarm.set(clockIn, alarmStartTime, alarmInterval);
  }

  clock = clockIn;

}

// Write a two dimensional lammps particle variable to a netcdf file.
template <typename kokkos_view_type_2d>
void ParticlesWrite::write_particle_variable_2d(kokkos_view_type_2d fieldOut, std::string varname, const int ncID, const int nParticlesTotal, const int* nParticlesPerProc, const int* particleDispls) {

  char mpiErrBuffer[MPI_MAX_ERROR_STRING];
  int mpiErrLen;

  int err;

  int nSend = *(particles->nParticles)*2;
  int nRecv[partition->nprocs()];
  int displacements[partition->nprocs()];
  for (int iProc=0 ; iProc < partition->nprocs() ; iProc++) {
    nRecv[iProc] = nParticlesPerProc[iProc]*2;
    displacements[iProc] = particleDispls[iProc]*2;
  }

  double *arraySend = new double[*(particles->nParticles)*2];
  double *arrayRecv;
  if (partition->on_master_proc()) {
    arrayRecv = new double[nParticlesTotal*2];
  }

  // move data from device to host
  auto fieldOut_mirror = Kokkos::create_mirror(fieldOut);
  Kokkos::deep_copy(fieldOut_mirror, fieldOut);

  for (int iParticle=0 ; iParticle<*(particles->nParticles) ; iParticle++) {
    arraySend[iParticle*2  ] = fieldOut_mirror(iParticle,0);
    arraySend[iParticle*2+1] = fieldOut_mirror(iParticle,1);
  }
  err = MPI_Gatherv(arraySend, nSend, MPI_DOUBLE, arrayRecv, nRecv, displacements, MPI_DOUBLE, partition->master_proc(), partition->comm());
  MPI_Error_string(err, mpiErrBuffer, &mpiErrLen);
  log->check(err == MPI_SUCCESS, "Problem gathering particle array: ", (std::string) mpiErrBuffer, ", for variable: ", varname);

  if (partition->on_master_proc()) {

    int varid;
    err = nc_inq_varid(ncID, varname.c_str(), &varid);
    log->check(err == NC_NOERR, "Problem getting varid: ", nc_strerror(err), ", for variable: ", varname);

    double arrayWrite[nParticlesTotal][2];
    for (int iParticle=0 ; iParticle<nParticlesTotal ; iParticle++) {
      arrayWrite[iParticle][0] = arrayRecv[iParticle*2  ];
      arrayWrite[iParticle][1] = arrayRecv[iParticle*2+1];
    }
    delete [] arrayRecv;

    size_t start[2] = {0, 0};
    size_t count[2] = {(size_t) nParticlesTotal, 2};
    err = nc_put_vara_double(ncID, varid, start, count, &arrayWrite[0][0]);
    log->check(err == NC_NOERR, "Problem writing variable: ", nc_strerror(err), ", for variable: ", varname);

  }

  delete [] arraySend;

}

// Write a one dimensional lammps particle variable to a netcdf file.
void ParticlesWrite::write_particle_variable(kokkos_view_type_1d_float fieldOut, std::string varname, const int ncID, const int nParticlesTotal, const int* nParticlesPerProc, const int* particleDispls) {

  char mpiErrBuffer[MPI_MAX_ERROR_STRING];
  int mpiErrLen;

  int nSend = *(particles->nParticles);
  int nRecv[partition->nprocs()];
  int displacements[partition->nprocs()];
  for (int iProc=0 ; iProc < partition->nprocs() ; iProc++) {
    nRecv[iProc] = nParticlesPerProc[iProc];
    displacements[iProc] = particleDispls[iProc];
  }

  double *arraySend = new double[*(particles->nParticles)];
  double *arrayRecv;
  if (partition->on_master_proc()) {
    arrayRecv = new double[nParticlesTotal];
  }

  // move data from device to host
  auto fieldOut_mirror = Kokkos::create_mirror(fieldOut);
  Kokkos::deep_copy(fieldOut_mirror, fieldOut);

  for (int iParticle=0 ; iParticle<*(particles->nParticles) ; iParticle++) {
    arraySend[iParticle] = fieldOut_mirror[iParticle];
  }

  int err = MPI_Gatherv(arraySend, nSend, MPI_DOUBLE, arrayRecv, nRecv, displacements, MPI_DOUBLE, partition->master_proc(), partition->comm());
  MPI_Error_string(err, mpiErrBuffer, &mpiErrLen);
  log->check(err == MPI_SUCCESS, "Problem gathering particle array: ", (std::string) mpiErrBuffer, ", for variable: ", varname);

  if (partition->on_master_proc()) {

    int varid;
    err = nc_inq_varid(ncID, varname.c_str(), &varid);
    log->check(err == NC_NOERR, "Problem geting varid: ", nc_strerror(err), ", for variable: ", varname);

    size_t start[1] = {0};
    size_t count[1] = {(size_t) nParticlesTotal};
    err = nc_put_vara_double(ncID, varid, start, count, &arrayRecv[0]);
    log->check(err == NC_NOERR, "Problem writing variable: ", nc_strerror(err), ", for variable: ", varname);

    delete [] arrayRecv;

  }

  delete [] arraySend;

}

// Write a one dimensional lammps globalID particle variable to a netcdf file.
void ParticlesWrite::write_particle_variable(kokkos_view_type_tagint_host fieldOut, std::string varname, const int ncID, const int nParticlesTotal, const int* nParticlesPerProc, const int* particleDispls) {

  char mpiErrBuffer[MPI_MAX_ERROR_STRING];
  int mpiErrLen;

  int nSend = *(particles->nParticles);
  int nRecv[partition->nprocs()];
  int displacements[partition->nprocs()];
  for (int iProc=0 ; iProc < partition->nprocs() ; iProc++) {
    nRecv[iProc] = nParticlesPerProc[iProc];
    displacements[iProc] = particleDispls[iProc];
  }

  int *arraySend = new int[*(particles->nParticles)];
  int *arrayRecv;
  if (partition->on_master_proc()) {
    arrayRecv = new int[nParticlesTotal];
  }

  for (int iParticle=0 ; iParticle<*(particles->nParticles) ; iParticle++) {
    arraySend[iParticle] = (int) fieldOut[iParticle];
  }

  int err = MPI_Gatherv(arraySend, nSend, MPI_INT, arrayRecv, nRecv, displacements, MPI_INT, partition->master_proc(), partition->comm());
  MPI_Error_string(err, mpiErrBuffer, &mpiErrLen);
  log->check(err == MPI_SUCCESS, "Problem gathering particle array: ", (std::string) mpiErrBuffer, ", for variable: ", varname);

  if (partition->on_master_proc()) {

    int varid;
    err = nc_inq_varid(ncID, varname.c_str(), &varid);
    log->check(err == NC_NOERR, "Problem geting varid: ", nc_strerror(err), ", for variable: ", varname);

    size_t start[1] = {0};
    size_t count[1] = {(size_t) nParticlesTotal};
    err = nc_put_vara_int(ncID, varid, start, count, &arrayRecv[0]);
    log->check(err == NC_NOERR, "Problem writing variable: ", nc_strerror(err), ", for variable: ", varname);

    delete [] arrayRecv;

  }

  delete [] arraySend;

}

// Write a double column variable to a netcdf file.
void ParticlesWrite::write_particle_variable(DEMSI::ColumnVariable<double>* columnVariable, const int ncID, const int nParticlesTotal, const int* nParticlesPerProc, const int* particleDispls) {

  char mpiErrBuffer[MPI_MAX_ERROR_STRING];
  int mpiErrLen;

  int nSend = *(particles->nParticles) * columnVariable->size_per_particle();
  int nRecv[partition->nprocs()];
  int displacements[partition->nprocs()];
  for (int iProc=0 ; iProc < partition->nprocs() ; iProc++) {
    nRecv[iProc] = nParticlesPerProc[iProc] * columnVariable->size_per_particle();
    displacements[iProc] = particleDispls[iProc] * columnVariable->size_per_particle();
  }

  double *arraySend = new double[*(particles->nParticles) * columnVariable->size_per_particle()];
  double *arrayRecv;
  if (partition->on_master_proc()) {
    arrayRecv = new double[nParticlesTotal * columnVariable->size_per_particle()];
  }

  double* columnVariablePtr = columnVariable->get(0);
  for (int iParticle=0 ; iParticle < *(particles->nParticles) ; iParticle++) {
    for (int j = 0 ; j < columnVariable->size_per_particle() ; j++) {
      arraySend[iParticle*columnVariable->size_per_particle() + j] = *columnVariablePtr;
      columnVariablePtr++;
    }
  }

  std::string varname = columnVariable->name();
  int err = MPI_Gatherv(arraySend, nSend, MPI_DOUBLE, arrayRecv, nRecv, displacements, MPI_DOUBLE, partition->master_proc(), partition->comm());
  MPI_Error_string(err, mpiErrBuffer, &mpiErrLen);
  log->check(err == MPI_SUCCESS, "Problem gathering particle array: ", (std::string) mpiErrBuffer, ", for variable: ", varname);

  if (partition->on_master_proc()) {

    int varid;
    err = nc_inq_varid(ncID, varname.c_str(), &varid);
    log->check(err == NC_NOERR, "Problem geting varid: ", nc_strerror(err), ", for variable: ", varname);

    std::vector<int> dimSizes = columnVariable->get_dimension_sizes();
    size_t start[dimSizes.size()+1];
    size_t count[dimSizes.size()+1];
    start[0] = 0;
    count[0] = (size_t) nParticlesTotal;
    for (int iDim = 0 ; iDim < dimSizes.size() ; iDim++) {
      start[iDim+1] = 0;
      count[iDim+1] = (size_t) dimSizes[iDim];
    }
    err = nc_put_vara_double(ncID, varid, start, count, &arrayRecv[0]);
    log->check(err == NC_NOERR, "Problem writing variable: ", nc_strerror(err), ", for variable: ", varname);

    delete [] arrayRecv;

  }

  delete [] arraySend;

}

// Write an int column variable to a netcdf file.
void ParticlesWrite::write_particle_variable(DEMSI::ColumnVariable<int>* columnVariable, const int ncID, const int nParticlesTotal, const int* nParticlesPerProc, const int* particleDispls) {

  char mpiErrBuffer[MPI_MAX_ERROR_STRING];
  int mpiErrLen;

  int nSend = *(particles->nParticles) * columnVariable->size_per_particle();
  int nRecv[partition->nprocs()];
  int displacements[partition->nprocs()];
  for (int iProc=0 ; iProc < partition->nprocs() ; iProc++) {
    nRecv[iProc] = nParticlesPerProc[iProc] * columnVariable->size_per_particle();
    displacements[iProc] = particleDispls[iProc] * columnVariable->size_per_particle();
  }

  int *arraySend = new int[*(particles->nParticles) * columnVariable->size_per_particle()];
  int *arrayRecv;
  if (partition->on_master_proc()) {
    arrayRecv = new int[nParticlesTotal * columnVariable->size_per_particle()];
  }

  int* columnVariablePtr = columnVariable->get(0);
  for (int iParticle=0 ; iParticle<*(particles->nParticles) ; iParticle++) {
    for (int j = 0 ; j < columnVariable->size_per_particle() ; j++) {
      arraySend[iParticle*columnVariable->size_per_particle() + j] = *columnVariablePtr;
      columnVariablePtr++;
    }
  }

  std::string varname = columnVariable->name();
  int err = MPI_Gatherv(arraySend, nSend, MPI_INT, arrayRecv, nRecv, displacements, MPI_INT, partition->master_proc(), partition->comm());
  MPI_Error_string(err, mpiErrBuffer, &mpiErrLen);
  log->check(err == MPI_SUCCESS, "Problem gathering particle array: ", (std::string) mpiErrBuffer, ", for variable: ", varname);

  if (partition->on_master_proc()) {

    int varid;
    err = nc_inq_varid(ncID, varname.c_str(), &varid);
    log->check(err == NC_NOERR, "Problem geting varid: ", nc_strerror(err), ", for variable: ", varname);

    size_t start[1] = {0};
    size_t count[1] = {(size_t) nParticlesTotal*columnVariable->size_per_particle()};
    err = nc_put_vara_int(ncID, varid, start, count, &arrayRecv[0]);
    log->check(err == NC_NOERR, "Problem writing variable: ", nc_strerror(err), ", for variable: ", varname);

    delete [] arrayRecv;

  }

  delete [] arraySend;

}

// Write bond variable for globalIDs
void ParticlesWrite::write_bond_variable(Kokkos::View<LAMMPS_NS::tagint*[2]>fieldOut, std::string varname, const int ncID, const int nBondsTotal, const int* nBondsPerProc, const int* bondDispls) {

  char mpiErrBuffer[MPI_MAX_ERROR_STRING];
  int mpiErrLen;

  int nBonds = fieldOut.extent(0);
  int nSend = nBonds * 2;

  // move data from device to host
  auto fieldOut_mirror = Kokkos::create_mirror(fieldOut);
  Kokkos::deep_copy(fieldOut_mirror, fieldOut);

  int* arraySend = new int[2*nBonds];

  for (int iBond=0 ; iBond<nBonds ; iBond++) {
    arraySend[iBond*2]   = fieldOut_mirror(iBond,0);
    arraySend[iBond*2+1] = fieldOut_mirror(iBond,1);
  }

  int* nRecv;
  int* displacements;
  int* arrayRecv;
  if (partition->on_master_proc() and nBondsTotal > 0) {
    nRecv = new int[partition->nprocs()];
    displacements = new int[partition->nprocs()];
    for (int iProc=0 ; iProc < partition->nprocs() ; iProc++) {
      nRecv        [iProc] = nBondsPerProc[iProc] * 2;
      displacements[iProc] = bondDispls   [iProc] * 2;
    } // iProc
    arrayRecv = new int[2*nBondsTotal];
  }

  int err = MPI_Gatherv(arraySend, nSend, MPI_INT, arrayRecv, nRecv, displacements, MPI_INT, partition->master_proc(), partition->comm());
  MPI_Error_string(err, mpiErrBuffer, &mpiErrLen);
  log->check(err == MPI_SUCCESS, "Problem gathering bond array: ", (std::string) mpiErrBuffer, ", for variable: ", varname);

  delete [] arraySend;

  if (partition->on_master_proc()) {

    int arrayWrite[nBondsTotal][2];

    for (int iBond=0 ; iBond<nBondsTotal ; iBond++) {
      arrayWrite[iBond][0] = arrayRecv[iBond*2];
      arrayWrite[iBond][1] = arrayRecv[iBond*2+1];
    } // iBond

    delete [] arrayRecv;
    delete [] nRecv;
    delete [] displacements;

    int varid;
    err = nc_inq_varid(ncID, varname.c_str(), &varid);
    log->check(err == NC_NOERR, "Problem geting varid: ", nc_strerror(err), ", for variable: ", varname);

    size_t start[2] = {0, 0};
    size_t count[2] = {(size_t) nBondsTotal, (size_t) 2};
    err = nc_put_vara_int(ncID, varid, start, count, &arrayWrite[0][0]);
    log->check(err == NC_NOERR, "Problem writing variable: ", nc_strerror(err), ", for variable: ", varname);

  }

}

// Write bond variable for double
void ParticlesWrite::write_bond_variable(Kokkos::View<double*> fieldOut, std::string varname, const int ncID, const int nBondsTotal, const int* nBondsPerProc, const int* bondDispls) {

  char mpiErrBuffer[MPI_MAX_ERROR_STRING];
  int mpiErrLen;

  int nBonds = fieldOut.extent(0);
  int nSend = nBonds;

  // move data from device to host
  auto fieldOut_mirror = Kokkos::create_mirror(fieldOut);
  Kokkos::deep_copy(fieldOut_mirror, fieldOut);

  double* arraySend = new double[nBonds];

  for (int iBond=0 ; iBond<nBonds ; iBond++) {
    arraySend[iBond] = fieldOut_mirror[iBond];
  } // iBond

  int* nRecv;
  int* displacements;
  double* arrayRecv;
  if (partition->on_master_proc() and nBondsTotal > 0) {
    nRecv = new int[partition->nprocs()];
    displacements = new int[partition->nprocs()];
    for (int iProc=0 ; iProc < partition->nprocs() ; iProc++) {
      nRecv        [iProc] = nBondsPerProc[iProc];
      displacements[iProc] = bondDispls   [iProc];
    } // iProc
    arrayRecv = new double[nBondsTotal];
  }

  int err = MPI_Gatherv(arraySend, nSend, MPI_DOUBLE, arrayRecv, nRecv, displacements, MPI_DOUBLE, partition->master_proc(), partition->comm());
  MPI_Error_string(err, mpiErrBuffer, &mpiErrLen);
  log->check(err == MPI_SUCCESS, "Problem gathering bond array: ", (std::string) mpiErrBuffer, ", for variable: ", varname);

  delete [] arraySend;

  if (partition->on_master_proc()) {

    int varid;
    err = nc_inq_varid(ncID, varname.c_str(), &varid);
    log->check(err == NC_NOERR, "Problem geting varid: ", nc_strerror(err), ", for variable: ", varname);

    size_t start[1] = {0};
    size_t count[1] = {(size_t) nBondsTotal};
    err = nc_put_vara_double(ncID, varid, start, count, &arrayRecv[0]);
    log->check(err == NC_NOERR, "Problem writing variable: ", nc_strerror(err), ", for variable: ", varname);

    delete [] arrayRecv;
    delete [] nRecv;
    delete [] displacements;

  }

}

// Write bond variable for double array size 2
void ParticlesWrite::write_bond_variable(Kokkos::View<double*[2]> fieldOut, std::string varname, const int ncID, const int nBondsTotal, const int* nBondsPerProc, const int* bondDispls) {

  char mpiErrBuffer[MPI_MAX_ERROR_STRING];
  int mpiErrLen;

  int nBonds = fieldOut.extent(0);
  int nSend = nBonds * 2;

  // move data from device to host
  auto fieldOut_mirror = Kokkos::create_mirror(fieldOut);
  Kokkos::deep_copy(fieldOut_mirror, fieldOut);

  double* arraySend = new double[2*nBonds];

  for (int iBond=0 ; iBond<nBonds ; iBond++) {
    arraySend[iBond*2]   = fieldOut_mirror(iBond,0);
    arraySend[iBond*2+1] = fieldOut_mirror(iBond,1);
  } // iBond

  int* nRecv;
  int* displacements;
  double* arrayRecv;
  if (partition->on_master_proc() and nBondsTotal > 0) {
    nRecv = new int[partition->nprocs()];
    displacements = new int[partition->nprocs()];
    for (int iProc=0 ; iProc < partition->nprocs() ; iProc++) {
      nRecv        [iProc] = nBondsPerProc[iProc] * 2;
      displacements[iProc] = bondDispls   [iProc] * 2;
    } // iProc
    arrayRecv = new double[2*nBondsTotal];
  }

  int err = MPI_Gatherv(arraySend, nSend, MPI_DOUBLE, arrayRecv, nRecv, displacements, MPI_DOUBLE, partition->master_proc(), partition->comm());
  MPI_Error_string(err, mpiErrBuffer, &mpiErrLen);
  log->check(err == MPI_SUCCESS, "Problem gathering bond array: ", (std::string) mpiErrBuffer, ", for variable: ", varname);

  delete [] arraySend;

  if (partition->on_master_proc()) {

    double arrayWrite[nBondsTotal][2];

    for (int iBond=0 ; iBond<nBondsTotal ; iBond++) {
      arrayWrite[iBond][0] = arrayRecv[iBond*2];
      arrayWrite[iBond][1] = arrayRecv[iBond*2+1];
    } // iBond

    delete [] arrayRecv;
    delete [] nRecv;
    delete [] displacements;

    int varid;
    err = nc_inq_varid(ncID, varname.c_str(), &varid);
    log->check(err == NC_NOERR, "Problem geting varid: ", nc_strerror(err), ", for variable: ", varname);

    size_t start[2] = {0, 0};
    size_t count[2] = {(size_t) nBondsTotal, (size_t) 2};
    err = nc_put_vara_double(ncID, varid, start, count, &arrayWrite[0][0]);
    log->check(err == NC_NOERR, "Problem writing variable: ", nc_strerror(err), ", for variable: ", varname);

  }

}

// Write the particles data to a netcdf file at the current time.
void ParticlesWrite::write(std::string postfix) {

  // get the current time from the simulation clock
  Time currentTime = clock->get_time();

  // check if we should write an output
  if (active and outputAlarm.is_ringing()) {

    (*log)(DEMSI::LOG::DEBUG) << "ParticlesWrite::write alarm ringing" << std::endl;

    int err, ncID;
    std::string filenameOut;

    char mpiErrBuffer[MPI_MAX_ERROR_STRING];
    int mpiErrLen;

    int dimidTWO;

    // get transfer info
    int nParticlesPerProc[partition->nprocs()];
    int particleDispls[partition->nprocs()];
    int nParticlesTotal = 0;

    int ierr;
    int nParticlesGather = *(particles->nParticles);
    ierr = MPI_Gather(&nParticlesGather, 1, MPI_INT, &nParticlesPerProc, 1, MPI_INT, partition->master_proc(), partition->comm());
    MPI_Error_string(ierr, mpiErrBuffer, &mpiErrLen);
    log->check(ierr == MPI_SUCCESS, "Problem gathering nParticles: ", (std::string) mpiErrBuffer);

    if (partition->on_master_proc()) {

      // total particles
      for (int iProc = 0 ; iProc<partition->nprocs() ; iProc++) {
        nParticlesTotal = nParticlesTotal + nParticlesPerProc[iProc];
      }

      // particle displacements
      particleDispls[0] = 0;
      for (int iProc = 1 ; iProc<partition->nprocs() ; iProc++) {
        particleDispls[iProc] = particleDispls[iProc-1] + nParticlesPerProc[iProc-1];
      }

      int dimidTime, dimidStrLen, dimidParticles, dimidCategories;
      int varidTime, varidGlobalID, varidType, varidX, varidV, varidForcing;
      int varidRadius, varidIceThickness, varidIceFraction;

      // create filename
      filenameOut = expand_filename_template(filenameTemplate, currentTime);
      if (postfix.length() > 0) {
        filenameOut = filenameOut + "_" + postfix;
      }

      // attempt to create
      int cmode;
      if (clobber) {
        cmode = NC_CLOBBER;
      } else {
        cmode = NC_NOCLOBBER;
      }
      err = nc_create(filenameOut.c_str(), cmode, &ncID);
      log->check(err == NC_NOERR, "Problem creating particle output file: ", nc_strerror(err), ", for file: ", filenameOut);

      // define dimensions
      err = nc_def_dim(ncID, "strLen", 64, &dimidStrLen);
      log->check(err == NC_NOERR, "Problem defining strLen dimension: ", nc_strerror(err), ", for file: ", filenameOut);

      err = nc_def_dim(ncID, "nParticles", nParticlesTotal, &dimidParticles);
      log->check(err == NC_NOERR, "Problem defining nParticles dimension: ", nc_strerror(err), ", for file: ", filenameOut);

      err = nc_def_dim(ncID, "TWO", 2, &dimidTWO);
      log->check(err == NC_NOERR, "Problem defining TWO dimension: ", nc_strerror(err), ", for file: ", filenameOut);

      int nCategories = column->columnDimensions->size("nCategories");
      err = nc_def_dim(ncID, "nCategories", nCategories, &dimidCategories);
      log->check(err == NC_NOERR, "Problem defining nCategories dimension: ", nc_strerror(err), ", for file: ", filenameOut);

      // column dimensions
      for (int iField = 0 ; iField < outputFieldNames.size() ; iField++) {
        if (column->columnVariables->exists(outputFieldNames[iField])) {
          std::vector<std::string> dimNames = column->columnVariables->dimension_names(outputFieldNames[iField]);
          std::vector<int>         dimSizes = column->columnVariables->dimension_sizes(outputFieldNames[iField]);
          for (int iDim = 0 ; iDim < dimNames.size() ; iDim++) {
            if (dimNames[iDim] != "nCategories") { // remove once iceareacategory pulled into column
              int dimid;
              err = nc_def_dim(ncID, dimNames[iDim].c_str(), dimSizes[iDim], &dimid);
              log->check(err == NC_NOERR or err == NC_ENAMEINUSE, "Problem defining nCategories dimension: ", dimNames[iDim], ", :", nc_strerror(err), ", for file: ", filenameOut);
            }
          } // variable dimensions
        } // field exists
      } // field loop

      // define variables
      int dimidsTime[1] = {dimidStrLen};
      err = nc_def_var(ncID, "Time", NC_CHAR, 1, dimidsTime, &varidTime);
      log->check(err == NC_NOERR, "Problem defining Time variable: ", nc_strerror(err), ", for file: ", filenameOut);

      int dimidsGlobalID[1] = {dimidParticles};
      err = nc_def_var(ncID, "globalID", NC_INT, 1, dimidsGlobalID, &varidGlobalID);
      log->check(err == NC_NOERR, "Problem defining globalID variable: ", nc_strerror(err), ", for file: ", filenameOut);

      int dimidsType[1] = {dimidParticles};
      err = nc_def_var(ncID, "type", NC_INT, 1, dimidsType, &varidType);
      log->check(err == NC_NOERR, "Problem defining type variable: ", nc_strerror(err), ", for file: ", filenameOut);

      int dimidsX[2] = {dimidParticles,dimidTWO};
      err = nc_def_var(ncID, "x", NC_DOUBLE, 2, dimidsX, &varidX);
      log->check(err == NC_NOERR, "Problem defining x variable: ", nc_strerror(err), ", for file: ", filenameOut);

      int dimidsV[2] = {dimidParticles,dimidTWO};
      err = nc_def_var(ncID, "v", NC_DOUBLE, 2, dimidsV, &varidV);
      log->check(err == NC_NOERR, "Problem defining v variable: ", nc_strerror(err), ", for file: ", filenameOut);

      int dimidsForcing[2] = {dimidParticles,dimidTWO};
      err = nc_def_var(ncID, "forcing", NC_DOUBLE, 2, dimidsForcing, &varidForcing);
      log->check(err == NC_NOERR, "Problem defining forcing variable: ", nc_strerror(err), ", for file: ", filenameOut);

      int dimidsRadius[1] = {dimidParticles};
      err = nc_def_var(ncID, "radius", NC_DOUBLE, 1, dimidsRadius, &varidRadius);
      log->check(err == NC_NOERR, "Problem defining radius variable: ", nc_strerror(err), ", for file: ", filenameOut);

      // double column variables
      for (int iField = 0 ; iField < outputFieldNames.size() ; iField++) {
        if (column->columnVariables->exists(outputFieldNames[iField])) {
          std::vector<std::string> dimNames = column->columnVariables->dimension_names(outputFieldNames[iField]);
          std::vector<int> dimids;
          int dimid;
          err = nc_inq_dimid(ncID, "nParticles", &dimid);
          log->check(err == NC_NOERR, "Problem getting dimid for nParticles dimension: ", nc_strerror(err), ", for file: ", filenameOut);
          dimids.push_back(dimid);
          for (int iDimension = 0 ; iDimension < dimNames.size() ; iDimension++) {

            err = nc_inq_dimid(ncID, dimNames[iDimension].c_str(), &dimid);
            log->check(err == NC_NOERR, "Problem getting dimid for ", dimNames[iDimension], " dimension: ", nc_strerror(err), ", for file: ", filenameOut);
            dimids.push_back(dimid);
          }
          int varid;
          if (column->columnVariables->type(outputFieldNames[iField]) == DEMSI::FIELDTYPE::DOUBLE_FIELD) {
            err = nc_def_var(ncID, outputFieldNames[iField].c_str(), NC_DOUBLE, dimids.size(), &dimids[0], &varid);
          } else {
            err = nc_def_var(ncID, outputFieldNames[iField].c_str(), NC_INT, dimids.size(), &dimids[0], &varid);
          }
          log->check(err == NC_NOERR, "Problem defining ", outputFieldNames[iField], " variable: ", nc_strerror(err), ", for file: ", filenameOut);
        }
      }

    } // on_master_proc

    // get bonds info
    Kokkos::View<LAMMPS_NS::tagint*[2]> bondGlobalIDs;
    Kokkos::View<double*> bondLength;
    Kokkos::View<double*> bondThickness;
    Kokkos::View<double*[2]> bondCrackFraction;
    Kokkos::View<double*[2]> bondEndPoint1Particle1;
    Kokkos::View<double*[2]> bondEndPoint2Particle1;
    Kokkos::View<double*[2]> bondEndPoint1Particle2;
    Kokkos::View<double*[2]> bondEndPoint2Particle2;

    // bond transfer info
    int nBondsPerProc[partition->nprocs()];
    int bondDispls[partition->nprocs()];
    int nBondsTotal = 0;

    if (contacts->bonds_active()) {

      if (contacts->get_contact_type() == DEMSI::CONTACTTYPE::HOPKINS) {
	    contacts->get_bonds_info_hopkins(bondGlobalIDs, bondLength, bondThickness, bondCrackFraction, 
                bondEndPoint1Particle1, bondEndPoint2Particle1, bondEndPoint1Particle2, bondEndPoint2Particle2);
      }
      int nBondsGather = bondGlobalIDs.extent(0);

      // get the number of bonds on each processor on the master task
      ierr = MPI_Gather(&nBondsGather, 1, MPI_INT, &nBondsPerProc, 1, MPI_INT, partition->master_proc(), partition->comm());
      MPI_Error_string(ierr, mpiErrBuffer, &mpiErrLen);
      log->check(ierr == MPI_SUCCESS, "Problem gathering nBonds: ", (std::string) mpiErrBuffer);

      int dimidBonds;
      int varidBondGlobalIDs, varidBondCrackFraction, varidBondLength, varidBondThickness;
      int varidBondEndPoint1Particle1, varidBondEndPoint2Particle1, varidBondEndPoint1Particle2, varidBondEndPoint2Particle2;

      // total bonds
      if (partition->on_master_proc()) {
        for (int iProc = 0 ; iProc<partition->nprocs() ; iProc++) {
          nBondsTotal += nBondsPerProc[iProc];
        }
      }
      ierr = MPI_Bcast(&nBondsTotal, 1, MPI_INT, partition->master_proc(), partition->comm());
      MPI_Error_string(ierr, mpiErrBuffer, &mpiErrLen);
      log->check(ierr == MPI_SUCCESS, "Problem broadcasting nBondsTotal: ", (std::string) mpiErrBuffer);

      if (partition->on_master_proc()) {

        // particle displacements
        bondDispls[0] = 0;
        for (int iProc = 1 ; iProc<partition->nprocs() ; iProc++) {
          bondDispls[iProc] = bondDispls[iProc-1] + nBondsPerProc[iProc-1];
        }

        // check if we need to output bonds
        if (nBondsTotal > 0) {

          // define bond dimensions
          err = nc_def_dim(ncID, "nBonds", nBondsTotal, &dimidBonds);
          log->check(err == NC_NOERR, "Problem defining nBonds dimension: ", nc_strerror(err), ", for file: ", filenameOut);

          // define variables
          int dimidsBondGlobalIDs[2] = {dimidBonds, dimidTWO};
          err = nc_def_var(ncID, "bondGlobalIDs", NC_INT, 2, dimidsBondGlobalIDs, &varidBondGlobalIDs);
          log->check(err == NC_NOERR, "Problem defining bondGlobalIDs variable: ", nc_strerror(err), ", for file: ", filenameOut);

          if (contacts->get_contact_type() == DEMSI::CONTACTTYPE::HOPKINS) {

            int dimidsBondCrackFraction[2] = {dimidBonds, dimidTWO};
            err = nc_def_var(ncID, "bondCrackFraction", NC_DOUBLE, 2, dimidsBondCrackFraction, &varidBondCrackFraction);
            log->check(err == NC_NOERR, "Problem defining bondCrackFraction variable: ", nc_strerror(err), ", for file: ", filenameOut);

            int dimidsBondLength[1] = {dimidBonds};
            err = nc_def_var(ncID, "bondLength", NC_DOUBLE, 1, dimidsBondLength, &varidBondLength);
            log->check(err == NC_NOERR, "Problem defining bondLength variable: ", nc_strerror(err), ", for file: ", filenameOut);

            int dimidsBondThickness[1] = {dimidBonds};
            err = nc_def_var(ncID, "bondThickness", NC_DOUBLE, 1, dimidsBondThickness, &varidBondThickness);
            log->check(err == NC_NOERR, "Problem defining bondThickness variable: ", nc_strerror(err), ", for file: ", filenameOut);

            int dimidsBondEndPoint1Particle1[2] = {dimidBonds, dimidTWO};
            err = nc_def_var(ncID, "bondEndPoint1Particle1", NC_DOUBLE, 2, dimidsBondEndPoint1Particle1, &varidBondEndPoint1Particle1);
            log->check(err == NC_NOERR, "Problem defining bondEndPoint1Particle1 variable: ", nc_strerror(err), ", for file: ", filenameOut);

            int dimidsBondEndPoint2Particle1[2] = {dimidBonds, dimidTWO};
            err = nc_def_var(ncID, "bondEndPoint2Particle1", NC_DOUBLE, 2, dimidsBondEndPoint2Particle1, &varidBondEndPoint2Particle1);
            log->check(err == NC_NOERR, "Problem defining bondEndPoint2Particle1 variable: ", nc_strerror(err), ", for file: ", filenameOut);

            int dimidsBondEndPoint1Particle2[2] = {dimidBonds, dimidTWO};
            err = nc_def_var(ncID, "bondEndPoint1Particle2", NC_DOUBLE, 2, dimidsBondEndPoint1Particle2, &varidBondEndPoint1Particle2);
            log->check(err == NC_NOERR, "Problem defining bondEndPoint1Particle2 variable: ", nc_strerror(err), ", for file: ", filenameOut);

            int dimidsBondEndPoint2Particle2[2] = {dimidBonds, dimidTWO};
            err = nc_def_var(ncID, "bondEndPoint2Particle2", NC_DOUBLE, 2, dimidsBondEndPoint2Particle2, &varidBondEndPoint2Particle2);
            log->check(err == NC_NOERR, "Problem defining bondEndPoint2Particle2 variable: ", nc_strerror(err), ", for file: ", filenameOut);

          }

        } // nBondsTotal > 0

      } // on_master_proc

    } // bonds active

    if (partition->on_master_proc()) {

      // end the definition phase
      err = nc_enddef(ncID);
      log->check(err == NC_NOERR, "Problem with enddef: ", nc_strerror(err), ", for file: ", filenameOut);

      // write the time
      char* strOut = currentTime.get_char64();

      size_t startTime[1] = {0};
      size_t countTime[1] = {strlen(strOut)+1};

      int varidTime;
      err = nc_inq_varid(ncID, "Time", &varidTime);
      log->check(err == NC_NOERR, "Problem geting varid: ", nc_strerror(err), ", for variable: Time");

      err = nc_put_vara_text(ncID, varidTime, startTime, countTime, &strOut[0]);
      log->check(err == NC_NOERR, "Problem writing time: ", nc_strerror(err), ", for file: ", filenameOut);
      delete [] strOut;

    }

    // write out particle fields
    write_particle_variable   (particles->globalID, (std::string) "globalID", ncID, nParticlesTotal, nParticlesPerProc, particleDispls);
    write_particle_variable   (particles->type,     (std::string) "type",     ncID, nParticlesTotal, nParticlesPerProc, particleDispls);
    write_particle_variable_2d(particles->x,        (std::string) "x",        ncID, nParticlesTotal, nParticlesPerProc, particleDispls);
    write_particle_variable_2d(particles->v,        (std::string) "v",        ncID, nParticlesTotal, nParticlesPerProc, particleDispls);
    write_particle_variable_2d(particles->forcing,  (std::string) "forcing",  ncID, nParticlesTotal, nParticlesPerProc, particleDispls);
    write_particle_variable   (particles->radius,   (std::string) "radius",   ncID, nParticlesTotal, nParticlesPerProc, particleDispls);

    // column variables
    for (int iField = 0 ; iField < outputFieldNames.size() ; iField++) {
      if (column->columnVariables->exists(outputFieldNames[iField])) {
        if (column->columnVariables->type(outputFieldNames[iField]) == DEMSI::FIELDTYPE::DOUBLE_FIELD) {
          write_particle_variable(column->columnVariables->double_field(outputFieldNames[iField]), ncID, nParticlesTotal, nParticlesPerProc, particleDispls);
        } else {
          write_particle_variable(column->columnVariables->int_field   (outputFieldNames[iField]), ncID, nParticlesTotal, nParticlesPerProc, particleDispls);
        }
      }
    }

    if (contacts->bonds_active() and nBondsTotal > 0) {

      // write out bond fields
      write_bond_variable(bondGlobalIDs,          (std::string) "bondGlobalIDs",          ncID, nBondsTotal, nBondsPerProc, bondDispls);

      if (contacts->get_contact_type() == DEMSI::CONTACTTYPE::HOPKINS) {
        write_bond_variable(bondLength,             (std::string) "bondLength",             ncID, nBondsTotal, nBondsPerProc, bondDispls);
        write_bond_variable(bondThickness,          (std::string) "bondThickness",          ncID, nBondsTotal, nBondsPerProc, bondDispls);
        write_bond_variable(bondCrackFraction,      (std::string) "bondCrackFraction",      ncID, nBondsTotal, nBondsPerProc, bondDispls);
        write_bond_variable(bondEndPoint1Particle1, (std::string) "bondEndPoint1Particle1", ncID, nBondsTotal, nBondsPerProc, bondDispls);
        write_bond_variable(bondEndPoint2Particle1, (std::string) "bondEndPoint2Particle1", ncID, nBondsTotal, nBondsPerProc, bondDispls);
        write_bond_variable(bondEndPoint1Particle2, (std::string) "bondEndPoint1Particle2", ncID, nBondsTotal, nBondsPerProc, bondDispls);
        write_bond_variable(bondEndPoint2Particle2, (std::string) "bondEndPoint2Particle2", ncID, nBondsTotal, nBondsPerProc, bondDispls);
      }

    } // bonds active

    if (partition->on_master_proc()) {
      err = nc_close(ncID);
      log->check(err == NC_NOERR, "Problem closing file: ", nc_strerror(err), ", for file: ", filenameOut);
    }

    // reset output alarm
    outputAlarm.reset();

  } // output alarm ringing

}

// stream name
std::string ParticlesWrite::name(void) {
  return streamName;
}

//------------------------------------------------------------------------------
// class ParticlesOutputStreams
//------------------------------------------------------------------------------

// Constructor for the ParticlesOutputStreams class.
ParticlesOutputStreams::ParticlesOutputStreams(DEMSI::Partition* partition, DEMSI::Log* logIn, DEMSI::Configs* configs, DEMSI::Contacts* contacts, DEMSI::Particles* particles, DEMSI::Column* column, DEMSI::Clock* clock, DEMSI::Time alarmStartTime) {

  log = logIn;

  // get particle output stream names
  std::vector<std::string> outputStreamNames = configs->get_array({"ConfigGroup:particleOutput","Array:particleOutputStreams"}, "ParticleOutputStream", "name");

  // iterate over output streams
  for (int iStream=0 ; iStream<outputStreamNames.size() ; iStream++) {

    std::string streamName = outputStreamNames[iStream];

    // output stream configs location
    std::list<std::string> configLocationRoot, configLocation;
    configLocationRoot.push_back("ConfigGroup:particleOutput");
    configLocationRoot.push_back("Array:particleOutputStreams");
    configLocationRoot.push_back("ParticleOutputStream:" + streamName);

    // get output stream configs
    std::string writeFilenameTemplate;
    configLocation = configLocationRoot;
    configLocation.push_back("Config:particleWriteFilenameTemplate");
    configs->get(configLocation, writeFilenameTemplate);

    std::string writeOutputDirectory;
    configLocation = configLocationRoot;
    configLocation.push_back("Config:particleWriteOutputDirectory");
    configs->get(configLocation, writeOutputDirectory);

    writeFilenameTemplate = writeOutputDirectory + writeFilenameTemplate;

    std::string writeIntervalStr;
    configLocation = configLocationRoot;
    configLocation.push_back("Config:particleWriteInterval");
    configs->get(configLocation, writeIntervalStr);

    DEMSI::TimeInterval writeInterval(writeIntervalStr, log);

    bool writeClobber;
    configLocation = configLocationRoot;
    configLocation.push_back("Config:particleWriteClobber");
    configs->get(configLocation, writeClobber);

    bool writeBonds;
    configLocation = configLocationRoot;
    configLocation.push_back("Config:particleWriteBonds");
    configs->get(configLocation, writeBonds);

    // get output field names
    std::list<std::string> outputFieldsLocation = configLocationRoot;
    outputFieldsLocation.push_back("Array:OutputFields");
    std::vector<std::string> outputFieldNames = configs->get_array(outputFieldsLocation, "OutputField");

    // create a new ParticlesWrite object
    DEMSI::ParticlesWrite* particlesWrite = new DEMSI::ParticlesWrite(partition, log, contacts, particles, column, streamName, writeFilenameTemplate, clock, alarmStartTime, writeInterval, writeClobber, outputFieldNames);
    streams.push_back(particlesWrite);

  } // iStream

}

// Write out all particle output streams
void ParticlesOutputStreams::write(void) {

  for (int iStream = 0 ; iStream<streams.size() ; iStream++) {
    (*log)(DEMSI::LOG::DEBUG) << "ParticlesOutputStreams::write: " << iStream << " " << streams[iStream]->name() << std::endl;
    streams[iStream]->write();
  }

}

} // DEMSI namespace
