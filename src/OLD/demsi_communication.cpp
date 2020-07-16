#include "demsi_communication.h"

namespace DEMSI {

// calculate the modulo and do it correctly for negative numbers
int modulo(int i, int n) {
  int mod = i % n;
  if (mod < 0) {
    mod = mod + n;
  }
  return mod;
}

// send a list to another processor
void send_recv_particle_list(std::vector<int>* particleList, int iProcSend, int iProcRecv, DEMSI::Partition* partition, DEMSI::Log* log) {

  int mpiErr;
  char mpiErrBuffer[MPI_MAX_ERROR_STRING];
  int mpiErrLen;

  MPI_Request requests[2];
  MPI_Status statuses[2];

  int* particleListSend = new int[particleList->size()];
  for (int i=0 ; i < particleList->size() ; i++) {
    particleListSend[i] = (*particleList)[i];
  }

  mpiErr = MPI_Isend(&particleListSend[0], particleList->size(), MPI_INT, iProcSend, 0, partition->comm(), &requests[0]);
  MPI_Error_string(mpiErr, mpiErrBuffer, &mpiErrLen);
  log->check(mpiErr == MPI_SUCCESS, ": ", (std::string) mpiErrBuffer);

  // get the status of the receive we are expecting
  MPI_Status probeStatus;
  mpiErr = MPI_Probe(iProcRecv, 0, partition->comm(), &probeStatus);
  MPI_Error_string(mpiErr, mpiErrBuffer, &mpiErrLen);
  log->check(mpiErr == MPI_SUCCESS, ": ", (std::string) mpiErrBuffer);

  // from the status determine the number of particles we can expect
  int nRecv;
  mpiErr = MPI_Get_count(&probeStatus, MPI_INT, &nRecv);
  MPI_Error_string(mpiErr, mpiErrBuffer, &mpiErrLen);
  log->check(mpiErr == MPI_SUCCESS, ": ", (std::string) mpiErrBuffer);

  // receive the send lists from previous processor
  int* particleListRecv = new int[nRecv];
  mpiErr = MPI_Irecv(&particleListRecv[0], nRecv, MPI_INT, iProcRecv, 0, partition->comm(), &requests[1]);
  MPI_Error_string(mpiErr, mpiErrBuffer, &mpiErrLen);
  log->check(mpiErr == MPI_SUCCESS, ": ", (std::string) mpiErrBuffer);

  MPI_Waitall(2, requests, statuses);

  delete [] particleListSend;

  particleList->clear();
  for (int i=0 ; i < nRecv ; i++) {
    particleList->push_back(particleListRecv[i]);
  }

  delete [] particleListRecv;

}

// send a list to another processor
void send_recv_particle_list(std::vector<double>* particleList, int iProcSend, int iProcRecv, DEMSI::Partition* partition, DEMSI::Log* log) {

  int mpiErr;
  char mpiErrBuffer[MPI_MAX_ERROR_STRING];
  int mpiErrLen;

  MPI_Request requests[2];
  MPI_Status statuses[2];

  double* particleListSend = new double[particleList->size()];
  for (int i=0 ; i < particleList->size() ; i++) {
    particleListSend[i] = (*particleList)[i];
  }

  mpiErr = MPI_Isend(&particleListSend[0], particleList->size(), MPI_DOUBLE, iProcSend, 0, partition->comm(), &requests[0]);
  MPI_Error_string(mpiErr, mpiErrBuffer, &mpiErrLen);
  log->check(mpiErr == MPI_SUCCESS, ": ", (std::string) mpiErrBuffer);

  // get the status of the receive we are expecting
  MPI_Status probeStatus;
  mpiErr = MPI_Probe(iProcRecv, 0, partition->comm(), &probeStatus);
  MPI_Error_string(mpiErr, mpiErrBuffer, &mpiErrLen);
  log->check(mpiErr == MPI_SUCCESS, ": ", (std::string) mpiErrBuffer);

  // from the status determine the number of particles we can expect
  int nRecv;
  mpiErr = MPI_Get_count(&probeStatus, MPI_DOUBLE, &nRecv);
  MPI_Error_string(mpiErr, mpiErrBuffer, &mpiErrLen);
  log->check(mpiErr == MPI_SUCCESS, ": ", (std::string) mpiErrBuffer);

  // receive the send lists from previous processor
  double* particleListRecv = new double[nRecv];
  mpiErr = MPI_Irecv(&particleListRecv[0], nRecv, MPI_DOUBLE, iProcRecv, 0, partition->comm(), &requests[1]);
  MPI_Error_string(mpiErr, mpiErrBuffer, &mpiErrLen);
  log->check(mpiErr == MPI_SUCCESS, ": ", (std::string) mpiErrBuffer);

  MPI_Waitall(2, requests, statuses);

  delete [] particleListSend;

  particleList->clear();
  for (int i=0 ; i < nRecv ; i++) {
    particleList->push_back(particleListRecv[i]);
  }

  delete [] particleListRecv;

}

// Constructor for ShareLists object
ShareLists::ShareLists(std::vector<int>* listIntegerIn, std::vector<double>* listDoubleIn, DEMSI::Partition* partitionIn, DEMSI::Log* logIn) {

  partition = partitionIn;
  log = logIn;

  if (listIntegerIn == NULL) {
    hasIntegerList = false;
  } else {
    hasIntegerList = true;
    listInteger = listIntegerIn;
  }

  if (listDoubleIn == NULL) {
    hasDoubleList = false;
  } else {
    hasDoubleList = true;
    listDouble = listDoubleIn;
  }

  iProcTransit = -1;
  iStep = 0;

}

// Move lists to next processor
bool ShareLists::iterate(void) {

  if (iStep >= partition->nprocs()-1) return false;

  // processor current list is from
  iProcTransit = modulo((partition->this_proc() - iStep - 1), partition->nprocs());

  // where we are sending the lists
  int iProcSendThisStep = modulo((partition->this_proc() + 1), partition->nprocs());

  // where we are getting the lists from
  int iProcRecvThisStep = modulo((partition->this_proc() - 1), partition->nprocs());

  // Transfer integer list
  if (hasIntegerList) {

    send_recv_particle_list(listInteger, iProcSendThisStep, iProcRecvThisStep, partition, log);

  } // hasDoubleList

  // Transfer double list
  if (hasDoubleList) {

    send_recv_particle_list(listDouble, iProcSendThisStep, iProcRecvThisStep, partition, log);

  } // hasDoubleList

  // increment step
  iStep++;

  return true;

}

// Return processor ID where the current list originally came from.
int ShareLists::iProc_origin(void) {

  return iProcTransit;

}

} // namespace DEMSI
