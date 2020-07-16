/*
 * demsi_column_transfer.cpp
 *
 *  Created on: Oct 17, 2017
 *      Author: dsbolin
 */

#include "demsi_column_variables.h"
#include "demsi_communication.h"
#include <algorithm>
#include <map>
#include <mpi.h>

namespace DEMSI {

// convert the proc ID maps from storing globalID info to particle index info.
void convert_globalID_to_index(std::map<int, std::vector<int>>* procIDs, std::map<int,int> globalIDsToIndices) {

  // iterator for procIDsItr
  std::map<int, std::vector<int>>::iterator procIDsItr;

  procIDsItr = procIDs->begin();
  while (procIDsItr != procIDs->end()) {

    std::sort(procIDsItr->second.begin(), procIDsItr->second.end());

    for (int i=0 ; i < procIDsItr->second.size() ; i++) {
      procIDsItr->second[i] = globalIDsToIndices.find(procIDsItr->second[i])->second;
    }

    procIDsItr++;
  }

}

// transfer a single particle column variable between processors
void transfer_variable(DEMSI::ColumnVariable<double>* columnVariable, std::map<int, std::vector<int>> procIDsForLostParticles, std::map<int, std::vector<int>> procIDsForGainedParticles, DEMSI::Partition* partition, DEMSI::Particles* particles, DEMSI::Log* log) {

  int mpiErr;
  char mpiErrBuffer[MPI_MAX_ERROR_STRING];
  int mpiErrLen;
  MPI_Status status;

  // tmp array for new values of variable
  double newVariable[std::max(*(particles->nParticles) * columnVariable->size_per_particle(),1)];

  // again loop around processors
  for (int iProc=0 ; iProc < partition->nprocs()-1 ; iProc++) {

    // which processor we are sending data to
    int iProcSend = modulo((partition->this_proc() + iProc + 1), partition->nprocs());

    // send the data
    std::map<int, std::vector<int>>::iterator sendItr = procIDsForLostParticles.find(iProcSend);
    // check if there are particles we need to send
    if (sendItr != procIDsForLostParticles.end() ) {

      // There are particles we need to send
      double dataSend[sendItr->second.size() * columnVariable->size_per_particle()];
      for (int iParticle = 0 ; iParticle < sendItr->second.size() ; iParticle++) {
	for (int j = 0 ; j < columnVariable->size_per_particle() ; j++) {
	  int ij = iParticle * columnVariable->size_per_particle() + j;
	  dataSend[ij] = (*columnVariable)(sendItr->second[iParticle],j);
	} // j
      } // iParticle
      mpiErr = MPI_Send(&dataSend[0], sendItr->second.size() * columnVariable->size_per_particle(), MPI_DOUBLE, iProcSend, 0, partition->comm());
      MPI_Error_string(mpiErr, mpiErrBuffer, &mpiErrLen);
      log->check(mpiErr == MPI_SUCCESS, ": ", (std::string) mpiErrBuffer);
    }

    // which processor are we going to receive data from
    int iProcRecv = modulo((partition->this_proc() - iProc - 1), partition->nprocs());

    // receive the data
    std::map<int, std::vector<int>>::iterator recvItr = procIDsForGainedParticles.find(iProcRecv);
    // check if there are particles we need to recv
    if (recvItr != procIDsForGainedParticles.end() ) {

      // There are particles we need to recv
      double dataRecv[recvItr->second.size() * columnVariable->size_per_particle()];
      mpiErr = MPI_Recv(&dataRecv[0], recvItr->second.size() * columnVariable->size_per_particle(), MPI_DOUBLE, iProcRecv, 0, partition->comm(), &status);
      MPI_Error_string(mpiErr, mpiErrBuffer, &mpiErrLen);
      log->check(mpiErr == MPI_SUCCESS, ": ", (std::string) mpiErrBuffer);

      for (int iParticle = 0 ; iParticle < recvItr->second.size() ; iParticle++) {
	for (int j = 0 ; j < columnVariable->size_per_particle() ; j++) {
	  int ij  = iParticle * columnVariable->size_per_particle() + j;
	  int ij2 = recvItr->second[iParticle] * columnVariable->size_per_particle() + j;
	  newVariable[ij2] = dataRecv[ij];
	} // j
      } // iParticle
    }

  } // processor loop

  // fill in particles we keep
  for (int i = 0 ; i < procIDsForGainedParticles[partition->this_proc()].size() ; i++) {
    for (int j = 0 ; j < columnVariable->size_per_particle() ; j++) {
      int ij = procIDsForGainedParticles[partition->this_proc()][i]  * columnVariable->size_per_particle() + j;
      newVariable[ij] = (*columnVariable)(procIDsForLostParticles[partition->this_proc()][i],j);
    } // j
  } // i

  // resize view
  columnVariable->resize(*(particles->nParticles));

  // copy new array to view
  for (int iParticle=0 ; iParticle < *(particles->nParticles) ; iParticle++) {
    for (int j = 0 ; j < columnVariable->size_per_particle() ; j++) {
      int ij = iParticle * columnVariable->size_per_particle() + j;
      (*columnVariable)(iParticle,j) = newVariable[ij];
    } // j
  } // iParticle

}

// transfer a single particle column variable between processors
void transfer_variable(DEMSI::ColumnVariable<int>* columnVariable, std::map<int, std::vector<int>> procIDsForLostParticles, std::map<int, std::vector<int>> procIDsForGainedParticles, DEMSI::Partition* partition, DEMSI::Particles* particles, DEMSI::Log* log) {

  int mpiErr;
  char mpiErrBuffer[MPI_MAX_ERROR_STRING];
  int mpiErrLen;
  MPI_Status status;

  // tmp array for new values of variable
  int newVariable[std::max(*(particles->nParticles) * columnVariable->size_per_particle(),1)];

  // again loop around processors
  for (int iProc=0 ; iProc < partition->nprocs()-1 ; iProc++) {

    // which processor we are sending data to
    int iProcSend = modulo((partition->this_proc() + iProc + 1), partition->nprocs());

    // send the data
    std::map<int, std::vector<int>>::iterator sendItr = procIDsForLostParticles.find(iProcSend);
    // check if there are particles we need to send
    if (sendItr != procIDsForLostParticles.end() ) {

      // There are particles we need to send
      int dataSend[sendItr->second.size() * columnVariable->size_per_particle()];
      for (int iParticle = 0 ; iParticle < sendItr->second.size() ; iParticle++) {
	for (int j = 0 ; j < columnVariable->size_per_particle() ; j++) {
	  int ij = iParticle * columnVariable->size_per_particle() + j;
	  dataSend[ij] = (*columnVariable)(sendItr->second[iParticle],j);
	} // j
      } // iParticle
      mpiErr = MPI_Send(&dataSend[0], sendItr->second.size() * columnVariable->size_per_particle(), MPI_INT, iProcSend, 0, partition->comm());
      MPI_Error_string(mpiErr, mpiErrBuffer, &mpiErrLen);
      log->check(mpiErr == MPI_SUCCESS, ": ", (std::string) mpiErrBuffer);
    }

    // which processor are we going to receive data from
    int iProcRecv = modulo((partition->this_proc() - iProc - 1), partition->nprocs());

    // receive the data
    std::map<int, std::vector<int>>::iterator recvItr = procIDsForGainedParticles.find(iProcRecv);
    // check if there are particles we need to recv
    if (recvItr != procIDsForGainedParticles.end() ) {

      // There are particles we need to recv
      int dataRecv[recvItr->second.size() * columnVariable->size_per_particle()];
      mpiErr = MPI_Recv(&dataRecv[0], recvItr->second.size() * columnVariable->size_per_particle(), MPI_INT, iProcRecv, 0, partition->comm(), &status);
      MPI_Error_string(mpiErr, mpiErrBuffer, &mpiErrLen);
      log->check(mpiErr == MPI_SUCCESS, ": ", (std::string) mpiErrBuffer);

      for (int iParticle = 0 ; iParticle < recvItr->second.size() ; iParticle++) {
	for (int j = 0 ; j < columnVariable->size_per_particle() ; j++) {
	  int ij  = iParticle * columnVariable->size_per_particle() + j;
	  int ij2 = recvItr->second[iParticle] * columnVariable->size_per_particle() + j;
	  newVariable[ij2] = dataRecv[ij];
	} // j
      } // iParticle
    }

  } // processor loop

  // fill in particles we keep
  for (int i = 0 ; i < procIDsForGainedParticles[partition->this_proc()].size() ; i++) {
    for (int j = 0 ; j < columnVariable->size_per_particle() ; j++) {
      int ij = procIDsForGainedParticles[partition->this_proc()][i]  * columnVariable->size_per_particle() + j;
      newVariable[ij] = (*columnVariable)(procIDsForLostParticles[partition->this_proc()][i],j);
    } // j
  } // i

  // resize view
  columnVariable->resize(*(particles->nParticles));

  // copy new array to view
  for (int iParticle=0 ; iParticle < *(particles->nParticles) ; iParticle++) {
    for (int j = 0 ; j < columnVariable->size_per_particle() ; j++) {
      int ij = iParticle * columnVariable->size_per_particle() + j;
      (*columnVariable)(iParticle,j) = newVariable[ij];
    } // j
  } // iParticle

}

// Initialize transfer column data between processors.
void ColumnVariables::init_processor_transfer(void) {

  // set the set of owned particle IDs
  previousOwnedParticles.clear();
  for (int iParticle = 0; iParticle < *(particles->nParticles); iParticle++) {
    previousOwnedParticles.insert(particles->globalID[iParticle]);
  }

  // map to easily get particle index from globalID
  globalIDsToIndicesPrevious.clear();
  for (int iParticle=0 ; iParticle < *(particles->nParticles) ; iParticle++) {
    globalIDsToIndicesPrevious[particles->globalID[iParticle]] = iParticle;
  }

}

// Initialize transfer column data between processors after remapping.
void ColumnVariables::init_processor_transfer_remap(const int nParticles) {

  char mpiErrBuffer[MPI_MAX_ERROR_STRING];
  int mpiErrLen;
  int mpiErr;

  // number of owned init particles on each processor
  int *nParticlesOnAllProcs = new int[partition->nprocs()];

  mpiErr = MPI_Allgather(&nParticles, 1, MPI_INT, nParticlesOnAllProcs, 1, MPI_INT, partition->comm());
  MPI_Error_string(mpiErr, mpiErrBuffer, &mpiErrLen);
  log->check(mpiErr == MPI_SUCCESS, "init_processor_transfer_remap: MPI_Allgather failed: ", (std::string) mpiErrBuffer);

  // total number of particles
  int nParticlesTotal = 0;
  for (int iProc = 0 ; iProc < partition->nprocs() ; iProc++) {
    nParticlesTotal += nParticlesOnAllProcs[iProc];
  } // nParticlesTotal

  // cumulative number of owned init particles on previous processors
  int *nParticlesOnPrevProc = new int[partition->nprocs()];

  nParticlesOnPrevProc[0] = 0;
  for (int iProc = 1 ; iProc < partition->nprocs() ; iProc++) {
    nParticlesOnPrevProc[iProc] = nParticlesOnPrevProc[iProc-1] + nParticlesOnAllProcs[iProc-1];
  } // iProc

  // set the set of owned particle IDs
  previousOwnedParticles.clear();
  for (int iParticle = 0; iParticle < nParticles; iParticle++) {
    int globalID = iParticle + 1 + nParticlesOnPrevProc[partition->this_proc()];
    previousOwnedParticles.insert(globalID);
  }

  // map to easily get particle index from globalID
  globalIDsToIndicesPrevious.clear();
  for (int iParticle=0 ; iParticle < nParticles ; iParticle++) {
    int globalID = iParticle + 1 + nParticlesOnPrevProc[partition->this_proc()];
    globalIDsToIndicesPrevious[globalID] = iParticle;
  }

}

// Transfer column data between processors.
void ColumnVariables::processor_transfer(void) {

  // mpi error variables
  int mpiErr;
  char mpiErrBuffer[MPI_MAX_ERROR_STRING];
  int mpiErrLen;
  MPI_Status status;

  // Current particles on this processor
  std::set <int> currentOwnedParticles;
  for (int iParticle = 0; iParticle < *(particles->nParticles); iParticle++) {
    currentOwnedParticles.insert(particles->globalID[iParticle]);
  }

  // Any particles that are in previous owned list but not in current must be sent
  std::vector <int> particlesLost;
  std::set_difference(previousOwnedParticles.begin(), previousOwnedParticles.end(), currentOwnedParticles.begin(), currentOwnedParticles.end(), std::back_inserter(particlesLost));
  std::sort(particlesLost.begin(), particlesLost.end());

  // Any particles that are not in my previous owned list but are in current list must be received
  std::vector <int> particlesGained;
  std::set_difference(currentOwnedParticles.begin(), currentOwnedParticles.end(), previousOwnedParticles.begin(), previousOwnedParticles.end(), std::back_inserter(particlesGained));
  std::sort(particlesGained.begin(), particlesGained.end());

  // Any particles we are going to keep
  std::vector <int> particlesKept;
  std::set_intersection(previousOwnedParticles.begin(), previousOwnedParticles.end(), currentOwnedParticles.begin(), currentOwnedParticles.end(), std::back_inserter(particlesKept));
  std::sort(particlesKept.begin(), particlesKept.end());

  // processor IDs and particle global IDs to send lost particles to
  std::map<int, std::vector<int>> procIDsForLostParticles; // map(procID, vector of globalIDs (converted to particle indices))

  // processor IDs and particle global IDs to receive gaimed particles from
  std::map<int, std::vector<int>> procIDsForGainedParticles; // map(procID, vector of globalIDs (converted to particle indices))

  // send and recv lists to send around processors
  std::vector<int> particlesLostTransit;
  for (int i=0 ; i < particlesLost.size() ; i++) particlesLostTransit.push_back(particlesLost[i]);
  std::vector<int> particlesGainedTransit;
  for (int i=0 ; i < particlesGained.size() ; i++) particlesGainedTransit.push_back(particlesGained[i]);

  //--------------------------------------------------------------
  // lists of particles lost from the original processor
  //--------------------------------------------------------------
  DEMSI::ShareLists* shareLostLists = new DEMSI::ShareLists(&particlesLostTransit,NULL,partition,log);
  while (shareLostLists->iterate()) {

    // find if any of the list we just received is in the owned receive list for this processor
    std::vector<int>::iterator particlesLostTransitItr;
    for (int i=0 ; i < particlesGained.size() ; i++) {

      particlesLostTransitItr = std::find(particlesLostTransit.begin(), particlesLostTransit.end(), particlesGained[i]);
      if (particlesLostTransitItr != particlesLostTransit.end()) {
	procIDsForGainedParticles[shareLostLists->iProc_origin()].push_back(particlesGained[i]); // add particle to this processor receive list
	particlesLostTransit.erase(particlesLostTransitItr); // remove from list we are sending around all processors
      }

    }

  } // share lost lists
  delete shareLostLists;

  //--------------------------------------------------------------
  // lists of particles gained by the original processor
  //--------------------------------------------------------------
  DEMSI::ShareLists* shareGainedLists = new DEMSI::ShareLists(&particlesGainedTransit,NULL,partition,log);
  while (shareGainedLists->iterate()) {

    // find if any of the sendList we just received is in the owned receive list for this processor
    std::vector<int>::iterator particlesGainedTransitItr;
    for (int i=0 ; i < particlesLost.size() ; i++) {

      particlesGainedTransitItr = std::find(particlesGainedTransit.begin(), particlesGainedTransit.end(), particlesLost[i]);
      if (particlesGainedTransitItr != particlesGainedTransit.end()) {
	procIDsForLostParticles[shareGainedLists->iProc_origin()].push_back(particlesLost[i]); // add particle to this processor send list
	particlesGainedTransit.erase(particlesGainedTransitItr); // remove from list we are sending around all processors
      }

    }

  } // share gained lists
  delete shareGainedLists;

  // map to easily get particle index from globalID
  std::map<int,int> globalIDsToIndicesCurrent;
  for (int iParticle=0 ; iParticle < *(particles->nParticles) ; iParticle++) {
    globalIDsToIndicesCurrent[particles->globalID[iParticle]] = iParticle;
  }

  // convert globalIDs to local particle indices
  convert_globalID_to_index(&procIDsForGainedParticles, globalIDsToIndicesCurrent);
  convert_globalID_to_index(&procIDsForLostParticles, globalIDsToIndicesPrevious);

  // find indices for particles we keep
  for (int i=0; i < particlesKept.size() ; i++) {
    procIDsForLostParticles[partition->this_proc()].push_back(globalIDsToIndicesPrevious.find(particlesKept[i])->second); // old position
    procIDsForGainedParticles[partition->this_proc()].push_back(globalIDsToIndicesCurrent.find(particlesKept[i])->second); // new position
  }

  // transfer double variables
  for (int iVar = 0 ; iVar < doubleVariables.size() ; iVar++) {
    if (doubleVariables[iVar]->active()) {
      transfer_variable(doubleVariables[iVar], procIDsForLostParticles, procIDsForGainedParticles, partition, particles, log);
    } // active
  } // iVar

  // transfer integer variables
  for (int iVar = 0 ; iVar < intVariables.size() ; iVar++) {
    if (intVariables[iVar]->active()) {
      transfer_variable(intVariables[iVar], procIDsForLostParticles, procIDsForGainedParticles, partition, particles, log);
    } // active
  } // iVar

  // transfer current to previous
  previousOwnedParticles = currentOwnedParticles;
  globalIDsToIndicesPrevious = globalIDsToIndicesCurrent;

}

} // namespace DEMSI
