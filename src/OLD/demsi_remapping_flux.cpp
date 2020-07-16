#include "demsi_remapping.h"
#include "demsi_communication.h"
#include "demsi_column.h"

#include <permonqps.h>

namespace DEMSI {

//#define DEBUG_REMAPPING

//------------------------------------------------------------------------------
// Remapping correction with flux based otimization
//------------------------------------------------------------------------------

// correct weights matrix using flux method with PermonQP
void Remapping::weights_correction(Mat &Weights, std::vector<double> effectiveElementArea) {

  // MPI error checking
  char mpiErrBuffer[MPI_MAX_ERROR_STRING];
  int mpiErrLen;
  int mpiErr;


  // effective area ice without correction
  double* effectiveAreaIce = new double[tessellation->nParticlesAllInit];
  for (int iParticleInit = 0 ; iParticleInit < tessellation->nParticlesAllInit ; iParticleInit++) {
    effectiveAreaIce[iParticleInit] = 0.0;
  } // iParticleInit
  remap_variable(1, effectiveAreaIce, &effectiveElementArea[0], Weights);
  tessellation->update_init_halo_particles(effectiveAreaIce);


  // particle area
  double* elementArea = new double[tessellation->nParticlesAllInit];
  for (int iParticleInit = 0 ; iParticleInit < tessellation->nParticlesInit ; iParticleInit++) {
    elementArea[iParticleInit] = tessellation->initPolygons[iParticleInit].area();
  } // iParticleInit


  // find new elements in use
  int *useParticle = new int[tessellation->nParticlesAllInit];
  for (int iParticleInit = 0 ; iParticleInit < tessellation->nParticlesInit ; iParticleInit++) {
    if (effectiveAreaIce[iParticleInit] > 0.0 and
	tessellation->typeInit[iParticleInit] == 0) {
      useParticle[iParticleInit] = 1;
    } else {
      useParticle[iParticleInit] = 0;
    }
  } // iParticleInit
  tessellation->update_init_halo_particles(useParticle);

  int* useParticleTmp = new int[tessellation->nParticlesAllInit];

  const int nExpand = 1;
  for (int iExpand = 0 ; iExpand < nExpand ; iExpand++) {

    for (int iParticleInit = 0 ; iParticleInit < tessellation->nParticlesAllInit ; iParticleInit++) {
      useParticleTmp[iParticleInit] = useParticle[iParticleInit];
    } // iParticleInit

    for (int iParticleInit = 0 ; iParticleInit < tessellation->nParticlesInit ; iParticleInit++) {

      for (int iCellOnCell = 0 ; iCellOnCell < tessellation->nVerticesOnCellInit[iParticleInit] ; iCellOnCell++) {
	int iParticleInitNeighbour = tessellation->cellsOnCellLocal[iParticleInit][iCellOnCell];
	if (iParticleInitNeighbour != -1 and
	    useParticleTmp[iParticleInit] == 0 and
	    useParticleTmp[iParticleInitNeighbour] == 1 and
	    tessellation->typeInit[iParticleInit] == 0) {
	  useParticle[iParticleInit] = 1;
	}
      } // iCellOnCell

    } // iParticleInit

    tessellation->update_init_halo_particles(useParticle);

  } // iExpand

  delete [] useParticleTmp;

#ifdef DEBUG_REMAPPING
  int nParticlesUseTotal;
  int nParticlesUseLocal = 0;
  int nParticlesUseHalo = 0;
  for (int iParticleInit = 0 ; iParticleInit < tessellation->nParticlesInit ; iParticleInit++) {
    if (useParticle[iParticleInit]) nParticlesUseLocal++;
  } // iParticleInit
  for (int iParticleInit = tessellation->nParticlesInit ; iParticleInit < tessellation->nParticlesAllInit ; iParticleInit++) {
    if (useParticle[iParticleInit]) nParticlesUseHalo++;
  } // iParticleInit
  mpiErr = MPI_Allreduce(&nParticlesUseLocal, &nParticlesUseTotal, 1, MPI_INT, MPI_SUM, partition->comm());
  MPI_Error_string(mpiErr, mpiErrBuffer, &mpiErrLen);
  log->check(mpiErr == MPI_SUCCESS, "weights_correction: MPI_Allreduce nParticlesUse failed: ", (std::string) mpiErrBuffer);
  (*log)(DEMSI::LOG::DEBUG) << "nParticlesUseLocal: " << nParticlesUseLocal << std::endl;
  (*log)(DEMSI::LOG::DEBUG) << "nParticlesUseHalo:  " << nParticlesUseHalo  << " (nParticlesHaloInit: " << tessellation->nParticlesHaloInit << ")" << std::endl;
  (*log)(DEMSI::LOG::DEBUG) << "nParticlesUseTotal: " << nParticlesUseTotal << std::endl;

  for (int iParticleInit = 0 ; iParticleInit < tessellation->nParticlesInit ; iParticleInit++) {
    useParticleDebug[iParticleInit] = (double) useParticle[iParticleInit];
    elementAreaDebug[iParticleInit] = elementArea[iParticleInit];
    effectiveAreaIceDebug[iParticleInit] = effectiveAreaIce[iParticleInit];
  } // iParticleInit

  tessellationOutputStreams->write("remappingDebugOutput", true);
#endif


#ifdef DEBUG_REMAPPING
  // initial excess and total area
  double excessAreaStartLocal = 0.0;
  double totalAreaStartLocal  = 0.0;
  double availableAreaLocal = 0.0;
  for (int iParticleInit = 0 ; iParticleInit < tessellation->nParticlesInit ; iParticleInit++) {
    excessAreaStartLocal += std::max(effectiveAreaIce[iParticleInit] - elementArea[iParticleInit],0.0);
    totalAreaStartLocal  += effectiveAreaIce[iParticleInit];
    if (useParticle[iParticleInit] == 1 and tessellation->typeInit[iParticleInit] != 2) {
      availableAreaLocal += elementArea[iParticleInit];
    }
  } // iParticleInit

  double excessAreaStart = 0.0;
  double totalAreaStart  = 0.0;
  double availableArea = 0.0;

  mpiErr = MPI_Allreduce(&excessAreaStartLocal, &excessAreaStart, 1, MPI_DOUBLE, MPI_SUM, partition->comm());
  MPI_Error_string(mpiErr, mpiErrBuffer, &mpiErrLen);
  log->check(mpiErr == MPI_SUCCESS, "weights_correction: MPI_Allreduce excessAreaStart failed: ", (std::string) mpiErrBuffer);

  mpiErr = MPI_Allreduce(&totalAreaStartLocal, &totalAreaStart, 1, MPI_DOUBLE, MPI_SUM, partition->comm());
  MPI_Error_string(mpiErr, mpiErrBuffer, &mpiErrLen);
  log->check(mpiErr == MPI_SUCCESS, "weights_correction: MPI_Allreduce totalAreaStart failed: ", (std::string) mpiErrBuffer);

  mpiErr = MPI_Allreduce(&availableAreaLocal, &availableArea, 1, MPI_DOUBLE, MPI_SUM, partition->comm());
  MPI_Error_string(mpiErr, mpiErrBuffer, &mpiErrLen);
  log->check(mpiErr == MPI_SUCCESS, "weights_correction: MPI_Allreduce availableArea failed: ", (std::string) mpiErrBuffer);

  (*log)(DEMSI::LOG::DEBUG) << std::endl;
  (*log)(DEMSI::LOG::DEBUG) << "Remapping areas:" << std::endl;
  (*log)(DEMSI::LOG::DEBUG) << "   Excess area:    " << excessAreaStart << std::endl;
  (*log)(DEMSI::LOG::DEBUG) << "   Total area:     " << totalAreaStart << std::endl;
  (*log)(DEMSI::LOG::DEBUG) << "   Ratio:          " << excessAreaStart / totalAreaStart << std::endl;
  (*log)(DEMSI::LOG::DEBUG) << "Remapping correction:" << std::endl;
  (*log)(DEMSI::LOG::DEBUG) << "   Available area: " << availableArea << std::endl;
  (*log)(DEMSI::LOG::DEBUG) << "   Required area:  " << totalAreaStart << std::endl;
  (*log)(DEMSI::LOG::DEBUG) << "   Ratio:          " << availableArea / totalAreaStart << std::endl;
  (*log)(DEMSI::LOG::DEBUG) << std::endl;
#endif


  // particles in use
  std::vector<int> iParticleInitUse;
  for (int iParticleInit = 0 ; iParticleInit < tessellation->nParticlesInit ; iParticleInit++) {
    if (useParticle[iParticleInit] == 1) {
      iParticleInitUse.push_back(iParticleInit);
    }
  } // iParticleInit
  int nParticles = iParticleInitUse.size(); // not halo

  // number of owned init particles on each processor
  int *nParticlesOnAllProcs = new int[partition->nprocs()];

  mpiErr = MPI_Allgather(&nParticles, 1, MPI_INT, nParticlesOnAllProcs, 1, MPI_INT, partition->comm());
  MPI_Error_string(mpiErr, mpiErrBuffer, &mpiErrLen);
  log->check(mpiErr == MPI_SUCCESS, "weights_correction: MPI_Allgather failed: ", (std::string) mpiErrBuffer);

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

  // mapping from all particles to reduced particles
  std::map<int,int> iParticleMap;
  for (int iParticle = 0 ; iParticle < nParticles ; iParticle++) {
    iParticleMap[iParticleInitUse[iParticle]] = iParticle;
  } // iParticle

  int *iParticleGlobal = new int[tessellation->nParticlesAllInit];
  for (int iParticleInit = 0 ; iParticleInit < tessellation->nParticlesAllInit ; iParticleInit++) {
    iParticleGlobal[iParticleInit] = 0;
  } // iParticleInit
  for (int iParticle = 0 ; iParticle < iParticleInitUse.size() ; iParticle++) {
    iParticleGlobal[iParticleInitUse[iParticle]] = iParticle + nParticlesOnPrevProc[partition->this_proc()];
  } // iParticle
  tessellation->update_init_halo_particles(iParticleGlobal);

  int *iParticleGlobal1 = new int[tessellation->nParticlesAllInit];
  for (int iParticleInit = 0 ; iParticleInit < tessellation->nParticlesAllInit ; iParticleInit++) {
    iParticleGlobal1[iParticleInit] = 0;
  } // iParticleInit
  for (int iParticle = 0 ; iParticle < iParticleInitUse.size() ; iParticle++) {
    iParticleGlobal1[iParticleInitUse[iParticle]] = iParticle + 2*nParticlesOnPrevProc[partition->this_proc()];
  } // iParticle
  tessellation->update_init_halo_particles(iParticleGlobal1);

  int *iParticleGlobal2 = new int[tessellation->nParticlesAllInit];
  for (int iParticleInit = 0 ; iParticleInit < tessellation->nParticlesAllInit ; iParticleInit++) {
    iParticleGlobal2[iParticleInit] = 0;
  } // iParticleInit
  for (int iParticle = 0 ; iParticle < iParticleInitUse.size() ; iParticle++) {
    iParticleGlobal2[iParticleInitUse[iParticle]] = iParticle + 2*nParticlesOnPrevProc[partition->this_proc()] + nParticles;
  } // iParticle
  tessellation->update_init_halo_particles(iParticleGlobal2);

#ifdef DEBUG_REMAPPING
  (*log)(DEMSI::LOG::DEBUG) << "nParticlesTotal: " << nParticlesTotal << std::endl;
  for (int iProc = 0 ; iProc < partition->nprocs() ; iProc++) {
    (*log)(DEMSI::LOG::DEBUG) << "nParticlesOnAllProcs nParticlesOnPrevProc: " << iProc << " ";
    (*log)(DEMSI::LOG::DEBUG) << nParticlesOnAllProcs[iProc] << " " << nParticlesOnPrevProc[iProc] << std::endl;
  } // nParticlesTotal
  (*log)(DEMSI::LOG::DEBUG) << std::endl;
#endif



  // edges in use
  std::vector<int> edges;
  for (int iEdgeInit = 0 ; iEdgeInit < tessellation->nEdgesInit ; iEdgeInit++) {
    int iParticleInit1 = tessellation->cellsOnEdgeLocal[iEdgeInit][0];
    int iParticleInit2 = tessellation->cellsOnEdgeLocal[iEdgeInit][1];
    if (iParticleInit1 != -1 and
	iParticleInit2 != -1 and
	useParticle[iParticleInit1] == 1 and
	useParticle[iParticleInit2] == 1) {
      edges.push_back(iEdgeInit);
    }
  } // iEdgeInit

  int nEdges = edges.size();

  // number of owned init particles on each processor
  int *nEdgesOnAllProcs = new int[partition->nprocs()];

  mpiErr = MPI_Allgather(&nEdges, 1, MPI_INT, nEdgesOnAllProcs, 1, MPI_INT, partition->comm());
  MPI_Error_string(mpiErr, mpiErrBuffer, &mpiErrLen);
  log->check(mpiErr == MPI_SUCCESS, "weights_correction: MPI_Allgather failed: ", (std::string) mpiErrBuffer);

  // cumulative number of owned init particles on previous processors
  int *nEdgesOnPrevProc = new int[partition->nprocs()];

  nEdgesOnPrevProc[0] = 0;
  for (int iProc = 1 ; iProc < partition->nprocs() ; iProc++) {
    nEdgesOnPrevProc[iProc] = nEdgesOnPrevProc[iProc-1] + nEdgesOnAllProcs[iProc-1];
  } // iProc

  int nEdgesTotal = 0;
  for (int iProc = 0 ; iProc < partition->nprocs() ; iProc++) {
    nEdgesTotal += nEdgesOnAllProcs[iProc];
  } // iProc

  // mapping from all edges to reduced edges
  std::map<int,int> iEdgeMap;
  for (int iEdge = 0 ; iEdge < nEdges ; iEdge++) {
    iEdgeMap[edges[iEdge]] = iEdge;
  }

  int* iEdgeGlobal = new int[tessellation->nEdgesAllInit];
  for (int iEdgeInit = 0 ; iEdgeInit < tessellation->nParticlesAllInit ; iEdgeInit++) {
    iEdgeGlobal[iEdgeInit] = 0;
  } // iParticleInit
  for (int iEdge = 0 ; iEdge < nEdges ; iEdge++) {
    iEdgeGlobal[edges[iEdge]] = iEdge + nEdgesOnPrevProc[partition->this_proc()];
  } // iParticle
  tessellation->update_init_halo_edges(iEdgeGlobal);

#ifdef DEBUG_REMAPPING
  (*log)(DEMSI::LOG::DEBUG) << "nEdgesTotal: " << nEdgesTotal << std::endl;
  for (int iProc = 0 ; iProc < partition->nprocs() ; iProc++) {
    (*log)(DEMSI::LOG::DEBUG) << "nEdgesOnAllProcs nEdgesOnPrevProc: " << iProc << " ";
    (*log)(DEMSI::LOG::DEBUG) << nEdgesOnAllProcs[iProc] << " " << nEdgesOnPrevProc[iProc] << std::endl;
  } // nParticlesTotal
  (*log)(DEMSI::LOG::DEBUG) << std::endl;
#endif

  double *fluxes = new double[tessellation->nEdgesAllInit];
  for (int iEdgeInit = 0 ; iEdgeInit < tessellation->nEdgesAllInit ; iEdgeInit++) {
    fluxes[iEdgeInit] = 0;
  } // iEdgeInit


    // PERMON QP
  int ierrPermon;

  // Hessian
  Mat P;
  (*log)(DEMSI::LOG::DEBUG) << "Create P matrix" << std::endl;
  ierrPermon = MatCreate(PETSC_COMM_WORLD, &P);
  log->check(ierrPermon == 0, "weights_correction: Permon QP: Problem with MatCreate for P: ", std::to_string(ierrPermon));
  ierrPermon = MatSetType(P, MATMPIAIJ);
  log->check(ierrPermon == 0, "weights_correction: Permon QP: Problem with MatSetType for P: ", std::to_string(ierrPermon));
  ierrPermon = MatSetSizes(P, nEdges, nEdges, PETSC_DETERMINE, PETSC_DETERMINE);
  log->check(ierrPermon == 0, "weights_correction: Permon QP: Problem with MatSetSizes for P: ", std::to_string(ierrPermon));
  ierrPermon = MatMPIAIJSetPreallocation(P, 1, NULL, 0, NULL);
  log->check(ierrPermon == 0, "weights_correction: Permon QP: Problem with MatMPIAIJSetPreallocation for P: ", std::to_string(ierrPermon));
  ierrPermon = MatSetFromOptions(P);
  log->check(ierrPermon == 0, "weights_correction: Permon QP: Problem with MatSetFromOptions for P: ", std::to_string(ierrPermon));
  ierrPermon = MatSetUp(P);
  log->check(ierrPermon == 0, "weights_correction: Permon QP: Problem with MatSetUp for P: ", std::to_string(ierrPermon));

  for (int iEdge = 0 ; iEdge < nEdges ; iEdge++) {
    int iEdgeInit = edges[iEdge];
    PetscInt iEdgeGlobal1 = (PetscInt) iEdgeGlobal[iEdgeInit];
    PetscScalar value = (PetscScalar) 1.0;
    ierrPermon = MatSetValue(P, iEdgeGlobal1, iEdgeGlobal1, value, INSERT_VALUES);
  } // iEdge

  ierrPermon = MatAssemblyBegin(P, MAT_FINAL_ASSEMBLY);
  log->check(ierrPermon == 0, "weights_correction: Permon QP: Problem with MatAssemblyBegin for P: ", std::to_string(ierrPermon));
  ierrPermon = MatAssemblyEnd(P, MAT_FINAL_ASSEMBLY);
  log->check(ierrPermon == 0, "weights_correction: Permon QP: Problem with MatAssemblyEnd for P: ", std::to_string(ierrPermon));

  /*PetscViewer viewer;
    ierrPermon = PetscViewerDrawOpen(PETSC_COMM_WORLD,NULL,NULL,0,0,1200,1200,&viewer);
    ierrPermon = PetscObjectSetName((PetscObject)viewer,"Line graph Plot");
    ierrPermon = PetscViewerPushFormat(viewer,PETSC_VIEWER_DRAW_LG);
    ierrPermon = PetscViewerDrawSetPause(viewer, 10.0);
    ierrPermon = MatView(P, viewer);*/


  // Right hand side
  Vec Q;
  (*log)(DEMSI::LOG::DEBUG) << "Create Q vector" << std::endl;
  ierrPermon = VecCreateMPI(PETSC_COMM_WORLD, nEdges, PETSC_DETERMINE, &Q);
  log->check(ierrPermon == 0, "weights_correction: Permon QP: Problem with VecCreateMPI for Q: ", std::to_string(ierrPermon));
  ierrPermon = VecSet(Q, 0.0);
  log->check(ierrPermon == 0, "weights_correction: Permon QP: Problem with VecSet for Q: ", std::to_string(ierrPermon));
  ierrPermon = VecAssemblyBegin(Q);
  log->check(ierrPermon == 0, "weights_correction: Permon QP: Problem with VecAssemblyBegin for Q: ", std::to_string(ierrPermon));
  ierrPermon = VecAssemblyEnd(Q);
  log->check(ierrPermon == 0, "weights_correction: Permon QP: Problem with VecAssemblyEnd for Q: ", std::to_string(ierrPermon));


  // Inequality constraint matrix
  Mat A;
  (*log)(DEMSI::LOG::DEBUG) << "Create A matrix" << std::endl;
  ierrPermon = MatCreate(PETSC_COMM_WORLD, &A);
  log->check(ierrPermon == 0, "weights_correction: Permon QP: Problem with MatCreate for A: ", std::to_string(ierrPermon));
  ierrPermon = MatSetType(A, MATMPIAIJ);
  log->check(ierrPermon == 0, "weights_correction: Permon QP: Problem with MatSetType for A: ", std::to_string(ierrPermon));
  ierrPermon = MatSetSizes(A, 2*nParticles, nEdges, PETSC_DETERMINE, PETSC_DETERMINE);
  log->check(ierrPermon == 0, "weights_correction: Permon QP: Problem with MatSetSizes for A: ", std::to_string(ierrPermon));
  ierrPermon = MatMPIAIJSetPreallocation(A, 8, NULL, 8, NULL);
  log->check(ierrPermon == 0, "weights_correction: Permon QP: Problem with MatMPIAIJSetPreallocation for A: ", std::to_string(ierrPermon));
  ierrPermon = MatSetFromOptions(A);
  log->check(ierrPermon == 0, "weights_correction: Permon QP: Problem with MatSetFromOptions for A: ", std::to_string(ierrPermon));
  ierrPermon = MatSetUp(A);
  log->check(ierrPermon == 0, "weights_correction: Permon QP: Problem with MatSetUp for A: ", std::to_string(ierrPermon));

  for (int iEdge = 0 ; iEdge < nEdges ; iEdge++) {
    int iEdgeInit = edges[iEdge];
    int iParticleInit1 = tessellation->cellsOnEdgeLocal[iEdgeInit][0];
    int iParticleInit2 = tessellation->cellsOnEdgeLocal[iEdgeInit][1];
    if (iParticleInit1 != -1 and
	iParticleInit2 != -1 and
	useParticle[iParticleInit1] == 1 and
	useParticle[iParticleInit2] == 1) {
      PetscInt iEdgeGlobal1 = (PetscInt) iEdgeGlobal[iEdgeInit];
      PetscInt iParticleGlobal11 = (PetscInt) iParticleGlobal1[iParticleInit1];
      PetscInt iParticleGlobal22 = (PetscInt) iParticleGlobal1[iParticleInit2];
      if (iParticleGlobal11 > iParticleGlobal22) {
	PetscScalar value1 = 1.0 ; PetscScalar value2 = -1.0 ;
	ierrPermon = MatSetValues(A, 1, &iParticleGlobal11, 1, &iEdgeGlobal1, &value1, INSERT_VALUES);
	ierrPermon = MatSetValues(A, 1, &iParticleGlobal22, 1, &iEdgeGlobal1, &value2, INSERT_VALUES);
      } else {
	PetscScalar value1 = -1.0 ; PetscScalar value2 = 1.0 ;
	ierrPermon = MatSetValues(A, 1, &iParticleGlobal11, 1, &iEdgeGlobal1, &value1, INSERT_VALUES);
	ierrPermon = MatSetValues(A, 1, &iParticleGlobal22, 1, &iEdgeGlobal1, &value2, INSERT_VALUES);
      }
      iParticleGlobal11 = (PetscInt) iParticleGlobal2[iParticleInit1];
      iParticleGlobal22 = (PetscInt) iParticleGlobal2[iParticleInit2];
      if (iParticleGlobal11 > iParticleGlobal22) {
	PetscScalar value1 = -1.0 ; PetscScalar value2 = 1.0 ;
	ierrPermon = MatSetValues(A, 1, &iParticleGlobal11, 1, &iEdgeGlobal1, &value1, INSERT_VALUES);
	ierrPermon = MatSetValues(A, 1, &iParticleGlobal22, 1, &iEdgeGlobal1, &value2, INSERT_VALUES);
      } else {
	PetscScalar value1 = 1.0 ; PetscScalar value2 = -1.0 ;
	ierrPermon = MatSetValues(A, 1, &iParticleGlobal11, 1, &iEdgeGlobal1, &value1, INSERT_VALUES);
	ierrPermon = MatSetValues(A, 1, &iParticleGlobal22, 1, &iEdgeGlobal1, &value2, INSERT_VALUES);
      }
    }
  } // iEdgeInit

  ierrPermon = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  log->check(ierrPermon == 0, "weights_correction: Permon QP: Problem with MatAssemblyBegin for A: ", std::to_string(ierrPermon));
  ierrPermon = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
  log->check(ierrPermon == 0, "weights_correction: Permon QP: Problem with MatAssemblyEnd for A: ", std::to_string(ierrPermon));

  /*ierrPermon = PetscViewerDrawOpen(PETSC_COMM_WORLD, NULL, NULL, 0, 0, 1200, 1200, &viewer);
    ierrPermon = PetscObjectSetName((PetscObject) viewer, "Line graph Plot");
    ierrPermon = PetscViewerPushFormat(viewer, PETSC_VIEWER_DRAW_LG);
    ierrPermon = PetscViewerDrawSetPause(viewer, 10.0);
    ierrPermon = MatView(A, viewer);*/


  // Inequality constraint vector
  Vec C;
  (*log)(DEMSI::LOG::DEBUG) << "Create C vector" << std::endl;
  ierrPermon = VecCreateMPI(PETSC_COMM_WORLD, 2*nParticles, PETSC_DETERMINE, &C);
  log->check(ierrPermon == 0, "weights_correction: Permon QP: Problem with VecCreateMPI for C: ", std::to_string(ierrPermon));

  for (int iParticle = 0 ; iParticle < nParticles ; iParticle++) {
    int iParticleInit = iParticleInitUse[iParticle];
    PetscInt iParticleGlobal11 = (PetscInt) iParticleGlobal1[iParticleInit];
    PetscInt iParticleGlobal22 = (PetscInt) iParticleGlobal2[iParticleInit];
    ierrPermon = VecSetValue(C, iParticleGlobal11, (PetscScalar) (elementArea[iParticleInit] - effectiveAreaIce[iParticleInit]), INSERT_VALUES);
    ierrPermon = VecSetValue(C, iParticleGlobal22, (PetscScalar) (                             effectiveAreaIce[iParticleInit]), INSERT_VALUES);
  } // iParticle

  ierrPermon = VecAssemblyBegin(C);
  log->check(ierrPermon == 0, "weights_correction: Permon QP: Problem with VecAssemblyBegin for C: ", std::to_string(ierrPermon));
  ierrPermon = VecAssemblyEnd(C);
  log->check(ierrPermon == 0, "weights_correction: Permon QP: Problem with VecAssemblyEnd for C: ", std::to_string(ierrPermon));

  //PetscViewer viewer;
  /*ierrPermon = PetscViewerDrawOpen(PETSC_COMM_WORLD, NULL, NULL, 0, 0, 1200, 1200, &viewer);
    ierrPermon = PetscObjectSetName((PetscObject) viewer, "Line graph Plot");
    ierrPermon = PetscViewerPushFormat(viewer, PETSC_VIEWER_DRAW_LG);
    ierrPermon = PetscViewerDrawSetPause(viewer, 10.0);
    ierrPermon = VecView(C, viewer);*/


  // Flux solution
  Vec Flux;
  (*log)(DEMSI::LOG::DEBUG) << "Create Flux vector" << std::endl;
  ierrPermon = VecCreateMPI(PETSC_COMM_WORLD, nEdges, PETSC_DETERMINE, &Flux);
  log->check(ierrPermon == 0, "weights_correction: Permon QP: Problem with VecCreateMPI for Flux: ", std::to_string(ierrPermon));
  ierrPermon = VecSet(Flux, 0.0);
  log->check(ierrPermon == 0, "weights_correction: Permon QP: Problem with VecSet for Flux: ", std::to_string(ierrPermon));
  ierrPermon = VecAssemblyBegin(Flux);
  log->check(ierrPermon == 0, "weights_correction: Permon QP: Problem with VecAssemblyBegin for Flux: ", std::to_string(ierrPermon));
  ierrPermon = VecAssemblyEnd(Flux);
  log->check(ierrPermon == 0, "weights_correction: Permon QP: Problem with VecAssemblyEnd for Flux: ", std::to_string(ierrPermon));

  /*PetscViewer viewerIneq;
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, "ineq.bin", FILE_MODE_WRITE, &viewerIneq);
    MatView(A, viewerIneq);
    VecView(C, viewerIneq);
    PetscViewerDestroy(&viewerIneq);*/


  // Empty null space matrix for dualization
  Mat R;
  (*log)(DEMSI::LOG::DEBUG) << "Create R matrix" << std::endl;
  ierrPermon = MatCreateAIJ(PETSC_COMM_WORLD, nEdges, 0, PETSC_DETERMINE, PETSC_DETERMINE, 0, NULL, 0, NULL, &R);
  log->check(ierrPermon == 0, "weights_correction: Permon QP: Problem with MatCreateAIJ for R: ", std::to_string(ierrPermon));
  ierrPermon = MatAssemblyBegin(R, MAT_FINAL_ASSEMBLY);
  log->check(ierrPermon == 0, "weights_correction: Permon QP: Problem with MatAssemblyBegin for R: ", std::to_string(ierrPermon));
  ierrPermon = MatAssemblyEnd(R, MAT_FINAL_ASSEMBLY);
  log->check(ierrPermon == 0, "weights_correction: Permon QP: Problem with MatAssemblyEnd for R: ", std::to_string(ierrPermon));


  // set up permon qp system
  QP qp;
  (*log)(DEMSI::LOG::DEBUG) << "Create QP object" << std::endl;
  ierrPermon = QPCreate(PETSC_COMM_WORLD, &qp);
  log->check(ierrPermon == 0, "weights_correction: Permon QP: Problem with QPCreate: ", std::to_string(ierrPermon));
  ierrPermon = QPSetOperator(qp, P);
  log->check(ierrPermon == 0, "weights_correction: Permon QP: Problem with QPSetOperator: ", std::to_string(ierrPermon));
  ierrPermon = QPSetRhs(qp, Q);
  log->check(ierrPermon == 0, "weights_correction: Permon QP: Problem with QPSetRhs: ", std::to_string(ierrPermon));
  ierrPermon = QPSetInitialVector(qp, Flux);
  log->check(ierrPermon == 0, "weights_correction: Permon QP: Problem with QPSetInitialVector: ", std::to_string(ierrPermon));
  ierrPermon = QPSetIneq(qp, A, C);
  log->check(ierrPermon == 0, "weights_correction: Permon QP: Problem with QPSetIneq: ", std::to_string(ierrPermon));
  ierrPermon = QPSetOperatorNullSpace(qp, R);
  log->check(ierrPermon == 0, "weights_correction: Permon QP: Problem with QPSetOperatorNullSpace: ", std::to_string(ierrPermon));
  ierrPermon = PetscOptionsInsertString(NULL, "-qpt_dualize_B_nest_extension 0 -qpt_dualize_G_explicit 0");
  log->check(ierrPermon == 0, "weights_correction: Permon QP: Problem with PetscOptionsInsertString: ", std::to_string(ierrPermon));
  ierrPermon = PetscOptionsInsertString(NULL, "-dual_mat_inv_ksp_type preonly -dual_mat_inv_pc_type none -dual_mat_inv_ksp_rtol 1e-4");
  log->check(ierrPermon == 0, "weights_correction: Permon QP: Problem with PetscOptionsInsertString: ", std::to_string(ierrPermon));
  ierrPermon = PetscOptionsInsertString(NULL, "-qps_smalxe_monitor -qp_chain_view_kkt");
  log->check(ierrPermon == 0, "weights_correction: Permon QP: Problem with PetscOptionsInsertString: ", std::to_string(ierrPermon));
  ierrPermon = QPTDualize(qp, MAT_INV_MONOLITHIC, MAT_REG_NONE);
  log->check(ierrPermon == 0, "weights_correction: Permon QP: Problem with QPTDualize: ", std::to_string(ierrPermon));
  ierrPermon = QPSetFromOptions(qp);
  log->check(ierrPermon == 0, "weights_correction: Permon QP: Problem with QPSetFromOptions: ", std::to_string(ierrPermon));


  // set up solver
  QPS qps;
  (*log)(DEMSI::LOG::DEBUG) << "Create QPS object" << std::endl;
  ierrPermon = QPSCreate(PetscObjectComm((PetscObject)qp), &qps);
  log->check(ierrPermon == 0, "weights_correction: Permon QP: Problem with QPSCreate: ", std::to_string(ierrPermon));
  ierrPermon = QPSSetQP(qps, qp);
  log->check(ierrPermon == 0, "weights_correction: Permon QP: Problem with QPSSetQP: ", std::to_string(ierrPermon));
  PetscReal rtol = 1e-50;
  PetscReal atol = 1e-5;
  PetscReal dtol = 1e+4;
  PetscInt maxit = PETSC_DEFAULT;
  ierrPermon = QPSSetTolerances(qps, rtol, atol, dtol, maxit);
  log->check(ierrPermon == 0, "weights_correction: Permon QP: Problem with QPSSetTolerances: ", std::to_string(ierrPermon));


  // Solve QP
  (*log)(DEMSI::LOG::DEBUG) << "Solve QPS object" << std::endl;
  ierrPermon = QPSSolve(qps);
  log->check(ierrPermon == 0, "weights_correction: Permon QP: Problem with QPSSolve: ", std::to_string(ierrPermon));


#ifdef DEBUG_REMAPPING
  ierrPermon = QPSView(qps, PETSC_VIEWER_STDOUT_(PETSC_COMM_WORLD));
  log->check(ierrPermon == 0, "weights_correction: Permon QP: Problem with QPSView: ", std::to_string(ierrPermon));
  ierrPermon = QPSViewConvergence(qps, PETSC_VIEWER_STDOUT_(PETSC_COMM_WORLD));
  log->check(ierrPermon == 0, "weights_correction: Permon QP: Problem with QPSViewConvergence: ", std::to_string(ierrPermon));
#endif


  // Check that QPS converged
  (*log)(DEMSI::LOG::DEBUG) << "Check convergence" << std::endl;
  PetscBool converged;
  ierrPermon = QPIsSolved(qp, &converged);
  log->check(ierrPermon == 0, "weights_correction: Permon QP: Problem with QPIsSolved: ", std::to_string(ierrPermon));
  log->check((bool) converged, "weights_correction: Permon QP did not converge.");


#ifdef DEBUG_REMAPPING
  // check feasibility A Flux <= C
  Vec AreaChange;
  (*log)(DEMSI::LOG::DEBUG) << "AreaChange" << std::endl;
  ierrPermon = VecCreateMPI(PETSC_COMM_WORLD, 2*nParticles, PETSC_DETERMINE, &AreaChange);
  log->check(ierrPermon == 0, "weights_correction: Permon QP: Problem with VecCreateMPI for AreaChange: ", std::to_string(ierrPermon));
  ierrPermon = VecAssemblyBegin(AreaChange);
  log->check(ierrPermon == 0, "weights_correction: Permon QP: Problem with VecAssemblyBegin for AreaChange: ", std::to_string(ierrPermon));
  ierrPermon = VecAssemblyEnd(AreaChange);
  log->check(ierrPermon == 0, "weights_correction: Permon QP: Problem with VecAssemblyEnd for AreaChange: ", std::to_string(ierrPermon));

  ierrPermon = MatMult(A, Flux, AreaChange);
  log->check(ierrPermon == 0, "weights_correction: Permon QP: Problem with MatMult for A, Flux and AreaChange: ", std::to_string(ierrPermon));

  int low, high;
  ierrPermon = VecGetOwnershipRange(AreaChange, &low, &high);
  (*log)(DEMSI::LOG::DEBUG) << "low, high: " << low << " " << high << std::endl;

  double totalExcessAreaUpperBeforeLocal = 0.0;
  double totalExcessAreaUpperAfterLocal = 0.0;
  for (int iParticle = 0 ; iParticle < nParticles ; iParticle++) {
    int iParticleInit = iParticleInitUse[iParticle];

    PetscScalar areaChangeUpper;
    PetscInt iParticleGlobal = (PetscInt) iParticleGlobal1[iParticleInit];
    ierrPermon = VecGetValues(AreaChange, 1, &iParticleGlobal, &areaChangeUpper);

    double excessAreaUpperBefore  = std::max(0.0,effectiveAreaIce[iParticleInit] - elementArea[iParticleInit]);
    totalExcessAreaUpperBeforeLocal  += excessAreaUpperBefore ;

    double excessAreaUpperAfter  = std::max(0.0,effectiveAreaIce[iParticleInit] + (double) areaChangeUpper - elementArea[iParticleInit]);
    totalExcessAreaUpperAfterLocal  += excessAreaUpperAfter ;
  } // iParticle

  double totalExcessAreaLowerBeforeLocal = 0.0;
  double totalExcessAreaLowerAfterLocal = 0.0;
  for (int iParticle = 0 ; iParticle < nParticles ; iParticle++) {
    int iParticleInit = iParticleInitUse[iParticle];

    PetscScalar areaChangeLower;
    PetscInt iParticleGlobal = (PetscInt) iParticleGlobal2[iParticleInit];
    ierrPermon = VecGetValues(AreaChange, 1, &iParticleGlobal, &areaChangeLower);
    areaChangeLower = -areaChangeLower;

    double excessAreaLowerBefore = std::min(0.0,effectiveAreaIce[iParticleInit]);
    totalExcessAreaLowerBeforeLocal += excessAreaLowerBefore;

    double excessAreaLowerAfter = std::min(0.0,effectiveAreaIce[iParticleInit] + (double) areaChangeLower);
    totalExcessAreaLowerAfterLocal += excessAreaLowerAfter;
  } // iParticle

  double totalExcessAreaUpperBefore;
  double totalExcessAreaUpperAfter;
  double totalExcessAreaLowerBefore;
  double totalExcessAreaLowerAfter;

  mpiErr = MPI_Allreduce(&totalExcessAreaUpperBeforeLocal, &totalExcessAreaUpperBefore, 1, MPI_DOUBLE, MPI_SUM, partition->comm());
  MPI_Error_string(mpiErr, mpiErrBuffer, &mpiErrLen);
  log->check(mpiErr == MPI_SUCCESS, "weights_correction: MPI_Allreduce totalExcessAreaUpperBefore failed: ", (std::string) mpiErrBuffer);

  mpiErr = MPI_Allreduce(&totalExcessAreaUpperAfterLocal, &totalExcessAreaUpperAfter, 1, MPI_DOUBLE, MPI_SUM, partition->comm());
  MPI_Error_string(mpiErr, mpiErrBuffer, &mpiErrLen);
  log->check(mpiErr == MPI_SUCCESS, "weights_correction: MPI_Allreduce totalExcessAreaUpperAfter failed: ", (std::string) mpiErrBuffer);

  mpiErr = MPI_Allreduce(&totalExcessAreaLowerBeforeLocal, &totalExcessAreaLowerBefore, 1, MPI_DOUBLE, MPI_SUM, partition->comm());
  MPI_Error_string(mpiErr, mpiErrBuffer, &mpiErrLen);
  log->check(mpiErr == MPI_SUCCESS, "weights_correction: MPI_Allreduce totalExcessAreaLowerBefore failed: ", (std::string) mpiErrBuffer);

  mpiErr = MPI_Allreduce(&totalExcessAreaLowerAfterLocal, &totalExcessAreaLowerAfter, 1, MPI_DOUBLE, MPI_SUM, partition->comm());
  MPI_Error_string(mpiErr, mpiErrBuffer, &mpiErrLen);
  log->check(mpiErr == MPI_SUCCESS, "weights_correction: MPI_Allreduce totalExcessAreaLowerAfter failed: ", (std::string) mpiErrBuffer);

  (*log)(DEMSI::LOG::DEBUG) << "Permon solution:" << std::endl;
  (*log)(DEMSI::LOG::DEBUG) << "   Excess areas Before: Lower: " << totalExcessAreaLowerBefore << ", Upper: " << totalExcessAreaUpperBefore << std::endl;
  (*log)(DEMSI::LOG::DEBUG) << "   Excess areas After:  Lower: " << totalExcessAreaLowerAfter  << ", Upper: " << totalExcessAreaUpperAfter  << std::endl;
  (*log)(DEMSI::LOG::DEBUG) << std::endl;

  ierrPermon = VecDestroy(&AreaChange);
  log->check(ierrPermon == 0, "weights_correction: Permon QP: Problem with VecDestroy for AreaChange: ", std::to_string(ierrPermon));
#endif


  // extract primal solution
  (*log)(DEMSI::LOG::DEBUG) << "Extract flux" << std::endl;
  for (int iEdge = 0 ; iEdge < nEdges ; iEdge++) {
    int iEdgeInit = edges[iEdge];
    PetscInt iEdgeGlobal1 = (PetscInt) iEdgeGlobal[iEdgeInit];
    PetscScalar flux;
    ierrPermon = VecGetValues(Flux, 1, &iEdgeGlobal1, &flux);
    fluxes[edges[iEdge]] = (double) flux;
  } // iEdge
  tessellation->update_init_halo_edges(fluxes);


  // clean up
  ierrPermon = MatDestroy(&P);
  log->check(ierrPermon == 0, "weights_correction: Permon QP: Problem with MatDestroy for P: ", std::to_string(ierrPermon));
  ierrPermon = VecDestroy(&Q);
  log->check(ierrPermon == 0, "weights_correction: Permon QP: Problem with VecDestroy for Q: ", std::to_string(ierrPermon));
  ierrPermon = MatDestroy(&A);
  log->check(ierrPermon == 0, "weights_correction: Permon QP: Problem with MatDestroy for A: ", std::to_string(ierrPermon));
  ierrPermon = VecDestroy(&C);
  log->check(ierrPermon == 0, "weights_correction: Permon QP: Problem with VecDestroy for C: ", std::to_string(ierrPermon));
  ierrPermon = MatDestroy(&R);
  log->check(ierrPermon == 0, "weights_correction: Permon QP: Problem with MAtDestroy for R: ", std::to_string(ierrPermon));
  ierrPermon = VecDestroy(&Flux);
  log->check(ierrPermon == 0, "weights_correction: Permon QP: Problem with VecDestroy for Flux: ", std::to_string(ierrPermon));

  ierrPermon = QPDestroy(&qp);
  log->check(ierrPermon == 0, "weights_correction: Permon QP: Problem with QPDestroy: ", std::to_string(ierrPermon));
  ierrPermon = QPSDestroy(&qps);
  log->check(ierrPermon == 0, "weights_correction: Permon QP: Problem with QPSDestroy: ", std::to_string(ierrPermon));


#ifdef DEBUG_REMAPPING

  // updated effective area
  double* effectiveAreaIceDebug = new double[tessellation->nParticlesAllInit];
  for (int iParticleInit = 0 ; iParticleInit < tessellation->nParticlesInit ; iParticleInit++) {
    effectiveAreaIceDebug[iParticleInit] = effectiveAreaIce[iParticleInit];
  } // iParticleInit

  // correct areas with fluxes
  for (int iEdgeInit = 0 ; iEdgeInit < tessellation->nEdgesAllInit ; iEdgeInit++) {
    int iParticleInit1 = tessellation->cellsOnEdgeLocal[iEdgeInit][0];
    int iParticleInit2 = tessellation->cellsOnEdgeLocal[iEdgeInit][1];
    if (iParticleInit1 != -1 and
	iParticleInit2 != -1 and
	useParticle[iParticleInit1] == 1 and
	useParticle[iParticleInit2] == 1) {
      int iParticleInitGlobal1 = (PetscInt) iParticleGlobal1[iParticleInit1];
      int iParticleInitGlobal2 = (PetscInt) iParticleGlobal1[iParticleInit2];
      if (iParticleInitGlobal1 > iParticleInitGlobal2) {
	effectiveAreaIceDebug[iParticleInit1] += fluxes[iEdgeInit];
	effectiveAreaIceDebug[iParticleInit2] -= fluxes[iEdgeInit];
      } else {
	effectiveAreaIceDebug[iParticleInit1] -= fluxes[iEdgeInit];
	effectiveAreaIceDebug[iParticleInit2] += fluxes[iEdgeInit];
      }
    }
  } // iEdge

  // final excess and total area
  double excessAreaEndLocal = 0.0;
  double totalAreaEndLocal  = 0.0;
  for (int iParticleInit = 0 ; iParticleInit < tessellation->nParticlesInit ; iParticleInit++) {
    excessAreaEndLocal += std::max(effectiveAreaIceDebug[iParticleInit] - elementArea[iParticleInit],0.0);
    totalAreaEndLocal  += effectiveAreaIceDebug[iParticleInit];
  } // iParticleInit

  double excessAreaEnd;
  double totalAreaEnd;

  mpiErr = MPI_Allreduce(&excessAreaEndLocal, &excessAreaEnd, 1, MPI_DOUBLE, MPI_SUM, partition->comm());
  MPI_Error_string(mpiErr, mpiErrBuffer, &mpiErrLen);
  log->check(mpiErr == MPI_SUCCESS, "weights_correction: MPI_Allreduce excessAreaEnd failed: ", (std::string) mpiErrBuffer);

  mpiErr = MPI_Allreduce(&totalAreaEndLocal, &totalAreaEnd, 1, MPI_DOUBLE, MPI_SUM, partition->comm());
  MPI_Error_string(mpiErr, mpiErrBuffer, &mpiErrLen);
  log->check(mpiErr == MPI_SUCCESS, "weights_correction: MPI_Allreduce totalAreaEnd failed: ", (std::string) mpiErrBuffer);

  (*log)(DEMSI::LOG::DEBUG) << "Correction solution: " << std::endl;

  (*log)(DEMSI::LOG::DEBUG) << "   Initial Excess Area: " << excessAreaStart;
  (*log)(DEMSI::LOG::DEBUG) << "   Final Excess Area:   " << excessAreaEnd;
  (*log)(DEMSI::LOG::DEBUG) << "   Ratio Excess Area:   " << excessAreaEnd / excessAreaStart << std::endl;

  (*log)(DEMSI::LOG::DEBUG) << "   Initial Total Area:  " << totalAreaStart;
  (*log)(DEMSI::LOG::DEBUG) << "   Final Total Area:    " << totalAreaEnd;
  (*log)(DEMSI::LOG::DEBUG) << "   Diff Total Area:     " << totalAreaEnd - totalAreaStart << std::endl;

  (*log)(DEMSI::LOG::DEBUG) << std::endl;

  delete [] effectiveAreaIceDebug;

#endif



  // Other fluxes
  double* fluxesLeft = new double[tessellation->nEdgesAllInit];
  double* fluxesPrev = new double[tessellation->nEdgesAllInit];
  double* fluxesUsed = new double[tessellation->nEdgesAllInit];
  for (int iEdgeInit = 0 ; iEdgeInit < tessellation->nEdgesAllInit ; iEdgeInit++) {
    fluxesLeft[iEdgeInit] = fluxes[iEdgeInit];
    fluxesPrev[iEdgeInit] = fluxes[iEdgeInit];
    fluxesUsed[iEdgeInit] = fluxes[iEdgeInit];
  } // iEdge



  double* effectiveAreaIcePrev = new double[tessellation->nParticlesAllInit];
  double* effectiveAreaIceNew = new double[tessellation->nParticlesAllInit];
  for (int iParticleInit = 0 ; iParticleInit < tessellation->nParticlesAllInit ; iParticleInit++) {
    effectiveAreaIcePrev[iParticleInit] = effectiveAreaIce[iParticleInit];
    effectiveAreaIceNew [iParticleInit] = effectiveAreaIce[iParticleInit];
  } // iParticleInit



  double totalFluxLeftLocal = 0.0;
  for (int iEdgeInit = 0 ; iEdgeInit < tessellation->nEdgesInit ; iEdgeInit++) {
    totalFluxLeftLocal += std::abs(fluxesLeft[iEdgeInit]);
  } // iEdgeInit
  double totalFluxLeft;
  mpiErr = MPI_Allreduce(&totalFluxLeftLocal, &totalFluxLeft, 1, MPI_DOUBLE, MPI_SUM, partition->comm());
  MPI_Error_string(mpiErr, mpiErrBuffer, &mpiErrLen);
  log->check(mpiErr == MPI_SUCCESS, "weights_correction: MPI_Allreduce failed for totalFluxLeft: ", (std::string) mpiErrBuffer);
#ifdef DEBUG_REMAPPING
  (*log)(DEMSI::LOG::DEBUG) << "Iteration:" << std::endl;
  (*log)(DEMSI::LOG::DEBUG) << "   Iteration number: " << 0 << ", total flux left: " << totalFluxLeft << std::endl;
#endif


  double *fluxRatioUsed = new double[tessellation->nParticlesAllInit];
  double *diagonalValue = new double[tessellation->nParticlesAllInit];

  int ierrPetsc;

  int it = 0;
  while (totalFluxLeft > 1.0) {
    it++;
    log->check(it<100, "weights_correction: flux correction did not converge");

    // diagonal
    for (int i = 0 ; i < tessellation->nParticlesAllInit ; i++) {
      diagonalValue[i] = 1.0;
    } // i

    // off-diagonal
    std::vector<int> iOffDiagonal;
    std::vector<int> jOffDiagonal;
    std::vector<double> offDiagonalValue;

    for (int iParticleInit = 0 ; iParticleInit < tessellation->nParticlesInit ; iParticleInit++) {

      // total outward flux from cell
      fluxRatioUsed[iParticleInit] = 0.0;
      for (int iCellOnCell = 0 ; iCellOnCell < tessellation->nVerticesOnCellInit[iParticleInit] ; iCellOnCell++) {
	int iParticleInitNeighbour = tessellation->cellsOnCellLocal[iParticleInit][iCellOnCell];
	if (iParticleInitNeighbour != -1 and
	    useParticle[iParticleInit] == 1 and
	    useParticle[iParticleInitNeighbour] == 1) {
	  int iEdgeInit = tessellation->edgesOnCellLocal[iParticleInit][iCellOnCell];
	  int iParticleInitGlobal          = tessellation->localToGlobalIndicesInit[iParticleInit];
	  int iParticleInitNeighbourGlobal = tessellation->localToGlobalIndicesInit[iParticleInitNeighbour];
	  if ((iParticleInitGlobal > iParticleInitNeighbourGlobal and fluxesLeft[iEdgeInit] < 0.0) or
	      (iParticleInitGlobal < iParticleInitNeighbourGlobal and fluxesLeft[iEdgeInit] > 0.0)) {
	    fluxRatioUsed[iParticleInit] += std::fabs(fluxesPrev[iEdgeInit]);
	  }
	}
      } // iCellOnCell
      if (fluxRatioUsed[iParticleInit] > 0.0) {
	fluxRatioUsed[iParticleInit] = std::min(effectiveAreaIcePrev[iParticleInit] / fluxRatioUsed[iParticleInit], 1.0);
      } else{
	fluxRatioUsed[iParticleInit] = 1.0;
      }

    } // iParticleInit
    tessellation->update_init_halo_particles(fluxRatioUsed);

    for (int iEdgeInit = 0 ; iEdgeInit < tessellation->nEdgesInit ; iEdgeInit++) {

      int iParticleInit1 = tessellation->cellsOnEdgeLocal[iEdgeInit][0];
      int iParticleInit2 = tessellation->cellsOnEdgeLocal[iEdgeInit][1];

      if (iParticleInit1 != -1 and
	  iParticleInit2 != -1) {

	if (useParticle[iParticleInit1] == 1 and
	    useParticle[iParticleInit2] == 1) {

	  double ratioUsed = std::min(fluxRatioUsed[iParticleInit1],fluxRatioUsed[iParticleInit2]);
	  fluxesUsed[iEdgeInit] *= ratioUsed;
	}

      }

    } // iEdgeInit
    tessellation->update_init_halo_edges(fluxesUsed);

    for (int iEdgeInit = 0 ; iEdgeInit < tessellation->nEdgesAllInit ; iEdgeInit++) {

      int iParticleInit1 = tessellation->cellsOnEdgeLocal[iEdgeInit][0];
      int iParticleInit2 = tessellation->cellsOnEdgeLocal[iEdgeInit][1];

      if (iParticleInit1 != -1 and
	  iParticleInit2 != -1 and
	  useParticle[iParticleInit1] == 1 and
	  useParticle[iParticleInit2] == 1) {

	int iParticleInitGlobal1 = tessellation->localToGlobalIndicesInit[iParticleInit1];
	int iParticleInitGlobal2 = tessellation->localToGlobalIndicesInit[iParticleInit2];

	double fluxAbs = std::fabs(fluxesUsed[iEdgeInit]);

	if ((iParticleInitGlobal1 > iParticleInitGlobal2 and fluxesUsed[iEdgeInit] < 0.0) or
	    (iParticleInitGlobal1 < iParticleInitGlobal2 and fluxesUsed[iEdgeInit] > 0.0)) {
	  // flux from iParticleInitGlobal1 to iParticleInitGlobal2

	  if (effectiveAreaIcePrev[iParticleInit1] > 0.0) {

	    // update area from restricted flux
	    effectiveAreaIceNew[iParticleInit1] -= fluxAbs;
	    effectiveAreaIceNew[iParticleInit2] += fluxAbs;

	    // diagonal of correction matrix
	    if (iParticleInit1 < tessellation->nParticlesInit) {

	      diagonalValue[iParticleInit1] -= fluxAbs / effectiveAreaIcePrev[iParticleInit1];

	      // off diagonal
	      iOffDiagonal.push_back(iParticleInitGlobal2);
	      jOffDiagonal.push_back(iParticleInitGlobal1);
	      offDiagonalValue.push_back(fluxAbs / effectiveAreaIcePrev[iParticleInit1]);

	    }

	  }

	} else {
	  // flux from iParticleInitGlobal2 to iParticleInitGlobal1

	  if (effectiveAreaIcePrev[iParticleInit2] > 0.0) {

	    // update area from restricted flux
	    effectiveAreaIceNew[iParticleInit2] -= fluxAbs;
	    effectiveAreaIceNew[iParticleInit1] += fluxAbs;

	    // diagonal of correction matrix
	    if (iParticleInit2 < tessellation->nParticlesInit) {

	      diagonalValue[iParticleInit2] -= fluxAbs / effectiveAreaIcePrev[iParticleInit2];

	      // off diagonal
	      iOffDiagonal.push_back(iParticleInitGlobal1);
	      jOffDiagonal.push_back(iParticleInitGlobal2);
	      offDiagonalValue.push_back(fluxAbs / effectiveAreaIcePrev[iParticleInit2]);

	    }

	  }

	} // useParticle

      } // not boundary

    } // iEdgeInit


    for (int iEdgeInit = 0 ; iEdgeInit < tessellation->nEdgesInit ; iEdgeInit++) {
      fluxesLeft[iEdgeInit] -= fluxesUsed[iEdgeInit];
    } // iEdgeInit

    totalFluxLeftLocal = 0.0;
    for (int iEdgeInit = 0 ; iEdgeInit < tessellation->nEdgesInit ; iEdgeInit++) {
      totalFluxLeftLocal += std::abs(fluxesLeft[iEdgeInit]);
    } // iEdge
    mpiErr = MPI_Allreduce(&totalFluxLeftLocal, &totalFluxLeft, 1, MPI_DOUBLE, MPI_SUM, partition->comm());
    MPI_Error_string(mpiErr, mpiErrBuffer, &mpiErrLen);
    log->check(mpiErr == MPI_SUCCESS, "weights_correction: MPI_Allreduce failed for totalFluxLeft: ", (std::string) mpiErrBuffer);
#ifdef DEBUG_REMAPPING
    (*log)(DEMSI::LOG::DEBUG) << "   Iteration number: " << it << ", total flux left: " << totalFluxLeft << std::endl;
#endif



    tessellation->update_init_halo_edges(fluxesLeft);
    for (int iEdgeInit = 0 ; iEdgeInit < tessellation->nEdgesAllInit ; iEdgeInit++) {
      fluxesPrev[iEdgeInit] = fluxesLeft[iEdgeInit];
      fluxesUsed[iEdgeInit] = fluxesLeft[iEdgeInit];
    } // iEdge

    tessellation->update_init_halo_particles(effectiveAreaIceNew);
    for (int iParticleInit = 0 ; iParticleInit < tessellation->nParticlesAllInit ; iParticleInit++) {
      effectiveAreaIcePrev[iParticleInit] = effectiveAreaIceNew[iParticleInit];
    } // iParticleInit



    // assemble correction matrix
    Mat U;
    (*log)(DEMSI::LOG::DEBUG) << "Create U matrix" << std::endl;
    ierrPetsc = MatCreate(PETSC_COMM_WORLD, &U);
    log->check(ierrPetsc == 0, "weights_correction: Petsc: Problem with MatCreate for U: ", std::to_string(ierrPetsc));
    ierrPetsc = MatSetType(U, MATMPIAIJ);
    log->check(ierrPetsc == 0, "weights_correction: Petsc: Problem with MatSetType for U: ", std::to_string(ierrPetsc));
    ierrPetsc = MatSetSizes(U, (PetscInt) tessellation->nParticlesInit, (PetscInt) tessellation->nParticlesInit, PETSC_DETERMINE, PETSC_DETERMINE);
    log->check(ierrPetsc == 0, "weights_correction: Petsc: Problem with MatSetSizes for U: ", std::to_string(ierrPetsc));
    ierrPetsc = MatMPIAIJSetPreallocation(U, 20, NULL, 20, NULL);
    log->check(ierrPetsc == 0, "weights_correction: Petsc: Problem with MatMPIAIJSetPreallocation for U: ", std::to_string(ierrPetsc));
    ierrPetsc = MatSetFromOptions(U);
    log->check(ierrPetsc == 0, "weights_correction: Petsc: Problem with MatSetFromOptions for U: ", std::to_string(ierrPetsc));
    ierrPetsc = MatSetUp(U);
    log->check(ierrPetsc == 0, "weights_correction: Petsc: Problem with MatSetUp for U: ", std::to_string(ierrPetsc));

    // diagonal
    for (int iParticleInit = 0 ; iParticleInit < tessellation->nParticlesInit ; iParticleInit++) {
      PetscInt iParticleInitGlobal = (PetscInt) tessellation->localToGlobalIndicesInit[iParticleInit];
      PetscScalar value = (PetscScalar) diagonalValue[iParticleInit];
      ierrPetsc = MatSetValues(U, 1, &iParticleInitGlobal, 1, &iParticleInitGlobal, &value, INSERT_VALUES);
    } // iParticleInit

    // off diagnonal
    for (int ij = 0 ; ij < offDiagonalValue.size() ; ij++) {
      PetscInt i = (PetscInt) iOffDiagonal[ij];
      PetscInt j = (PetscInt) jOffDiagonal[ij];
      PetscScalar value = (PetscScalar) offDiagonalValue[ij];
      ierrPetsc = MatSetValues(U, 1, &i, 1, &j, &value, INSERT_VALUES);
    }

    ierrPetsc = MatAssemblyBegin(U, MAT_FINAL_ASSEMBLY);
    log->check(ierrPetsc == 0, "weights_correction: Petsc: Problem with MatAssemblyBegin for U: ", std::to_string(ierrPetsc));
    ierrPetsc = MatAssemblyEnd(U, MAT_FINAL_ASSEMBLY);
    log->check(ierrPetsc == 0, "weights_correction: Petsc: Problem with MatAssemblyEnd for U: ", std::to_string(ierrPetsc));

    /*PetscViewer viewer;
      ierrPetsc = PetscViewerDrawOpen(PETSC_COMM_WORLD, NULL, NULL, 0, 0, 1200, 1200, &viewer);
      ierrPetsc = PetscObjectSetName((PetscObject) viewer, "Line graph Plot");
      ierrPetsc = PetscViewerPushFormat(viewer, PETSC_VIEWER_DRAW_LG);
      ierrPetsc = PetscViewerDrawSetPause(viewer, 20.0);
      ierrPetsc = MatView(Weights, viewer);

      ierrPetsc = PetscViewerDrawOpen(PETSC_COMM_WORLD, NULL, NULL, 0, 0, 1200, 1200, &viewer);
      ierrPetsc = PetscObjectSetName((PetscObject) viewer, "Line graph Plot");
      ierrPetsc = PetscViewerPushFormat(viewer, PETSC_VIEWER_DRAW_LG);
      ierrPetsc = PetscViewerDrawSetPause(viewer, 20.0);
      ierrPetsc = MatView(U, viewer);*/

    // matrix multiply W * U
    Mat UW;
    (*log)(DEMSI::LOG::DEBUG) << "MatMatMult UW" << std::endl;
    ierrPetsc = MatMatMult(U, Weights, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &UW);
    log->check(ierrPetsc == 0, "weights_correction: Petsc: Problem with MatMatMult for U and W: ", std::to_string(ierrPetsc));

    /*ierrPetsc = PetscViewerDrawOpen(PETSC_COMM_WORLD, NULL, NULL, 0, 0, 1200, 1200, &viewer);
      ierrPetsc = PetscObjectSetName((PetscObject) viewer, "Line graph Plot");
      ierrPetsc = PetscViewerPushFormat(viewer, PETSC_VIEWER_DRAW_LG);
      ierrPetsc = PetscViewerDrawSetPause(viewer, 20.0);
      ierrPetsc = MatView(UW, viewer);*/

    (*log)(DEMSI::LOG::DEBUG) << "MatCopy UW to Weights" << std::endl;
    ierrPetsc = MatCopy(UW, Weights, DIFFERENT_NONZERO_PATTERN);
    log->check(ierrPetsc == 0, "weights_correction: Petsc: Problem with MatCopy for UW and W: ", std::to_string(ierrPetsc));

    /*ierrPetsc = PetscViewerDrawOpen(PETSC_COMM_WORLD, NULL, NULL, 0, 0, 1200, 1200, &viewer);
      ierrPetsc = PetscObjectSetName((PetscObject) viewer, "Line graph Plot");
      ierrPetsc = PetscViewerPushFormat(viewer, PETSC_VIEWER_DRAW_LG);
      ierrPetsc = PetscViewerDrawSetPause(viewer, 20.0);
      ierrPetsc = MatView(Weights, viewer);*/

    // destroy matrices
    ierrPetsc = MatDestroy(&U);
    log->check(ierrPetsc == 0, "weights_correction: Petsc: Problem with MatDestroy for U: ", std::to_string(ierrPetsc));
    ierrPetsc = MatDestroy(&UW);
    log->check(ierrPetsc == 0, "weights_correction: Petsc: Problem with MatDestroy for UW: ", std::to_string(ierrPetsc));



  } // it
#ifdef DEBUG_REMAPPING
  (*log)(DEMSI::LOG::DEBUG) << std::endl;
#endif

#ifdef DEBUG_REMAPPING
  double *effectiveAreaIceDebugRemap = new double[tessellation->nParticlesAllInit];
  remap_variable(1, effectiveAreaIceDebugRemap, &effectiveElementArea[0], Weights);

  // final excess and total area
  excessAreaEndLocal = 0.0;
  totalAreaEndLocal  = 0.0;
  for (int iParticleInit = 0 ; iParticleInit < tessellation->nParticlesInit ; iParticleInit++) {
    excessAreaEndLocal += std::max(effectiveAreaIceDebugRemap[iParticleInit] - elementArea[iParticleInit],0.0);
    totalAreaEndLocal  += effectiveAreaIceDebugRemap[iParticleInit];
  } // iParticleInit

  mpiErr = MPI_Allreduce(&excessAreaEndLocal, &excessAreaEnd, 1, MPI_DOUBLE, MPI_SUM, partition->comm());
  MPI_Error_string(mpiErr, mpiErrBuffer, &mpiErrLen);
  log->check(mpiErr == MPI_SUCCESS, "weights_correction: MPI_Allreduce excessAreaEnd failed: ", (std::string) mpiErrBuffer);

  mpiErr = MPI_Allreduce(&totalAreaEndLocal, &totalAreaEnd, 1, MPI_DOUBLE, MPI_SUM, partition->comm());
  MPI_Error_string(mpiErr, mpiErrBuffer, &mpiErrLen);
  log->check(mpiErr == MPI_SUCCESS, "weights_correction: MPI_Allreduce totalAreaEnd failed: ", (std::string) mpiErrBuffer);

  (*log)(DEMSI::LOG::DEBUG) << "Weight matrix solution after correction:" << std::endl;

  (*log)(DEMSI::LOG::DEBUG) << "   Initial Excess Area: " << excessAreaStart;
  (*log)(DEMSI::LOG::DEBUG) << "   Final Excess Area:   " << excessAreaEnd;
  (*log)(DEMSI::LOG::DEBUG) << "   Ratio Excess Area:   " << excessAreaEnd / excessAreaStart << std::endl;

  (*log)(DEMSI::LOG::DEBUG) << "   Initial Total Area:  " << totalAreaStart;
  (*log)(DEMSI::LOG::DEBUG) << "   Final Total Area:    " << totalAreaEnd;
  (*log)(DEMSI::LOG::DEBUG) << "   Diff Total Area:     " << totalAreaEnd - totalAreaStart << std::endl;

  (*log)(DEMSI::LOG::DEBUG) << std::endl;

  delete [] effectiveAreaIceDebugRemap;
#endif

  // cleanup
  delete [] effectiveAreaIce;
  delete [] elementArea;
  delete [] useParticle;
  delete [] nParticlesOnAllProcs;
  delete [] nParticlesOnPrevProc;
  delete [] iParticleGlobal;
  delete [] nEdgesOnAllProcs;
  delete [] nEdgesOnPrevProc;
  delete [] iEdgeGlobal;
  delete [] fluxes;
  delete [] fluxRatioUsed;
  delete [] fluxesPrev;
  delete [] fluxesUsed;
  delete [] fluxesLeft;
  delete [] diagonalValue;

} // Remapping::weights_correction

} // namespace DEMSI
