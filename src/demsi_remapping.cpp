#include "demsi_remapping.h"
#include "demsi_communication.h"
#include "demsi_column.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <set>
#include <numeric>
#include <iomanip>
#include <algorithm>
#include <random>

#include <lammps/atom.h>
#include <lammps/math_const.h>
#include <lammps/domain.h>
#include <lammps/update.h>
#include <lammps/neighbor.h>

namespace DEMSI {

//------------------------------------------------------------------------------
// Remapping tracer class
//------------------------------------------------------------------------------

RemappingTracer::RemappingTracer(DEMSI::Tessellation* tessellationIn, DEMSI::ColumnVariable<double>* tracerIn) {

  tessellation = tessellationIn;
  tracer = tracerIn;

  int sizePerParticle = tracer->size_per_particle();

  // arrays for old particles
  int nSizeOld = tracer->num_particles()*sizePerParticle;
  extensiveTracer.resize(nSizeOld);

  for (int i = 0 ; i < nSizeOld ; i++) {
    extensiveTracer[i] = 0.0;
  } // i

  // arrays for potential new particles
  int nSizeNew = tessellation->nParticlesAllInit*sizePerParticle;
  extensiveTracerRemap.resize(nSizeNew);
  intensiveTracerRemap.resize(nSizeNew);

  for (int i = 0 ; i < nSizeNew ; i++) {
    extensiveTracerRemap[i] = 0.0;
    intensiveTracerRemap[i] = 0.0;
  } // i

};

//------------------------------------------------------------------------------
// Initialization
//------------------------------------------------------------------------------

// Remapping object constructor
Remapping::Remapping(DEMSI::Partition* partitionIn, DEMSI::Log* logIn, DEMSI::Configs* configsIn, DEMSI::LammpsInstance* lammpsInstanceIn, DEMSI::Grid* gridIn, DEMSI::Contacts* contactsIn, DEMSI::Particles* particlesIn, DEMSI::Tessellation* tessellationIn, DEMSI::TessellationOutputStreams* tessellationOutputStreamsIn, DEMSI::Column* columnIn, DEMSI::Ocean* oceanIn, DEMSI::Clock* simulationClock) {

  partition = partitionIn;
  log = logIn;
  configs = configsIn;
  lammpsInstance = lammpsInstanceIn;
  grid = gridIn;
  contacts = contactsIn;
  particles = particlesIn;
  tessellation = tessellationIn;
  tessellationOutputStreams = tessellationOutputStreamsIn;
  column = columnIn;
  ocean = oceanIn;

  // perform remapping
  if (configs->exists({"ConfigGroup:remapping"})) {
    configs->get({"ConfigGroup:remapping","Config:useRemapping"}, useRemapping);
  } // remapping group exists

  if (useRemapping) {

    // check conservation
    configs->get({"ConfigGroup:remapping","Config:checkConservation"}, checkConservation);

    // check weights
    configs->get({"ConfigGroup:remapping","Config:checkRemappingWeights"}, checkRemappingWeights);
    log->check(not checkRemappingWeights or (checkRemappingWeights and partition->nprocs() == 1), "Cannot check remapping weights with more than one processor");

    // check effective areas
    configs->get({"ConfigGroup:remapping","Config:checkEffectiveArea"}, checkEffectiveArea);

    // check ice concentration bounds
    configs->get({"ConfigGroup:remapping","Config:checkIceConcentrationBounds"}, checkIceConcentrationBounds);

    // max error allowed on ice concentration bounds
    configs->get({"ConfigGroup:remapping","Config:maxIceConcBoundError"}, maxIceConcBoundError);

    // remapping correction type
    configs->get({"ConfigGroup:remapping","Config:remappingCorrectionType"}, remappingCorrectionType);

    // check for optimization-based remap
    if (remappingCorrectionType == "obr") {
       // for now turn off check of effective area if obr is selected
       // effective area is not corrected, iceArea is corrected directly
       checkEffectiveArea = 0.0;
    }

    // set the forcing group alarm
    std::string remappingIntervalStr;
    configs->get({"ConfigGroup:remapping","Config:remappingInterval"}, remappingIntervalStr);

    DEMSI::Time currentTime = simulationClock->get_time();
    DEMSI::TimeInterval remapInterval(remappingIntervalStr, log);

    remappingAlarm.set(simulationClock, currentTime, remapInterval);
    remappingAlarm.reset();

    // initialize the neighbour lists for the initial element distribution
    init_neighbour_lists();

    column->columnDimensions->add("TWO", "The number two", 2);
    column->columnDimensions->add("maxVerticesInit", "max number of vertices", 10);

    column->nVerticesInit = new DEMSI::ColumnVariable<int>(particles, column->columnDimensions,
	"nVerticesInit",
	{},
	true);
    column->columnVariables->add(column->nVerticesInit);
    column->verticesInit = new DEMSI::ColumnVariable<double>(particles, column->columnDimensions,
	"verticesInit",
	{"maxVerticesInit","TWO"},
	true);
    column->columnVariables->add(column->verticesInit);
    column->areaInit = new DEMSI::ColumnVariable<double>(particles, column->columnDimensions,
	"areaInit",
	{},
	true);
    column->columnVariables->add(column->areaInit);
    column->areaInitRemap = new DEMSI::ColumnVariable<double>(particles, column->columnDimensions,
	"areaInitRemap",
	{},
	true);
    column->columnVariables->add(column->areaInitRemap);

    for (int iParticle = 0 ; iParticle < *(particles->nParticles) ; iParticle++) {

      int iParticleInitClosest;
      double minDistance = 1e30;

      for (int iParticleInit = 0 ; iParticleInit < tessellation->nParticlesInit ; iParticleInit++) {

	double dx = particles->x(iParticle,0) - tessellation->xInit[iParticleInit];
	double dy = particles->x(iParticle,1) - tessellation->yInit[iParticleInit];

	double distance = std::sqrt(dx*dx + dy*dy);

	if (distance < minDistance) {
	  minDistance = distance;
	  iParticleInitClosest = iParticleInit;
	}

      } // iParticleInit
      int iParticleInit = iParticleInitClosest;

      log->check(minDistance < 10.0, "Did not find new particle amongst init", std::to_string(minDistance));

      (*column->nVerticesInit)(iParticle) = tessellation->initPolygons[iParticleInit].nvertices();
      (*column->areaInit)(iParticle) = tessellation->initPolygons[iParticleInit].area();
      (*column->areaInitRemap)(iParticle) = tessellation->initPolygons[iParticleInit].area();
      for (int iVertex = 0 ; iVertex < tessellation->initPolygons[iParticleInit].nvertices() ; iVertex++) {
	(*column->verticesInit)(iParticle,iVertex,0) = tessellation->initPolygons[iParticleInit].get_x(iVertex) - tessellation->xInit[iParticleInit];
	(*column->verticesInit)(iParticle,iVertex,1) = tessellation->initPolygons[iParticleInit].get_y(iVertex) - tessellation->yInit[iParticleInit];
      }

    } // iParticle

#ifdef DEBUG_REMAPPING
    // debugging fields
    useParticleDebug = new double[tessellation->nParticlesInit];
    elementAreaDebug = new double[tessellation->nParticlesInit];
    effectiveAreaIceDebug = new double[tessellation->nParticlesInit];

    tessellation->register_output_field(useParticleDebug, "useParticleDebug");
    tessellation->register_output_field(elementAreaDebug, "elementAreaDebug");
    tessellation->register_output_field(effectiveAreaIceDebug, "effectiveAreaIceDebug");
#endif

  } // useRemapping

}

// Initialize neighbour lists for fast finding of new elements
void Remapping::init_neighbour_lists(void) {

  char mpiErrBuffer[MPI_MAX_ERROR_STRING];
  int mpiErrLen;
  int err;

  // find the maximum polygon size - this will define the neighbour grid size
  double maxPolygonSizeLocal = 0.0;
  for (int iParticleInit = 0 ; iParticleInit < tessellation->nParticlesInit ; iParticleInit++) {
    maxPolygonSizeLocal = std::max(maxPolygonSizeLocal, tessellation->initPolygons[iParticleInit].max_size());
  } // iParticleInit

  err = MPI_Reduce(&maxPolygonSizeLocal, &maxPolygonSize, 1, MPI_DOUBLE, MPI_MAX, partition->master_proc(), partition->comm());
  MPI_Error_string(err, mpiErrBuffer, &mpiErrLen);
  log->check(err == MPI_SUCCESS, "Problem getting max polygon size of init particles: ", (std::string) mpiErrBuffer);

  err = MPI_Bcast(&maxPolygonSize, 1, MPI_DOUBLE, partition->master_proc(), partition->comm());
  MPI_Error_string(err, mpiErrBuffer, &mpiErrLen);
  log->check(err == MPI_SUCCESS, "Problem broadcasting max polygon size of init particles: ", (std::string) mpiErrBuffer);

  maxPolygonSize = maxPolygonSize * 1.01;

  // define the neighbour grid
  nxNeighbours = std::ceil(grid->lx() / maxPolygonSize);
  nyNeighbours = std::ceil(grid->ly() / maxPolygonSize);
  dxNeighbours = grid->lx() / (double) nxNeighbours;
  dyNeighbours = grid->ly() / (double) nyNeighbours;

  // set the neighbour lists
  for (int iParticleInit = 0 ; iParticleInit < tessellation->nParticlesAllInit ; iParticleInit++) {
    for (int i = 0 ; i < nxNeighbours ; i++) {
      for (int j = 0 ; j < nyNeighbours ; j++) {
	double xmin = ((double) i + 0.5) * dxNeighbours - 1.5 * dxNeighbours;
	double xmax = ((double) i + 0.5) * dxNeighbours + 1.5 * dxNeighbours;
	double ymin = ((double) j + 0.5) * dyNeighbours - 1.5 * dyNeighbours;
	double ymax = ((double) j + 0.5) * dyNeighbours + 1.5 * dyNeighbours;
	if (tessellation->xInit[iParticleInit] >= xmin and tessellation->xInit[iParticleInit] <= xmax and
	    tessellation->yInit[iParticleInit] >= ymin and tessellation->yInit[iParticleInit] <= ymax) {
	  indexNeighbours[std::make_pair(i,j)].push_back(iParticleInit);
	}
      } // j
    } // i
  } // iParticleInit

}

//------------------------------------------------------------------------------
// Run time
//------------------------------------------------------------------------------

// Calculate the remapping weights for the current particle distribution
  void Remapping::remapping_weights(Mat &Weights, std::vector<double> &weights, std::vector<int> &iParticlesOldWeights, std::vector<int> &iParticlesInitWeights) {

  // loop over old particles that are being remapped
  for (int iParticleOld = 0 ; iParticleOld < *(particles->nParticles) ; iParticleOld++) {

    DEMSI::Polygon oldPolygon;
    for (int iVertex = 0 ; iVertex < (*column->nVerticesInit)(iParticleOld); iVertex++) {
      double xVertex = (*column->verticesInit)(iParticleOld,iVertex,0) + particles->x(iParticleOld,0);
      double yVertex = (*column->verticesInit)(iParticleOld,iVertex,1) + particles->x(iParticleOld,1);
      oldPolygon.add_vertex(xVertex,yVertex);
    } // iVertex

    double polygonAreaOld = oldPolygon.area();

    //PetscInt iParticleOldGlobal = iParticleOldStart + (PetscInt) iParticleOld;

    // loop over the full initial particle distribution
    //for (int iParticleInit = 0 ; iParticleInit < tessellation->nParticlesAllInit ; iParticleInit++) {
    int ix = std::floor(particles->x(iParticleOld,0) / dxNeighbours);
    int iy = std::floor(particles->x(iParticleOld,1) / dyNeighbours);

    for (int iNeighbour = 0 ; iNeighbour < indexNeighbours[std::make_pair(ix,iy)].size() ; iNeighbour++) {

      int iParticleInit = indexNeighbours[std::make_pair(ix,iy)][iNeighbour];

      // coastal particles
      if (particles->type[iParticleOld] == 2 and tessellation->typeInit[iParticleInit] == 2 and
	  std::fabs(particles->x(iParticleOld,0) - tessellation->xInit[iParticleInit]) < 10.0 and
	  std::fabs(particles->x(iParticleOld,1) - tessellation->yInit[iParticleInit]) < 10.0) {

	weights.push_back(1.0);
	iParticlesInitWeights.push_back(iParticleInit);
	iParticlesOldWeights.push_back(iParticleOld);

      // Non coastal particles
      } else if (particles->type[iParticleOld] != 2 and tessellation->typeInit[iParticleInit] != 2) {

	DEMSI::Polygon intersectionPolygon = oldPolygon.intersection(tessellation->initPolygons[iParticleInit]);

	double polygonAreaIntersection = intersectionPolygon.area();

	if (polygonAreaIntersection > 0.0) {

	  double weight = polygonAreaIntersection / polygonAreaOld;

	  weights.push_back(weight);
	  iParticlesInitWeights.push_back(iParticleInit);
	  iParticlesOldWeights.push_back(iParticleOld);

	} // positive intersection

      } // element type

    } // iParticleInit

  } // iParticle

  // rescale for particles overlapping coastline
  // total weight per old particle
  double *totalWeights = new double[*(particles->nParticles)];
  for (int iParticleOld = 0 ; iParticleOld < *(particles->nParticles) ; iParticleOld++) {
    totalWeights[iParticleOld] = 0.0;
  } // iParticleOld
  for (int iWeight = 0 ; iWeight < weights.size() ; iWeight++) {
    totalWeights[iParticlesOldWeights[iWeight]] += weights[iWeight];
  } // iWeight

  // rescale weights so total weight per old particle is one
  for (int iWeight = 0 ; iWeight < weights.size() ; iWeight++) {
    weights[iWeight] /= totalWeights[iParticlesOldWeights[iWeight]];
  } // iWeight

  // create weights matrix
  int ierrPetsc;

  ierrPetsc = MatCreate(PETSC_COMM_WORLD, &Weights);
  log->check(ierrPetsc == 0, "remapping_weights: Petsc: Problem with MatCreate for Weights: ", std::to_string(ierrPetsc));
  ierrPetsc = MatSetType(Weights, MATMPIAIJ);
  log->check(ierrPetsc == 0, "remapping_weights: Petsc: Problem with MatSetType for Weights: ", std::to_string(ierrPetsc));
  ierrPetsc = MatSetSizes(Weights, (PetscInt) tessellation->nParticlesInit, (PetscInt) *(particles->nParticles), PETSC_DETERMINE, PETSC_DETERMINE);
  log->check(ierrPetsc == 0, "remapping_weights: Petsc: Problem with MatSetSizes for Weights: ", std::to_string(ierrPetsc));
  ierrPetsc = MatMPIAIJSetPreallocation(Weights, 20, NULL, 20, NULL);
  log->check(ierrPetsc == 0, "remapping_weights: Petsc: Problem with MatMPIAIJSetPreallocation for Weights: ", std::to_string(ierrPetsc));
  ierrPetsc = MatSetFromOptions(Weights);
  log->check(ierrPetsc == 0, "remapping_weights: Petsc: Problem with MatSetFromOptions for Weights: ", std::to_string(ierrPetsc));
  ierrPetsc = MatSetUp(Weights);
  log->check(ierrPetsc == 0, "remapping_weights: Petsc: Problem with MatSetUp for Weights: ", std::to_string(ierrPetsc));

  PetscInt iRowStart, iRowEnd;
  ierrPetsc = MatGetOwnershipRange(Weights, &iRowStart, &iRowEnd);
  log->check(ierrPetsc == 0, "remapping_weights: Petsc: Problem with MatGetOwnershipRange for Weights: ", std::to_string(ierrPetsc));

  PetscInt iParticleOldStart, iParticleOldEnd;
  ierrPetsc = MatGetOwnershipRangeColumn(Weights, &iParticleOldStart, &iParticleOldEnd);
  log->check(ierrPetsc == 0, "remapping_weights: Petsc: Problem with MatGetOwnershipRangeColumn for Weights: ", std::to_string(ierrPetsc));

  for (int iWeight = 0 ; iWeight < weights.size() ; iWeight++) {

    int iParticleOld  = iParticlesOldWeights [iWeight];
    int iParticleInit = iParticlesInitWeights[iWeight];

    PetscInt iParticleOldGlobal  = iParticleOldStart + (PetscInt) iParticleOld;;
    PetscInt iParticleInitGlobal = (PetscInt) tessellation->localToGlobalIndicesInit[iParticleInit];;

    PetscScalar weight = (PetscScalar) weights[iWeight];

    ierrPetsc = MatSetValues(Weights, 1, &iParticleInitGlobal, 1, &iParticleOldGlobal, &weight, INSERT_VALUES);

  } // iWeight

  ierrPetsc = MatAssemblyBegin(Weights, MAT_FINAL_ASSEMBLY);
  log->check(ierrPetsc == 0, "remapping_weights: Petsc: Problem with MatAssemblyBegin for Weights: ", std::to_string(ierrPetsc));
  ierrPetsc = MatAssemblyEnd(Weights, MAT_FINAL_ASSEMBLY);
  log->check(ierrPetsc == 0, "remapping_weights: Petsc: Problem with MatAssemblyEnd for Weights: ", std::to_string(ierrPetsc));

} // Remapping::remapping_weights

// remap a column variable back to the initial distribution
void Remapping::remap_variable(const int sizePerParticle, double *extensiveTracerRemap, double *extensiveTracer, const Mat Weights) {

  PetscInt nParticlesOld = (PetscInt) *(particles->nParticles);
  PetscInt nParticlesInit = (PetscInt) tessellation->nParticlesInit;

  PetscInt ierrPetsc;

  for (PetscInt iSizePerParticle = 0 ; iSizePerParticle < sizePerParticle ; iSizePerParticle++) {

    // ExtensiveTracer
    Vec ExtensiveTracer;
    ierrPetsc = VecCreateMPI(PETSC_COMM_WORLD, nParticlesOld, PETSC_DETERMINE, &ExtensiveTracer);

    PetscInt iRowStartExtensiveTracer, iRowEndExtensiveTracer;
    ierrPetsc = VecGetOwnershipRange(ExtensiveTracer, &iRowStartExtensiveTracer, &iRowEndExtensiveTracer);

    for (PetscInt iParticleOld = 0 ; iParticleOld < nParticlesOld ; iParticleOld++) {
      PetscInt jOld  = iParticleOld * sizePerParticle + iSizePerParticle;
      PetscInt j = iRowStartExtensiveTracer + iParticleOld;
      ierrPetsc = VecSetValue(ExtensiveTracer, j, (PetscScalar) extensiveTracer[jOld], INSERT_VALUES);
    } // iParticle

    ierrPetsc = VecAssemblyBegin(ExtensiveTracer);
    ierrPetsc = VecAssemblyEnd(ExtensiveTracer);

    // ExtensiveTracerRemap
    Vec ExtensiveTracerRemap;
    ierrPetsc = VecCreateMPI(PETSC_COMM_WORLD, nParticlesInit, PETSC_DETERMINE, &ExtensiveTracerRemap);
    ierrPetsc = VecSet(ExtensiveTracerRemap, 0.0);
    ierrPetsc = VecAssemblyBegin(ExtensiveTracerRemap);
    ierrPetsc = VecAssemblyEnd(ExtensiveTracerRemap);

    PetscInt iRowStartExtensiveTracerRemap, iRowEndExtensiveTracerRemap;
    ierrPetsc = VecGetOwnershipRange(ExtensiveTracerRemap, &iRowStartExtensiveTracerRemap, &iRowEndExtensiveTracerRemap);

    // remap
    ierrPetsc = MatMult(Weights, ExtensiveTracer, ExtensiveTracerRemap);

    // ExtensiveTracerRemap
    for (PetscInt iParticleInit = 0 ; iParticleInit < nParticlesInit ; iParticleInit++) {
      PetscInt jInit = iParticleInit * sizePerParticle + iSizePerParticle;
      PetscScalar value;
      PetscInt j = iRowStartExtensiveTracerRemap + iParticleInit;
      ierrPetsc = VecGetValues(ExtensiveTracerRemap, 1, &j, &value);
      extensiveTracerRemap[jInit] = (double) value;
      if (std::isnan(value)) {
	log->abort("remap_variable is nan");
      }
    } // iParticle

    // destroy vectors
    ierrPetsc = VecDestroy(&ExtensiveTracer);
    ierrPetsc = VecDestroy(&ExtensiveTracerRemap);

  } // iSizePerParticle

} // Remapping::remap_variable

// remap a column variable back to the initial distribution
void Remapping::remap_variable(const int sizePerParticle, std::vector<double> &extensiveTracerRemap, std::vector<double> extensiveTracer, const Mat Weights) {

  PetscInt nParticlesOld = (PetscInt) *(particles->nParticles);
  PetscInt nParticlesInit = (PetscInt) tessellation->nParticlesInit;

  PetscInt ierrPetsc;

  for (PetscInt iSizePerParticle = 0 ; iSizePerParticle < sizePerParticle ; iSizePerParticle++) {

    // ExtensiveTracer
    Vec ExtensiveTracer;
    ierrPetsc = VecCreateMPI(PETSC_COMM_WORLD, nParticlesOld, PETSC_DETERMINE, &ExtensiveTracer);

    PetscInt iRowStartExtensiveTracer, iRowEndExtensiveTracer;
    ierrPetsc = VecGetOwnershipRange(ExtensiveTracer, &iRowStartExtensiveTracer, &iRowEndExtensiveTracer);

    for (PetscInt iParticleOld = 0 ; iParticleOld < nParticlesOld ; iParticleOld++) {
      PetscInt jOld  = iParticleOld * sizePerParticle + iSizePerParticle;
      PetscInt j = iRowStartExtensiveTracer + iParticleOld;
      ierrPetsc = VecSetValue(ExtensiveTracer, j, (PetscScalar) extensiveTracer[jOld], INSERT_VALUES);
    } // iParticle

    ierrPetsc = VecAssemblyBegin(ExtensiveTracer);
    ierrPetsc = VecAssemblyEnd(ExtensiveTracer);

    // ExtensiveTracerRemap
    Vec ExtensiveTracerRemap;
    ierrPetsc = VecCreateMPI(PETSC_COMM_WORLD, nParticlesInit, PETSC_DETERMINE, &ExtensiveTracerRemap);
    ierrPetsc = VecSet(ExtensiveTracerRemap, 0.0);
    ierrPetsc = VecAssemblyBegin(ExtensiveTracerRemap);
    ierrPetsc = VecAssemblyEnd(ExtensiveTracerRemap);

    PetscInt iRowStartExtensiveTracerRemap, iRowEndExtensiveTracerRemap;
    ierrPetsc = VecGetOwnershipRange(ExtensiveTracerRemap, &iRowStartExtensiveTracerRemap, &iRowEndExtensiveTracerRemap);

    // remap
    ierrPetsc = MatMult(Weights, ExtensiveTracer, ExtensiveTracerRemap);

    // ExtensiveTracerRemap
    for (PetscInt iParticleInit = 0 ; iParticleInit < nParticlesInit ; iParticleInit++) {
      PetscInt jInit = iParticleInit * sizePerParticle + iSizePerParticle;
      PetscScalar value;
      PetscInt j = iRowStartExtensiveTracerRemap + iParticleInit;
      ierrPetsc = VecGetValues(ExtensiveTracerRemap, 1, &j, &value);
      extensiveTracerRemap[jInit] = (double) value;
    } // iParticle

    // destroy vectors
    ierrPetsc = VecDestroy(&ExtensiveTracer);
    ierrPetsc = VecDestroy(&ExtensiveTracerRemap);

  } // iSizePerParticle

} // Remapping::remap_variable

// determine the unique set of init indices that form the new set of elements
void Remapping::new_element_indices(const int nCategories, const std::vector<double> effectiveElementArea, const std::vector<double> iceConcentration, const std::vector<int> newElements, std::vector<int> &iParticleInitFromNew, int &nParticlesNew, std::vector<std::pair<int,int>> &bonds) {

  double iceAreaLimit = 1.0;

  bool* particleHasIce = new bool[tessellation->nParticlesAllInit];
  for (int iParticleInit = 0 ; iParticleInit < tessellation->nParticlesAllInit ; iParticleInit++) {
    if (effectiveElementArea[iParticleInit] > iceAreaLimit) {
      particleHasIce[iParticleInit] = true;
    } else {
      particleHasIce[iParticleInit] = false;
    }
  } // iParticleInit

  std::set<int> iParticleInitUnique;
  for (int iParticleInit = 0 ; iParticleInit < tessellation->nParticlesAllInit ; iParticleInit++) {
    if (particleHasIce[iParticleInit]) {
      iParticleInitUnique.insert(iParticleInit);
    }
  } // iParticleInit

  // get the unique halo particle indices that are active going forward
  /*std::vector<int> globalIndexInitHalo;
  std::copy(globalIndexInitHaloUnique.begin(), globalIndexInitHaloUnique.end(), std::back_inserter(globalIndexInitHalo));

  // send halo indices to other processors
  DEMSI::ShareLists* haloLists = new DEMSI::ShareLists(&globalIndexInitHalo,NULL,partition,log);
  while (haloLists->iterate()) {

    std::vector<int>::iterator globalIndexInitHaloItr;
    for (int iParticleInit = 0 ; iParticleInit < tessellation->nParticlesInit ; iParticleInit++) {

      globalIndexInitHaloItr = std::find(globalIndexInitHalo.begin(), globalIndexInitHalo.end(), tessellation->globalIndexInit[iParticleInit]);
      if (globalIndexInitHaloItr != globalIndexInitHalo.end()) {

	iParticleInitUnique.insert(iParticleInit);
	globalIndexInitHalo.erase(globalIndexInitHaloItr);

      }

    }

  }
  delete haloLists;*/

  // add new elements
  std::copy(newElements.begin(), newElements.end(), std::inserter(iParticleInitUnique, iParticleInitUnique.end()));

  // get the unique init particle indices that are active going forward
  std::copy(iParticleInitUnique.begin(), iParticleInitUnique.end(), std::back_inserter(iParticleInitFromNew));

  // sort the indices
  std::sort(iParticleInitFromNew.begin(), iParticleInitFromNew.end());

  nParticlesNew = iParticleInitFromNew.size();

  // iParticleNewFromInit
  std::map<int,int> iParticleNewFromInit;
  for (int iParticleNew = 0 ; iParticleNew < nParticlesNew ; iParticleNew++) {
    int iParticleInit = iParticleInitFromNew[iParticleNew];
    iParticleNewFromInit[iParticleInit] = iParticleNew;
  }

  // bonds
  for (int iParticleInit = 0 ; iParticleInit < tessellation->nParticlesAllInit ; iParticleInit++) {
    for (int iCellOnCell = 0 ; iCellOnCell < tessellation->nVerticesOnCellInit[iParticleInit] ; iCellOnCell++) {
      int iParticleInitNeighbour = tessellation->cellsOnCellLocal[iParticleInit][iCellOnCell];
      if (iParticleInitNeighbour > iParticleInit and particleHasIce[iParticleInit] and particleHasIce[iParticleInitNeighbour]) {
	int iParticleNew =          iParticleNewFromInit[iParticleInit];
	int iParticleNewNeighbour = iParticleNewFromInit[iParticleInitNeighbour];
	bonds.push_back(std::pair<int,int>(iParticleNew+1,iParticleNewNeighbour+1));
      }
    } // iCellOnCell
  } // iParticleInit

  delete [] particleHasIce;

} // Remapping::new_element_indices

// write a lammps particle input file for the new remaped particle distribution
void Remapping::write_lammps_file(int nParticlesNew, double *xNew, double *yNew, double *rNew, int *typeNew, double maxRadius) {

  char mpiErrBuffer[MPI_MAX_ERROR_STRING];
  int mpiErrLen;
  int err;

  int nParticlesNewAll[partition->nprocs()];
  err = MPI_Gather(&nParticlesNew, 1, MPI_INT, &nParticlesNewAll, 1, MPI_INT, partition->master_proc(), partition->comm());
  MPI_Error_string(err, mpiErrBuffer, &mpiErrLen);
  log->check(err == MPI_SUCCESS, "Problem getting total number of active particles for each proc: ", (std::string) mpiErrBuffer);

  int nParticlesNewTotal = 0;
  for (int iProc = 0 ; iProc < partition->nprocs() ; iProc++) {
    nParticlesNewTotal += nParticlesNewAll[iProc];
  } // iProc

  double maxRadiusAll;
  err = MPI_Reduce(&maxRadius, &maxRadiusAll, 1, MPI_DOUBLE, MPI_MAX, partition->master_proc(), partition->comm());
  MPI_Error_string(err, mpiErrBuffer, &mpiErrLen);
  log->check(err == MPI_SUCCESS, "Problem getting max radius of active particles: ", (std::string) mpiErrBuffer);

  std::ofstream lammpsInputFile;

  if (partition->on_master_proc()) {

    // open lammps input file
    std::string filenameOutLammps = "lammps_particle_input_remap.dat";
    lammpsInputFile.open(filenameOutLammps);
    lammpsInputFile << std::setprecision(12);

    // lammps file header
    lammpsInputFile << "# particle input file" << std::endl << std::endl;
    lammpsInputFile << nParticlesNewTotal << " atoms" << std::endl;
    lammpsInputFile << std::endl;
    lammpsInputFile << 2 << " atom types" << std::endl;
    lammpsInputFile << std::endl;
    lammpsInputFile << 0.0 << " " << grid->lx() << " xlo xhi" << std::endl;
    lammpsInputFile << 0.0 << " " << grid->ly() << " ylo yhi" << std::endl;
    lammpsInputFile << -maxRadiusAll << " " << maxRadiusAll << " zlo zhi" << std::endl << std::endl;

    // write by particle data
    lammpsInputFile << "Atoms # sphere" << std::endl << std::endl;

  } // on master proc

  int displacements[partition->nprocs()];
  displacements[0] = 0;
  for (int iProc = 1 ; iProc < partition->nprocs() ; iProc++) {
    displacements[iProc] = displacements[iProc-1] + nParticlesNewAll[iProc-1];
  } // iProc

  double *xAll, *yAll, *rAll;
  int *typeAll;

  if (partition->on_master_proc()) {
    xAll = new double[nParticlesNewTotal];
    yAll = new double[nParticlesNewTotal];
    rAll = new double[nParticlesNewTotal];
    typeAll = new int[nParticlesNewTotal];
  } // on master proc

  err = MPI_Gatherv(&xNew[0], nParticlesNew, MPI_DOUBLE, &xAll[0], nParticlesNewAll, displacements, MPI_DOUBLE, partition->master_proc(), partition->comm());
  MPI_Error_string(err, mpiErrBuffer, &mpiErrLen);
  log->check(err == MPI_SUCCESS, "Problem gatherving xNew: ", (std::string) mpiErrBuffer);

  err = MPI_Gatherv(&yNew[0], nParticlesNew, MPI_DOUBLE, &yAll[0], nParticlesNewAll, displacements, MPI_DOUBLE, partition->master_proc(), partition->comm());
  MPI_Error_string(err, mpiErrBuffer, &mpiErrLen);
  log->check(err == MPI_SUCCESS, "Problem gatherving xNew: ", (std::string) mpiErrBuffer);

  err = MPI_Gatherv(&rNew[0], nParticlesNew, MPI_DOUBLE, &rAll[0], nParticlesNewAll, displacements, MPI_DOUBLE, partition->master_proc(), partition->comm());
  MPI_Error_string(err, mpiErrBuffer, &mpiErrLen);
  log->check(err == MPI_SUCCESS, "Problem gatherving xNew: ", (std::string) mpiErrBuffer);

  err = MPI_Gatherv(&typeNew[0], nParticlesNew, MPI_INT, &typeAll[0], nParticlesNewAll, displacements, MPI_INT, partition->master_proc(), partition->comm());
  MPI_Error_string(err, mpiErrBuffer, &mpiErrLen);
  log->check(err == MPI_SUCCESS, "Problem gatherving typeNew: ", (std::string) mpiErrBuffer);

  if (partition->on_master_proc()) {

    for (int iParticleNew = 0 ; iParticleNew < nParticlesNewTotal ; iParticleNew++) {

      double iceThickness = 1.0;
      double iceFraction = 1.0;
      double density = iceThickness * iceFraction * 900;

      lammpsInputFile << iParticleNew+1 << " "; // atom id, start at 1
      lammpsInputFile << std::max(typeAll[iParticleNew],1) << " "; // atom type
      lammpsInputFile << rAll[iParticleNew]*2.0 << " "; // diameter
      lammpsInputFile << density << " "; // density placeholder
      lammpsInputFile << xAll[iParticleNew] << " "; // x position
      lammpsInputFile << yAll[iParticleNew] << " "; // y position
      lammpsInputFile << 0.0 << std::endl; // z position

    } // iParticleNew

    lammpsInputFile << std::endl << std::endl;
    lammpsInputFile.close();

    delete [] xAll;
    delete [] yAll;
    delete [] rAll;
    delete [] typeAll;

  } // on master proc

}

// performing the remapping
void Remapping::remap(bool const remapNow, bool const init) {

  // check if the alarm is ringing for remapping
  if (useRemapping and (remappingAlarm.is_ringing() or remapNow)) {

    Kokkos::Profiling::pushRegion("Remapping");

    int nCategories = column->columnDimensions->size("nCategories");

    //--------------------------------------------------------------------------
    // Remapping objects
    //--------------------------------------------------------------------------

    Kokkos::Profiling::pushRegion("Remapping variables init");
    (*log)(DEMSI::LOG::DEBUG) << "   ...Remapping variables init..." << std::endl;

    std::map<DEMSI::ColumnVariable<double>*,DEMSI::RemappingTracer*> tracerMap;

    // effective element area
    DEMSI::RemappingTracer* remappingTracer = new DEMSI::RemappingTracer(tessellation, column->effectiveElementArea);
    tracerMap.insert(std::pair<DEMSI::ColumnVariable<double>*,DEMSI::RemappingTracer*>(column->effectiveElementArea,remappingTracer));

    // tracers
    for (int iPtr = 0 ; iPtr < column->columnTracerPtrs.size() ; iPtr++) {
      DEMSI::RemappingTracer* remappingTracer = new DEMSI::RemappingTracer(tessellation, column->columnTracerPtrs[iPtr]);
      tracerMap.insert(std::pair<DEMSI::ColumnVariable<double>*,DEMSI::RemappingTracer*>(column->columnTracerPtrs[iPtr],remappingTracer));
    }

    // other variables
    std::vector<DEMSI::ColumnVariable<double>*> fluxVariables;
    if (column->useColumn) {
      fluxVariables.push_back(column->oceanHeatFlux);
      fluxVariables.push_back(column->oceanShortwaveFlux);
      fluxVariables.push_back(column->pondSnowDepthDifference); // required why?
      fluxVariables.push_back(column->pondLidMeltFluxFraction); // required why?
    } // useColumn

    std::map<DEMSI::ColumnVariable<double>*,DEMSI::RemappingTracer*> fluxesMap;

    for (int iFlux = 0 ; iFlux < fluxVariables.size() ; iFlux++) {
      DEMSI::RemappingTracer* remappingFlux = new DEMSI::RemappingTracer(tessellation, fluxVariables[iFlux]);
      fluxesMap.insert(std::pair<DEMSI::ColumnVariable<double>*,DEMSI::RemappingTracer*>(fluxVariables[iFlux],remappingFlux));
    } // iFlux

    Kokkos::Profiling::popRegion();

    //--------------------------------------------------------------------------
    // Remapping weights
    //--------------------------------------------------------------------------

    Kokkos::Profiling::pushRegion("Remapping weights");
    (*log)(DEMSI::LOG::DEBUG) << "   ......" << std::endl;

    // get geometric remapping weights
    Mat Weights;
    std::vector<double> weights; // weights vector of length nWeights
    std::vector<int> iParticlesOldWeights; // indices of old particles, length nWeights
    std::vector<int> iParticlesInitWeights; // indices of new particles, length nWeights

    remapping_weights(Weights, weights, iParticlesOldWeights, iParticlesInitWeights);

    Kokkos::Profiling::popRegion();

    //--------------------------------------------------------------------------
    // Particles in use
    //--------------------------------------------------------------------------

    std::vector<int> useParticle;

    if (remappingCorrectionType == "obr") {
      particles_in_use(Weights, column->effectiveElementArea, useParticle);
    }

    //--------------------------------------------------------------------------
    // Preparation before remapping
    //--------------------------------------------------------------------------

    Kokkos::Profiling::pushRegion("Remapping pre-remap");
    (*log)(DEMSI::LOG::DEBUG) << "   ......" << std::endl;

    // firstly convert volumes to thicknesses
    column->to_thickness();

    // effective element area
    DEMSI::RemappingTracer* effectiveAreaRemap = tracerMap[column->effectiveElementArea];

    for (int iParticle = 0 ; iParticle < column->effectiveElementArea->num_particles() ; iParticle++) {
      effectiveAreaRemap->extensiveTracer[iParticle] = (*column->effectiveElementArea)(iParticle);
    } // iParticle

    // iterate over tracers parents before children to get extensive tracer values
    DEMSI::TracerTree* tracerTree = column->tracerTree;
    while (tracerTree != NULL) {

      DEMSI::ColumnVariable<double>* tracer = tracerTree->tracer;
      DEMSI::ColumnVariable<double>* parent = tracerTree->parent;
      DEMSI::RemappingTracer* tracerRemap = tracerMap[tracer];
      DEMSI::RemappingTracer* parentRemap = tracerMap[parent];

      (*log)(DEMSI::LOG::DEBUG) << "      ...Calc. extensive tracer " << tracer->name() << std::endl;

      if (tracer->name() == "iceAreaCategory") {

	// extensive tracer and total tracer before
	for (int iParticle = 0 ; iParticle < tracer->num_particles() ; iParticle++) {
	  for (int iChild = 0 ; iChild < tracer->size_per_particle() ; iChild++) {
	    int ijc = iParticle * tracer->size_per_particle() + iChild;
	    int iTracer = iChild;
	    tracerRemap->extensiveTracer[ijc] = (*tracer)(ijc) * (*column->effectiveElementArea)(iParticle);
	  } // iChild
	} // iParticle

      } else {

	int sizeChildPerParent = tracer->size_per_particle() / parent->size_per_particle();

	// other tracers
	for (int iParticle = 0 ; iParticle < tracer->num_particles() ; iParticle++) {
	  for (int iParent = 0 ; iParent < parent->size_per_particle() ; iParent++) {
	    int ijp = iParticle * parent->size_per_particle() + iParent;
	    for (int iChild = 0 ; iChild < sizeChildPerParent ; iChild++) {
	      int ijp = iParticle * parent->size_per_particle() + iParent;
	      int ijc = iParticle * tracer->size_per_particle() + iParent * sizeChildPerParent + iChild;
	      int iTracer = iParent * sizeChildPerParent + iChild;
	      tracerRemap->extensiveTracer[ijc] = (*tracer)(ijc) * parentRemap->extensiveTracer[ijp];
	    } // iChild
	  } // iParent
	} // iParticle

      }

      tracerTree = tracerTree->next;
    }

    // iterate over other variables
    for (int iFlux = 0 ; iFlux < fluxVariables.size() ; iFlux++) {

      DEMSI::ColumnVariable<double>* fluxVariable = fluxVariables[iFlux];
      DEMSI::RemappingTracer* fluxRemap = fluxesMap[fluxVariable];

      (*log)(DEMSI::LOG::DEBUG) << "      ...Calc. extensive tracer " << fluxVariable->name() << std::endl;

      for (int iParticle = 0 ; iParticle < fluxVariable->num_particles() ; iParticle++) {
	for (int iChild = 0 ; iChild < fluxVariable->size_per_particle() ; iChild++) {
	  int ijc = iParticle * fluxVariable->size_per_particle() + iChild;
	  fluxRemap->extensiveTracer[ijc] = (*fluxVariable)(ijc) * (*column->effectiveElementArea)(iParticle);
	} // iChild
      } // iParticle

    } // iFlux

    Kokkos::Profiling::popRegion();

    //--------------------------------------------------------------------------
    // Remapping correction fluxes
    //--------------------------------------------------------------------------

    if (remappingCorrectionType == "flux") {

      Kokkos::Profiling::pushRegion("Remapping correction fluxes");
      (*log)(DEMSI::LOG::DEBUG) << "   ......" << std::endl;

      weights_correction(Weights, effectiveAreaRemap->extensiveTracer);

      Kokkos::Profiling::popRegion();

    } // useFluxCorrection

    //--------------------------------------------------------------------------
    // Conservation before
    //--------------------------------------------------------------------------

    std::vector<double> totalTracerBefore;

    if (checkConservation) {

      (*log)(DEMSI::LOG::DEBUG) << "   ...Remapping conservation before" << std::endl;

      // effectiveElementArea
      double total = 0.0;
      for (int i = 0 ; i < effectiveAreaRemap->extensiveTracer.size() ; i++) {
	total += effectiveAreaRemap->extensiveTracer[i];
      } // i
      totalTracerBefore.push_back(total);

      // tracers
      tracerTree = column->tracerTree;
      while (tracerTree != NULL) {

	DEMSI::ColumnVariable<double>* tracer = tracerTree->tracer;
	DEMSI::RemappingTracer* tracerRemap = tracerMap[tracer];

	double total = 0.0;
	for (int i = 0 ; i < tracerRemap->extensiveTracer.size() ; i++) {
	  total += tracerRemap->extensiveTracer[i];
	} // i
	totalTracerBefore.push_back(total);

	tracerTree = tracerTree->next;
      }

    } // checkConservation

    //--------------------------------------------------------------------------
    // Perform geometric remapping
    //--------------------------------------------------------------------------

    Kokkos::Profiling::pushRegion("Geometric remapping");

    // remap the effective element area
    (*log)(DEMSI::LOG::DEBUG) << "      ...Geometric remap for effectiveElementArea" << std::endl;

    remap_variable(column->effectiveElementArea->size_per_particle(), effectiveAreaRemap->extensiveTracerRemap, effectiveAreaRemap->extensiveTracer, Weights);

    // iterate over tracers parents before children
    tracerTree = column->tracerTree;
    while (tracerTree != NULL) {

      DEMSI::ColumnVariable<double>* tracer = tracerTree->tracer;
      DEMSI::RemappingTracer* tracerRemap = tracerMap[tracer];

      (*log)(DEMSI::LOG::DEBUG) << "      ...Geometric remap for " << tracer->name() << std::endl;

      remap_variable(tracer->size_per_particle(), tracerRemap->extensiveTracerRemap, tracerRemap->extensiveTracer, Weights);

      tracerTree = tracerTree->next;
    }

    // iterate over other variables
    for (int iFlux = 0 ; iFlux < fluxVariables.size() ; iFlux++) {

      DEMSI::ColumnVariable<double>* fluxVariable = fluxVariables[iFlux];
      DEMSI::RemappingTracer* fluxRemap = fluxesMap[fluxVariable];

      (*log)(DEMSI::LOG::DEBUG) << "      ...Geometric remap for " << fluxVariable->name() << std::endl;

      remap_variable(fluxVariable->size_per_particle(), fluxRemap->extensiveTracerRemap, fluxRemap->extensiveTracer, Weights);

    } // iFlux

    // clean up Weights matrix
    int ierrPetsc = MatDestroy(&Weights);
    log->check(ierrPetsc == 0, "Petsc: Problem with MatDestroy for Weights: ", std::to_string(ierrPetsc));

    Kokkos::Profiling::popRegion();

    //--------------------------------------------------------------------------
    // OBR correct tracers
    //--------------------------------------------------------------------------

    if (remappingCorrectionType == "obr") {

      Kokkos::Profiling::pushRegion("Remapping OBR");

      // compute value of intensive tracer - used for target in optimization
      tracerTree = column->tracerTree;
      while (tracerTree != NULL) {

	DEMSI::ColumnVariable<double>* tracer = tracerTree->tracer;
	DEMSI::ColumnVariable<double>* parent = tracerTree->parent;
	DEMSI::RemappingTracer* tracerRemap = tracerMap[tracer];
	DEMSI::RemappingTracer* parentRemap = tracerMap[parent];

	(*log)(DEMSI::LOG::DEBUG) << "      ...Convert to tracer (not extensive) " << tracer->name() << std::endl;

	if (tracer->name() == "iceAreaCategory") {

	  for (int iParticleInit = 0 ; iParticleInit < tessellation->nParticlesAllInit ; iParticleInit++) {
	    for (int iChild = 0 ; iChild < tracer->size_per_particle() ; iChild++) {
	      int ijc    = iParticleInit * tracer->size_per_particle() + iChild;
	      if (effectiveAreaRemap->extensiveTracerRemap[iParticleInit] > 0) {
		tracerRemap->intensiveTracerRemap[ijc] = tracerRemap->extensiveTracerRemap[ijc] / effectiveAreaRemap->extensiveTracerRemap[iParticleInit];
	      } else {
		tracerRemap->intensiveTracerRemap[ijc] = 0.0;
	      }
	    } // iChild
	  } // iParticle

	} else {

	  int sizeChildPerParent = tracer->size_per_particle() / parent->size_per_particle();

	  // other tracers
	  for (int iParticleInit = 0 ; iParticleInit < tessellation->nParticlesAllInit ; iParticleInit++) {
	    for (int iParent = 0 ; iParent < parent->size_per_particle() ; iParent++) {
	      int ijp = iParticleInit * parent->size_per_particle() + iParent;
	      for (int iChild = 0 ; iChild < sizeChildPerParent ; iChild++) {
		int ijc    = iParticleInit * tracer->size_per_particle() + iParent * sizeChildPerParent + iChild;
		if (parentRemap->extensiveTracerRemap[ijp] != 0) {
		  tracerRemap->intensiveTracerRemap[ijc] = tracerRemap->extensiveTracerRemap[ijc] / parentRemap->extensiveTracerRemap[ijp];
		} else {
		  tracerRemap->intensiveTracerRemap[ijc] = 0.0;
		}
	      } // iChild
	    } // iParent
	  } // iParticle
	}

	tracerTree = tracerTree->next;
      }

      for (int iFlux = 0 ; iFlux < fluxVariables.size() ; iFlux++) {

	DEMSI::ColumnVariable<double>* fluxVariable = fluxVariables[iFlux];
	DEMSI::RemappingTracer* fluxRemap = fluxesMap[fluxVariable];

	for (int iParticleInit = 0 ; iParticleInit < tessellation->nParticlesAllInit ; iParticleInit++) {
	  for (int iChild = 0 ; iChild < fluxVariable->size_per_particle() ; iChild++) {
	    int ijc    = iParticleInit * fluxVariable->size_per_particle() + iChild;
	    if (effectiveAreaRemap->extensiveTracerRemap[iParticleInit] != 0) {
	      fluxRemap->intensiveTracerRemap[ijc] = fluxRemap->extensiveTracerRemap[ijc] / effectiveAreaRemap->extensiveTracerRemap[iParticleInit];
	    } else {
	      fluxRemap->intensiveTracerRemap[ijc] = 0.0;
	    }
	  }
	}
      }

      // correct area testing
      //correct_area_obr(column->effectiveElementArea, effectiveAreaRemap->extensiveTracerRemap, effectiveAreaRemap->extensiveTracer, useParticle);

     // perform optimization step to enforce conservation and bounds
     tracerTree = column->tracerTree;
     while (tracerTree != NULL) {

       DEMSI::ColumnVariable<double>* tracer = tracerTree->tracer;
       DEMSI::ColumnVariable<double>* parent = tracerTree->parent;
       DEMSI::RemappingTracer* tracerRemap = tracerMap[tracer];
       DEMSI::RemappingTracer* parentRemap = tracerMap[parent];

       (*log)(DEMSI::LOG::DEBUG) << "      ...OBR Correction for " << tracer->name() << std::endl;

       if (tracer->name() == "iceAreaCategory") {
	 correct_tracer_obr(tracer, tracerRemap->intensiveTracerRemap, tracerRemap->extensiveTracerRemap, tracerRemap->extensiveTracer, column->effectiveElementArea, effectiveAreaRemap->extensiveTracerRemap, effectiveAreaRemap->extensiveTracer, useParticle, weights, iParticlesOldWeights, iParticlesInitWeights);
       } else {
	 correct_tracer_obr(tracer, tracerRemap->intensiveTracerRemap, tracerRemap->extensiveTracerRemap, tracerRemap->extensiveTracer, parent, parentRemap->extensiveTracerRemap, parentRemap->extensiveTracer, useParticle, weights, iParticlesOldWeights, iParticlesInitWeights);
       }

       tracerTree = tracerTree->next;
     }

     // iterate over other variables
     for (int iFlux = 0 ; iFlux < fluxVariables.size() ; iFlux++) {

	DEMSI::ColumnVariable<double>* fluxVariable = fluxVariables[iFlux];
	DEMSI::RemappingTracer* fluxRemap = fluxesMap[fluxVariable];

	(*log)(DEMSI::LOG::DEBUG) << "      ...OBR correction for " << fluxVariable->name() << std::endl;

	correct_tracer_obr(fluxVariable, fluxRemap->intensiveTracerRemap, fluxRemap->extensiveTracerRemap, fluxRemap->extensiveTracer, column->effectiveElementArea, effectiveAreaRemap->extensiveTracerRemap, effectiveAreaRemap->extensiveTracer, useParticle, weights, iParticlesOldWeights, iParticlesInitWeights);

     } // iFlux

     Kokkos::Profiling::popRegion();

   } // remappingCorrectionType

    //--------------------------------------------------------------------------
    // Conservation after
    //--------------------------------------------------------------------------

    std::vector<double> totalTracerAfter;
    std::vector<std::string> variableNames;

    if (checkConservation) {

      (*log)(DEMSI::LOG::DEBUG) << "   ...Remapping conservation after" << std::endl;

      // effectiveElementArea
      double total = 0.0;
      for (int i = 0 ; i < effectiveAreaRemap->extensiveTracerRemap.size() ; i++) {
	total += effectiveAreaRemap->extensiveTracerRemap[i];
      } // i
      totalTracerAfter.push_back(total);
      variableNames.push_back("effectiveElementArea");

      // tracers
      tracerTree = column->tracerTree;
      while (tracerTree != NULL) {

	DEMSI::ColumnVariable<double>* tracer = tracerTree->tracer;
	DEMSI::RemappingTracer* tracerRemap = tracerMap[tracer];

	double total = 0.0;
	for (int i = 0 ; i < tracerRemap->extensiveTracerRemap.size() ; i++) {
	  total += tracerRemap->extensiveTracerRemap[i];
	} // i
	totalTracerAfter.push_back(total);
	variableNames.push_back(tracer->name());

	tracerTree = tracerTree->next;
      }

      // check tracer conservation
      for (int iCheck = 0 ; iCheck < variableNames.size() ; iCheck++) {
	double totalBefore = totalTracerBefore[iCheck];
	double totalAfter  = totalTracerAfter [iCheck];
	double relativeError = std::abs(totalAfter - totalBefore) / totalBefore;
	if (relativeError > 1e-11) {
	  (*log)(DEMSI::LOG::CRITICAL) << "Conservation error in remapping for " << variableNames[iCheck] << ", rel error: " << relativeError << std::endl;
	}
      } // iCheck

    } // checkConservation

    //--------------------------------------------------------------------------
    // Convert from extensive to intensive
    //--------------------------------------------------------------------------

    // check effective element area is within geometric bounds
    if (checkEffectiveArea) {

      (*log)(DEMSI::LOG::DEBUG) << "   ...Check effective area limits..." << std::endl;

      for (int iParticleInit = 0 ; iParticleInit < tessellation->nParticlesAllInit ; iParticleInit++) {
	double effectiveAreaExcess = effectiveAreaRemap->extensiveTracerRemap[iParticleInit] - tessellation->initPolygons[iParticleInit].area();
	log->check(effectiveAreaExcess < 1.0, "Excess effective area in remapping: ", std::to_string(effectiveAreaRemap->extensiveTracerRemap[iParticleInit]), " ", std::to_string(tessellation->initPolygons[iParticleInit].area())," for particle: ", std::to_string(iParticleInit));
      } // iParticle

    } // checkEffectiveAreaLimit

    // tracers
    tracerTree = column->tracerTree;
    while (tracerTree != NULL) {

      DEMSI::ColumnVariable<double>* tracer = tracerTree->tracer;
      DEMSI::ColumnVariable<double>* parent = tracerTree->parent;
      DEMSI::RemappingTracer* tracerRemap = tracerMap[tracer];
      DEMSI::RemappingTracer* parentRemap = tracerMap[parent];

      (*log)(DEMSI::LOG::DEBUG) << "      ...Calc. intensive tracer for " << tracer->name() << std::endl;

      if (tracer->name() == "iceAreaCategory") {

	// extensive tracer and total tracer before
	for (int iParticleInit = 0 ; iParticleInit < tessellation->nParticlesAllInit ; iParticleInit++) {
	  for (int iChild = 0 ; iChild < tracer->size_per_particle() ; iChild++) {
	    int ijc = iParticleInit * tracer->size_per_particle() + iChild;
	    double iceArea = tracerRemap->extensiveTracerRemap[ijc] / tessellation->initPolygons[iParticleInit].area();
	    if (iceArea > 1e-6) {
	      tracerRemap->intensiveTracerRemap[ijc] = iceArea;
	    } else {
	      tracerRemap->intensiveTracerRemap[ijc] = 0.0;
	    }
	  } // iChild
	} // iParticle

      } else {

	int sizeChildPerParent = tracer->size_per_particle() / parent->size_per_particle();

	// other tracers
	for (int iParticleInit = 0 ; iParticleInit < tessellation->nParticlesAllInit ; iParticleInit++) {
	  for (int iParent = 0 ; iParent < parent->size_per_particle() ; iParent++) {
	    int ijp = iParticleInit * parent->size_per_particle() + iParent;
	    for (int iChild = 0 ; iChild < sizeChildPerParent ; iChild++) {
	      int ijc = iParticleInit * tracer->size_per_particle() + iParent * sizeChildPerParent + iChild;
	      if (parentRemap->extensiveTracerRemap[ijp] != 0.0) {
		tracerRemap->intensiveTracerRemap[ijc] = tracerRemap->extensiveTracerRemap[ijc] / parentRemap->extensiveTracerRemap[ijp];
	      } else {
		tracerRemap->intensiveTracerRemap[ijc] = 0.0;
	      }
	    } // iChild
	  } // iParent
	} // iParticle

      }

      tracerTree = tracerTree->next;
    }

    // iterate over other variables
    for (int iFlux = 0 ; iFlux < fluxVariables.size() ; iFlux++) {

      DEMSI::ColumnVariable<double>* fluxVariable = fluxVariables[iFlux];
      DEMSI::RemappingTracer* fluxRemap = fluxesMap[fluxVariable];

      for (int iParticleInit = 0 ; iParticleInit < tessellation->nParticlesAllInit ; iParticleInit++) {
	for (int iChild = 0 ; iChild < fluxVariable->size_per_particle() ; iChild++) {
	  int ijc = iParticleInit * fluxVariable->size_per_particle() + iChild;
	  fluxRemap->intensiveTracerRemap[ijc] = fluxRemap->extensiveTracerRemap[ijc] / tessellation->initPolygons[iParticleInit].area();
	} // iChild
      } // iParticle

    } // iFlux

    //--------------------------------------------------------------------------
    // Remove thin ice categories
    //--------------------------------------------------------------------------

    DEMSI::RemappingTracer* iceConcRemap = tracerMap[column->iceAreaCategory];
    for (int i = 0 ; i < iceConcRemap->intensiveTracerRemap.size() ; i++) {
      if (iceConcRemap->intensiveTracerRemap[i] < 1.0e-11) {
	iceConcRemap->intensiveTracerRemap[i] = 0.0;
      }
    } // i

    //--------------------------------------------------------------------------
    // New elements
    //--------------------------------------------------------------------------

    Kokkos::Profiling::pushRegion("New elements");
    (*log)(DEMSI::LOG::DEBUG) << "   ...New elements..." << std::endl;

    // get new elements from ocean
    std::vector<int> newElements = ocean->get_new_ice_elements();

    // get new element indices
    std::vector<int> iParticleInitFromNew; // init index for new particles, length nParticlesNew
    int nParticlesNew; // total number of new particles generated

    // bonds
    std::vector<std::pair<int,int>> bonds;

    new_element_indices(nCategories, effectiveAreaRemap->extensiveTracerRemap, iceConcRemap->intensiveTracerRemap, newElements, iParticleInitFromNew, nParticlesNew, bonds);

    // reinitialize the processor transfer
    column->reset_processor_transfer_remap(nParticlesNew);

    Kokkos::Profiling::popRegion();

    //--------------------------------------------------------------------------
    // Resize DEMSI variables
    //--------------------------------------------------------------------------

    Kokkos::Profiling::pushRegion("Resize variables");
    (*log)(DEMSI::LOG::DEBUG) << "   ...Resize variables..." << std::endl;

    // double variables
    for (int iPtr = 0 ; iPtr < column->columnVariables->doubleVariables.size() ; iPtr++) {

      int sizePerParticle = column->columnVariables->doubleVariables[iPtr]->size_per_particle();
      int nSizeNew = nParticlesNew * sizePerParticle;

      column->columnVariables->doubleVariables[iPtr]->resize(nSizeNew);

      column->columnVariables->doubleVariables[iPtr]->set(0.0);

    } // iPtr

    // int variables
    for (int iPtr = 0 ; iPtr < column->columnVariables->intVariables.size() ; iPtr++) {

      int sizePerParticle = column->columnVariables->intVariables[iPtr]->size_per_particle();
      int nSizeNew = nParticlesNew * sizePerParticle;

      column->columnVariables->intVariables[iPtr]->resize(nSizeNew);

      column->columnVariables->intVariables[iPtr]->set(0);

    } // iPtr

    Kokkos::Profiling::popRegion();

    //--------------------------------------------------------------------------
    // Set DEMSI column variables
    //--------------------------------------------------------------------------

    // iterate over column tracers and remap
    Kokkos::Profiling::pushRegion("Set column variables");
    (*log)(DEMSI::LOG::DEBUG) << "   ...Set column variables..." << std::endl;

    // reset effectiveElementArea back to init area
    for (int iParticle = 0 ; iParticle < nParticlesNew ; iParticle++) {
      (*column->effectiveElementArea)(iParticle) = tessellation->initPolygons[iParticleInitFromNew[iParticle]].area();
    } // iParticle

    // tracers
    tracerTree = column->tracerTree;
    while (tracerTree != NULL) {

      DEMSI::ColumnVariable<double>* tracer = tracerTree->tracer;
      DEMSI::RemappingTracer* tracerRemap = tracerMap[tracer];

      (*log)(DEMSI::LOG::DEBUG) << "      ...Set column variable for " << tracer->name() << std::endl;

      // extensive tracer and total tracer before
      for (int iParticle = 0 ; iParticle < nParticlesNew ; iParticle++) {
	int iParticleInit = iParticleInitFromNew[iParticle];
	for (int iChild = 0 ; iChild < tracer->size_per_particle() ; iChild++) {
	  int ijc    = iParticleInit * tracer->size_per_particle() + iChild;
	  int ijcNew = iParticle     * tracer->size_per_particle() + iChild;
	  (*tracer)(ijcNew) = tracerRemap->intensiveTracerRemap[ijc];
	} // iChild
      } // iParticle

      tracerTree = tracerTree->next;
    }

    // iterate over other variables
    for (int iFlux = 0 ; iFlux < fluxVariables.size() ; iFlux++) {

      DEMSI::ColumnVariable<double>* fluxVariable = fluxVariables[iFlux];
      DEMSI::RemappingTracer* fluxRemap = fluxesMap[fluxVariable];

      (*log)(DEMSI::LOG::DEBUG) << "      ... Set column variable for " << fluxVariable->name() << std::endl;

      for (int iParticle = 0 ; iParticle < nParticlesNew ; iParticle++) {
	int iParticleInit = iParticleInitFromNew[iParticle];
	for (int iChild = 0 ; iChild < fluxVariable->size_per_particle() ; iChild++) {
	  int ijc    = iParticleInit * fluxVariable->size_per_particle() + iChild;
	  int ijcNew = iParticle     * fluxVariable->size_per_particle() + iChild;
	  (*fluxVariable)(ijcNew) = fluxRemap->intensiveTracerRemap[ijc];
	} // iChild
      } // iParticle

    } // iFlux

    Kokkos::Profiling::popRegion();

    //--------------------------------------------------------------------------
    // Check ice concentrations bounds
    //--------------------------------------------------------------------------

    if (checkIceConcentrationBounds) {

      (*log)(DEMSI::LOG::DEBUG) << "      ...check ice concentration bounds" << std::endl;

      for (int iParticle = 0 ; iParticle < nParticlesNew ; iParticle++) {

	std::stringstream ssError;
	double iceConcentration = 0.0;
	for (int iCategory = 0 ; iCategory < nCategories ; iCategory++) {

	  iceConcentration += (*column->iceAreaCategory)(iParticle,iCategory);

	  ssError << std::scientific << std::setprecision(17) << (*column->iceAreaCategory)(iParticle,iCategory);
	  log->check((*column->iceAreaCategory)(iParticle,iCategory) >= -maxIceConcBoundError,
		     "Remapping: iceAreaCategory < 0.0 (", ssError.str(),
		     ") for iParticle: ", std::to_string(iParticle), " and iCategory: ", std::to_string(iCategory));
	  ssError.str(std::string());

	  ssError << std::scientific << std::setprecision(17) << (1.0-(*column->iceAreaCategory)(iParticle,iCategory));
	  log->check((*column->iceAreaCategory)(iParticle,iCategory) <= (1.0 + maxIceConcBoundError),
		     "Remapping: iceAreaCategory > 1.0 (", ssError.str(),
		     ") for iParticle: ", std::to_string(iParticle), " and iCategory: ", std::to_string(iCategory));
	  ssError.str(std::string());

	} // iCategory
	ssError << std::scientific << std::setprecision(17) << iceConcentration;
	log->check(iceConcentration >= -maxIceConcBoundError,
		   "Remapping: iceConcentration < 0.0 (", ssError.str(), ") for iParticle: ", std::to_string(iParticle));
	ssError.str(std::string());

	ssError << std::scientific << std::setprecision(17) << 1.0-iceConcentration;
	log->check(iceConcentration <= (1.0 + maxIceConcBoundError),
		   "Remapping: iceConcentration > 1.0 (", ssError.str(), ") for iParticle: ", std::to_string(iParticle));
	ssError.str(std::string());
      } // iParticle

    } // checkIceConcentrationBounds

    //--------------------------------------------------------------------------
    // Update LAMMPS particle distribution
    //--------------------------------------------------------------------------

    Kokkos::Profiling::pushRegion("Remapping lammps");

    // get new element positions
    double* xNew = new double[nParticlesNew];
    double* yNew = new double[nParticlesNew];
    double* rNew = new double[nParticlesNew];
    int* typeNew = new int[nParticlesNew];

    double maxRadius = 0.0;
    for (int iParticleNew = 0 ; iParticleNew < nParticlesNew ; iParticleNew++) {
      xNew[iParticleNew] = tessellation->xInit[iParticleInitFromNew[iParticleNew]];
      yNew[iParticleNew] = tessellation->yInit[iParticleInitFromNew[iParticleNew]];
      rNew[iParticleNew] = tessellation->radiusInit[iParticleInitFromNew[iParticleNew]];
      typeNew[iParticleNew] = tessellation->typeInit[iParticleInitFromNew[iParticleNew]];
      maxRadius = std::max(maxRadius,rNew[iParticleNew]);
    } // iParticleNew

    // now we recreate the lammps distribution
    write_lammps_file(nParticlesNew, xNew, yNew, rNew, typeNew, maxRadius);

    delete [] xNew;
    delete [] yNew;
    delete [] rNew;
    delete [] typeNew;

    // destroy old atoms
    lammpsInstance->one("delete_atoms group all");

    // create new atoms
    lammpsInstance->one("read_data lammps_particle_input_remap.dat add merge");

    // set the groups
    lammpsInstance->one("group mobile type 1");
    lammpsInstance->one("group coastline type 2");

    lammpsInstance->one("run 0 pre yes");

    // reset the demsi pointers
    particles->set_demsi_pointers();

    Kokkos::Profiling::popRegion();

    //--------------------------------------------------------------------------
    // Set polygons for new particles
    //--------------------------------------------------------------------------

    column->nVerticesInit->resize(nParticlesNew);
    column->verticesInit->resize(nParticlesNew);
    column->areaInit->resize(nParticlesNew);

    for (int iParticleNew = 0 ; iParticleNew < nParticlesNew ; iParticleNew++) {

      int iParticleInit = iParticleInitFromNew[iParticleNew];

      (*column->nVerticesInit)(iParticleNew) = tessellation->initPolygons[iParticleInit].nvertices();
      (*column->areaInit)(iParticleNew) = tessellation->initPolygons[iParticleInit].area();
      for (int iVertex = 0 ; iVertex < tessellation->initPolygons[iParticleInit].nvertices() ; iVertex++) {
	(*column->verticesInit)(iParticleNew,iVertex,0) = tessellation->initPolygons[iParticleInit].get_x(iVertex) - tessellation->xInit[iParticleInit];
	(*column->verticesInit)(iParticleNew,iVertex,1) = tessellation->initPolygons[iParticleInit].get_y(iVertex) - tessellation->yInit[iParticleInit];
      }
    } // iParticle

    //--------------------------------------------------------------------------
    // Post remap tasks
    //--------------------------------------------------------------------------

    Kokkos::Profiling::pushRegion("Remapping post-remap");

    // convert thicknesses back to volumes
    column->from_thickness();

    // update element mass
    column->update_mass();

    // reset velocity
    for (int iParticle = 0 ; iParticle < *(particles->nParticles) ; iParticle++) {
      particles->v(iParticle,0) = 0.0;
      particles->v(iParticle,1) = 0.0;
    } // iParticle

    // ocean coupling/update
    ocean->ice_to_ocean(nParticlesNew, iParticleInitFromNew);
    //if (not init) ocean->update_remap(nParticlesNew, iParticleInitFromNew); // comparison with CICE
    ocean->ocean_to_ice(nParticlesNew, iParticleInitFromNew);

    // reinitialize the contact
    if (contacts->get_contact_type() == DEMSI::CONTACTTYPE::HOPKINS) {
      contacts->initialize_bonded_contacts_hopkins();
    }
    column->compute_mean_min_thickness();
    lammpsInstance->comm_ghosts();

    // reset column transfer
    column->processor_transfer();

    // aggregate tracers
    bool useColumn;
    configs->get({"ConfigGroup:useSections","Config:useColumn"}, useColumn);
    if (useColumn) {
      column->aggregate();
    } else {
      int nCategories = column->columnDimensions->size("nCategories");
      column->iceAreaCell->resize(*(particles->nParticles));
      column->iceVolumeCell->resize(*(particles->nParticles));
      column->snowVolumeCell->resize(*(particles->nParticles));
      for (int iParticle = 0 ; iParticle < *(particles->nParticles) ; iParticle++) {
	(*column->iceAreaCell)(iParticle) = 0.0;
	(*column->iceVolumeCell)(iParticle) = 0.0;
	(*column->snowVolumeCell)(iParticle) = 0.0;
	for (int iCategory = 0 ; iCategory < nCategories ; iCategory++) {
	  (*column->iceAreaCell)(iParticle) += (*column->iceAreaCategory)(iParticle,iCategory);
	  (*column->iceVolumeCell)(iParticle) += (*column->iceVolumeCategory)(iParticle,iCategory);
	  (*column->snowVolumeCell)(iParticle) += (*column->snowVolumeCategory)(iParticle,iCategory);
	} // iCategory
      } // iParticle
    } // useColumn

    for (int iParticle = 0 ; iParticle < *(particles->nParticles) ; iParticle++) {
      for (int iCategory = 0 ; iCategory < nCategories ; iCategory++) {
	log->check(not std::isnan((*column->iceAreaCategory)(iParticle,iCategory)),
		   "NaN iceAreaCategory after remapping for iParticle: ", std::to_string(iParticle), " and iCategory: ", std::to_string(iCategory));
	log->check(not std::isnan((*column->iceVolumeCategory)(iParticle,iCategory)),
		   "NaN iceVolumeCategory after remapping for iParticle: ", std::to_string(iParticle), " and iCategory: ", std::to_string(iCategory));
      } // iCategory
    } // iParticle

    for (int iParticle = 0 ; iParticle < *(particles->nParticles) ; iParticle++) {
      for (int iCategory = 0 ; iCategory < nCategories ; iCategory++) {
	if ((*column->iceAreaCategory)(iParticle,iCategory) < 0.0) {
	  (*log)() << "Negative iceAreaCategory (" << std::scientific << std::setprecision(17) << (*column->iceAreaCategory)(iParticle,iCategory) << ") after remapping for iParticle: " << iParticle << " and iCategory: " << iCategory << std::endl;
	  log->abort();
	}
      } // iCategory
    } // iParticle

    // initialize column physics after remapping
    column->init_remap();

    // reset remap alarm
    remappingAlarm.reset();

    Kokkos::Profiling::popRegion();

    Kokkos::Profiling::popRegion();

    //output_areas();

    //exit(0);

  } // alarm ringing

}

} // namespace DEMSI
