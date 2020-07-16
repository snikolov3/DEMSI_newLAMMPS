#include "demsi_initial_tessellation.h"

#include "demsi_communication.h"

#include <iomanip>
#include <netcdf.h>
#include <lammps/neighbor.h>

namespace DEMSI {

// download an init particle distribution variable on the master proc and send to owning processors
void get_init_variable_and_scatter(int *&arrayLocal,
				   const int ncID,
				   const std::string varname,
				   const int nParticlesGlobal,
				   const int sizePerParticle,
				   const int* sendCountsParticle,
				   const int* displacementsParticle,
				   DEMSI::Partition *partition,
				   DEMSI::Log *log) {

  char mpiErrBuffer[MPI_MAX_ERROR_STRING];
  int mpiErrLen;
  int err;

  int *arrayGlobal;

  if (partition->on_master_proc()) {

    // get variable id
    int varid;
    err = nc_inq_varid(ncID, varname.c_str(), &varid);
    log->check(err == NC_NOERR, "Problem getting varid for :", varname, ": ", nc_strerror(err));

    arrayGlobal = new int[nParticlesGlobal*sizePerParticle];

    err = nc_get_var_int(ncID, varid, &arrayGlobal[0]);
    log->check(err == NC_NOERR, "Problem reading data for :", varname, ": ", nc_strerror(err));

  } // on master proc

  int nParticlesLocal = sendCountsParticle[partition->this_proc()];
  arrayLocal = new int[nParticlesLocal*sizePerParticle];

  int* sendCounts = new int[partition->nprocs()];
  int* displacements = new int[partition->nprocs()];
  for (int iProc = 0 ; iProc < partition->nprocs() ; iProc++) {
    sendCounts[iProc] = sendCountsParticle[iProc] * sizePerParticle;
    displacements[iProc] = displacementsParticle[iProc] * sizePerParticle;
  } // iProc

  err = MPI_Scatterv(&arrayGlobal[0], sendCounts, displacements, MPI_INT, &arrayLocal[0], nParticlesLocal*sizePerParticle, MPI_INT, partition->master_proc(), partition->comm());
  MPI_Error_string(err, mpiErrBuffer, &mpiErrLen);
  log->check(err == MPI_SUCCESS, "Problem scatterving array: ", (std::string) mpiErrBuffer, ", for variable: ", varname);

  if (partition->on_master_proc()) delete [] arrayGlobal;

  delete [] sendCounts;
  delete [] displacements;

}

// download an init particle distribution variable on the master proc and send to owning processors
void get_init_variable_and_scatter(double *&arrayLocal,
				   const int ncID,
				   const std::string varname,
				   const int nParticlesGlobal,
				   const int sizePerParticle,
				   const int* sendCountsParticle,
				   const int* displacementsParticle,
				   DEMSI::Partition *partition,
				   DEMSI::Log *log) {

  char mpiErrBuffer[MPI_MAX_ERROR_STRING];
  int mpiErrLen;
  int err;

  double *arrayGlobal;

  if (partition->on_master_proc()) {

    // get variable id
    int varid;
    err = nc_inq_varid(ncID, varname.c_str(), &varid);
    log->check(err == NC_NOERR, "Problem getting varid for :", varname, ": ", nc_strerror(err));

    arrayGlobal = new double[nParticlesGlobal*sizePerParticle];

    err = nc_get_var_double(ncID, varid, &arrayGlobal[0]);
    log->check(err == NC_NOERR, "Problem reading data for :", varname, ": ", nc_strerror(err));

  } // on master proc

  int nParticlesLocal = sendCountsParticle[partition->this_proc()];
  arrayLocal = new double[nParticlesLocal*sizePerParticle];

  int* sendCounts = new int[partition->nprocs()];
  int* displacements = new int[partition->nprocs()];
  for (int iProc = 0 ; iProc < partition->nprocs() ; iProc++) {
    sendCounts[iProc] = sendCountsParticle[iProc] * sizePerParticle;
    displacements[iProc] = displacementsParticle[iProc] * sizePerParticle;
  } // iProc

  err = MPI_Scatterv(&arrayGlobal[0], sendCounts, displacements, MPI_DOUBLE, &arrayLocal[0], nParticlesLocal*sizePerParticle, MPI_DOUBLE, partition->master_proc(), partition->comm());
  MPI_Error_string(err, mpiErrBuffer, &mpiErrLen);
  log->check(err == MPI_SUCCESS, "Problem scatterving array: ", (std::string) mpiErrBuffer, ", for variable: ", varname);

  if (partition->on_master_proc()) delete [] arrayGlobal;

  delete [] sendCounts;
  delete [] displacements;

}

// download an init particle distribution variable on the master proc and send to owning processors
void get_init_variable_and_scatter(double *&xArrayLocal,
				   double *&yArrayLocal,
				   const int ncID,
				   const std::string varname,
				   const int nParticlesGlobal,
				   const int* sendCountsParticle,
				   const int* displacementsParticle,
				   DEMSI::Partition *partition,
				   DEMSI::Log *log) {

  char mpiErrBuffer[MPI_MAX_ERROR_STRING];
  int mpiErrLen;
  int err;

  double *xArrayGlobal;
  double *yArrayGlobal;

  if (partition->on_master_proc()) {

    // get variable id
    int varid;
    err = nc_inq_varid(ncID, varname.c_str(), &varid);
    log->check(err == NC_NOERR, "Problem getting varid for :", varname, ": ", nc_strerror(err));

    xArrayGlobal = new double[nParticlesGlobal];
    yArrayGlobal = new double[nParticlesGlobal];

    size_t xStart[2] = {0, 0};
    size_t xCount[2] = {(size_t) nParticlesGlobal, 1};
    err = nc_get_vara_double(ncID, varid, xStart, xCount, &xArrayGlobal[0]);
    log->check(err == NC_NOERR, "Problem reading data for x:", varname, ": ", nc_strerror(err));

    size_t yStart[2] = {0, 1};
    size_t yCount[2] = {(size_t) nParticlesGlobal, 1};
    err = nc_get_vara_double(ncID, varid, yStart, yCount, &yArrayGlobal[0]);
    log->check(err == NC_NOERR, "Problem reading data for y:", varname, ": ", nc_strerror(err));

  } // on master proc

  int nParticlesLocal = sendCountsParticle[partition->this_proc()];
  xArrayLocal = new double[nParticlesLocal];
  yArrayLocal = new double[nParticlesLocal];

  int* sendCounts = new int[partition->nprocs()];
  int* displacements = new int[partition->nprocs()];
  for (int iProc = 0 ; iProc < partition->nprocs() ; iProc++) {
    sendCounts[iProc] = sendCountsParticle[iProc];
    displacements[iProc] = displacementsParticle[iProc];
  } // iProc

  err = MPI_Scatterv(&xArrayGlobal[0], sendCounts, displacements, MPI_DOUBLE, &xArrayLocal[0], nParticlesLocal, MPI_DOUBLE, partition->master_proc(), partition->comm());
  MPI_Error_string(err, mpiErrBuffer, &mpiErrLen);
  log->check(err == MPI_SUCCESS, "Problem scatterving x array: ", (std::string) mpiErrBuffer, ", for variable: ", varname);

  err = MPI_Scatterv(&yArrayGlobal[0], sendCounts, displacements, MPI_DOUBLE, &yArrayLocal[0], nParticlesLocal, MPI_DOUBLE, partition->master_proc(), partition->comm());
  MPI_Error_string(err, mpiErrBuffer, &mpiErrLen);
  log->check(err == MPI_SUCCESS, "Problem scatterving y array: ", (std::string) mpiErrBuffer, ", for variable: ", varname);

  if (partition->on_master_proc()) {
    delete [] xArrayGlobal;
    delete [] yArrayGlobal;
  }

  delete [] sendCounts;
  delete [] displacements;

}

// class constructor
Tessellation::Tessellation(DEMSI::Partition* partitionIn, DEMSI::Log* logIn, DEMSI::Configs* configsIn, DEMSI::LammpsInstance* lammpsInstanceIn, DEMSI::Grid* gridIn, DEMSI::Particles* particlesIn) {

  partition = partitionIn;
  log = logIn;
  configs = configsIn;
  lammpsInstance = lammpsInstanceIn;
  grid = gridIn;
  particles = particlesIn;

  init_initial_distribution();
  
}

// setup initial particle distribution
void Tessellation::init_initial_distribution(void) {

  // check if perform remapping
  bool useRemapping = false;
  if (configs->exists({"ConfigGroup:remapping"})) {
    configs->get({"ConfigGroup:remapping","Config:useRemapping"}, useRemapping);
  } // remapping group exists

  if (useRemapping) {

    // list of global indices for halo elements
    std::vector<int> initHaloIndices;

    // read in the initial particle distribution
    read_initial_distribution(initHaloIndices);

    // initialize the init particle exchange
    init_particles_halo(initHaloIndices);

    // set the local to global particle index mapping
    set_local_to_global_index_map_init();

    // set init particle polygons
    set_init_particles_polygons();

    // convert global indexing of cellsOnCell to local
    set_local_cells_on_cell();

    // set up the polygon edge connectivity
    set_edge_connectivity();

    // set up edge halo exchange
    init_edge_halo();

  } // useRemapping

} // Tessellation::init_initial_distribution

// read in the initial particle distribution
void Tessellation::read_initial_distribution(std::vector<int> &initHaloIndices) {

  // mpi/netcdf errors
  char mpiErrBuffer[MPI_MAX_ERROR_STRING];
  int mpiErrLen;
  int err;

  // input configs
  std::string particleInputFilename;
  configs->get({"ConfigGroup:particleInput","Config:particleInputFile"}, particleInputFilename);

  // netcdf
  int ncID;

  // broadcast variables
  int nParticlesInitTotal;
  int maxVerticesInit;

  // master task reads in particles
  if (partition->on_master_proc()) {

    // open file
    err = nc_open(particleInputFilename.c_str(), 0, &ncID);
    log->check(err == NC_NOERR, "Problem opening file: ", nc_strerror(err), ", for file: ", particleInputFilename);

    // get particle number
    int dimidParticles;
    err = nc_inq_dimid(ncID, "nParticles", &dimidParticles);
    log->check(err == NC_NOERR, "Problem getting dimid for nParticles: ", nc_strerror(err), ", for file: ", particleInputFilename);

    size_t nParticlesTmp;
    err = nc_inq_dimlen(ncID, dimidParticles, &nParticlesTmp);
    log->check(err == NC_NOERR, "Problem getting nParticles dimension: ", nc_strerror(err), ", for file: ", particleInputFilename);
    nParticlesInitTotal = (int) nParticlesTmp;

    // get max vertices number
    int dimidMaxVertices;
    err = nc_inq_dimid(ncID, "maxVertices", &dimidMaxVertices);
    log->check(err == NC_NOERR, "Problem getting dimid for maxVertices: ", nc_strerror(err), ", for file: ", particleInputFilename);

    size_t maxVerticesTmp;
    err = nc_inq_dimlen(ncID, dimidMaxVertices, &maxVerticesTmp);
    log->check(err == NC_NOERR, "Problem getting maxVertices dimension: ", nc_strerror(err), ", for file: ", particleInputFilename);
    maxVerticesInit = (int) maxVerticesTmp;

  } // on master proc

  // send broadcast variables to all processors
  err = MPI_Bcast(&nParticlesInitTotal, 1, MPI_INT, partition->master_proc(), partition->comm());
  MPI_Error_string(err, mpiErrBuffer, &mpiErrLen);
  log->check(err == MPI_SUCCESS, "Problem broadcasting nParticlesInitTotal: ", (std::string) mpiErrBuffer);

  err = MPI_Bcast(&maxVerticesInit, 1, MPI_INT, partition->master_proc(), partition->comm());
  MPI_Error_string(err, mpiErrBuffer, &mpiErrLen);
  log->check(err == MPI_SUCCESS, "Problem broadcasting nParticlesInitTotal: ", (std::string) mpiErrBuffer);

  // how to distribute the input data amongst all processors
  int *sendCounts = new int[partition->nprocs()];
  int *displacements = new int[partition->nprocs()];
  int evenDivision = nParticlesInitTotal / partition->nprocs();
  int remainder = nParticlesInitTotal - evenDivision*partition->nprocs();

  for (int iProc = 0 ; iProc < partition->nprocs() ; iProc++) {
    sendCounts[iProc] = evenDivision;
  } // iProc
  for (int iProc = 0 ; iProc < remainder ; iProc++) {
    sendCounts[iProc] += 1;
  } // iProc

  // local nParticlesInit
  int nParticlesInitList = sendCounts[partition->this_proc()];

  displacements[0] = 0;
  for (int iProc = 1 ; iProc < partition->nprocs() ; iProc++) {
    displacements[iProc] = displacements[iProc-1] + sendCounts[iProc-1];
  } // iProc

  // local arrays for receiving variables
  double *xLocal, *yLocal, *radiusLocal, *xVertexLocal, *yVertexLocal;
  int *typeLocal, *nVerticesOnCellLocal, *cellsOnCellLocal;

  get_init_variable_and_scatter(typeLocal, ncID, (std::string) "type", nParticlesInitTotal, 1, sendCounts, displacements, partition, log);
  get_init_variable_and_scatter(xLocal, yLocal, ncID, (std::string) "x", nParticlesInitTotal, sendCounts, displacements, partition, log);
  get_init_variable_and_scatter(radiusLocal, ncID, (std::string) "radius", nParticlesInitTotal, 1, sendCounts, displacements, partition, log);
  get_init_variable_and_scatter(nVerticesOnCellLocal, ncID, (std::string) "nVerticesOnCell", nParticlesInitTotal, 1, sendCounts, displacements, partition, log);
  get_init_variable_and_scatter(xVertexLocal, ncID, (std::string) "xVertex", nParticlesInitTotal, maxVerticesInit, sendCounts, displacements, partition, log);
  get_init_variable_and_scatter(yVertexLocal, ncID, (std::string) "yVertex", nParticlesInitTotal, maxVerticesInit, sendCounts, displacements, partition, log);
  get_init_variable_and_scatter(cellsOnCellLocal, ncID, (std::string) "cellsOnCell", nParticlesInitTotal, maxVerticesInit, sendCounts, displacements, partition, log);

  delete [] sendCounts;

  // maximum cell size
  double maxCellSizeInitLocal = 0.0;
  for (int iParticle = 0 ; iParticle < nParticlesInitList ; iParticle++) {
    for (int iVertex1 = 0 ; iVertex1 < nVerticesOnCellLocal[iParticle] ; iVertex1++) {
      int ij1 = iVertex1 + maxVerticesInit*iParticle;
      double xVertex1 = xVertexLocal[ij1];
      double yVertex1 = yVertexLocal[ij1];
      for (int iVertex2 = iVertex1 + 1 ; iVertex2 < nVerticesOnCellLocal[iParticle] ; iVertex2++) {
	int ij2 = iVertex2 + maxVerticesInit*iParticle;
	double xVertex2 = xVertexLocal[ij2];
	double yVertex2 = yVertexLocal[ij2];
	maxCellSizeInitLocal = std::max(std::sqrt(std::pow((xVertex2 - xVertex1), 2) + std::pow((yVertex2 - yVertex1), 2)),maxCellSizeInitLocal);
      } // iVertex2
    } // iVertex1
  } // iParticle
  err = MPI_Allreduce(&maxCellSizeInitLocal, &maxCellSizeInit, 1, MPI_DOUBLE, MPI_MAX, partition->comm());
  MPI_Error_string(err, mpiErrBuffer, &mpiErrLen);
  log->check(err == MPI_SUCCESS, "read_initial_distribution: MPI_Allreduce maxCellSizeInit failed: ", (std::string) mpiErrBuffer);
  (*log)() << "maxCellSizeInit: " << maxCellSizeInit << std::endl;

  // maximum half cell size
  double maxHalfCellSizeInitLocal = 0.0;
  for (int iParticle = 0 ; iParticle < nParticlesInitList ; iParticle++) {
    double x = xLocal[iParticle];
    double y = yLocal[iParticle];
    for (int iVertex = 0 ; iVertex < nVerticesOnCellLocal[iParticle] ; iVertex++) {
      int ij = iVertex + maxVerticesInit*iParticle;
      double xVertex = xVertexLocal[ij];
      double yVertex = yVertexLocal[ij];
      maxHalfCellSizeInitLocal = std::max(std::sqrt(std::pow((x - xVertex), 2) + std::pow((y - yVertex), 2)),maxHalfCellSizeInitLocal);
    } // iVertex
  } // iParticle
  err = MPI_Allreduce(&maxHalfCellSizeInitLocal, &maxHalfCellSizeInit, 1, MPI_DOUBLE, MPI_MAX, partition->comm());
  MPI_Error_string(err, mpiErrBuffer, &mpiErrLen);
  log->check(err == MPI_SUCCESS, "read_initial_distribution: MPI_Allreduce maxHalfCellSizeInit failed: ", (std::string) mpiErrBuffer);
  (*log)() << "maxHalfCellSizeInit: " << maxHalfCellSizeInit << std::endl;

  int nDoublesPerParticle = 3 + 2 * maxVerticesInit;
  std::vector<double> initParticleListDouble;
  initParticleListDouble.reserve(nParticlesInitList*nDoublesPerParticle);
  for (int iParticle = 0 ; iParticle < nParticlesInitList ; iParticle++) {
    initParticleListDouble.push_back(xLocal[iParticle]);
    initParticleListDouble.push_back(yLocal[iParticle]);
    initParticleListDouble.push_back(radiusLocal[iParticle]);
    for (int iVertex = 0 ; iVertex < maxVerticesInit ; iVertex++) {
      int ij = iVertex + maxVerticesInit*iParticle;
      initParticleListDouble.push_back(xVertexLocal[ij]);
      initParticleListDouble.push_back(yVertexLocal[ij]);
    } // iVertex
  } // iParticle

  int nIntegersPerParticle = 3 + 1 * maxVerticesInit;
  std::vector<int> initParticleListInt;
  initParticleListInt.reserve(nParticlesInitList*nIntegersPerParticle);
  for (int iParticle = 0 ; iParticle < nParticlesInitList ; iParticle++) {
    int globalIndexInit = displacements[partition->this_proc()] + iParticle;
    initParticleListInt.push_back(globalIndexInit);
    initParticleListInt.push_back(nVerticesOnCellLocal[iParticle]);
    initParticleListInt.push_back(typeLocal[iParticle]);
    for (int iVertex = 0 ; iVertex < maxVerticesInit ; iVertex++) {
      int ij = iVertex + maxVerticesInit*iParticle;
      initParticleListInt.push_back(cellsOnCellLocal[ij]);
    } // iVertex
  } // iParticle

  delete [] displacements;

  // clean
  delete [] xLocal;
  delete [] yLocal;
  delete [] radiusLocal;
  delete [] xVertexLocal;
  delete [] yVertexLocal;
  delete [] nVerticesOnCellLocal;
  delete [] cellsOnCellLocal;

  // copy the lists for second pass
  std::vector<double> initParticleListDouble2;
  std::vector<int> initParticleListInt2;
  std::copy(initParticleListDouble.begin(), initParticleListDouble.end(), std::back_inserter(initParticleListDouble2));
  std::copy(initParticleListInt.begin(), initParticleListInt.end(), std::back_inserter(initParticleListInt2));

  // list that accumulates elements belonging to this processor
  std::vector<double> initParticleListDoubleThisProc;
  std::vector<int> initParticleListIntThisProc;

  std::vector<int> initOwnedIndices;

  // now send around all processors
  DEMSI::ShareLists* particleLists = new DEMSI::ShareLists(&initParticleListInt,&initParticleListDouble,partition,log);
  do {

    int nParticlesThisList = initParticleListDouble.size() / nDoublesPerParticle;

    for (int iParticle = 0 ; iParticle < nParticlesThisList ; iParticle++) {

      // position of the element
      double x = initParticleListDouble[iParticle*nDoublesPerParticle];
      double y = initParticleListDouble[iParticle*nDoublesPerParticle+1];
      int globalIndexInit = initParticleListInt[iParticle*nIntegersPerParticle];

      // check if element is in our local domain - then it belongs to this processor
      if (grid->in_domain(x, y)) {

	initOwnedIndices.push_back(globalIndexInit);

	for (int i = 0 ; i < nDoublesPerParticle ; i++) {
	  initParticleListDoubleThisProc.push_back(initParticleListDouble[iParticle*nDoublesPerParticle+i]);
	}
	for (int i = 0 ; i < nIntegersPerParticle ; i++) {
	  initParticleListIntThisProc.push_back(initParticleListInt[iParticle*nIntegersPerParticle+i]);
	}

      }

    } // iParticle

  } while (particleLists->iterate());
  delete particleLists;

  // width of the halo region around domain
  double skin = lammpsInstance->lmp->neighbor->skin;
  double haloWidth = 2.0*maxHalfCellSizeInit + skin;
  (*log)() << "haloWidth: " << haloWidth << " " << this->maxCellSizeInit << " " << skin << std::endl;

  std::vector<double> initParticleListDoubleHalo;
  std::vector<int> initParticleListIntHalo;

  // now send around all processors again for halo
  particleLists = new DEMSI::ShareLists(&initParticleListInt2,&initParticleListDouble2,partition,log);
  do {

    int nParticlesThisList = initParticleListDouble2.size() / nDoublesPerParticle;

    for (int iParticle = 0 ; iParticle < nParticlesThisList ; iParticle++) {

      double x = initParticleListDouble2[iParticle*nDoublesPerParticle];
      double y = initParticleListDouble2[iParticle*nDoublesPerParticle+1];
      int globalIndexInit = initParticleListInt2[iParticle*nIntegersPerParticle];
      int nVerticesOnCell = initParticleListInt2[iParticle*nIntegersPerParticle+1];

      bool lNeighbour = false;
      if (grid->in_halo_region(x, y, 3.0*maxCellSizeInit)) {
	for (int iCellOnCell = 0 ; iCellOnCell < nVerticesOnCell ; iCellOnCell++) {
	  int globalIndexInitNeighbour = initParticleListInt2[iParticle*nIntegersPerParticle + 3 + iCellOnCell];
	  std::vector<int>::iterator initOwnedIndicesItr = std::find(initOwnedIndices.begin(), initOwnedIndices.end(), globalIndexInitNeighbour);
	  if (initOwnedIndicesItr != initOwnedIndices.end()) lNeighbour = true;
	} // iVertex
      }

      if (lNeighbour or grid->in_halo_region(x, y, haloWidth)) {

	initHaloIndices.push_back(globalIndexInit);

	for (int i = 0 ; i < nDoublesPerParticle ; i++) {
	  initParticleListDoubleHalo.push_back(initParticleListDouble2[iParticle*nDoublesPerParticle+i]);
	}
	for (int i = 0 ; i < nIntegersPerParticle ; i++) {
	  initParticleListIntHalo.push_back(initParticleListInt2[iParticle*nIntegersPerParticle+i]);
	}

      } // in halo region

    } // iParticle

  } while (particleLists->iterate());
  delete particleLists;

  // final size of initial distribution owned by this processor
  this->nParticlesInit = initParticleListIntThisProc.size() / nIntegersPerParticle;

  // final size of halo owned by this processor
  this->nParticlesHaloInit = initParticleListIntHalo.size() / nIntegersPerParticle;

  // allocate initial distribution
  this->nParticlesAllInit = this->nParticlesInit + this->nParticlesHaloInit;

  (*log)() << "nParticlesInit:    " << this->nParticlesInit << std::endl;
  (*log)() << "nParticlesHaloInit:" << this->nParticlesHaloInit << std::endl;
  (*log)() << "nParticlesAllInit: " << this->nParticlesAllInit << std::endl;

  this->maxVerticesInit = maxVerticesInit;

  this->globalIndexInit = new int[this->nParticlesAllInit];

  this->typeInit = new int[this->nParticlesAllInit];
  this->xInit = new double[this->nParticlesAllInit];
  this->yInit = new double[this->nParticlesAllInit];
  this->radiusInit = new double[this->nParticlesAllInit];

  this->nVerticesOnCellInit = new int[this->nParticlesAllInit];

  this->xVertexInit = new double*[this->nParticlesAllInit];
  this->yVertexInit = new double*[this->nParticlesAllInit];
  for (int i = 0 ; i < this->nParticlesAllInit ; i++) {
    this->xVertexInit[i] = new double[maxVerticesInit];
    this->yVertexInit[i] = new double[maxVerticesInit];
  } // i

  this->cellsOnCellGlobal = new int*[this->nParticlesAllInit];
  this->cellsOnCellLocal = new int*[this->nParticlesAllInit];
  for (int i = 0 ; i < this->nParticlesAllInit ; i++) {
    this->cellsOnCellGlobal[i] = new int[maxVerticesInit];
    this->cellsOnCellLocal[i] = new int[maxVerticesInit];
  } // i

  // owned
  for (int iParticle = 0 ; iParticle < this->nParticlesInit ; iParticle++) {

    // integer quantities
    int indexPrevInt = iParticle * nIntegersPerParticle;
    this->globalIndexInit[iParticle] = initParticleListIntThisProc[indexPrevInt];
    this->nVerticesOnCellInit[iParticle] = initParticleListIntThisProc[indexPrevInt + 1];
    this->typeInit[iParticle] = initParticleListIntThisProc[indexPrevInt + 2];
    for (int iVertex = 0 ; iVertex < maxVerticesInit ; iVertex++) {
      this->cellsOnCellGlobal[iParticle][iVertex] = initParticleListIntThisProc[indexPrevInt + 3 + iVertex];
    } // iVertex

    // double quantities
    int indexPrevDouble = iParticle * nDoublesPerParticle;
    this->xInit     [iParticle] = initParticleListDoubleThisProc[indexPrevDouble];
    this->yInit     [iParticle] = initParticleListDoubleThisProc[indexPrevDouble + 1];
    this->radiusInit[iParticle] = initParticleListDoubleThisProc[indexPrevDouble + 2];
    for (int iVertex = 0 ; iVertex < maxVerticesInit ; iVertex++) {
      this->xVertexInit[iParticle][iVertex] = initParticleListDoubleThisProc[indexPrevDouble + 3 + iVertex * 2];
      this->yVertexInit[iParticle][iVertex] = initParticleListDoubleThisProc[indexPrevDouble + 3 + iVertex * 2 + 1];
    } // iVertex

  } // iParticle

  // halo
  for (int iParticleHalo = 0 ; iParticleHalo < this->nParticlesHaloInit ; iParticleHalo++) {

    int iParticle = iParticleHalo + this->nParticlesInit;

    // integer quantities
    int indexPrevInt = iParticleHalo * nIntegersPerParticle;
    this->globalIndexInit[iParticle] = initParticleListIntHalo[indexPrevInt];
    this->nVerticesOnCellInit[iParticle] = initParticleListIntHalo[indexPrevInt + 1];
    this->typeInit[iParticle] = initParticleListIntHalo[indexPrevInt + 2];
    for (int iVertex = 0 ; iVertex < maxVerticesInit ; iVertex++) {
      this->cellsOnCellGlobal[iParticle][iVertex] = initParticleListIntThisProc[indexPrevInt + 3 + iVertex];
    } // iVertex

    // double quantities
    int indexPrevDouble = iParticleHalo * nDoublesPerParticle;
    this->xInit     [iParticle] = initParticleListDoubleHalo[indexPrevDouble];
    this->yInit     [iParticle] = initParticleListDoubleHalo[indexPrevDouble + 1];
    this->radiusInit[iParticle] = initParticleListDoubleHalo[indexPrevDouble + 2];
    for (int iVertex = 0 ; iVertex < maxVerticesInit ; iVertex++) {
      this->xVertexInit[iParticle][iVertex] = initParticleListDoubleHalo[indexPrevDouble + 3 + iVertex * 2];
      this->yVertexInit[iParticle][iVertex] = initParticleListDoubleHalo[indexPrevDouble + 3 + iVertex * 2 + 1];
    }

  } // iParticle

  // correct type to be 0 or 2
  for (int iParticleInit = 0 ; iParticleInit < this->nParticlesAllInit ; iParticleInit++) {
    if (typeInit[iParticleInit] == 1) {
      typeInit[iParticleInit] = 0;
    }
  } // iParticleInit

  // max numver of vertices
  nVerticesOnCellInitMax = *std::max_element(nVerticesOnCellInit, nVerticesOnCellInit+nParticlesAllInit);

  // Debugging init distribution
  //init_distribution_diagnostics();

} // Tessellation::read_initial_distribution

// Initialize halo exchange for initial tessellation particles
void Tessellation::init_particles_halo(std::vector<int> initHaloIndices) {

  (*log)() << "initHaloIndices:" << initHaloIndices.size() << std::endl;

  // global ids of halo particles on owned processor
  std::vector<int> initOwnedIndices;

  // now send halo global indices for this processor to all other processors
  std::vector<int> initHaloIndicesTransit(initHaloIndices);
  DEMSI::ShareLists* haloLists = new DEMSI::ShareLists(&initHaloIndicesTransit,NULL,partition,log);
  while (haloLists->iterate()) {

    // find local index in halo list
    std::vector<int>::iterator initHaloIndicesItr;
    for (int i=0 ; i < this->nParticlesInit ; i++) {

      initHaloIndicesItr = std::find(initHaloIndicesTransit.begin(), initHaloIndicesTransit.end(), this->globalIndexInit[i]);
      if (initHaloIndicesItr != initHaloIndicesTransit.end()) {
	localOwnedParticlesIndicesExchInit[haloLists->iProc_origin()].push_back(this->globalIndexInit[i]);
	initOwnedIndices.push_back(this->globalIndexInit[i]);
      }

    }

  }
  delete haloLists;

  // now push the owned list around all the processors
  DEMSI::ShareLists* ownedLists = new DEMSI::ShareLists(&initOwnedIndices,NULL,partition,log);
  while (ownedLists->iterate()) {

    std::vector<int>::iterator initOwnedIndicesItr;
    for (int i=0 ; i < initHaloIndices.size() ; i++) {

      initOwnedIndicesItr = std::find(initOwnedIndices.begin(), initOwnedIndices.end(), initHaloIndices[i]);
      if (initOwnedIndicesItr != initOwnedIndices.end()) {
	localHaloParticlesIndicesExchInit[ownedLists->iProc_origin()].push_back(initHaloIndices[i]);
      }

    }

  }
  delete ownedLists;

  // now each processor has a map of processor id to global indices for both sending and receiving
  // Make sure global indices are in ascending order
  for (int iProc = 0 ; iProc < partition->nprocs() ; iProc++) {
    std::sort(localOwnedParticlesIndicesExchInit[iProc].begin(), localOwnedParticlesIndicesExchInit[iProc].end());
    std::sort(localHaloParticlesIndicesExchInit[iProc].begin(), localHaloParticlesIndicesExchInit[iProc].end());
  } // iProc

  // now convert these maps from global indices to local indices
  // first make a mapping from global index to local index
  std::map<int,int> globalToLocalIndex;
  for (int iParticle = 0 ; iParticle < this->nParticlesInit + this->nParticlesHaloInit ; iParticle++) {
    globalToLocalIndex[this->globalIndexInit[iParticle]] = iParticle;
  } // iParticle

  // now convert the indices
  for (int iProc = 0 ; iProc < partition->nprocs() ; iProc++) {
    for (int i = 0 ; i < localOwnedParticlesIndicesExchInit[iProc].size() ; i++) {
      localOwnedParticlesIndicesExchInit[iProc][i] = globalToLocalIndex[localOwnedParticlesIndicesExchInit[iProc][i]];
    } // i
    for (int i = 0 ; i < localHaloParticlesIndicesExchInit[iProc].size() ; i++) {
      localHaloParticlesIndicesExchInit[iProc][i] = globalToLocalIndex[localHaloParticlesIndicesExchInit[iProc][i]];
    } // i
  } // iProc

  // global id test
  int *globalIndexInitHalo = new int[nParticlesAllInit];
  for (int iParticleInit = 0 ; iParticleInit < nParticlesInit ; iParticleInit++) {
    globalIndexInitHalo[iParticleInit] = globalIndexInit[iParticleInit];
  } // iParticleInit

  update_init_halo_particles(globalIndexInitHalo);

  bool fail = false;
  for (int iParticleInit = nParticlesInit ; iParticleInit < nParticlesAllInit ; iParticleInit++) {
    if (globalIndexInitHalo[iParticleInit] != globalIndexInit[iParticleInit]) {
      (*log)() << "init_particles_halo: diff global id: " << iParticleInit << " " << globalIndexInitHalo[iParticleInit] << " " << globalIndexInit[iParticleInit] << std::endl;
      fail = true;
    }
  } // iParticleInit
  log->check(not fail, "init_particles_halo test failed.");

  delete [] globalIndexInitHalo;

} // Tessellation::init_particles_halo

// create the polygon objects of the initial tessellation
void Tessellation::set_init_particles_polygons(void) {

  // fill polygon vector
  for (int iParticleInit = 0 ; iParticleInit < this->nParticlesAllInit ; iParticleInit++) {
    std::vector<double> xVertices, yVertices;
    for (int iVertexOnCell = 0 ; iVertexOnCell < this->nVerticesOnCellInit[iParticleInit] ; iVertexOnCell++) {
      xVertices.push_back(this->xVertexInit[iParticleInit][iVertexOnCell]);
      yVertices.push_back(this->yVertexInit[iParticleInit][iVertexOnCell]);
    } // iVertex
    this->initPolygons.push_back({xVertices,yVertices});
  } // iParticleInit

} // Tessellation::set_init_particles_polygons

// Determine local cellsOnCell index from global
void Tessellation::set_local_cells_on_cell(void) {

  // map from global to local index
  std::map<int,int> globalToLocalIndex;
  for (int iParticleInit = 0 ; iParticleInit < nParticlesAllInit ; iParticleInit++) {
    globalToLocalIndex[globalIndexInit[iParticleInit]] = iParticleInit;
  } // iParticleInit

  for (int iParticleInit = 0 ; iParticleInit < nParticlesAllInit ; iParticleInit++) {
    for (int iVertexOnCell = 0 ; iVertexOnCell < nVerticesOnCellInit[iParticleInit] ; iVertexOnCell++) {
      int iParticleInitNeighbour = cellsOnCellGlobal[iParticleInit][iVertexOnCell];
      if (iParticleInitNeighbour == -1) {
	cellsOnCellLocal[iParticleInit][iVertexOnCell] = -1;
      } else {
	if (globalToLocalIndex.find(iParticleInitNeighbour) != globalToLocalIndex.end()) {
	  cellsOnCellLocal[iParticleInit][iVertexOnCell] = globalToLocalIndex[iParticleInitNeighbour];
	} else {
	  cellsOnCellLocal[iParticleInit][iVertexOnCell] = -2;
	}
      }
    } // iVertexOnCell
  } // iParticleInit

} // Tessellation::set_local_cells_on_cell

// set the edge to cell connectivity
void Tessellation::set_edge_connectivity(void) {

  // count edge types
  nEdgesInit = 0;
  nEdgesHaloInit = 0;

  for (int iParticleInit = 0 ; iParticleInit < nParticlesInit ; iParticleInit++) {
    int iParticleInitGlobal = globalIndexInit[iParticleInit];

    for (int iCellOnCell = 0 ; iCellOnCell < nVerticesOnCellInit[iParticleInit] ; iCellOnCell++) {
      int iParticleInitNeighbour = cellsOnCellLocal[iParticleInit][iCellOnCell];

      if (iParticleInitNeighbour == -1) { // On boundary
	nEdgesInit++;
      } else {

	int iParticleInitNeighbourGlobal = globalIndexInit[iParticleInitNeighbour];

	if (iParticleInitNeighbour < nParticlesInit) { // both particles owned
	  if (iParticleInit < iParticleInitNeighbour) { // prevent double counting
	    nEdgesInit++;
	  }
	} else { // neighbour particle not owned
	  if (iParticleInitGlobal < iParticleInitNeighbourGlobal) { // smaller global index owns
	    nEdgesInit++;
	  } else {
	    nEdgesHaloInit++;
	  }
	}

      }

    } // iVertexOnCell

  } // iParticleInit
  nEdgesAllInit = nEdgesInit + nEdgesHaloInit;

  // nEdges on procs
  char mpiErrBuffer[MPI_MAX_ERROR_STRING];
  int mpiErrLen;
  int mpiErr;

  // number of owned init particles on each processor
  nEdgesInitOnAllProcs = new int[partition->nprocs()];

  mpiErr = MPI_Allgather(&nEdgesInit, 1, MPI_INT, nEdgesInitOnAllProcs, 1, MPI_INT, partition->comm());
  MPI_Error_string(mpiErr, mpiErrBuffer, &mpiErrLen);
  log->check(mpiErr == MPI_SUCCESS, "set_edge_connectivity: MPI_Allgather failed: ", (std::string) mpiErrBuffer);

  // cumulative number of owned init particles on previous processors
  nEdgesInitOnPrevProc = new int[partition->nprocs()];

  nEdgesInitOnPrevProc[0] = 0;
  for (int iProc = 1 ; iProc < partition->nprocs() ; iProc++) {
    nEdgesInitOnPrevProc[iProc] = nEdgesInitOnPrevProc[iProc-1] + nEdgesInitOnAllProcs[iProc-1];
  } // iProc

  int nEdgesTotal = 0;
  for (int iProc = 0 ; iProc < partition->nprocs() ; iProc++) {
    nEdgesTotal += nEdgesInitOnAllProcs[iProc];
  } // iProc

  (*log)(DEMSI::LOG::DEBUG) << "nEdgesInit:     " << nEdgesInit << std::endl;
  (*log)(DEMSI::LOG::DEBUG) << "nEdgesHaloInit: " << nEdgesHaloInit << std::endl;
  (*log)(DEMSI::LOG::DEBUG) << "nEdgesAllInit:  " << nEdgesAllInit << std::endl;
  (*log)(DEMSI::LOG::DEBUG) << "nEdgesTotal:    " << nEdgesTotal << std::endl;
  for (int iProc = 0 ; iProc < partition->nprocs() ; iProc++) {
    (*log)(DEMSI::LOG::DEBUG) << "nEdgesInitOnAllProcs nEdgesInitOnPrevProc: " << nEdgesInitOnAllProcs[iProc] << " " << nEdgesInitOnPrevProc[iProc] << std::endl;
  } // iProc

  // cellsOnEdgeLocal
  cellsOnEdgeGlobal = new int*[nEdgesAllInit];
  cellsOnEdgeLocal = new int*[nEdgesAllInit];
  for (int i = 0 ; i < nEdgesAllInit ; i++) {
    cellsOnEdgeGlobal[i] = new int[2];
    cellsOnEdgeLocal[i] = new int[2];
  } // i

  // uniquely loop through edges
  int iEdgeInit = 0;
  int iEdgeHaloInit = 0;

  for (int iParticleInit = 0 ; iParticleInit < nParticlesInit ; iParticleInit++) {
    int iParticleInitGlobal = globalIndexInit[iParticleInit];

    for (int iCellOnCell = 0 ; iCellOnCell < nVerticesOnCellInit[iParticleInit] ; iCellOnCell++) {
      int iParticleInitNeighbour = cellsOnCellLocal[iParticleInit][iCellOnCell];

      if (iParticleInitNeighbour == -1) { // On boundary
	cellsOnEdgeLocal [iEdgeInit][0] = iParticleInit;
	cellsOnEdgeLocal [iEdgeInit][1] = -1;
	cellsOnEdgeGlobal[iEdgeInit][0] = iParticleInitGlobal;
	cellsOnEdgeGlobal[iEdgeInit][1] = -1;
	iEdgeInit++;
      } else {

	int iParticleInitNeighbourGlobal = globalIndexInit[iParticleInitNeighbour];

	if (iParticleInitNeighbour < nParticlesInit) { // both particles owned
	  if (iParticleInit < iParticleInitNeighbour) { // prevent double counting
	    cellsOnEdgeLocal [iEdgeInit][0] = iParticleInit;
	    cellsOnEdgeLocal [iEdgeInit][1] = iParticleInitNeighbour;
	    cellsOnEdgeGlobal[iEdgeInit][0] = iParticleInitGlobal;
	    cellsOnEdgeGlobal[iEdgeInit][1] = iParticleInitNeighbourGlobal;
	    iEdgeInit++;
	  }
	} else { // neighbour particle not owned
	  if (iParticleInitGlobal < iParticleInitNeighbourGlobal) { // smaller global index owns
	    cellsOnEdgeLocal [iEdgeInit][0] = iParticleInit;
	    cellsOnEdgeLocal [iEdgeInit][1] = iParticleInitNeighbour;
	    cellsOnEdgeGlobal[iEdgeInit][0] = iParticleInitGlobal;
	    cellsOnEdgeGlobal[iEdgeInit][1] = iParticleInitNeighbourGlobal;
	    iEdgeInit++;
	  } else {
	    cellsOnEdgeLocal [nEdgesInit + iEdgeHaloInit][0] = iParticleInit;
	    cellsOnEdgeLocal [nEdgesInit + iEdgeHaloInit][1] = iParticleInitNeighbour;
	    cellsOnEdgeGlobal[nEdgesInit + iEdgeHaloInit][0] = iParticleInitGlobal;
	    cellsOnEdgeGlobal[nEdgesInit + iEdgeHaloInit][1] = iParticleInitNeighbourGlobal;
	    iEdgeHaloInit++;
	  }
	}

      }

    } // iVertexOnCell

  } // iParticleInit

  // edgesOnCell
  edgesOnCellGlobal = new int*[nParticlesAllInit];
  edgesOnCellLocal = new int*[nParticlesAllInit];
  for (int i = 0 ; i < nParticlesAllInit ; i++) {
    edgesOnCellGlobal[i] = new int[maxVerticesInit];
    edgesOnCellLocal[i] = new int[maxVerticesInit];
  } // i
  for (int iParticleInit = 0 ; iParticleInit < nParticlesAllInit ; iParticleInit++) {
    for (int iCellOnCell = 0 ; iCellOnCell < maxVerticesInit ; iCellOnCell++) {
      edgesOnCellLocal[iParticleInit][iCellOnCell] = -1;
    } // iCellOnCell
  } // iParticleInit

  int *iMissing = new int[nParticlesAllInit];
  for (int i = 0 ; i < nParticlesAllInit ; i++) {
    iMissing[i] = 0;
  } // i

  for (int iEdgeInit = 0 ; iEdgeInit < nEdgesAllInit ; iEdgeInit++) {

    int iParticleInit1 = cellsOnEdgeLocal[iEdgeInit][0];
    int iParticleInit2 = cellsOnEdgeLocal[iEdgeInit][1];

    // particle 1
    if (iParticleInit1 != -1) {
      if (iParticleInit2 != -1) {
	for (int iCellOnCell = 0 ; iCellOnCell < nVerticesOnCellInit[iParticleInit1] ; iCellOnCell++) {
	  int iParticleInitNeighbour = cellsOnCellLocal[iParticleInit1][iCellOnCell];
	  if (iParticleInit2 == iParticleInitNeighbour) {
	    edgesOnCellLocal [iParticleInit1][iCellOnCell] = iEdgeInit;
	    edgesOnCellGlobal[iParticleInit1][iCellOnCell] = iEdgeInit + nEdgesInitOnPrevProc[partition->this_proc()];
	  }
	} // iCellOnCell
      } else {
	int nMissing = 0;
	for (int iCellOnCell = 0 ; iCellOnCell < nVerticesOnCellInit[iParticleInit1] ; iCellOnCell++) {
	  int iParticleInitNeighbour = cellsOnCellLocal[iParticleInit1][iCellOnCell];
	  if (iParticleInitNeighbour == -1) {
	    if (nMissing == iMissing[iParticleInit1]) {
	      edgesOnCellLocal [iParticleInit1][iCellOnCell] = iEdgeInit;
	      edgesOnCellGlobal[iParticleInit1][iCellOnCell] = iEdgeInit + nEdgesInitOnPrevProc[partition->this_proc()];
	      iMissing[iParticleInit1]++;
	      break;
	    }
	    nMissing++;
	  }
	} // iCellOnCell
      }
    }

    // particle 2
    if (iParticleInit2 != -1) {
      if (iParticleInit1 != -1) {
	for (int iCellOnCell = 0 ; iCellOnCell < nVerticesOnCellInit[iParticleInit2] ; iCellOnCell++) {
	  int iParticleInitNeighbour = cellsOnCellLocal[iParticleInit2][iCellOnCell];
	  if (iParticleInit1 == iParticleInitNeighbour) {
	    edgesOnCellLocal [iParticleInit2][iCellOnCell] = iEdgeInit;
	    edgesOnCellGlobal[iParticleInit2][iCellOnCell] = iEdgeInit + nEdgesInitOnPrevProc[partition->this_proc()];
	  }
	} // iCellOnCell
      } else {
	int nMissing = 0;
	for (int iCellOnCell = 0 ; iCellOnCell < nVerticesOnCellInit[iParticleInit2] ; iCellOnCell++) {
	  int iParticleInitNeighbour = cellsOnCellLocal[iParticleInit2][iCellOnCell];
	  if (iParticleInitNeighbour == -1) {
	    if (nMissing == iMissing[iParticleInit2]) {
	      edgesOnCellLocal [iParticleInit2][iCellOnCell] = iEdgeInit;
	      edgesOnCellGlobal[iParticleInit2][iCellOnCell] = iEdgeInit + nEdgesInitOnPrevProc[partition->this_proc()];
	      iMissing[iParticleInit2]++;
	      break;
	    }
	    nMissing++;
	  }
	} // iCellOnCell
      }
    }

  } // iEdge

  // clean up
  delete [] iMissing;

} // Tessellation::set_edge_connectivity

// Initialize halo exchange for initial tessellation edges
void Tessellation::init_edge_halo(void) {

  // local indices of edges that are halo elsewhere
  std::vector<int> iEdgesHaloOwned;

  // local indices of edges that are owned elsewhere
  std::vector<int> iEdgesHaloDistant;

  // uniquely loop through edges
  int iEdgeInit = 0;
  int iEdgeHaloInit = 0;

  for (int iParticleInit = 0 ; iParticleInit < nParticlesInit ; iParticleInit++) {
    int iParticleInitGlobal = globalIndexInit[iParticleInit];

    for (int iCellOnCell = 0 ; iCellOnCell < nVerticesOnCellInit[iParticleInit] ; iCellOnCell++) {
      int iParticleInitNeighbour = cellsOnCellLocal[iParticleInit][iCellOnCell];

      if (iParticleInitNeighbour == -1) { // On boundary
	iEdgeInit++;
      } else {

	int iParticleInitNeighbourGlobal = globalIndexInit[iParticleInitNeighbour];

	if (iParticleInitNeighbour < nParticlesInit) { // both particles owned
	  if (iParticleInit < iParticleInitNeighbour) { // prevent double counting
	    iEdgeInit++;
	  }
	} else { // neighbour particle not owned
	  if (iParticleInitGlobal < iParticleInitNeighbourGlobal) { // smaller global index owns
	    iEdgesHaloOwned.push_back(iEdgeInit);
	    iEdgeInit++;
	  } else {
	    iEdgesHaloDistant.push_back(nEdgesInit + iEdgeHaloInit);
	    iEdgeHaloInit++;
	  }
	}

      }

    } // iVertexOnCell

  } // iParticleInit

  // push round procs local owned halo edges
  std::vector<int> ownedHaloEdges;
  for (int i = 0 ; i < iEdgesHaloOwned.size() ; i++) {
    ownedHaloEdges.push_back(iEdgesHaloOwned[i]);
    ownedHaloEdges.push_back(iEdgesHaloOwned[i]+nEdgesInitOnPrevProc[partition->this_proc()]);
    ownedHaloEdges.push_back(cellsOnEdgeGlobal[iEdgesHaloOwned[i]][0]);
    ownedHaloEdges.push_back(cellsOnEdgeGlobal[iEdgesHaloOwned[i]][1]);
  } // i

  std::vector<int> iEdgeGlobalDistant;

  DEMSI::ShareLists* haloLists = new DEMSI::ShareLists(&ownedHaloEdges,NULL,partition,log);
  while (haloLists->iterate()) {

    std::vector<std::pair<int,int>> localHaloEdgeIndices;

    for (int iEdgeDistant = 0 ; iEdgeDistant < (ownedHaloEdges.size() / 4) ; iEdgeDistant++) {

      // distant halo edge
      int iEdgeDistantHalo            = ownedHaloEdges[iEdgeDistant*4];
      int iEdgeDistantHaloGlobal      = ownedHaloEdges[iEdgeDistant*4 + 1];
      int iParticleInitDistantOnEdge1 = ownedHaloEdges[iEdgeDistant*4 + 2];
      int iParticleInitDistantOnEdge2 = ownedHaloEdges[iEdgeDistant*4 + 3];

      // check if we need this edge info
      for (int iEdgeLocal = 0 ; iEdgeLocal < iEdgesHaloDistant.size() ; iEdgeLocal++) {

	int iEdge = iEdgesHaloDistant[iEdgeLocal];
	if ((cellsOnEdgeGlobal[iEdge][0] == iParticleInitDistantOnEdge1 and cellsOnEdgeGlobal[iEdge][1] == iParticleInitDistantOnEdge2) or
	    (cellsOnEdgeGlobal[iEdge][1] == iParticleInitDistantOnEdge1 and cellsOnEdgeGlobal[iEdge][0] == iParticleInitDistantOnEdge2)) {

	  localHaloEdgeIndices.push_back(std::make_pair(iEdgeDistantHaloGlobal, iEdge));
	  iEdgeGlobalDistant.push_back(iEdgeDistantHaloGlobal);

	}

      } // iEdgeLocal

    } // iEdge

    // sort by global index
    std::sort(localHaloEdgeIndices.begin(), localHaloEdgeIndices.end());

    // add to indices to receive
    for (int i = 0 ; i < localHaloEdgeIndices.size() ; i++) {
      localHaloEdgesIndicesExchInit[haloLists->iProc_origin()].push_back(localHaloEdgeIndices[i].second);
    } // i

  }
  delete haloLists;

  // now push the owned list around all the processors
  DEMSI::ShareLists* ownedLists = new DEMSI::ShareLists(&iEdgeGlobalDistant,NULL,partition,log);
  while (ownedLists->iterate()) {

    std::vector<std::pair<int,int>> localHaloEdgeIndices;

    for (int i = 0 ; i < iEdgeGlobalDistant.size() ; i++) {

      int iEdgeGlobal = iEdgeGlobalDistant[i];

      for (int j = 0 ; j < iEdgesHaloOwned.size() ; j++) {
	if (iEdgeGlobal == iEdgesHaloOwned[j]+nEdgesInitOnPrevProc[partition->this_proc()]) {
	  localHaloEdgeIndices.push_back(std::make_pair(iEdgeGlobal, iEdgesHaloOwned[j]));
	}
      } // j

    } // i

    // sort by global index
    std::sort(localHaloEdgeIndices.begin(), localHaloEdgeIndices.end());

    // add to indices to receive
    for (int i = 0 ; i < localHaloEdgeIndices.size() ; i++) {
      localOwnedEdgesIndicesExchInit[ownedLists->iProc_origin()].push_back(localHaloEdgeIndices[i].second);
    } // i

  }
  delete ownedLists;

  for (int iProc = 0 ; iProc < partition->nprocs() ; iProc++) {
    (*log)(DEMSI::LOG::DEBUG) << "localHaloEdgesIndicesExchInit[iProc].size():  " << iProc << " " << localHaloEdgesIndicesExchInit[iProc].size() << std::endl;
    (*log)(DEMSI::LOG::DEBUG) << "localOwnedEdgesIndicesExchInit[iProc].size(): " << iProc << " " << localOwnedEdgesIndicesExchInit[iProc].size() << std::endl;
  } // iProc

  // test
  int* testVar1 = new int[nEdgesAllInit];
  int* testVar2 = new int[nEdgesAllInit];
  for (int iEdgeInit = 0 ; iEdgeInit < nEdgesInit ; iEdgeInit++) {
    int iParticleInit1 = cellsOnEdgeLocal[iEdgeInit][0];
    int iParticleInit2 = cellsOnEdgeLocal[iEdgeInit][1];
    int globalIndexInit1;
    if (iParticleInit1 != -1) {
      globalIndexInit1 = globalIndexInit[iParticleInit1];
    } else {
      globalIndexInit1 = -1;
    }
    int globalIndexInit2;
    if (iParticleInit2 != -1) {
      globalIndexInit2 = globalIndexInit[iParticleInit2];
    } else {
      globalIndexInit2 = -1;
    }
    if (globalIndexInit1 < globalIndexInit2) {
      testVar1[iEdgeInit] = globalIndexInit1;
      testVar2[iEdgeInit] = globalIndexInit2;
    } else {
      testVar1[iEdgeInit] = globalIndexInit2;
      testVar2[iEdgeInit] = globalIndexInit1;
    }
  } // iEdgeInit

  update_init_halo_edges(testVar1);
  update_init_halo_edges(testVar2);

  for (int iEdgeInit = 0 ; iEdgeInit < nEdgesAllInit ; iEdgeInit++) {
    int iParticleInit1 = cellsOnEdgeLocal[iEdgeInit][0];
    int iParticleInit2 = cellsOnEdgeLocal[iEdgeInit][1];
    int globalIndexInit1;
    if (iParticleInit1 != -1) {
      globalIndexInit1 = globalIndexInit[iParticleInit1];
    } else {
      globalIndexInit1 = -1;
    }
    int globalIndexInit2;
    if (iParticleInit2 != -1) {
      globalIndexInit2 = globalIndexInit[iParticleInit2];
    } else {
      globalIndexInit2 = -1;
    }
    int testVarCheck1;
    int testVarCheck2;
    if (globalIndexInit1 < globalIndexInit2) {
      testVarCheck1 = globalIndexInit1;
      testVarCheck2 = globalIndexInit2;
    } else {
      testVarCheck1 = globalIndexInit2;
      testVarCheck2 = globalIndexInit1;
    }
    log->check(testVar1[iEdgeInit] == testVarCheck1 and testVar2[iEdgeInit] == testVarCheck2, "init_edge_halo: test failed");
  } // iEdgeInit

  delete [] testVar1;
  delete [] testVar2;

} // Tessellation::init_edge_halo

// Write out diagnostics for the initial particle distribution
void Tessellation::init_distribution_diagnostics(void) {

  std::ostringstream filename;
  std::ofstream diagFile;

  (*log)(DEMSI::LOG::DEBUG) << "Remapping init: this->nParticlesInit:     " << this->nParticlesInit << std::endl;
  (*log)(DEMSI::LOG::DEBUG) << "Remapping init: this->nParticlesHaloInit: " << this->nParticlesHaloInit << std::endl;
  (*log)(DEMSI::LOG::DEBUG) << "Remapping init: this->nParticlesAllInit:  " << this->nParticlesAllInit << std::endl;

  filename << "part_init_" << partition->this_proc() << "_" << partition->nprocs() << ".txt";
  diagFile.open(filename.str());
  for (int iParticle = 0 ; iParticle < this->nParticlesInit ; iParticle++) {
    diagFile << this->xInit[iParticle] << " " << this->yInit[iParticle] << " " << this->radiusInit[iParticle] << std::endl;
  } // iParticle
  diagFile.close();

  filename.str(std::string());
  filename << "part_init_halo_" << partition->this_proc() << "_" << partition->nprocs() << ".txt";
  diagFile.open(filename.str());
  for (int iParticle = this->nParticlesInit ; iParticle < this->nParticlesInit+this->nParticlesHaloInit ; iParticle++) {
    diagFile << this->xInit[iParticle] << " " << this->yInit[iParticle] << " " << this->radiusInit[iParticle] << std::endl;
  } // iParticle
  diagFile.close();

  filename.str(std::string());
  filename << "part_init_poly_" << partition->this_proc() << "_" << partition->nprocs() << ".txt";
  diagFile.open(filename.str());
  for (int iParticle = 0 ; iParticle < this->nParticlesInit ; iParticle++) {
    diagFile << "set object " << iParticle+1 << " polygon from " << this->initPolygons[iParticle].get_x(0) << "," << this->initPolygons[iParticle].get_y(0);
    for (int iVertex = 1 ; iVertex < this->initPolygons[iParticle].nvertices() ; iVertex++) {
      diagFile << " to " << this->initPolygons[iParticle].get_x(iVertex) << "," << this->initPolygons[iParticle].get_y(iVertex);
    }
    diagFile << std::endl;
    diagFile << "set object " << iParticle+1 << " fillstyle empty" << std::endl;
  } // iParticle
  diagFile.close();

  for (int iProc = 0 ; iProc < partition->nprocs() ; iProc++) {
    (*log)(DEMSI::LOG::DEBUG) << "nHaloSend: " << iProc << " " << localHaloParticlesIndicesExchInit[iProc].size() << std::endl;
    (*log)(DEMSI::LOG::DEBUG) << "nHaloRecv: " << iProc << " " << localOwnedParticlesIndicesExchInit[iProc].size() << std::endl;
  } // iProc

  filename.str(std::string());
  filename << "halo_" << partition->this_proc() << "_" << partition->nprocs() << ".txt";
  diagFile.open(filename.str());
  for (int iParticle = this->nParticlesInit ; iParticle < this->nParticlesAllInit ; iParticle++) {
    diagFile << iParticle << " " << this->globalIndexInit[iParticle] << std::endl;
  }
  diagFile << std::endl;
  for (int iProc = 0 ; iProc < partition->nprocs() ; iProc++) {
    for (int i = 0 ; i < localHaloParticlesIndicesExchInit[iProc].size() ; i++) {
      diagFile << iProc << " " << i << " " << localHaloParticlesIndicesExchInit[iProc][i] << std::endl;
    } // i
  } // iProc
  diagFile << std::endl;
  for (int iProc = 0 ; iProc < partition->nprocs() ; iProc++) {
    for (int i = 0 ; i < localOwnedParticlesIndicesExchInit[iProc].size() ; i++) {
      diagFile << iProc << " " << i << " " << localOwnedParticlesIndicesExchInit[iProc][i] << std::endl;
    } // i
  } // iProc
  diagFile.close();

}

// Set mapping from local to global indices for initial tessellation particles
void Tessellation::set_local_to_global_index_map_init(void) {

  char mpiErrBuffer[MPI_MAX_ERROR_STRING];
  int mpiErrLen;
  int mpiErr;

  // number of owned init particles on each processor
  int *nParticlesInitOnAllProcs = new int[partition->nprocs()];

  mpiErr = MPI_Allgather(&nParticlesInit, 1, MPI_INT, nParticlesInitOnAllProcs, 1, MPI_INT, partition->comm());
  MPI_Error_string(mpiErr, mpiErrBuffer, &mpiErrLen);
  log->check(mpiErr == MPI_SUCCESS, "set_local_to_global_index_map_init MPI_Allgather failed: ", (std::string) mpiErrBuffer);

  // cumulative number of owned init particles on previous processors
  int *nParticlesInitOnPrevProc = new int[partition->nprocs()];

  nParticlesInitOnPrevProc[0] = 0;
  for (int iProc = 1 ; iProc < partition->nprocs() ; iProc++) {
    nParticlesInitOnPrevProc[iProc] = nParticlesInitOnPrevProc[iProc-1] + nParticlesInitOnAllProcs[iProc-1];
  } // iProc

  int nParticlesTotal = 0;
  for (int iProc = 0 ; iProc < partition->nprocs() ; iProc++) {
    nParticlesTotal += nParticlesInitOnAllProcs[iProc];
  } // iProc

  (*log)(DEMSI::LOG::DEBUG) << "nParticlesTotal: " << nParticlesTotal << std::endl;
  for (int iProc = 0 ; iProc < partition->nprocs() ; iProc++) {
    (*log)(DEMSI::LOG::DEBUG) << "nParticlesInitOnAllProcs nParticlesInitOnPrevProc: " << nParticlesInitOnAllProcs[iProc] << " " << nParticlesInitOnPrevProc[iProc] << std::endl;
  } // iProc

  std::vector<int> localInitIndices(nParticlesAllInit);
  std::vector<int> localInitProcID (nParticlesAllInit);

  for (int iParticleInit = 0 ; iParticleInit < nParticlesInit ; iParticleInit++) {
    localInitIndices[iParticleInit] = iParticleInit;
    localInitProcID [iParticleInit] = partition->this_proc();
  } // iParticleInit

  update_init_halo_particles(localInitIndices);
  update_init_halo_particles(localInitProcID);

  localToGlobalIndicesInit.resize(nParticlesAllInit);
  for (int iParticleInit = 0 ; iParticleInit < nParticlesAllInit ; iParticleInit++) {
    localToGlobalIndicesInit[iParticleInit] = nParticlesInitOnPrevProc[localInitProcID[iParticleInit]] + localInitIndices[iParticleInit];
  } // iParticleInit

  delete [] nParticlesInitOnAllProcs;
  delete [] nParticlesInitOnPrevProc;

} // Tessellation::set_local_to_global_index_map_init

// Output particle array to test file for debugging
void Tessellation::output_init_particles(const std::string filenamePrefix, double* array) {

  std::stringstream filenameOut;
  filenameOut << std::setfill('0') << std::setw(3) << filenamePrefix << "_" << partition->nprocs() << "_" <<  partition->this_proc() << ".txt";
  std::ofstream file;
  file.open(filenameOut.str());
  file << std::scientific << std::setprecision(17);
  for (int iParticleInit = 0 ; iParticleInit < nParticlesInit ; iParticleInit++) {
    file << iParticleInit << " " << globalIndexInit[iParticleInit] << " " << array[iParticleInit] << std::endl;
  } // iParticleInit
  file.close();

} // Tessellation::output_init_particles

// Output edge array to test file for debugging
void Tessellation::output_init_edges(const std::string filenamePrefix, double* array) {

  std::stringstream filenameOut;
  filenameOut << std::setfill('0') << std::setw(3) << filenamePrefix << "_" << partition->nprocs() << "_" <<  partition->this_proc() << ".txt";
  std::ofstream file;
  file.open(filenameOut.str());
  file << std::scientific << std::setprecision(17);
  for (int iEdgeInit = 0 ; iEdgeInit < nEdgesInit ; iEdgeInit++) {
    int iParticleInit1 = cellsOnEdgeGlobal[iEdgeInit][0];
    int iParticleInit2 = cellsOnEdgeGlobal[iEdgeInit][1];
    if (iParticleInit1 != -1 and
	iParticleInit2 != -1) {
      int iParticleInit11 = std::min(iParticleInit1,iParticleInit2);
      int iParticleInit22 = std::max(iParticleInit1,iParticleInit2);
      file << iEdgeInit << " " << iParticleInit11 << " " << iParticleInit22 << " " << array[iEdgeInit] << std::endl;
    }
  } // iEdgeInit
  file.close();

} // Tessellation::output_init_particles


//------------------------------------------------------------------------------
// Halo update
//------------------------------------------------------------------------------

#define DEBUG_HALO
void Tessellation::update_init_halo_particles(std::vector<int> &array) {

  char mpiErrBuffer[MPI_MAX_ERROR_STRING];
  int mpiErrLen;
  int mpiErr;

  // transfer halo element accumulation
  for (int iProc=0 ; iProc < partition->nprocs()-1 ; iProc++) {

    // which processor we are sending data to
    int iProcSend = modulo((partition->this_proc() + iProc + 1), partition->nprocs());

    int *arraySend;
    MPI_Request request;

#ifdef DEBUG_HALO
    (*log)(DEMSI::LOG::DEBUG) << "update_init_halo_particles: localOwnedParticlesIndicesExchInit[iProcSend].size(): " <<
      iProcSend << " " << localOwnedParticlesIndicesExchInit[iProcSend].size() << std::endl;
#endif

    // check we have elements to send to the send processor
    if (localOwnedParticlesIndicesExchInit[iProcSend].size() > 0) {

      arraySend = new int[localOwnedParticlesIndicesExchInit[iProcSend].size()];

      for (int iParticleSend = 0 ; iParticleSend < localOwnedParticlesIndicesExchInit[iProcSend].size() ; iParticleSend++) {
	arraySend[iParticleSend] = array[localOwnedParticlesIndicesExchInit[iProcSend][iParticleSend]];
      } // iParticleSend

      mpiErr = MPI_Isend(&arraySend[0], localOwnedParticlesIndicesExchInit[iProcSend].size(), MPI_INT, iProcSend, 0, partition->comm(), &request);
      MPI_Error_string(mpiErr, mpiErrBuffer, &mpiErrLen);
      log->check(mpiErr == MPI_SUCCESS, "update_init_halo MPI_Isend failed: ", (std::string) mpiErrBuffer);

    }

    // which processor are we going to receive data from
    int iProcRecv = modulo((partition->this_proc() - iProc - 1), partition->nprocs());

#ifdef DEBUG_HALO
    (*log)(DEMSI::LOG::DEBUG) << "update_init_halo_particles: localHaloParticlesIndicesExchInit[iProcRecv].size(): " <<
      iProcRecv << " " << localHaloParticlesIndicesExchInit[iProcRecv].size() << std::endl;
#endif

    // check we have elements to receive from the send processor
    if (localHaloParticlesIndicesExchInit[iProcRecv].size() > 0) {

      int *arrayRecv = new int[localHaloParticlesIndicesExchInit[iProcRecv].size()];

      MPI_Status status;
      mpiErr = MPI_Recv(&arrayRecv[0], localHaloParticlesIndicesExchInit[iProcRecv].size(), MPI_INT, iProcRecv, 0, partition->comm(), &status);
      MPI_Error_string(mpiErr, mpiErrBuffer, &mpiErrLen);
      log->check(mpiErr == MPI_SUCCESS, "update_init_halo MPI_Recv failed: ", (std::string) mpiErrBuffer);

      for (int iParticleRecv = 0 ; iParticleRecv < localHaloParticlesIndicesExchInit[iProcRecv].size() ; iParticleRecv++) {
	array[localHaloParticlesIndicesExchInit[iProcRecv][iParticleRecv]] = arrayRecv[iParticleRecv];
      } // iParticleSend

      delete [] arrayRecv;

    }

    if (localOwnedParticlesIndicesExchInit[iProcSend].size() > 0) {

      MPI_Status status;
      mpiErr = MPI_Wait(&request, &status);
      MPI_Error_string(mpiErr, mpiErrBuffer, &mpiErrLen);
      log->check(mpiErr == MPI_SUCCESS, "update_init_halo MPI_Wait failed: ", (std::string) mpiErrBuffer);

      delete [] arraySend;

    }

  } // iProc

} // Tessellation::update_init_halo_particles

void Tessellation::update_init_halo_particles(std::vector<double> &array) {

  char mpiErrBuffer[MPI_MAX_ERROR_STRING];
  int mpiErrLen;
  int mpiErr;

  // transfer halo element accumulation
  for (int iProc=0 ; iProc < partition->nprocs()-1 ; iProc++) {

    // which processor we are sending data to
    int iProcSend = modulo((partition->this_proc() + iProc + 1), partition->nprocs());

    double *arraySend;
    MPI_Request request;

#ifdef DEBUG_HALO
    (*log)(DEMSI::LOG::DEBUG) << "update_init_halo_particles: localOwnedParticlesIndicesExchInit[iProcSend].size(): " <<
      iProcSend << " " << localOwnedParticlesIndicesExchInit[iProcSend].size() << std::endl;
#endif

    // check we have elements to send to the send processor
    if (localOwnedParticlesIndicesExchInit[iProcSend].size() > 0) {

      arraySend = new double[localOwnedParticlesIndicesExchInit[iProcSend].size()];

      for (int iParticleSend = 0 ; iParticleSend < localOwnedParticlesIndicesExchInit[iProcSend].size() ; iParticleSend++) {
	arraySend[iParticleSend] = array[localOwnedParticlesIndicesExchInit[iProcSend][iParticleSend]];
      } // iParticleSend

      mpiErr = MPI_Isend(&arraySend[0], localOwnedParticlesIndicesExchInit[iProcSend].size(), MPI_DOUBLE, iProcSend, 0, partition->comm(), &request);
      MPI_Error_string(mpiErr, mpiErrBuffer, &mpiErrLen);
      log->check(mpiErr == MPI_SUCCESS, "update_init_halo MPI_Isend failed: ", (std::string) mpiErrBuffer);

    }

    // which processor are we going to receive data from
    int iProcRecv = modulo((partition->this_proc() - iProc - 1), partition->nprocs());

#ifdef DEBUG_HALO
    (*log)(DEMSI::LOG::DEBUG) << "update_init_halo_particles: localHaloParticlesIndicesExchInit[iProcRecv].size(): " <<
      iProcRecv << " " << localHaloParticlesIndicesExchInit[iProcRecv].size() << std::endl;
#endif

    // check we have elements to receive from the send processor
    if (localHaloParticlesIndicesExchInit[iProcRecv].size() > 0) {

      double *arrayRecv = new double[localHaloParticlesIndicesExchInit[iProcRecv].size()];

      MPI_Status status;
      mpiErr = MPI_Recv(&arrayRecv[0], localHaloParticlesIndicesExchInit[iProcRecv].size(), MPI_DOUBLE, iProcRecv, 0, partition->comm(), &status);
      MPI_Error_string(mpiErr, mpiErrBuffer, &mpiErrLen);
      log->check(mpiErr == MPI_SUCCESS, "update_init_halo MPI_Recv failed: ", (std::string) mpiErrBuffer);

      for (int iParticleRecv = 0 ; iParticleRecv < localHaloParticlesIndicesExchInit[iProcRecv].size() ; iParticleRecv++) {
	array[localHaloParticlesIndicesExchInit[iProcRecv][iParticleRecv]] = arrayRecv[iParticleRecv];
      } // iParticleSend

      delete [] arrayRecv;

    }

    if (localOwnedParticlesIndicesExchInit[iProcSend].size() > 0) {

      MPI_Status status;
      mpiErr = MPI_Wait(&request, &status);
      MPI_Error_string(mpiErr, mpiErrBuffer, &mpiErrLen);
      log->check(mpiErr == MPI_SUCCESS, "update_init_halo MPI_Wait failed: ", (std::string) mpiErrBuffer);

      delete [] arraySend;

    }

  } // iProc

} // Tessellation::update_init_halo_particles

void Tessellation::update_init_halo_particles(int *array) {

  char mpiErrBuffer[MPI_MAX_ERROR_STRING];
  int mpiErrLen;
  int mpiErr;

  // transfer halo element accumulation
  for (int iProc=0 ; iProc < partition->nprocs()-1 ; iProc++) {

    // which processor we are sending data to
    int iProcSend = modulo((partition->this_proc() + iProc + 1), partition->nprocs());

    int *arraySend;
    MPI_Request request;

#ifdef DEBUG_HALO
    (*log)(DEMSI::LOG::DEBUG) << "update_init_halo_particles: localOwnedParticlesIndicesExchInit[iProcSend].size(): " <<
      iProcSend << " " << localOwnedParticlesIndicesExchInit[iProcSend].size() << std::endl;
#endif

    // check we have elements to send to the send processor
    if (localOwnedParticlesIndicesExchInit[iProcSend].size() > 0) {

      arraySend = new int[localOwnedParticlesIndicesExchInit[iProcSend].size()];

      for (int iParticleSend = 0 ; iParticleSend < localOwnedParticlesIndicesExchInit[iProcSend].size() ; iParticleSend++) {
	arraySend[iParticleSend] = array[localOwnedParticlesIndicesExchInit[iProcSend][iParticleSend]];
      } // iParticleSend

      mpiErr = MPI_Isend(&arraySend[0], localOwnedParticlesIndicesExchInit[iProcSend].size(), MPI_INT, iProcSend, 0, partition->comm(), &request);
      MPI_Error_string(mpiErr, mpiErrBuffer, &mpiErrLen);
      log->check(mpiErr == MPI_SUCCESS, "update_init_halo MPI_Isend failed: ", (std::string) mpiErrBuffer);

    }

    // which processor are we going to receive data from
    int iProcRecv = modulo((partition->this_proc() - iProc - 1), partition->nprocs());

#ifdef DEBUG_HALO
    (*log)(DEMSI::LOG::DEBUG) << "update_init_halo_particles: localHaloParticlesIndicesExchInit[iProcRecv].size(): " <<
      iProcRecv << " " << localHaloParticlesIndicesExchInit[iProcRecv].size() << std::endl;
#endif

    // check we have elements to receive from the send processor
    if (localHaloParticlesIndicesExchInit[iProcRecv].size() > 0) {

      int *arrayRecv = new int[localHaloParticlesIndicesExchInit[iProcRecv].size()];

      MPI_Status status;
      mpiErr = MPI_Recv(&arrayRecv[0], localHaloParticlesIndicesExchInit[iProcRecv].size(), MPI_INT, iProcRecv, 0, partition->comm(), &status);
      MPI_Error_string(mpiErr, mpiErrBuffer, &mpiErrLen);
      log->check(mpiErr == MPI_SUCCESS, "update_init_halo MPI_Recv failed: ", (std::string) mpiErrBuffer);

      for (int iParticleRecv = 0 ; iParticleRecv < localHaloParticlesIndicesExchInit[iProcRecv].size() ; iParticleRecv++) {
	array[localHaloParticlesIndicesExchInit[iProcRecv][iParticleRecv]] = arrayRecv[iParticleRecv];
      } // iParticleSend

      delete [] arrayRecv;

    }

    if (localOwnedParticlesIndicesExchInit[iProcSend].size() > 0) {

      MPI_Status status;
      mpiErr = MPI_Wait(&request, &status);
      MPI_Error_string(mpiErr, mpiErrBuffer, &mpiErrLen);
      log->check(mpiErr == MPI_SUCCESS, "update_init_halo MPI_Wait failed: ", (std::string) mpiErrBuffer);

      delete [] arraySend;

    }

  } // iProc

} // Tessellation::update_init_halo_particles

void Tessellation::update_init_halo_particles(double *array) {

  char mpiErrBuffer[MPI_MAX_ERROR_STRING];
  int mpiErrLen;
  int mpiErr;

  // transfer halo element accumulation
  for (int iProc=0 ; iProc < partition->nprocs()-1 ; iProc++) {

    // which processor we are sending data to
    int iProcSend = modulo((partition->this_proc() + iProc + 1), partition->nprocs());

    double *arraySend;
    MPI_Request request;

#ifdef DEBUG_HALO
    (*log)(DEMSI::LOG::DEBUG) << "update_init_halo_particles: localOwnedParticlesIndicesExchInit[iProcSend].size(): " <<
      iProcSend << " " << localOwnedParticlesIndicesExchInit[iProcSend].size() << std::endl;
#endif

    // check we have elements to send to the send processor
    if (localOwnedParticlesIndicesExchInit[iProcSend].size() > 0) {

      arraySend = new double[localOwnedParticlesIndicesExchInit[iProcSend].size()];

      for (int iParticleSend = 0 ; iParticleSend < localOwnedParticlesIndicesExchInit[iProcSend].size() ; iParticleSend++) {
	arraySend[iParticleSend] = array[localOwnedParticlesIndicesExchInit[iProcSend][iParticleSend]];
      } // iParticleSend

      mpiErr = MPI_Isend(&arraySend[0], localOwnedParticlesIndicesExchInit[iProcSend].size(), MPI_DOUBLE, iProcSend, 0, partition->comm(), &request);
      MPI_Error_string(mpiErr, mpiErrBuffer, &mpiErrLen);
      log->check(mpiErr == MPI_SUCCESS, "update_init_halo MPI_Isend failed: ", (std::string) mpiErrBuffer);

    }

    // which processor are we going to receive data from
    int iProcRecv = modulo((partition->this_proc() - iProc - 1), partition->nprocs());

#ifdef DEBUG_HALO
    (*log)(DEMSI::LOG::DEBUG) << "update_init_halo_particles: localHaloParticlesIndicesExchInit[iProcRecv].size(): " <<
      iProcRecv << " " << localHaloParticlesIndicesExchInit[iProcRecv].size() << std::endl;
#endif

    // check we have elements to receive from the send processor
    if (localHaloParticlesIndicesExchInit[iProcRecv].size() > 0) {

      double *arrayRecv = new double[localHaloParticlesIndicesExchInit[iProcRecv].size()];

      MPI_Status status;
      mpiErr = MPI_Recv(&arrayRecv[0], localHaloParticlesIndicesExchInit[iProcRecv].size(), MPI_DOUBLE, iProcRecv, 0, partition->comm(), &status);
      MPI_Error_string(mpiErr, mpiErrBuffer, &mpiErrLen);
      log->check(mpiErr == MPI_SUCCESS, "update_init_halo MPI_Recv failed: ", (std::string) mpiErrBuffer);

      for (int iParticleRecv = 0 ; iParticleRecv < localHaloParticlesIndicesExchInit[iProcRecv].size() ; iParticleRecv++) {
	array[localHaloParticlesIndicesExchInit[iProcRecv][iParticleRecv]] = arrayRecv[iParticleRecv];
      } // iParticleSend

      delete [] arrayRecv;

    }

    if (localOwnedParticlesIndicesExchInit[iProcSend].size() > 0) {

      MPI_Status status;
      mpiErr = MPI_Wait(&request, &status);
      MPI_Error_string(mpiErr, mpiErrBuffer, &mpiErrLen);
      log->check(mpiErr == MPI_SUCCESS, "update_init_halo MPI_Wait failed: ", (std::string) mpiErrBuffer);

      delete [] arraySend;

    }

  } // iProc

} // Tessellation::update_init_halo_particles

void Tessellation::update_init_halo_edges(std::vector<int> &array) {

  char mpiErrBuffer[MPI_MAX_ERROR_STRING];
  int mpiErrLen;
  int mpiErr;

  // transfer halo element accumulation
  for (int iProc=0 ; iProc < partition->nprocs()-1 ; iProc++) {

    // which processor we are sending data to
    int iProcSend = modulo((partition->this_proc() + iProc + 1), partition->nprocs());

    int *arraySend;
    MPI_Request request;

#ifdef DEBUG_HALO
    (*log)(DEMSI::LOG::DEBUG) << "update_init_halo_edges: localOwnedEdgesIndicesExchInit[iProcSend].size(): " <<
      iProcSend << " " << localOwnedEdgesIndicesExchInit[iProcSend].size() << std::endl;
#endif

    // check we have elements to send to the send processor
    if (localOwnedEdgesIndicesExchInit[iProcSend].size() > 0) {

      arraySend = new int[localOwnedEdgesIndicesExchInit[iProcSend].size()];

      for (int iParticleSend = 0 ; iParticleSend < localOwnedEdgesIndicesExchInit[iProcSend].size() ; iParticleSend++) {
	arraySend[iParticleSend] = array[localOwnedEdgesIndicesExchInit[iProcSend][iParticleSend]];
      } // iParticleSend

      mpiErr = MPI_Isend(&arraySend[0], localOwnedEdgesIndicesExchInit[iProcSend].size(), MPI_INT, iProcSend, 0, partition->comm(), &request);
      MPI_Error_string(mpiErr, mpiErrBuffer, &mpiErrLen);
      log->check(mpiErr == MPI_SUCCESS, "update_init_halo MPI_Isend failed: ", (std::string) mpiErrBuffer);

    }

    // which processor are we going to receive data from
    int iProcRecv = modulo((partition->this_proc() - iProc - 1), partition->nprocs());

#ifdef DEBUG_HALO
    (*log)(DEMSI::LOG::DEBUG) << "update_init_halo_edges: localHaloEdgesIndicesExchInit[iProcRecv].size(): " <<
      iProcRecv << " " << localHaloEdgesIndicesExchInit[iProcRecv].size() << std::endl;
#endif

    // check we have elements to receive from the send processor
    if (localHaloEdgesIndicesExchInit[iProcRecv].size() > 0) {

      int *arrayRecv = new int[localHaloEdgesIndicesExchInit[iProcRecv].size()];

      MPI_Status status;
      mpiErr = MPI_Recv(&arrayRecv[0], localHaloEdgesIndicesExchInit[iProcRecv].size(), MPI_INT, iProcRecv, 0, partition->comm(), &status);
      MPI_Error_string(mpiErr, mpiErrBuffer, &mpiErrLen);
      log->check(mpiErr == MPI_SUCCESS, "update_init_halo MPI_Recv failed: ", (std::string) mpiErrBuffer);

      for (int iParticleRecv = 0 ; iParticleRecv < localHaloEdgesIndicesExchInit[iProcRecv].size() ; iParticleRecv++) {
	array[localHaloEdgesIndicesExchInit[iProcRecv][iParticleRecv]] = arrayRecv[iParticleRecv];
      } // iParticleSend

      delete [] arrayRecv;

    }

    if (localOwnedEdgesIndicesExchInit[iProcSend].size() > 0) {

      MPI_Status status;
      mpiErr = MPI_Wait(&request, &status);
      MPI_Error_string(mpiErr, mpiErrBuffer, &mpiErrLen);
      log->check(mpiErr == MPI_SUCCESS, "update_init_halo MPI_Wait failed: ", (std::string) mpiErrBuffer);

      delete [] arraySend;

    }

  } // iProc

} // Tessellation::update_init_halo_edges

void Tessellation::update_init_halo_edges(std::vector<double> &array) {

  char mpiErrBuffer[MPI_MAX_ERROR_STRING];
  int mpiErrLen;
  int mpiErr;

  // transfer halo element accumulation
  for (int iProc=0 ; iProc < partition->nprocs()-1 ; iProc++) {

    // which processor we are sending data to
    int iProcSend = modulo((partition->this_proc() + iProc + 1), partition->nprocs());

    double *arraySend;
    MPI_Request request;

#ifdef DEBUG_HALO
    (*log)(DEMSI::LOG::DEBUG) << "update_init_halo_edges: localOwnedEdgesIndicesExchInit[iProcSend].size(): " <<
      iProcSend << " " << localOwnedEdgesIndicesExchInit[iProcSend].size() << std::endl;
#endif

    // check we have elements to send to the send processor
    if (localOwnedEdgesIndicesExchInit[iProcSend].size() > 0) {

      arraySend = new double[localOwnedEdgesIndicesExchInit[iProcSend].size()];

      for (int iParticleSend = 0 ; iParticleSend < localOwnedEdgesIndicesExchInit[iProcSend].size() ; iParticleSend++) {
	arraySend[iParticleSend] = array[localOwnedEdgesIndicesExchInit[iProcSend][iParticleSend]];
      } // iParticleSend

      mpiErr = MPI_Isend(&arraySend[0], localOwnedEdgesIndicesExchInit[iProcSend].size(), MPI_DOUBLE, iProcSend, 0, partition->comm(), &request);
      MPI_Error_string(mpiErr, mpiErrBuffer, &mpiErrLen);
      log->check(mpiErr == MPI_SUCCESS, "update_init_halo MPI_Isend failed: ", (std::string) mpiErrBuffer);

    }

    // which processor are we going to receive data from
    int iProcRecv = modulo((partition->this_proc() - iProc - 1), partition->nprocs());

#ifdef DEBUG_HALO
    (*log)(DEMSI::LOG::DEBUG) << "update_init_halo_edges: localHaloEdgesIndicesExchInit[iProcRecv].size(): " <<
      iProcRecv << " " << localHaloEdgesIndicesExchInit[iProcRecv].size() << std::endl;
#endif

    // check we have elements to receive from the send processor
    if (localHaloEdgesIndicesExchInit[iProcRecv].size() > 0) {

      double *arrayRecv = new double[localHaloEdgesIndicesExchInit[iProcRecv].size()];

      MPI_Status status;
      mpiErr = MPI_Recv(&arrayRecv[0], localHaloEdgesIndicesExchInit[iProcRecv].size(), MPI_DOUBLE, iProcRecv, 0, partition->comm(), &status);
      MPI_Error_string(mpiErr, mpiErrBuffer, &mpiErrLen);
      log->check(mpiErr == MPI_SUCCESS, "update_init_halo MPI_Recv failed: ", (std::string) mpiErrBuffer);

      for (int iParticleRecv = 0 ; iParticleRecv < localHaloEdgesIndicesExchInit[iProcRecv].size() ; iParticleRecv++) {
	array[localHaloEdgesIndicesExchInit[iProcRecv][iParticleRecv]] = arrayRecv[iParticleRecv];
      } // iParticleSend

      delete [] arrayRecv;

    }

    if (localOwnedEdgesIndicesExchInit[iProcSend].size() > 0) {

      MPI_Status status;
      mpiErr = MPI_Wait(&request, &status);
      MPI_Error_string(mpiErr, mpiErrBuffer, &mpiErrLen);
      log->check(mpiErr == MPI_SUCCESS, "update_init_halo MPI_Wait failed: ", (std::string) mpiErrBuffer);

      delete [] arraySend;

    }

  } // iProc

} // Tessellation::update_init_halo_edges

void Tessellation::update_init_halo_edges(int *array) {

  char mpiErrBuffer[MPI_MAX_ERROR_STRING];
  int mpiErrLen;
  int mpiErr;

  // transfer halo element accumulation
  for (int iProc=0 ; iProc < partition->nprocs()-1 ; iProc++) {

    // which processor we are sending data to
    int iProcSend = modulo((partition->this_proc() + iProc + 1), partition->nprocs());

    int *arraySend;
    MPI_Request request;

#ifdef DEBUG_HALO
    (*log)(DEMSI::LOG::DEBUG) << "update_init_halo_edges: localOwnedEdgesIndicesExchInit[iProcSend].size(): " <<
      iProcSend << " " << localOwnedEdgesIndicesExchInit[iProcSend].size() << std::endl;
#endif

    // check we have elements to send to the send processor
    if (localOwnedEdgesIndicesExchInit[iProcSend].size() > 0) {

      arraySend = new int[localOwnedEdgesIndicesExchInit[iProcSend].size()];

      for (int iParticleSend = 0 ; iParticleSend < localOwnedEdgesIndicesExchInit[iProcSend].size() ; iParticleSend++) {
	arraySend[iParticleSend] = array[localOwnedEdgesIndicesExchInit[iProcSend][iParticleSend]];
      } // iParticleSend

      mpiErr = MPI_Isend(&arraySend[0], localOwnedEdgesIndicesExchInit[iProcSend].size(), MPI_INT, iProcSend, 0, partition->comm(), &request);
      MPI_Error_string(mpiErr, mpiErrBuffer, &mpiErrLen);
      log->check(mpiErr == MPI_SUCCESS, "update_init_halo MPI_Isend failed: ", (std::string) mpiErrBuffer);

    }

    // which processor are we going to receive data from
    int iProcRecv = modulo((partition->this_proc() - iProc - 1), partition->nprocs());

#ifdef DEBUG_HALO
    (*log)(DEMSI::LOG::DEBUG) << "update_init_halo_edges: localHaloEdgesIndicesExchInit[iProcRecv].size(): " <<
      iProcRecv << " " << localHaloEdgesIndicesExchInit[iProcRecv].size() << std::endl;
#endif

    // check we have elements to receive from the send processor
    if (localHaloEdgesIndicesExchInit[iProcRecv].size() > 0) {

      int *arrayRecv = new int[localHaloEdgesIndicesExchInit[iProcRecv].size()];

      MPI_Status status;
      mpiErr = MPI_Recv(&arrayRecv[0], localHaloEdgesIndicesExchInit[iProcRecv].size(), MPI_INT, iProcRecv, 0, partition->comm(), &status);
      MPI_Error_string(mpiErr, mpiErrBuffer, &mpiErrLen);
      log->check(mpiErr == MPI_SUCCESS, "update_init_halo MPI_Recv failed: ", (std::string) mpiErrBuffer);

      for (int iParticleRecv = 0 ; iParticleRecv < localHaloEdgesIndicesExchInit[iProcRecv].size() ; iParticleRecv++) {
	array[localHaloEdgesIndicesExchInit[iProcRecv][iParticleRecv]] = arrayRecv[iParticleRecv];
      } // iParticleSend

      delete [] arrayRecv;

    }

    if (localOwnedEdgesIndicesExchInit[iProcSend].size() > 0) {

      MPI_Status status;
      mpiErr = MPI_Wait(&request, &status);
      MPI_Error_string(mpiErr, mpiErrBuffer, &mpiErrLen);
      log->check(mpiErr == MPI_SUCCESS, "update_init_halo MPI_Wait failed: ", (std::string) mpiErrBuffer);

      delete [] arraySend;

    }

  } // iProc

} // Tessellation::update_init_halo_edges

void Tessellation::update_init_halo_edges(double *array) {

  char mpiErrBuffer[MPI_MAX_ERROR_STRING];
  int mpiErrLen;
  int mpiErr;

  // transfer halo element accumulation
  for (int iProc=0 ; iProc < partition->nprocs()-1 ; iProc++) {

    // which processor we are sending data to
    int iProcSend = modulo((partition->this_proc() + iProc + 1), partition->nprocs());

    double *arraySend;
    MPI_Request request;

#ifdef DEBUG_HALO
    (*log)(DEMSI::LOG::DEBUG) << "update_init_halo_edges: localOwnedEdgesIndicesExchInit[iProcSend].size(): " <<
      iProcSend << " " << localOwnedEdgesIndicesExchInit[iProcSend].size() << std::endl;
#endif

    // check we have elements to send to the send processor
    if (localOwnedEdgesIndicesExchInit[iProcSend].size() > 0) {

      arraySend = new double[localOwnedEdgesIndicesExchInit[iProcSend].size()];

      for (int iParticleSend = 0 ; iParticleSend < localOwnedEdgesIndicesExchInit[iProcSend].size() ; iParticleSend++) {
	arraySend[iParticleSend] = array[localOwnedEdgesIndicesExchInit[iProcSend][iParticleSend]];
      } // iParticleSend

      mpiErr = MPI_Isend(&arraySend[0], localOwnedEdgesIndicesExchInit[iProcSend].size(), MPI_DOUBLE, iProcSend, 0, partition->comm(), &request);
      MPI_Error_string(mpiErr, mpiErrBuffer, &mpiErrLen);
      log->check(mpiErr == MPI_SUCCESS, "update_init_halo MPI_Isend failed: ", (std::string) mpiErrBuffer);

    }

    // which processor are we going to receive data from
    int iProcRecv = modulo((partition->this_proc() - iProc - 1), partition->nprocs());

#ifdef DEBUG_HALO
    (*log)(DEMSI::LOG::DEBUG) << "update_init_halo_edges: localHaloEdgesIndicesExchInit[iProcRecv].size(): " <<
      iProcRecv << " " << localHaloEdgesIndicesExchInit[iProcRecv].size() << std::endl;
#endif

    // check we have elements to receive from the send processor
    if (localHaloEdgesIndicesExchInit[iProcRecv].size() > 0) {

      double *arrayRecv = new double[localHaloEdgesIndicesExchInit[iProcRecv].size()];

      MPI_Status status;
      mpiErr = MPI_Recv(&arrayRecv[0], localHaloEdgesIndicesExchInit[iProcRecv].size(), MPI_DOUBLE, iProcRecv, 0, partition->comm(), &status);
      MPI_Error_string(mpiErr, mpiErrBuffer, &mpiErrLen);
      log->check(mpiErr == MPI_SUCCESS, "update_init_halo MPI_Recv failed: ", (std::string) mpiErrBuffer);

      for (int iParticleRecv = 0 ; iParticleRecv < localHaloEdgesIndicesExchInit[iProcRecv].size() ; iParticleRecv++) {
	array[localHaloEdgesIndicesExchInit[iProcRecv][iParticleRecv]] = arrayRecv[iParticleRecv];
      } // iParticleSend

      delete [] arrayRecv;

    }

    if (localOwnedEdgesIndicesExchInit[iProcSend].size() > 0) {

      MPI_Status status;
      mpiErr = MPI_Wait(&request, &status);
      MPI_Error_string(mpiErr, mpiErrBuffer, &mpiErrLen);
      log->check(mpiErr == MPI_SUCCESS, "update_init_halo MPI_Wait failed: ", (std::string) mpiErrBuffer);

      delete [] arraySend;

    }

  } // iProc

} // Tessellation::update_init_halo_edges

//----------------------------------------------------------------------------
// IO
//----------------------------------------------------------------------------

// Register fields for output
void Tessellation::register_output_field(double* field, const std::string fieldName) {

  fieldsWrite.push_back(std::pair<double*,std::string>(field, fieldName));

} // Tessellation::register_output_field

} // DEMSI namespace
