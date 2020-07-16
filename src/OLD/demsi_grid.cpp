#include "demsi_grid.h"
#include "demsi_logging.h"
#include "demsi_file_utils.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <netcdf.h>
#include <vector>
#include <iomanip>

#include <lammps/atom.h>
#include <lammps/math_const.h>
#include <lammps/domain.h>
#include <lammps/update.h>
#include <lammps/neighbor.h>
#include <lammps/comm.h>

namespace DEMSI {

  // Default constructor for the GridRegion class.
  GridRegion::GridRegion() {};

  // Constructor for the GridRegion class.
  GridRegion::GridRegion(const int i0In, const int j0In, const int i1In, const int j1In) {
    i0 = i0In;
    j0 = j0In;
    i1 = i1In;
    j1 = j1In;
  };

  // Return the lower left x index of the grid region.
  int GridRegion::get_i0(void) const {
    return i0;
  }

  // Return the lower left y index of the grid region.
  int GridRegion::get_j0(void) const {
    return j0;
  }

  // Return the upper right x index of the grid region.
  int GridRegion::get_i1(void) const {
    return i1;
  }

  // Return the upper right y index of the grid region.
  int GridRegion::get_j1(void) const {
    return j1;
  }

  // Set the limits of the grid region.
  void GridRegion::set(const int i0In, const int j0In, const int i1In, const int j1In) {
    i0 = i0In;
    j0 = j0In;
    i1 = i1In;
    j1 = j1In;
  }

  // Return the size of the grid region in one dimension.
  int GridRegion::size(int dim) {
    if (dim == 0) {
      return std::max(i1 - i0 + 1, 0);
    } else {
      return std::max(j1 - j0 + 1, 0);
    }
  }

  // Return the total number of cells in a grid region.
  int GridRegion::size(void) {
    return std::max(i1 - i0 + 1, 0) * std::max(j1 - j0 + 1, 0);
  }

  // Constructor for the GridPartition class.
  GridPartition::GridPartition(const double x0, const double y0, const double x1, const double y1, const double skin, const double dx, const double dy, const int nx, const int ny, Log* logIn) {

    log = logIn;

    int i0Owned = owned_lower(x0 / dx);
    int j0Owned = owned_lower(y0 / dy);
    int i1Owned = owned_upper(x1 / dx, nx);
    int j1Owned = owned_upper(y1 / dy, ny);

    int i0Total = total_lower((x0 - skin*1.1) / dx);
    int j0Total = total_lower((y0 - skin*1.1) / dy);
    int i1Total = total_upper((x1 + skin*1.1) / dx, nx);
    int j1Total = total_upper((y1 + skin*1.1) / dy, ny);

    if (i1Owned == nx - 2) i1Owned = nx - 1;
    if (j1Owned == ny - 2) j1Owned = ny - 1;

    owned = new DEMSI::GridRegion(i0Owned, j0Owned, i1Owned, j1Owned);
    total = new DEMSI::GridRegion(i0Total, j0Total, i1Total, j1Total);

  }

  // Destructor for the GridPartition class.
  GridPartition::~GridPartition() {
    delete owned;
    delete total;
  }

  // Return lower index (x or y) of a owned region given the lower domain limit of the lammps partition.
  int GridPartition::owned_lower(const double x) {
    if (std::floor(x) == std::ceil(x)) {
      return std::max((int) std::floor(x),0);
    } else {
      return std::max((int) std::ceil(x),0);
    }
  }

  // Return upper index (x or y) of a owned region given the upper domain limit of the lammps partition.
  int GridPartition::owned_upper(const double x, const int n) {
    if (std::floor(x) == std::ceil(x)) {
      return std::min((int) std::floor(x) - 1, n-1);
    } else {
      return std::min((int) std::floor(x), n-1);
    }
  }

  // Return lower index (x or y) of a total region given the lower domain limit of the lammps partition.
  int GridPartition::total_lower(const double x) {
    if (std::floor(x) == std::ceil(x)) {
      return std::max((int) std::floor(x),0);
    } else {
      return std::max((int) std::floor(x),0);
    }
  }

  // Return upper index (x or y) of a total region given the upper domain limit of the lammps partition.
  int GridPartition::total_upper(const double x, const int n) {
    if (std::floor(x) == std::ceil(x)) {
      return std::min((int) std::floor(x), n-1);;
    } else {
      return std::min((int) std::ceil(x), n-1);;
    }
  }

  // Return the total grid region of a grid partition.
  GridRegion* GridPartition::get_total(void) {
    GridRegion* totalOut = new DEMSI::GridRegion();
    (*totalOut) = (*total);
    return totalOut;
  }

  // Return the owned grid region of a grid partition.
  GridRegion* GridPartition::get_owned(void) {
    GridRegion* ownedOut = new DEMSI::GridRegion();
    (*ownedOut) = (*total);
    return ownedOut;
  }

  // Return the size in one dimension of the total grid region.
  int GridPartition::size_total(int dim) const {
    if (dim == 0) {
      return total->size(0);
    } else {
      return total->size(1);
    }
  }

  // Return the size in one dimension of the owned grid region.
  int GridPartition::size_owned(int dim) const {
    if (dim == 0) {
      return owned->size(0);
    } else {
      return owned->size(1);
    }
  }

  // Return the total number of cells in the total region.
  int GridPartition::size_total(void) const {
    return total->size(0) * total->size(1);
  }

  // Return the total number of cells in the owned region.
  int GridPartition::size_owned(void) const {
    return owned->size(0) * owned->size(1);
  }

  // Given a global x index in the grid get the local x index in the owned region.
  int GridPartition::ilocal_owned(int iGlobal) const {
    return iGlobal - owned->get_i0();
  }

  // Given a global y index in the grid get the local y index in the owned region.
  int GridPartition::jlocal_owned(int jGlobal) const {
    return jGlobal - owned->get_j0();
  }

  // Given a local x index in the owned region get the global x index in the grid.
  int GridPartition::iglobal_owned(int iLocal) const {
    return iLocal + owned->get_i0();
  }

  // Given a local y index in the owned region get the global y index in the grid.
  int GridPartition::jglobal_owned(int jLocal) const {
    return jLocal + owned->get_j0();
  }

  // Given a global x index in the grid get the local x index in the total region.
  int GridPartition::ilocal_total(int iGlobal) const {
    return iGlobal - total->get_i0();
  }

  // Given a global y index in the grid get the local y index in the total region.
  int GridPartition::jlocal_total(int jGlobal) const {
    return jGlobal - total->get_j0();
  }

  // Given a local x index in the total region get the global x index in the grid.
  int GridPartition::iglobal_total(int iLocal) const {
    return iLocal + total->get_i0();
  }

  // Given a local y index in the total region get the global y index in the grid.
  int GridPartition::jglobal_total(int jLocal) const {
    return jLocal + total->get_j0();
  }

  // Constructor for Grid class.
  Grid::Grid(const std::string gridFilenameIn, Partition* partitionIn, Log* logIn, LammpsInstance* lammpsInstanceIn, Clock* simulationClockIn) {

    gridFilename = gridFilenameIn;
    partition = partitionIn;
    log = logIn;
    lammpsInstance = lammpsInstanceIn;
    simulationClock = simulationClockIn;

    // create views
    latitude  = Kokkos::View<double **>("latitude",0,0);
    longitude = Kokkos::View<double **>("longitude",0,0);

    // register fields for reading
    fieldsGridRead.push_back(&latitude);
    fieldsGridRead.push_back(&longitude);

    // register fields for writing
    fieldsGridWrite.push_back(&latitude);
    fieldsGridWrite.push_back(&longitude);

    // Read grid and initialize geometric information and forcing
    read_grid_size(gridFilename);

  }

  // register gridded field for possible write
  void Grid::register_field_for_write(Kokkos::View<double**>* fieldPtr) {
    fieldsGridWrite.push_back(fieldPtr);
  }

  // Return the number of cells in the x direction of the partition (total region).
  int Grid::nx(void) const {
    return dimensions[0];
  }

  // Return the number of cells in the y direction of the partition (total region).
  int Grid::ny(void) const {
    return dimensions[1];
  }

  // Return the size of the global domain in the x direction.
  double Grid::lx(void) const {
    return ((double) (globalDimensions[0] - 1)) * resolution[0];
  }

  // Return the size of the global domain in the y direction.
  double Grid::ly(void) const {
    return ((double) (globalDimensions[1] - 1)) * resolution[1];
  }

  // Return the grid resolution in the x direction
  int Grid::dx(void) const {
     return resolution[0];
  }

  // Return the grid resolution in the y direction
  int Grid::dy(void) const {
     return resolution[1];
  }

  // Return the grid partition
  std::vector<GridPartition*> Grid::getGridPartitions(void) {
     return gridPartitions;
  }

  // Read in the size and resolution of the domain.
  void Grid::read_grid_size(const std::string filenameIn) {

    int err;
    char mpiErrBuffer[MPI_MAX_ERROR_STRING];
    int mpiErrLen;

    int nc_id;

    if (partition->on_master_proc()) {

      // open
      err = nc_open(filenameIn.c_str(), NC_NOWRITE, &nc_id);
      log->check(err == NC_NOERR, "Couldn't open grid file: ", filenameIn, ", :", nc_strerror(err));

      // dimensions
      int nx_id;
      err = nc_inq_dimid(nc_id, "nx", &nx_id);
      log->check(err == NC_NOERR, "Problem getting dimid for nx: ", nc_strerror(err));

      int ny_id;
      err = nc_inq_dimid(nc_id, "ny", &ny_id);
      log->check(err == NC_NOERR, "Problem getting dimid for ny: ", nc_strerror(err));

      size_t tmp_len;
      err = nc_inq_dimlen(nc_id, nx_id, &tmp_len);
      log->check(err == NC_NOERR, "Problem getting dimension nx: ", nc_strerror(err));
      globalDimensions[0] = static_cast<int>(tmp_len);

      err = nc_inq_dimlen(nc_id, ny_id, &tmp_len);
      log->check(err == NC_NOERR, "Problem getting dimension ny: ", nc_strerror(err));
      globalDimensions[1] = static_cast<int>(tmp_len);

      // resolution
      int resolution_id;
      err = nc_inq_varid(nc_id, "resolution", &resolution_id);
      log->check(err == NC_NOERR, "Problem getting resolution varid: ", nc_strerror(err));

      size_t start[1] = {0};
      size_t count[1] = {2};
      err = nc_get_vara_double(nc_id, resolution_id, start, count, &resolution[0]);
      log->check(err == NC_NOERR, "Problem getting resolution: ", nc_strerror(err));

      // close the file
      err = nc_close(nc_id);
      log->check(err == NC_NOERR);

    } // master proc

    err = MPI_Bcast(&globalDimensions[0], 2, MPI_INT, partition->master_proc(), partition->comm());
    MPI_Error_string(err, mpiErrBuffer, &mpiErrLen);
    log->check(err == MPI_SUCCESS, "Problem broadcasting global dimensions: ", (std::string) mpiErrBuffer);

    err = MPI_Bcast(&resolution[0], 2, MPI_DOUBLE, partition->master_proc(), partition->comm());
    MPI_Error_string(err, mpiErrBuffer, &mpiErrLen);
    log->check(err == MPI_SUCCESS, "Problem broadcasting resolution: ", (std::string) mpiErrBuffer);

  }

  // Read in the grid fields from the grid netcdf input file.
  void Grid::read_grid(void) {

    int err;
    int nc_id;

    if (partition->on_master_proc()) {

      // open
      err = nc_open(gridFilename.c_str(), NC_NOWRITE, &nc_id);
      log->check(err == NC_NOERR, "Could not open grid file: ", gridFilename, ", :", nc_strerror(err));

    }

    get_grid_partitions();

    check_grid_partitions();

    dimensions[0] = gridPartitions[partition->this_proc()]->size_total(0);
    dimensions[1] = gridPartitions[partition->this_proc()]->size_total(1);

    ownedDimensions[0] = gridPartitions[partition->this_proc()]->size_owned(0);
    ownedDimensions[1] = gridPartitions[partition->this_proc()]->size_owned(1);

    for (int iField = 0; iField < fieldsGridRead.size() ; iField++) {

      std::string fieldName = fieldsGridRead[iField]->label();
      Grid::read_field_fixed(fieldsGridRead[iField], nc_id, fieldName);

    } // field loop

    if (partition->on_master_proc()) {

      err = nc_close(nc_id);
      log->check(err == NC_NOERR, "Problem closing file: ", gridFilename);

    } // master proc

  }

  // Read a Eulerian field with no time index from a netcdf file.
  void Grid::read_field_fixed(Kokkos::View<double**>* field, const int nc_id, const std::string varname) {

    size_t start[2] = {0, 0};
    size_t count[2] = {(size_t) globalDimensions[0], (size_t) globalDimensions[1]};

    read_field(field, nc_id, varname, start, count);

  }

  // Read a Eulerian field with a time index from a netcdf file.
  void Grid::read_field_varying(Kokkos::View<double**>* field, const int nc_id, const std::string varname, const int timeIndex) {

    size_t start[3] = {(size_t) timeIndex, 0, 0};
    size_t count[3] = {1, (size_t) globalDimensions[0], (size_t) globalDimensions[1]};

    read_field(field, nc_id, varname, start, count);

  }

  // Read a single Eulerian field from a netcdf file.
  void Grid::read_field(Kokkos::View<double**>* field, const int nc_id, const std::string varname, const size_t* start, const size_t* count) {

    int err;
    char mpiErrBuffer[MPI_MAX_ERROR_STRING];
    int mpiErrLen;

    int* sendCount;
    int* displacements;
    double* sendBuffer;

    if (partition->on_master_proc()) {

      sendCount = new int [partition->nprocs()];
      displacements = new int [partition->nprocs()];

      int sendCountTotal = 0;
      for (int iProc = 0 ; iProc<partition->nprocs() ; iProc++) {
	sendCount[iProc] = gridPartitions[iProc]->size_total();
	sendCountTotal = sendCountTotal + sendCount[iProc];
      }
      displacements[0] = 0;
      for (int iProc = 1 ; iProc<partition->nprocs() ; iProc++) {
	displacements[iProc] = displacements[iProc-1] + sendCount[iProc-1];
      }

      sendBuffer = new double [sendCountTotal];

      int varid;

      err = nc_inq_varid(nc_id, varname.c_str(), &varid);
      log->check(err == NC_NOERR, "Problem getting varid for field: ", varname, ": ", nc_strerror(err));

      double* arrayIn = new double [globalDimensions[0]*globalDimensions[1]];

      err = nc_get_vara_double(nc_id, varid, start, count, &arrayIn[0]);
      log->check(err == NC_NOERR, "Problem getting field: ", nc_strerror(err));

      int ij = 0;
      for (int iProc = 0 ; iProc < partition->nprocs() ; iProc++) {

	for (int i=0 ; i<gridPartitions[iProc]->size_total(0) ; i++) {
	  for (int j=0 ; j<gridPartitions[iProc]->size_total(1) ; j++) {

	    sendBuffer[ij] = arrayIn[gridPartitions[iProc]->iglobal_total(i) * globalDimensions[1] + gridPartitions[iProc]->jglobal_total(j)];
	    ij++;

	  } // j
	} // i

      } // iProc

      delete [] arrayIn;

    } // master proc

    int recvCount = dimensions[0] * dimensions[1];
    double* recvBuffer = new double [dimensions[0]*dimensions[1]];

    err = MPI_Scatterv(sendBuffer, sendCount, displacements, MPI_DOUBLE, recvBuffer, recvCount, MPI_DOUBLE, partition->master_proc(), partition->comm());
    MPI_Error_string(err, mpiErrBuffer, &mpiErrLen);
    log->check(err == MPI_SUCCESS, "Problem scattering field: ", (std::string) mpiErrBuffer);

    // resize device view for field
    Kokkos::resize((*field), dimensions[0], dimensions[1]);

    // read field into mirror on host
    Kokkos::View<double**>::HostMirror field_mirror = Kokkos::create_mirror(*field);
    Kokkos::parallel_for("Grid:read_field 2d", Kokkos::RangePolicy<host_execution_space>
        (0,dimensions[0]), [=] (const int i) {

      for (int j=0 ; j<dimensions[1] ; j++) {
        const int ij = i*dimensions[1] + j;
        field_mirror(i,j) = recvBuffer[ij];
      }

    });
    // copy field back to device
    Kokkos::deep_copy(*field, field_mirror);

    if (partition->on_master_proc()) {
      delete [] sendCount;
      delete [] displacements;
      delete [] sendBuffer;
    }
    delete [] recvBuffer;

  }

  // Determine the grid partitions on all processors.
  void Grid::get_grid_partitions() {

    // empty current vector
    for (int i = 0; i < gridPartitions.size(); ++i) {
      delete gridPartitions[i];
    }
    gridPartitions.clear();

    int err;
    char mpiErrBuffer[MPI_MAX_ERROR_STRING];
    int mpiErrLen;

    double x0Local = lammpsInstance->lmp->domain->sublo[0];
    double y0Local = lammpsInstance->lmp->domain->sublo[1];
    double x1Local = lammpsInstance->lmp->domain->subhi[0];
    double y1Local = lammpsInstance->lmp->domain->subhi[1];

    double x0[partition->nprocs()];
    err = MPI_Allgather(&x0Local, 1, MPI_DOUBLE, x0, 1, MPI_DOUBLE, partition->comm());
    MPI_Error_string(err, mpiErrBuffer, &mpiErrLen);
    log->check(err == MPI_SUCCESS, "All gather x0: ", (std::string) mpiErrBuffer);

    double y0[partition->nprocs()];
    err = MPI_Allgather(&y0Local, 1, MPI_DOUBLE, y0, 1, MPI_DOUBLE, partition->comm());
    MPI_Error_string(err, mpiErrBuffer, &mpiErrLen);
    log->check(err == MPI_SUCCESS, "All gather y0: ", (std::string) mpiErrBuffer);

    double x1[partition->nprocs()];
    err = MPI_Allgather(&x1Local, 1, MPI_DOUBLE, x1, 1, MPI_DOUBLE, partition->comm());
    MPI_Error_string(err, mpiErrBuffer, &mpiErrLen);
    log->check(err == MPI_SUCCESS, "All gather x1: ", (std::string) mpiErrBuffer);

    double y1[partition->nprocs()];
    err = MPI_Allgather(&y1Local, 1, MPI_DOUBLE, y1, 1, MPI_DOUBLE, partition->comm());
    MPI_Error_string(err, mpiErrBuffer, &mpiErrLen);
    log->check(err == MPI_SUCCESS, "All gather y1: ", (std::string) mpiErrBuffer);

    double skin = lammpsInstance->lmp->neighbor->skin;

    for (int iProc=0 ; iProc<partition->nprocs() ; iProc++) {

      GridPartition *gridPartition = new DEMSI::GridPartition(x0[iProc], y0[iProc], x1[iProc], y1[iProc], skin, resolution[0], resolution[1], globalDimensions[0], globalDimensions[1], log);
      gridPartitions.push_back(gridPartition);

    }

  }

  // Check that every grid cell is owned by exactly one processor.
  void Grid::check_grid_partitions() {

    // check owned partition first
    int* cellsInPartition = new int [globalDimensions[0] * globalDimensions[1]];

    for (int i = 0 ; i < globalDimensions[0] ; i++) {
      for (int j = 0 ; j < globalDimensions[1] ; j++) {
	cellsInPartition[i*globalDimensions[1]+j] = 0;
      }
    }

    for (int iPartition = 0 ; iPartition < gridPartitions.size() ; iPartition++) {

      for (int i = 0 ; i < gridPartitions[iPartition]->size_owned(0) ; i++) {
	for (int j = 0 ; j < gridPartitions[iPartition]->size_owned(1) ; j++) {

	  int iGlobal = gridPartitions[iPartition]->iglobal_owned(i);
	  int jGlobal = gridPartitions[iPartition]->jglobal_owned(j);

	  int ijGlobal = iGlobal * globalDimensions[1] + jGlobal;

	  cellsInPartition[ijGlobal] = cellsInPartition[ijGlobal] + 1;

	}
      }

    }

    bool lGridComplete = true;
    for (int ij = 0 ; ij < globalDimensions[0] * globalDimensions[1] ; ij++) {
      if (cellsInPartition[ij] != 1) {
	(*log)() << ij << " " << cellsInPartition[ij] << std::endl;
	lGridComplete = false;
      }
    }

    MPI_Barrier(partition->comm());
    log->check(lGridComplete, "Grid partition check failed.");

    delete [] cellsInPartition;

    // now check total partition for this processor
    GridRegion *total = gridPartitions[partition->this_proc()]->get_total();

    double skin = lammpsInstance->lmp->neighbor->skin;

    double x0Local = lammpsInstance->lmp->domain->sublo[0];
    double y0Local = lammpsInstance->lmp->domain->sublo[1];
    double x1Local = lammpsInstance->lmp->domain->subhi[0];
    double y1Local = lammpsInstance->lmp->domain->subhi[1];

    log->check(total->get_i0() == 0 or
	       (total->get_i0()      * resolution[0] <= x0Local - skin*1.1 and
	       (total->get_i0() + 1) * resolution[0] >  x0Local - skin*1.1),
		"Total grid partition incorrect: i0:",
		std::to_string(total->get_i0()), std::to_string(resolution[0]), std::to_string(x0Local), std::to_string(skin*1.1));

    log->check(total->get_j0() == 0 or
	       (total->get_j0()      * resolution[1] <= y0Local - skin*1.1 and
	       (total->get_j0() + 1) * resolution[1] >  y0Local - skin*1.1),
		"Total grid partition incorrect: j0:",
		std::to_string(total->get_j0()), std::to_string(resolution[1]), std::to_string(y0Local), std::to_string(skin*1.1));

    log->check(total->get_i1() == globalDimensions[0]-1 or
	       (total->get_i1()      * resolution[0] >= x1Local + skin*1.1 and
	       (total->get_i1() - 1) * resolution[0] <  x1Local + skin*1.1),
		"Total grid partition incorrect: i1:",
		std::to_string(total->get_i1()), std::to_string(resolution[0]), std::to_string(x1Local), std::to_string(skin*1.1));

    log->check(total->get_j1() == globalDimensions[1]-1 or
	       (total->get_j1()      * resolution[1] >= y1Local + skin*1.1 and
	       (total->get_j1() - 1) * resolution[1] <  y1Local + skin*1.1),
		"Total grid partition incorrect: j1:",
		std::to_string(total->get_j1()), std::to_string(resolution[1]), std::to_string(y1Local), std::to_string(skin*1.1));


    // check that owned is subset of total
    GridRegion *owned = gridPartitions[partition->this_proc()]->get_owned();

    log->check(owned->get_i0() >= total->get_i0(),
		"Grid partition: owned not subset of total for i0: Owned: ",
		std::to_string(owned->get_i0()), ", total: ", std::to_string(total->get_i0()));

    log->check(owned->get_j0() >= total->get_j0(),
		"Grid partition: owned not subset of total for j0: Owned: ",
		std::to_string(owned->get_j0()), ", total: ", std::to_string(total->get_j0()));

    log->check(owned->get_i1() <= total->get_i1(),
		"Grid partition: owned not subset of total for i1: Owned: ",
		std::to_string(owned->get_i1()), ", total: ", std::to_string(total->get_i1()));

    log->check(owned->get_j1() <= total->get_j1(),
		"Grid partition: owned not subset of total for j1: Owned: ",
		std::to_string(owned->get_j1()), ", total: ", std::to_string(total->get_j1()));

    // clean-up
    delete total;
    delete owned;

  }

  // For an arbitary position (x,y) determine the bilinear interpolation weights from the Eulerian mesh.
  void Grid::interpolation_weights(const double x, const double y, int (&i)[4], int (&j)[4], double (&weight)[4]) {

    i[0] = gridPartitions[partition->this_proc()]->ilocal_total(std::floor(x / resolution[0]));
    i[1] = gridPartitions[partition->this_proc()]->ilocal_total(std::ceil(x / resolution[0]));
    i[2] = i[1];
    i[3] = i[0];

    j[0] = gridPartitions[partition->this_proc()]->jlocal_total(std::floor(y / resolution[1]));
    j[2] = gridPartitions[partition->this_proc()]->jlocal_total(std::ceil(y / resolution[1]));
    j[1] = j[0];
    j[3] = j[2];

    double integral;

    double wx1 = std::modf(x / resolution[0], &integral);
    double wx0 = 1.0 - wx1;

    double wy1 = std::modf(y / resolution[1], &integral);
    double wy0 = 1.0 - wy1;

    weight[0] = wx0 * wy0;
    weight[1] = wx0 * wy1;
    weight[2] = wx1 * wy1;
    weight[3] = wx1 * wy0;

  }

  // Batch call for bilinear interpolation weights for arbitary positions (x,y) from the Eulerian mesh.
  void Grid::interpolation_weights(kokkos_view_type_x particles_x, Kokkos::View<int*[4]> i, Kokkos::View<int*[4]> j, Kokkos::View<double*[4]> weight) {

    Kokkos::View<double[2]> resolution_device("resolution device");
    Kokkos::View<double[2], Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> > resolution_host(&resolution[0]);
    Kokkos::deep_copy(resolution_device, resolution_host);

    auto total_get_i0 = -gridPartitions[partition->this_proc()]->ilocal_total(0);
    auto total_get_j0 = -gridPartitions[partition->this_proc()]->jlocal_total(0);

    Kokkos::parallel_for("Grid::interpolation_weights", Kokkos::RangePolicy<device_execution_space>(0, particles_x.dimension_0()), KOKKOS_LAMBDA(const int iParticle) {

      double x = particles_x(iParticle,0);
      double y = particles_x(iParticle,1);

      // need to get what this ilocal_total and jlocal_total are doing and 
      // then just do the calculation on the vector
      i(iParticle, 0) = std::floor(x / resolution_device[0]) - total_get_i0;
      i(iParticle, 1) = std::ceil(x / resolution_device[0]) - total_get_i0;
      i(iParticle, 2) = i(iParticle, 1);
      i(iParticle, 3) = i(iParticle, 0);

      j(iParticle, 0) = std::floor(y / resolution_device[1]) - total_get_j0;
      j(iParticle, 2) = std::ceil(y / resolution_device[1]) - total_get_j0;
      j(iParticle, 1) = j(iParticle, 0);
      j(iParticle, 3) = j(iParticle, 2);

      double integral;

      double wx1 = std::modf(x / resolution_device[0], &integral);
      double wx0 = 1.0 - wx1;

      double wy1 = std::modf(y / resolution_device[1], &integral);
      double wy0 = 1.0 - wy1;

      weight(iParticle, 0) = wx0 * wy0;
      weight(iParticle, 1) = wx0 * wy1;
      weight(iParticle, 2) = wx1 * wy1;
      weight(iParticle, 3) = wx1 * wy0;
    });
  }

  // check if a point is in the local domain
  bool Grid::in_domain(const double x, const double y) {

    double x0Local = lammpsInstance->lmp->domain->sublo[0];
    double y0Local = lammpsInstance->lmp->domain->sublo[1];
    double x1Local = lammpsInstance->lmp->domain->subhi[0];
    double y1Local = lammpsInstance->lmp->domain->subhi[1];

    if (x > x0Local and x <= x1Local and y > y0Local and y <= y1Local) {
      return true;
    } else {
      return false;
    }

  }

  // Check if a point is in the local halo region.
  bool Grid::in_halo_region(const double x, const double y, const double haloWidth) {

    bool haloLeft   = (x >= xlow() - haloWidth and x <  xlow()              and y >= ylow() - haloWidth and y <= yhigh() + haloWidth);
    bool haloRight  = (x >= xhigh()            and x <= xhigh() + haloWidth and y >= ylow() - haloWidth and y <= yhigh() + haloWidth);
    bool haloBottom = (y >= ylow() - haloWidth and y <  ylow()              and x >= xlow() - haloWidth and x <= xhigh() + haloWidth);
    bool haloTop    = (y >= yhigh()            and y <= yhigh() + haloWidth and x >= xlow() - haloWidth and x <= xhigh() + haloWidth);

    return (haloLeft or haloRight or haloBottom or haloTop);

  }

  // x position of lower x boundary of local domain
  double Grid::xlow(void) {
    return lammpsInstance->lmp->domain->sublo[0];
  }

  // x position of upper x boundary of local domain
  double Grid::xhigh(void) {
    return lammpsInstance->lmp->domain->subhi[0];
  }

  // y position of lower y boundary of local domain
  double Grid::ylow(void) {
    return lammpsInstance->lmp->domain->sublo[1];
  }

  // y position of upper y boundary of local domain
  double Grid::yhigh(void) {
    return lammpsInstance->lmp->domain->subhi[1];
  }

  void Grid::check_in_domain(const double x, const double y) {

    double x0Local = lammpsInstance->lmp->domain->sublo[0];
    double y0Local = lammpsInstance->lmp->domain->sublo[1];
    double x1Local = lammpsInstance->lmp->domain->subhi[0];
    double y1Local = lammpsInstance->lmp->domain->subhi[1];

    double skin = lammpsInstance->lmp->neighbor->skin;

    if (x < x0Local-skin*1.1 or x > x1Local+skin*1.1 or y < y0Local-skin*1.1 or y > y1Local+skin*1.1) {

      (*log)() << std::scientific << std::setprecision(10);
      (*log)() << "skin: " << skin*1.1 << std::endl;
      (*log)() << "x: " << x << " " << x0Local << " " << x1Local << std::endl;
      (*log)() << "y: " << y << " " << y0Local << " " << y1Local << std::endl;
      log->abort("Outside domain!");
    }

  }

  // Return the neighbor proc to the west of this grid region
  int Grid::neighbor_west(void) const {
    return lammpsInstance->lmp->comm->procneigh[0][0];
  }

  // Return the neighbor proc to the east of this grid region
  int Grid::neighbor_east(void) const {
    return lammpsInstance->lmp->comm->procneigh[0][1];
  }

  // Return the neighbor proc to the south of this grid region
  int Grid::neighbor_south(void) const {
    return lammpsInstance->lmp->comm->procneigh[1][0];
  }

  // Return the neighbor proc to the north of this grid region
  int Grid::neighbor_north(void) const {
    return lammpsInstance->lmp->comm->procneigh[1][1];
  }

} // namespace DEMSI
