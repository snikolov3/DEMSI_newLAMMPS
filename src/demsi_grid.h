#ifndef DEMSI_GRID_H_
#define DEMSI_GRID_H_

/*! \file   demsi_grid.h
    \brief  Header file for the DEMSI::Grid class and related classes.
*/

#include "demsi_typedefs.h"
#include "demsi_logging.h"
#include "demsi_time.h"
#include "demsi_partition.h"
#include "demsi_lmp_instance.h"
#include "demsi_configs.h"

#include <Kokkos_Core.hpp>
#include <vector>


namespace DEMSI {

/*! \class GridRegion
 \brief Class defining a region of the grid.

 This class defines a rectangular region of the grid, starting from some lower
 left grid cell with global coordinates (x0,y0) and extending to some upper
 right grid cell at global coordinates (x1,y1).

*/
class GridRegion {

public:

  /*! \brief Default constructor for the GridRegion class. */
  GridRegion();

  /*! \brief Constructor for the GridRegion class.
      \param i0In x index of the lower left corner of the grid region.
      \param j0In y index of the lower left corner of the grid region.
      \param i1In x index of the upper right corner of the grid region.
      \param j1In y index of the upper right corner of the grid region.
   */
  GridRegion(const int i0In, const int j0In, const int i1In, const int j1In);

  /*! \brief Default destructor for the GridRegion class. */
  ~GridRegion() = default;

  /*! \brief Return the lower left x index of the grid region.
      \return The lower left x index of the grid region.
  */
  int get_i0(void) const;

  /*! \brief Return the lower left y index of the grid region.
      \return The lower left y index of the grid region.
  */
  int get_j0(void) const;

  /*! \brief Return the upper right x index of the grid region.
      \return The upper right x index of the grid region.
  */
  int get_i1(void) const;

  /*! \brief Return the upper right y index of the grid region.
      \return The upper right y index of the grid region.
  */
  int get_j1(void) const;

  /*! \brief Set the limits of the grid region.
      \param i0In The lower left x index of the grid region.
      \param j0In The lower left y index of the grid region.
      \param i1In The upper right x index of the grid region.
      \param j1In The upper right y index of the grid region.
  */
  void set(const int i0In, const int j0In, const int i1In, const int j1In);

  /*! \brief Return the size of the grid region in one dimension.
      \param dim The dimension (0 or 1) to return the grid size of.
      \return The size of the grid region in one dimension.
  */
  int size(int dim);

  /*! \brief Return the total number of cells in a grid region.
      \return The total number of cells in a grid region.
  */
  int size(void);

private:

  /*! Indices of lower left grid cell. */
  int i0, j0;

  /*! Indices of the upper right grid cell. */
  int i1, j1;

};

/*! \class GridPartition
 \brief This class defines a grid partition on a particular processor.

 This class defines a grid partition on a particular processor. It consists of
 two grid regions, an owned region and a total region. The total region
 defines all cells needed on a processor to cover the lammps domain partition,
 whereas the owned region contains cells uniquely owned by a processor. Every
 cell is in exactly one owned region.
*/
class GridPartition {

public:

  /*! \brief Constructor for the GridPartition class.
      \param x0 Lower x limit of the local lammps domain partition.
      \param y0 Lower y limit of the local lammps domain partition.
      \param x1 Upper x limit of the local lammps domain partition.
      \param y1 Upper y limit of the local lammps domain partition.
      \param skin Distance LAMMPS particles can leave LAMMPS domain.
      \param dx Grid resolution in x direction.
      \param dy Grid resolution in y direction.
      \param nx Total grid size in x direction.
      \param ny Total grid size in y direction.
      \param logIn Pointer to the log object.
  */
  GridPartition(const double x0, const double y0, const double x1, const double y1, const double skin, const double dx, const double dy, const int nx, const int ny, Log* logIn);

  /*! \brief Destructor for the GridPartition class. */
  ~GridPartition();

  /*! \brief Return lower index (x or y) of a owned region given the lower domain limit of the lammps partition.
      \param x The lower domain limit of the lammps partition.
      \return Lower index (x or y) of a owned region.
  */
  int owned_lower(const double x);

  /*! \brief Return upper index (x or y) of a owned region given the upper domain limit of the lammps partition.
      \param x The upper domain limit of the lammps partition.
      \param n The upper size limit of the domain in desired direction.
      \return Upper index (x or y) of a owned region.
  */
  int owned_upper(const double x, const int n);

  /*! \brief Return lower index (x or y) of a total region given the lower domain limit of the lammps partition.
      \param x The lower domain limit of the lammps partition.
      \return Lower index (x or y) of a total region.
  */
  int total_lower(const double x);

  /*! \brief Return upper index (x or y) of a total region given the upper domain limit of the lammps partition.
      \param x The upper domain limit of the lammps partition.
      \param n The upper size limit of the domain in desired direction.
      \return Upper index (x or y) of a total region.
  */
  int total_upper(const double x, const int n);

  /*! \brief Return the total grid region of a grid partition.
      \return The total grid region of a grid partition.
  */
  GridRegion* get_total(void);

  /*! \brief Return the owned grid region of a grid partition.
      \return The owned grid region of a grid partition.
  */
  GridRegion* get_owned(void);

  /*! \brief Return the size in one dimension of the total grid region.
      \param dim The dimension to get the size of.
      \return The size in one dimension of the total grid region.
  */
  int size_total(int dim) const;

  /*! \brief Return the size in one dimension of the owned grid region.
      \param dim The dimension to get the size of.
      \return The size in one dimension of the owned grid region.
  */
  int size_owned(int dim) const;

  /*! \brief Return the total number of cells in the total region.
      \return The total number of cells in the total region.
  */
  int size_total(void) const;

  /*! \brief Return the total number of cells in the owned region.
      \return The total number of cells in the owned region.
  */
  int size_owned(void) const;

  /*! \brief Given a global x index in the grid get the local x index in the owned region.
      \param iGlobal Global x index in the grid.
      \return The local x index in the owned region.
  */
  int ilocal_owned(int iGlobal) const;

  /*! \brief Given a global y index in the grid get the local y index in the owned region.
      \param jGlobal Global y index in the grid.
      \return The local y index in the owned region.
  */
  int jlocal_owned(int jGlobal) const;

  /*! \brief Given a local x index in the owned region get the global x index in the grid.
      \param iLocal Local x index in the owned region.
      \return Global x index in the grid.
  */
  int iglobal_owned(int iLocal) const;

  /*! \brief Given a local y index in the owned region get the global y index in the grid.
      \param jLocal Local y index in the owned region.
      \return Global y index in the grid.
  */
  int jglobal_owned(int jLocal) const;

  /*! \brief Given a global x index in the grid get the local x index in the total region.
      \param iGlobal Global x index in the grid.
      \return The local x index in the total region.
  */
  int ilocal_total(int iGlobal) const;

  /*! \brief Given a global y index in the grid get the local y index in the total region.
      \param jGlobal Global y index in the grid.
      \return The local y index in the total region.
  */
  int jlocal_total(int jGlobal) const;

  /*! \brief Given a local x index in the total region get the global x index in the grid.
      \param iLocal Local x index in the total region.
      \return Global x index in the grid.
  */
  int iglobal_total(int iLocal) const;

  /*! \brief Given a local y index in the total region get the global y index in the grid.
      \param jLocal Local y index in the total region.
      \return Global y index in the grid.
  */
  int jglobal_total(int jLocal) const;

private:

  /*! Owned grid region of the grid partition. */
  GridRegion* owned;

  /*! Total grid region of the grid partition. */
  GridRegion* total;

  /*! Pointer to the log object. */
  Log* log;

};

/*! \class Grid
 \brief Class defining fixed grid used for applying external forcing.

  The Grid class initializes a fixed grid for each hemisphere, and
  reads in forcing data from netcdf grid files.

*/
class GridWrite;
class Grid {
  friend class DEMSI::GridWrite;
public:

  /*! Constructor for Grid class.
    \param gridFilenameIn Name of the input grid file.
    \param partitionIn Pointer to the partition object.
    \param logIn Pointer to the log object
    \param lammpsInstanceIn Pointer to the lammps instance object.
    \param simulationClockIn Pointer to the simulation clock.
   */
  Grid(const std::string gridFilenameIn, Partition* partitionIn, Log* logIn, LammpsInstance* lammpsInstanceIn, Clock* simulationClockIn);

  /*! Default destructor */
  ~Grid() = default;

  /*! Latitude of grid points */
  Kokkos::View<double**> latitude;

  /*! Longitude of grid points */
  Kokkos::View<double**> longitude;

  /*! \brief Register gridded field for possible write
      \param fieldPtr Pointer to the field to register for possible output
  */
  void register_field_for_write(Kokkos::View<double**>* fieldPtr);

  /*! \brief Return the number of cells in the x direction of the partition (total region).
      \return The number of cells in the x direction of the partition (total region).
  */
  int nx(void) const;

  /*! \brief Return the number of cells in the y direction of the partition (total region).
      \return The number of cells in the y direction of the partition (total region).
  */
  int ny(void) const;

  /*! \brief Return the size of the global domain in the x direction.
      \return The size of the global domain in the x direction.
  */
  double lx(void) const;

  /*! \brief Return the size of the global domain in the y direction.
      \return The size of the global domain in the y direction.
  */
  double ly(void) const;

  /*! \brief Return the grid resolution in the x direction.
      \return The grid resolution in the x direction.
  */
  int dx(void) const;

  /*! \brief Return the grid resolution in the y direction.
      \return The grid resolution in the y direction.
  */
  int dy(void) const;

  /*! \brief Return the grid partitions.
      \return The grid partition object.
  */
  std::vector<GridPartition*> getGridPartitions(void);

  /*! \brief Read in the grid fields from the grid netcdf input file.
  */
  void read_grid(void);

  /*! \brief Read a Eulerian field with no time index from a netcdf file.
      \param field Pointer to the field to be read in.
      \param nc_id Netcdf file ID of the file to read from.
      \param varname Variable name in the netcdf file to read in.
  */
  void read_field_fixed(Kokkos::View<double**>* field, const int nc_id, const std::string varname);

  /*! \brief Read a Eulerian field with a time index from a netcdf file.
      \param field Pointer to the field to be read in.
      \param nc_id Netcdf file ID of the file to read from.
      \param varname Variable name in the netcdf file to read in.
      \param timeIndex Time index to read in.
  */
  void read_field_varying(Kokkos::View<double**>* field, const int nc_id, const std::string varname, const int timeIndex);

  /*! \brief For an arbitary position (x,y) determine the bilinear interpolation weights from the Eulerian mesh.
      \param x x coordinate of the position to interpolate to.
      \param y y coordinate of the position to interpolate to.
      \param i Array of x local indices on the Eulerian grid to interpolate from.
      \param j Array of y local indices on the Eulerian grid to interpolate from.
      \param weight Array of weights for the interpolation points.
  */
  void interpolation_weights(const double x, const double y, int (&i)[4], int (&j)[4], double (&weight)[4]);

  /*! \brief Batch call for bilinear interpolation weights for arbitary positions (x,y) from the Eulerian mesh.
      \param x x coordinates of the position to interpolate to.
      \param y y coordinates of the position to interpolate to.
      \param i Array of x local indices on the Eulerian grid to interpolate from.
      \param j Array of y local indices on the Eulerian grid to interpolate from.
      \param weight Array of weights for the interpolation points.
  */
  void interpolation_weights(kokkos_view_type_x particles_x, Kokkos::View<int*[4]> i, Kokkos::View<int*[4]> j, Kokkos::View<double*[4]> weight);

  /*! \brief check if a point is in the local domain
      \param x x position of point to test
      \param y y position of point to test
      \return is point in domain
  */
  bool in_domain(const double x, const double y);

  /*! \brief check if a point is in the local halo region.
      \param x x position of point to test
      \param y y position of point to test
      \param haloWidth width of the halo region.
      \return is point in halo region
  */
  bool in_halo_region(const double x, const double y, const double haloWidth);

  /*! \brief Return x position of lower x boundary of local domain
      \return x position of lower x boundary of local domain
  */
  double xlow(void);

  /*! \brief Return x position of upper x boundary of local domain
      \return x position of upper x boundary of local domain
  */
  double xhigh(void);

  /*! \brief Return y position of lower y boundary of local domain
      \return y position of lower y boundary of local domain
  */
  double ylow(void);

  /*! \brief Return y position of upper y boundary of local domain
      \return y position of upper y boundary of local domain
  */
  double yhigh(void);

  /*! \brief check a position is within the local total domain.
      \param x x position to check.
      \param y y position to check.
  */
  void check_in_domain(const double x, const double y);

  /*! \brief Get processor to the west of current region.
      \return Processor id of region to the west.
  */
  int neighbor_west(void) const;

  /*! \brief Get processor to the east of current region.
      \return Processor id of region to the east.
  */
  int neighbor_east(void) const;
 
  /*! \brief Get processor to the south of current region.
      \return Processor id of region to the south.
  */
  int neighbor_south(void) const;

  /*! \brief Get processor to the north of current region.
      \return Processor id of region to the north.
  */
  int neighbor_north(void) const;

private:

  /*! Grid dimensions */
  int dimensions[2];

  /*! Grid dimensions */
  int globalDimensions[2];

  /*! Grid dimensions */
  int ownedDimensions[2];

  /*! Grid resolution (dx,dy) */
  double resolution[2];

  /*! Filename of the input grid file. */
  std::string gridFilename;

  /*! Pointer to the logging object. */
  Log* log;

  /*! Pointer to the lammps instance object. */
  LammpsInstance* lammpsInstance;

  /*! Pointer to the partition object. */
  Partition* partition;

  /*! Pointer to the configs object. */
  Configs* configs;

  /*! Vector of grid fields to read from the grid input file. */
  std::vector<Kokkos::View<double**>*> fieldsGridRead;

  /*! Vector of fields registered for possible output */
  std::vector<Kokkos::View<double**>*> fieldsGridWrite;

  /*! Vector of grid partitions for all processors. */
  std::vector<GridPartition*> gridPartitions;

  /*! Pointer to the simulation clock */
  Clock* simulationClock;

  /*! \brief Read in the size and resolution of the domain.
      \param filenameIn Netcdf file with grid data
   */
  void read_grid_size(const std::string filenameIn);

  /*! \brief Read a single Eulerian field from a netcdf file.
      \param field Pointer to the field to read.
      \param nc_id Netcdf file ID to read from.
      \param varname Netcdf variable name to read.
      \param start Array of starting indices for read.
      \param count Array of sizes for read.
  */
  void read_field(Kokkos::View<double**>* field, const int nc_id, const std::string varname, const size_t* start, const size_t* count);

  /*! \brief Determine the grid partitions on all processors.
  */
  void get_grid_partitions(void);

  /*! \brief Check that every grid cell is owned by exactly one processor.
  */
  void check_grid_partitions(void);

};

} // namespace DEMSI

#endif /* GRID_H_ */
