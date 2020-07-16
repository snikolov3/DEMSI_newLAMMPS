#ifndef DEMSI_INTERPOLATION_H_
#define DEMSI_INTERPOLATION_H_

/*! \file   demsi_interpolation.h
    \brief  Header file for the DEMSI::Interpolation class.
*/

#include "demsi_logging.h"
#include "demsi_lmp_instance.h"
#include "demsi_particles.h"
#include "demsi_partition.h"
#include "demsi_grid.h"
#include "demsi_configs.h"
#include "demsi_column.h"
#include "demsi_column_variables.h"

#include  <Kokkos_Core.hpp>


namespace DEMSI {
/*! \class Interpolation
 \brief Class for interpolating data between particles and structured grid.

 This class includes linear mapping routines from particle to structured grid
 and structured grid to particle location.
*/
class Interpolation{
public:

  /*! \brief Interpolation class Constructor
      \param partitionIn   [in] Pointer to partition object.
      \param logIn         [in] Pointer to log object.
      \param configsIn     [in] Pointer to configs object.
      \param gridIn        [in] Pointer to grid object.
      \param particlesIn   [in] Pointer to particles object.
      \param columnIn      [in] Pointer to column object.
   */
  Interpolation(DEMSI::Partition* partitionIn, DEMSI::Log* logIn, DEMSI::Configs* configsIn,
                DEMSI::Grid* gridIn, DEMSI::Particles* particlesIn, DEMSI::Column* columnIn);

  /*! \brief Default destructor */
  ~Interpolation() = default;

  /*! \brief Boolean indicating whether there are selected output fields that need interpolation. */
  bool interpToGrid;

  /*! \brief Interpolate vector of particle column variable fields to grid.
             Fields to interpolate are determined by grid fields selected for
             output in xml input file.
   */
  void particle_fields_to_grid();

  /*! \brief Interpolate scalar particle values in Kokkos view to grid
      \param partField            [in] - Kokkos view of particle field 
      \param gridField           [out] - Kokkos view of interpolated grid field
   */
  void particles_to_grid(Kokkos::View<double *> &partField, Kokkos::View<double **> &gridField);

  /*! \brief Interpolate vector particle values in Kokkos view to grid
      \param partField            [in] - Kokkos view of particle field 
      \param gridField           [out] - Kokkos view of interpolated grid field
   */
  void particles_to_grid(Kokkos::View<double **> &partField, Kokkos::View<double ***> &gridField);

  /*! \brief Interpolate column variable to grid
      \param partField            [in] - Pointer to column variable
      \param gridField           [out] - Kokkos view of interpolated grid field
   */
  void particles_to_grid(DEMSI::ColumnVariable<double>* partField, Kokkos::View<double **> &gridField);

  // Interpolated grid fields
  /*! Interpolated ice thickness */
  Kokkos::View<double**> iceThicknessGrid;
  /*! Interpolated ice concentration */
  Kokkos::View<double**> iceConcentrationGrid;

private:

  /*! \brief Interpolate particle column variable to grid field. 
      \param momentMat      [in] - moment matrix for least squares interpolation
      \param momentMatInv   [in] - first row of inverse of moment matrix
      \param partField      [in] - Kokkos view of particle field
      \param gridField     [out] - Kokkos view of interpolated grid field
      \param gridFieldLO   [out] - Kokkos view of first-order interpolated grid field
  */
   void interpolate_fields(Kokkos::View<double***> & momentMat, Kokkos::View<double***> & momentMatInv, 
                           DEMSI::ColumnVariable<double>* partField, Kokkos::View<double**> gridField,
                           Kokkos::View<double**> gridFieldLO);

  /*! \brief Compute moment matrix for least squares interpolation 

      Compute moment matrix before interpolating particle fields to grid
      whenever particles have moved since the last interpolation time.
  */
  void get_moment_matrix(Kokkos::View<double***> & momentMat);

  /*! \brief Compute inverse of moment matrix for least squares interpolation 
  */
  void get_moment_matrix_inv(Kokkos::View<double***> & momentMat,
                             Kokkos::View<double***> & momentMatInv);

  /*! \brief Communicate moment matrix edge values to neighbors and sum
  */
  void exchange_moment_matrix(Kokkos::View<double***> & momentMat);

  /*! \brief Communicate grid field edge values to neighbors and sum
  */
  void exchange_grid_field(Kokkos::View<double**> & gridField);

  /*! \brief Communicate grid field edge values to neighbors and take minimum
  */
  void exchange_grid_field_min(Kokkos::View<double**> & gridField);

  /*! \brief Communicate grid field edge values to neighbors and take maximum
  */
  void exchange_grid_field_max(Kokkos::View<double**> & gridField);

  /*! \brief Limit grid field values based on local min/max from particles
      \param gridField  [in/out] - Kokkos view of interpolated grid field
      \param gridFieldLO    [in] - Kokkos view of first-order interpolated grid field
      \param minField       [in] - Kokkos view of local minima on grid
      \param maxField       [in] - Kokkos view of local maxima on grid
  */
   void limit_grid_field(Kokkos::View<double**>  gridField,Kokkos::View<double**>  gridFieldLO,
                         Kokkos::View<double**> & minField, Kokkos::View<double**> & maxField);

  /*! \brief Compute local bounds on grid from particle field values
      \param partField      [in] - Kokkos view of particle field
      \param minField       [in] - Kokkos view of local minima on grid
      \param maxField       [in] - Kokkos view of local maxima on grid
  */
   void get_field_bounds(DEMSI::ColumnVariable<double>* partField,
                         Kokkos::View<double**> & minField, Kokkos::View<double**> & maxField);

  /*! \brief Locate particle on structured grid 
      \param x              [in] - x coordinate of particle
      \param y              [in] - y coordinate of particle
      \param grid_dx        [in] - grid resolution in x direction
      \param grid_dy        [in] - grid resolution in y direction
      \param iglobal       [out] - global x index of cell
      \param jglobal       [out] - global y index of cell
      \param xi            [out] - reference coordinate in cell
      \param eta           [out] - reference coordinate in cell
   */
  KOKKOS_INLINE_FUNCTION
  void locate_particle(double x, double y, double grid_dx, double grid_dy, 
                       int &iglobal, int &jglobal, double &xi, double &eta);

  /*! \brief Compute basis function values at particle location
      \param xi             [in] - reference coordinate in cell
      \param eta            [in] - reference coordinate in cell
      \param basis         [out] - basis function values indexed by cell node
   */
  KOKKOS_INLINE_FUNCTION
  void grid_basis(double xi, double eta, double basis[4]);


  /*! Vector of interpolated grid fields */
  //std::vector<Kokkos::View<double**>*> gridFieldsInterp;
  std::vector<Kokkos::View<double**>> gridFieldsInterp;
 
  /*! Vector of particles fields for interpolation */
  std::vector<DEMSI::ColumnVariable<double>*> partFieldsInterp;

  /*! Pointer to the particles object */
  DEMSI::Particles* particles;

  /*! Pointer to the partitions object */
  DEMSI::Partition* partition;

  /*! Pointer to the grid object */
  DEMSI::Grid* grid;

  /*! Pointer to the log object */
  DEMSI::Log* log;

  /*! Pointer to the configs object */
  DEMSI::Configs* configs;

  /*! Pointer to the column object */
  DEMSI::Column* column;

  /*! Grid partitions */
  std::vector<DEMSI::GridPartition*> gridPartitions;

};

} // namespace DEMSI

#endif /* INTERPOLATION_H_ */
