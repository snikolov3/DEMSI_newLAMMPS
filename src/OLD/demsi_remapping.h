#ifndef DEMSI_REMAPPING_H
#define DEMSI_REMAPPING_H

#include "demsi_configs.h"
#include "demsi_particles.h"
#include "demsi_initial_tessellation.h"
#include "demsi_initial_tessellation_io.h"
#include "demsi_partition.h"
#include "demsi_logging.h"
#include "demsi_grid.h"
#include "demsi_column.h"
#include "demsi_column_variables.h"
#include "demsi_ocean.h"

#include <petscvec.h>
#include <petscmat.h>

#define DEBUG_REMAPPING

namespace DEMSI {

//------------------------------------------------------------------------------
/*! \class RemappingTracer
    \brief

*/
//------------------------------------------------------------------------------
class RemappingTracer {

public:

  /*! \brief Constructor for RemappingVariable class */
  RemappingTracer(DEMSI::Tessellation* tessellationIn, DEMSI::ColumnVariable<double>* tracerIn);

  /*! \brief Default destructor for RemappingVariable class */
  ~RemappingTracer() = default;

  /*! Pointer to the tracer column variable associated with the tracer */
  DEMSI::ColumnVariable<double>* tracer;

  /*! extensive tracer values before remapping (size nParticlesOld*sizePerParticle) */
  std::vector<double> extensiveTracer;

  /*! Tracer concentration after remap and correction (size nParticlesInit*sizePerParticle) */
  std::vector<double> extensiveTracerRemap;

  /*! Intensive tracer after remap and correction (size nParticlesInit*sizePerParticle) */
  std::vector<double> intensiveTracerRemap;

private:

  /*! Pointer to particles object */
  DEMSI::Tessellation* tessellation;

}; // RemappingVariable class

//------------------------------------------------------------------------------
/*! \class Remapping
    \brief Remap the particle distribution back to the initial distribution


*/
//------------------------------------------------------------------------------
class Remapping {

public:

  /*! \brief Constructor for Remapping class */
  Remapping(DEMSI::Partition* partitionIn, DEMSI::Log* logIn, DEMSI::Configs* configsIn, DEMSI::LammpsInstance* lammpsInstanceIn, DEMSI::Grid* gridIn, DEMSI::Contacts* contactsIn, DEMSI::Particles* particlesIn, DEMSI::Tessellation* tessellationIn, DEMSI::TessellationOutputStreams* tessellationOutputStreamsIn, DEMSI::Column* columnIn, DEMSI::Ocean* oceanIn, DEMSI::Clock* simulationClock);

  /*! \brief Default destructor for Remapping class */
  ~Remapping() = default;

  /*! \brief perform the remapping
      \param remapNow Perform the remapping now even if remapping alarm is not ringing
  */
  void remap(bool const remapNow = false, bool const init = false);

  /*! Vector of fields registered for possible output */
  std::vector<std::pair<double*,std::string>> fieldsRemappingWrite;

#ifdef DEBUG_REMAPPING
  double* useParticleDebug;
  double* elementAreaDebug;
  double* effectiveAreaIceDebug;
#endif

private:

  /*! Alarm for performing remapping */
  DEMSI::Alarm remappingAlarm;

  /*! Flag whether remapping to be performed */
  bool useRemapping = false;

  /*! Flag whether doing tracer conservation test */
  bool checkConservation = false;

  /*! Flag whether to check weights */
  bool checkRemappingWeights = false;

  /*! Flag whether to check effective areas */
  bool checkEffectiveArea = false;

  /*! Flag whether to check ice concentration bounds */
  bool checkIceConcentrationBounds = false;

  /*! Type of remapping correction */
  std::string remappingCorrectionType = "flux";

  /*! Maximum allowable ice concentration bounds error */
  double maxIceConcBoundError;

  // neighbour lists for fast search of neighbour init elements

  /*! maximum polygon size for determining neighbour grid size */
  double maxPolygonSize;

  /*! element number of the neighbour grid */
  int nxNeighbours, nyNeighbours;

  /*! element size of the neighbour grid */
  double dxNeighbours, dyNeighbours;

  /*! indices of init elements in the neighbour grid */
  std::map<std::pair<int, int>, std::vector<int>> indexNeighbours;

  /*! Pointer to the partiton object */
  DEMSI::Partition* partition;

  /*! Pointer to the log object */
  DEMSI::Log* log;

  /*! Pointer to the configs object */
  DEMSI::Configs* configs;

  /*! Pointer to the lammps instance object */
  DEMSI::LammpsInstance* lammpsInstance;

  /*! Pointer to the grid object */
  DEMSI::Grid* grid;

  /*! Pointer to the contacts object */
  DEMSI::Contacts* contacts;

  /*! Pointer to the particles object */
  DEMSI::Particles* particles;

  /*! Pointer to the tessellation object */
  DEMSI::Tessellation* tessellation;

  /*! Pointer to the tessellation streams object */
  DEMSI::TessellationOutputStreams* tessellationOutputStreams;

  /*! Pointer to the column object */
  DEMSI::Column* column;

  /*! Pointer to the ocean object */
  DEMSI::Ocean* ocean;

  //------------------------------------------------------------------------------
  // geometric remap
  //------------------------------------------------------------------------------

  /*! Initialize neighbour lists for fast finding of new elements */
  void init_neighbour_lists(void);

  /*! Calculate the remapping weights for the current particle distribution */
  void remapping_weights(Mat &Weights, std::vector<double> &weights, std::vector<int> &iParticlesOldWeights, std::vector<int> &iParticlesInitWeights);

  /*! Remap an individual column variable */
  void remap_variable(const int sizePerParticle, std::vector<double> &extensiveTracerRemap, std::vector<double> extensiveTracer, const Mat Weights);
  void remap_variable(const int sizePerParticle, double *extensiveTracerRemap, double *extensiveTracer, const Mat Weights);

  /*! determine the unique set of init indices that form the new set of elements */
  void new_element_indices(const int nCategories, const std::vector<double> effectiveElementArea, const std::vector<double> iceConcentration, const std::vector<int> newElements, std::vector<int> &iParticleInitFromNew, int &nParticlesNew, std::vector<std::pair<int,int>> &bonds);

  /*! Write a lammps particle init file from the new particle distribution */
  void write_lammps_file(int nParticlesNewActive, double *xNew, double *yNew, double *rNew, int *typeAll, double maxRadius);

  //------------------------------------------------------------------------------
  // flux based correction to remapping
  //------------------------------------------------------------------------------

  /*! Correct geometric remapping weights for monotonicity */
  void weights_correction(Mat &Weights, std::vector<double> effectiveElementArea);

  //------------------------------------------------------------------------------
  // OBR correction to remapping
  //------------------------------------------------------------------------------

  /*! See if parts in use */
  void particles_in_use(const Mat Weights, DEMSI::ColumnVariable<double>* effectiveElementArea, std::vector<int> &useParticle);

  /*! Optimization-based tracer remap */
  void correct_tracer_obr(DEMSI::ColumnVariable<double> *tracer, std::vector<double> &intensiveTracerRemap, std::vector<double> &extensiveTracerRemap, std::vector<double> &extensiveTracer, DEMSI::ColumnVariable<double> *parent, std::vector<double> &extensiveTracerParentRemap, std::vector<double> &extensiveTracerParent, std::vector<int> useParticle, std::vector<double> weights, std::vector<int> iParticlesOldWeights, std::vector<int> iParticlesInitWeights);

  /*! Optimization-based area remap */
  void correct_area_obr(DEMSI::ColumnVariable<double> *effectiveArea, std::vector<double> &extensiveAreaRemap, std::vector<double> &extensiveArea, std::vector<int> useParticle);

  /*! secant algorithm for optimization-based remap */
  void limit_obr(std::vector<double> &tracerVal, std::vector<double> &parentVal, std::vector<double> &tracerMin, std::vector<double> &tracerMax, double Mass);

   /*! secant algorithm for optimization-based remap */
  void limit_obr_simple(std::vector<double> &tracerVal, std::vector<double> &parentVal, std::vector<double> &tracerMin, std::vector<double> &tracerMax, double Mass);

   /*! solve QP for OBR*/
  void solve_qp(std::vector<double> &tracerArray, std::vector<double> &tracerVal, std::vector<double> &parentVal, std::vector<double> &tracerMin, std::vector<double> &tracerMax, double Mass, double lambda, double &residual);

  //------------------------------------------------------------------------------
  // Debugging
  //------------------------------------------------------------------------------

  /*! Output remapped areas for debugging */
  void output_areas(void);

  /*! Output an init double array to a netcdf file for debugging */
  void output_init_particle_array(const std::string filenameOut, double array[]);

  /*! Output an init int array to a netcdf file for debugging */
  void output_init_particle_array(const std::string filenameOut, int array[]);

  /*! output column variable mid-remapping */
  void output_post_remapping(const std::string filenameOut, DEMSI::ColumnVariable<double>* columnTracer, int nParticlesNew, const std::vector<int> iParticleInitFromNew);

}; // Remapping class

} // namespace DEMSI

#endif /* DEMSI_REMAPPING_H_ */
