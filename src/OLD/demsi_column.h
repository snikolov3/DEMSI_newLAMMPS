#ifndef DEMSI_COLUMN_H_
#define DEMSI_COLUMN_H_

/*! \file   demsi_column.h
    \brief  Header file for the DEMSI::Column class
*/

#include "demsi_logging.h"
#include "demsi_configs.h"
#include "demsi_particles.h"
#include "demsi_initial_tessellation.h"
#include "demsi_partition.h"
#include "demsi_lmp_instance.h"
#include "demsi_column_tracers.h"
#include "demsi_column_variables.h"
#include "demsi_time.h"
#include "demsi_grid.h"
#include "demsi_forcing.h"

#include <Kokkos_Core.hpp>

namespace DEMSI {

/*! \class Column
 \brief Class providing wrapper for the column icepack library.
*/
class Column{
public:

/*! \brief constructor for Column class. */
Column(DEMSI::Partition* partitionIn, DEMSI::Log *logIn, DEMSI::Configs *configs, DEMSI::LammpsInstance* lammpsInstanceIn, DEMSI::Grid* gridIn, DEMSI::Particles *particlesIn, DEMSI::Tessellation* tessellationIn, DEMSI::Forcing* forcingIn, DEMSI::Clock* simulationClockIn, DEMSI::TimeInterval* timeStepIntervalIn);

  /*! \brief default destructor for Column class. */
  ~Column() = default;

  /*! column variable dimensions */
  DEMSI::ColumnDimensions* columnDimensions;

  /*! pointers to all column variables */
  DEMSI::ColumnVariables* columnVariables;

  /*! init the column object */
  void init(void);

  /*! init the column object after remapping */
  void init_remap(void);

  /*! \brief column time stepping before dynamics */
  void pre_dynamics(void);

  /*! \brief column ridging */
  void ridging(void);

  /*! \brief column time stepping during dynamics */
  void dynamics(void);

  /*! \brief column time stepping after dynamics */
  void post_dynamics(void);

  /*! \brief init the shortwave calculation */
  void init_shortwave(void);

  /*! \brief init the slab ocean sst */
  void init_sst(void);

  /*! \brief set ridging variables needed for dynamics */
  void set_ridging_variables(void);

  /*! \brief perform column ridging */
  void ridge(void);

  /*! \brief Interpolate column forcing to particles and post-process */
  void interpolate_forcing(void);

  /*! \brief Aggregate column by category tracers over categories */
  void aggregate(void);

  /*! \brief reset fluxes */
  void finalize_timestep(void);

  /*! \brief Reset oceanic coupler fluxes */
  void reset_ocean_coupler_fluxes(void);

  /*! \brief update particle mass */
  void update_mass(void);

  /*! \brief compute mean and minimum element thickness */
  void compute_mean_min_thickness(void);

  /*! \brief convert ice/snow volumes to thicknesses */
  void to_thickness(void);

  /*! \brief convert ice/snow thicknesses to volumes */
  void from_thickness(void);

  /*! \brief delete particles that have too little sea ice */
  void delete_particles(void);

  /*! \brief Transfer column variables between processors */
  void processor_transfer(void);

  /*! \brief Reset the transfer across processors for column variables */
  void reset_processor_transfer(void);

  /*! \brief Reset transfer column data between procqessors after remapping.
      \param nParticles Number of new particles
  */
  void reset_processor_transfer_remap(const int nParticles);

  // column variables
#include "demsi_column_variables_declaration.inc"

  /*! Tracer tree */
  DEMSI::TracerTree* tracerTree;

  /*! Vector of column tracer pointers for active ones */
  std::vector<DEMSI::ColumnVariable<double>*> columnTracerPtrs;

  /*! Are we using column physics? */
  bool useColumn;

private:

  /*! Pointer to the partition object */
  DEMSI::Partition *partition;

  /*! Pointer to the log object */
  DEMSI::Log *log;

  /*! Pointer to the configs object */
  DEMSI::Configs *configs;

  /*! Pointer to the lammps instance object */
  DEMSI::LammpsInstance *lammpsInstance;

  /*! Pointer to the grid object */
  DEMSI::Grid* grid;

  /*! Pointer to the particles object */
  DEMSI::Particles *particles;

  /*! Pointer to the tessellation object */
  DEMSI::Tessellation *tessellation;

  /*! Pointer to the forcing object */
  DEMSI::Forcing* forcing;

  /*! Pointer to the simulation clock */
  DEMSI::Clock* simulationClock;

  /*! Pointer to the thermodynamic timestep time interval object. */
  DEMSI::TimeInterval* timeStepInterval;

  /*! Are we using shortwave physics? */
  bool useColumnShortwave = false;

  /*! Are we using vertical thermodynamics? */
  bool useColumnVerticalThermodynamics = false;

  /*! Are we using itd thermodynamics? */
  bool useColumnItdThermodynamics = false;

  /*! Are we using ridging? */
  bool useColumnRidging = false;

  /*! What ocean type are we using */
  std::string oceanType;

  /*! icepack string length */
  int columnStrLen;

  /*! ice thickness category limits */
  double* categoryThicknessLimits;

  /*! Pointer to the column tracer object */
  DEMSI::ColumnTracers* columnTracers;

  /*! \brief Initialize DEMSI column dimensions */
  void init_demsi_column_dimensions(void);

  /*! \brief Allocate DEMSI mass related variables */
  void allocate_demsi_mass_variables(void);

  /*! \brief Initialize area and volume variables from input file */
  void init_demsi_mass_variables_file(void);

  /*! \brief Initialize area and volume variables from constant config values */
  void init_demsi_mass_variables_const(void);

  /*! \brief Initialize area and volume variables from config thickness distribution */
  void init_demsi_mass_variables_thickness_dist(void);

  /*! \brief Initialize area and volume variables from cice method */
  void init_demsi_mass_variables_cice(void);

  /*! \brief Initialize DEMSI mass related variables */
  void init_demsi_mass_variables(void);

  /*! \brief init element area */
  void effective_element_area(void);

  /*! \brief Initialize DEMSI column variables */
  void init_demsi_column_variables(void);

  /*! \brief Initialize the vertical thermo profiles: init_thermo_vertical */
  void init_thermo_vertical(void);

  /*! \brief Initialize the ice thickness distribution */
  void init_ice_thickness_distribution(void);

  /*! \brief Initialize the ice thickness distribution for no column */
  void init_ice_thickness_distribution_no_column(void);

  /*! \brief Initialize the column state by particle */
  void init_column_state(void);

  /*! \brief Create a tree of the tracer hirearchy */
  void create_tracer_tree(void);

  /*! \brief Initialize the atmospheric coupler fields */
  void initialize_atmos_coupler_fields(void);

  /*! \brief Initialize the oceanic coupler fields */
  void initialize_ocean_coupler_fields(void);

  /*! \brief Reset atmospheric coupler fluxes */
  void reset_atmospheric_coupler_fluxes(void);

  /*! \brief Reinitialize diagnostic variables */
  void reinitialize_diagnostics(void);

  /*! \brief Print icepack warning messages */
  bool check_warnings(const std::string message);

  /*! \brief Print icepack warning messages in particle loop */
  bool check_warnings(const std::string message, const int iParticle);

  /*! \brief Set column configs in Icepack */
  void column_configs(void);

  /*! \brief Set an individual column parameter in icepack
      \param configGroupName Config group name for the config
      \param configName Name of the config
  */
  void set_column_parameter(const std::string configGroupName, const std::string configName);

  /*! \brief Time stepping of the icepack radiation preparation subroutine */
  void run_prep_radiation(void);

  /*! \brief Time stepping of the icepack radiation subroutine
      \param lInitialization True if at start of simulation
  */
  void run_step_radiation(const bool lInitialization);

  /*! \brief Time stepping of the icepack vertical thermodyanmic subroutine */
  void run_step_therm1(void);

  /*! \brief Time stepping of the icepack ridging subroutine */
  void run_step_ridge(void);

  /*! \brief Time stepping of the icepack vertical itd subroutine */
  void run_step_therm2(void);

  /*! \brief Coupling prep */
  void coupling_prep(void);

  /*! \brief Ocean mixed layer calculation */
  void ocean_mixed_layer(void);

  /*! \brief Limit air temperature over ice */
  void limit_temperature(void);

  /*! \brief Limit specific humidity */
  void limit_specific_humidity(void);

  /*! \brief Get shortwave radiation from cloudiness and location/time */
  void shortwave_from_cloud_fraction(void);

  /*! \brief Split shortwave between bands */
  void split_shortwave(void);

  /*! \brief Limit forcing variables to physically realistic values */
  void apply_physical_limits(void);

  /*! \brief Calculate the downwelling longwave radiation */
  void calculate_longwave_down(void);

  /*! \brief Calculate the downwelling longwave radiation from Rosati and Miyakoda */
  void longwave_rosati_miyakoda(void);

  /*! \brief Calculate the downwelling longwave radiation from Parkinson and Washington */
  void longwave_parkinson_and_washington(void);

  /*! \brief Scale precipitation according to units */
  void scale_precipitation(void);

  /*! \brief Split precipitation between snow and rain by air temperature */
  void split_precipitation(void);

  /*! Column variable to test transfer between processors */
  DEMSI::ColumnVariable<int>* columnVariablesTestTransfer = NULL;

  /*! \brief Initialize transfer of column variables between processors */
  void init_processor_transfer(void);

}; // Column class

} // namespace DEMSI

#endif /* COLUMN_H_ */
