#ifndef DEMSI_OCEAN_H_
#define DEMSI_OCEAN_H_

/*! \file   demsi_ocean.h
    \brief  Header file for the DEMSI::Ocean classes
*/

#include "demsi_particles.h"
#include "demsi_initial_tessellation.h"
#include "demsi_forcing.h"
#include "demsi_time.h"
#include "demsi_grid.h"
#include "demsi_column.h"
#include "demsi_column_variables.h"
#include <vector>

namespace DEMSI {

//------------------------------------------------------------------------------
/*! \class Ocean
    \brief Slab ocean model fixed with the initial element distribution.

    This class defines a slab ocean model with fixed elements. Grid cells are
    the initial element distribution, and are coupled with the mobile sea ice
    elements during geometric remapping.
*/
//------------------------------------------------------------------------------
class Ocean {

public:

  /*! \brief Ocean class constructor
      \param configsIn Pointer to the configs object
      \param logIn Pointer to the log object
      \param tessellationIn Pointer to tessellation object
      \param columnIn Pointer to column object
      \param forcingIn Pointer to forcing object
      \param gridIn Pointer to grid object
      \param simulationClockIn Pointer to the simulation clock
      \param timeStepIntervalIn Pointer to the thermodynamic timestep time interval object
   */
  Ocean(DEMSI::Log* logIn, DEMSI::Configs* configsIn, DEMSI::Grid* gridIn, DEMSI::Tessellation* tessellationIn, DEMSI::Forcing* forcingIn, DEMSI::Column* columnIn, DEMSI::Clock* simulationClockIn, DEMSI::TimeInterval* timeStepIntervalIn);

  /*! \brief Default destructor.*/
  ~Ocean() = default;

  /*! \brief Initialize the ocean SST */
  void init_sst(void);

  /*! \brief Interpolate forcing to ocean. */
  void interpolate_forcing(void);

  /*! \brief Update the state of the slab ocean model */
  void update(void);

  /*! \brief Update the state of the slab ocean model */
  void update_remap(const int nParticlesNew, const std::vector<int> iParticleInitFromNew);

  /*! \brief Get the indices of elements from the initial distribution that have frazil
      \return vector of indices of elements from the initial distribution that have frazil
  */
  std::vector<int> get_new_ice_elements(void);

  /*! \brief Copy variables from ice remapped elements to ocean init elements
      \param nParticlesNew Number of new remapped ice elements
      \param iParticleInitFromNew Indices of new remapped ice elements
   */
  void ice_to_ocean(const int nParticlesNew, const std::vector<int> iParticleInitFromNew);

  /*! \brief Copy variables from ocean init elements to ice remapped elements
      \param nParticlesNew Number of new remapped ice elements
      \param iParticleInitFromNew Indices of new remapped ice elements
  */
  void ocean_to_ice(const int nParticlesNew, const std::vector<int> iParticleInitFromNew);

  /*! Sea surface temperature */
  double* seaSurfaceTemperature;

  /*! Vector of fields registered for possible output */
  std::vector<std::pair<double*,std::string>> fieldsOceanWrite;

  /*! What ocean type are we using */
  std::string oceanType;

private:

  /*! Pointer to the log object */
  DEMSI::Log* log;

  /*! Pointer to the configs object */
  DEMSI::Configs* configs;

  /*! Pointer to the grid object */
  DEMSI::Grid* grid;

  /*! Pointer to the tessellation object */
  DEMSI::Tessellation* tessellation;

  /*! Pointer to the forcing object */
  DEMSI::Forcing* forcing;

  /*! Pointer to the column object */
  DEMSI::Column* column;

  /*! Pointer to the simulation clock object */
  DEMSI::Clock* simulationClock;

  /*! Pointer to the thermodynamic timestep time interval object. */
  DEMSI::TimeInterval* timeStepInterval;

  // ocean state
  /*! Sea surface salinity */
  double* seaSurfaceSalinity;

  /*! Sea freezing temperature */
  double* seaFreezingTemperature;

  // ice state on ocean grid
  /*! ice concentration */
  double* iceConcentration;

  /*! surface temperature of ice */
  double* iceSurfaceTemperature;

  // sea ice to ocean fluxes
  /*! heat flux from ice to ocean */
  double* oceanHeatFlux;

  /*! shortwave flux into ocean */
  double* oceanShortwaveFlux;

  // ocean to sea ice fluxes
  /*! freezing melting temperature accumlated over coupling interval */
  double* freezingMeltingPotential;

  // latitude/longitude
  double* latitude;
  double* longitude;

  // atmospheric forcing
  double* airTemperature;
  double* airSpecificHumidity;
  double* shortwaveDown;
  double* longwaveDown;
  double* uAirVelocity;
  double* vAirVelocity;
  double* snowfallRate;
  double* rainfallRate;
  double* precipitationRate;
  double* cloudFraction;

  // shortwave
  double* shortwaveVisibleDirectDown;
  double* shortwaveVisibleDiffuseDown;
  double* shortwaveIRDirectDown;
  double* shortwaveIRDiffuseDown;

  // other atmospheric forcing
  double* airDensity;
  double* windSpeed;
  double* airPotentialTemperature;
  double* sensibleTransferCoefficient;
  double* latentTransferCoefficient;
  double* uAirStress;
  double* vAirStress;
  double* airLevelHeight;
  double* atmosReferenceTemperature2mOcean;
  double* atmosReferenceHumidity2mOcean;
  double* airDragCoefficient;
  double* airOceanDragCoefficientRatio;
  double* albedoVisibleDirectOcean;
  double* albedoIRDirectOcean;
  double* albedoVisibleDiffuseOcean;
  double* albedoIRDiffuseOcean;
  double* longwaveUpOcean;
  double* sensibleHeatFluxOcean;
  double* latentHeatFluxOcean;
  double* evaporativeWaterFluxOcean;

  // ocean forcing
  double* uOceanVelocity;
  double* vOceanVelocity;
  double* seaSurfaceTiltU;
  double* seaSurfaceTiltV;
  double* oceanMixedLayerDepth;
  double* oceanHeatFluxConvergence;
  double* airStressOceanU;
  double* airStressOceanV;

  // testing
  double* testVar;

};

} // namespace DEMSI

#endif
