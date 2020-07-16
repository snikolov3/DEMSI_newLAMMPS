#include "demsi_ocean.h"

#include "demsi_time.h"
#include <Kokkos_Core.hpp>
//#include <netcdf.h>
//#include "demsi_file_utils.h"
#include <iomanip>
#include <fstream>
#include <cmath>

extern "C" {double column_get_constant(char*);}
extern "C" {void column_query_parameters_logical(char*, int*);}
extern "C" {double column_liquidus_temperature(double*);}
extern "C" {void column_limit_temperature(double*, double*);}
extern "C" {void column_limit_specific_humidity(double*, double*);}
extern "C" {void column_longwave_rosati_miyakoda(double*, double*, double*, double*, double*, double*, double*);}
extern "C" {void column_longwave_parkinson_and_washington(double*, double*, double*);}
extern "C" {void column_shortwave_from_cloud_fraction(double*, double*, double*, double*, double*, double*, int*);}
extern "C" {void column_split_precipitation(double*, double*, double*, double*);}
extern "C" {void column_ocean_mixed_layer(double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*);}
extern "C" {void column_postprocess_ocean_forcing(double*, double*, double*);}
extern "C" {void column_initialize_ocean_coupler_fields(int*, double*, double*, double*);}

namespace DEMSI {

// interpolate a forcing field for one particle
void interpolate_field(const Kokkos::View<double**>* fieldEulerian, double &fieldParticle, const double weight[], const int i[], const int j[]) {

  if (fieldEulerian->use_count() > 0) {
    fieldParticle =
      (*fieldEulerian)(i[0],j[0]) * weight[0] +
      (*fieldEulerian)(i[1],j[1]) * weight[1] +
      (*fieldEulerian)(i[2],j[2]) * weight[2] +
      (*fieldEulerian)(i[3],j[3]) * weight[3];
  }

}

// constructor for Ocean class
Ocean::Ocean(DEMSI::Log* logIn, DEMSI::Configs* configsIn, DEMSI::Grid* gridIn, DEMSI::Tessellation* tessellationIn, DEMSI::Forcing* forcingIn, DEMSI::Column* columnIn, DEMSI::Clock* simulationClockIn, DEMSI::TimeInterval* timeStepIntervalIn) {

  log = logIn;
  configs = configsIn;
  grid = gridIn;
  tessellation = tessellationIn;
  forcing = forcingIn;
  column = columnIn;

  simulationClock = simulationClockIn;
  timeStepInterval = timeStepIntervalIn;

  // what is the ocean type?
  oceanType = "embedded";
  if (configs->exists({"ConfigGroup:ocean","Config:oceanType"})) {
    configs->get({"ConfigGroup:ocean","Config:oceanType"}, oceanType);
  } // oceanType option exists

  if (oceanType == "init") {

    // allocate ocean state
    seaSurfaceTemperature  = new double[tessellation->nParticlesInit];
    seaSurfaceSalinity     = new double[tessellation->nParticlesInit];
    seaFreezingTemperature = new double[tessellation->nParticlesInit];

    // allocate sea ice state
    iceConcentration      = new double[tessellation->nParticlesInit];
    iceSurfaceTemperature = new double[tessellation->nParticlesInit];

    // allocate sea ice to ocean fluxes
    oceanHeatFlux      = new double[tessellation->nParticlesInit];
    oceanShortwaveFlux = new double[tessellation->nParticlesInit];

    // allocate ocean to sea ice fluxes
    freezingMeltingPotential = new double[tessellation->nParticlesInit];

    // longitude/latitude
    tessellation->latitudeInit  = new double[tessellation->nParticlesInit];
    tessellation->longitudeInit = new double[tessellation->nParticlesInit];

    // latitude/longitude
    latitude  = new double[tessellation->nParticlesInit];
    longitude = new double[tessellation->nParticlesInit];

    // atmospheric forcing
    airTemperature      = new double[tessellation->nParticlesInit];
    airSpecificHumidity = new double[tessellation->nParticlesInit];
    shortwaveDown       = new double[tessellation->nParticlesInit];
    longwaveDown        = new double[tessellation->nParticlesInit];
    uAirVelocity        = new double[tessellation->nParticlesInit];
    vAirVelocity        = new double[tessellation->nParticlesInit];
    snowfallRate        = new double[tessellation->nParticlesInit];
    rainfallRate        = new double[tessellation->nParticlesInit];
    precipitationRate   = new double[tessellation->nParticlesInit];
    cloudFraction       = new double[tessellation->nParticlesInit];

    // shortwave
    shortwaveVisibleDirectDown  = new double[tessellation->nParticlesInit];
    shortwaveVisibleDiffuseDown = new double[tessellation->nParticlesInit];
    shortwaveIRDirectDown       = new double[tessellation->nParticlesInit];
    shortwaveIRDiffuseDown      = new double[tessellation->nParticlesInit];

    // derived atmospheric forcing
    airDensity                       = new double[tessellation->nParticlesInit];
    windSpeed                        = new double[tessellation->nParticlesInit];
    airPotentialTemperature          = new double[tessellation->nParticlesInit];
    sensibleTransferCoefficient      = new double[tessellation->nParticlesInit];
    latentTransferCoefficient        = new double[tessellation->nParticlesInit];
    uAirStress                       = new double[tessellation->nParticlesInit];
    vAirStress                       = new double[tessellation->nParticlesInit];
    airLevelHeight                   = new double[tessellation->nParticlesInit];
    atmosReferenceTemperature2mOcean = new double[tessellation->nParticlesInit];
    atmosReferenceHumidity2mOcean    = new double[tessellation->nParticlesInit];
    airDragCoefficient               = new double[tessellation->nParticlesInit];
    airOceanDragCoefficientRatio     = new double[tessellation->nParticlesInit];
    albedoVisibleDirectOcean         = new double[tessellation->nParticlesInit];
    albedoIRDirectOcean              = new double[tessellation->nParticlesInit];
    albedoVisibleDiffuseOcean        = new double[tessellation->nParticlesInit];
    albedoIRDiffuseOcean             = new double[tessellation->nParticlesInit];
    longwaveUpOcean                  = new double[tessellation->nParticlesInit];
    sensibleHeatFluxOcean            = new double[tessellation->nParticlesInit];
    latentHeatFluxOcean              = new double[tessellation->nParticlesInit];
    evaporativeWaterFluxOcean        = new double[tessellation->nParticlesInit];

    // ocean forcing
    uOceanVelocity           = new double[tessellation->nParticlesInit];
    vOceanVelocity           = new double[tessellation->nParticlesInit];
    seaSurfaceTiltU          = new double[tessellation->nParticlesInit];
    seaSurfaceTiltV          = new double[tessellation->nParticlesInit];
    oceanMixedLayerDepth     = new double[tessellation->nParticlesInit];
    oceanHeatFluxConvergence = new double[tessellation->nParticlesInit];
    airStressOceanU          = new double[tessellation->nParticlesInit];
    airStressOceanV          = new double[tessellation->nParticlesInit];

    // testing
    testVar = new double[tessellation->nParticlesInit];

    // register ocean particles for output
    tessellation->register_output_field(seaSurfaceTemperature, "seaSurfaceTemperature");
    tessellation->register_output_field(seaSurfaceSalinity, "seaSurfaceSalinity");
    tessellation->register_output_field(seaFreezingTemperature, "seaFreezingTemperature");
    tessellation->register_output_field(iceConcentration, "iceConcentration");
    tessellation->register_output_field(iceSurfaceTemperature, "iceSurfaceTemperature");
    tessellation->register_output_field(oceanHeatFlux, "oceanHeatFlux");
    tessellation->register_output_field(oceanShortwaveFlux, "oceanShortwaveFlux");
    tessellation->register_output_field(freezingMeltingPotential, "freezingMeltingPotential");
    tessellation->register_output_field(latitude, "latitude");
    tessellation->register_output_field(longitude, "longitude");
    tessellation->register_output_field(airTemperature, "airTemperature");
    tessellation->register_output_field(airSpecificHumidity, "airSpecificHumidity");
    tessellation->register_output_field(shortwaveDown, "shortwaveDown");
    tessellation->register_output_field(longwaveDown, "longwaveDown");
    tessellation->register_output_field(uAirVelocity, "uAirVelocity");
    tessellation->register_output_field(vAirVelocity, "vAirVelocity");
    tessellation->register_output_field(snowfallRate, "snowfallRate");
    tessellation->register_output_field(rainfallRate, "rainfallRate");
    tessellation->register_output_field(precipitationRate, "precipitationRate");
    tessellation->register_output_field(cloudFraction, "cloudFraction");
    tessellation->register_output_field(shortwaveVisibleDirectDown, "shortwaveVisibleDirectDown");
    tessellation->register_output_field(shortwaveVisibleDiffuseDown, "shortwaveVisibleDiffuseDown");
    tessellation->register_output_field(shortwaveIRDirectDown, "shortwaveIRDirectDown");
    tessellation->register_output_field(shortwaveIRDiffuseDown, "shortwaveIRDiffuseDown");
    tessellation->register_output_field(airDensity, "airDensity");
    tessellation->register_output_field(windSpeed, "windSpeed");
    tessellation->register_output_field(airPotentialTemperature, "airPotentialTemperature");
    tessellation->register_output_field(sensibleTransferCoefficient, "sensibleTransferCoefficient");
    tessellation->register_output_field(latentTransferCoefficient, "latentTransferCoefficient");
    tessellation->register_output_field(uAirStress, "uAirStress");
    tessellation->register_output_field(vAirStress, "vAirStress");
    tessellation->register_output_field(airLevelHeight, "airLevelHeight");
    tessellation->register_output_field(atmosReferenceTemperature2mOcean, "atmosReferenceTemperature2mOcean");
    tessellation->register_output_field(atmosReferenceHumidity2mOcean, "atmosReferenceHumidity2mOcean");
    tessellation->register_output_field(airDragCoefficient, "airDragCoefficient");
    tessellation->register_output_field(airOceanDragCoefficientRatio, "airOceanDragCoefficientRatio");
    tessellation->register_output_field(albedoVisibleDirectOcean, "albedoVisibleDirectOcean");
    tessellation->register_output_field(albedoIRDirectOcean, "albedoIRDirectOcean");
    tessellation->register_output_field(albedoVisibleDiffuseOcean, "albedoVisibleDiffuseOcean");
    tessellation->register_output_field(albedoIRDiffuseOcean, "albedoIRDiffuseOcean");
    tessellation->register_output_field(longwaveUpOcean, "longwaveUpOcean");
    tessellation->register_output_field(sensibleHeatFluxOcean, "sensibleHeatFluxOcean");
    tessellation->register_output_field(latentHeatFluxOcean, "latentHeatFluxOcean");
    tessellation->register_output_field(evaporativeWaterFluxOcean, "evaporativeWaterFluxOcean");
    tessellation->register_output_field(uOceanVelocity, "uOceanVelocity");
    tessellation->register_output_field(vOceanVelocity, "vOceanVelocity");
    tessellation->register_output_field(seaSurfaceTiltU, "seaSurfaceTiltU");
    tessellation->register_output_field(seaSurfaceTiltV, "seaSurfaceTiltV");
    tessellation->register_output_field(oceanMixedLayerDepth, "oceanMixedLayerDepth");
    tessellation->register_output_field(oceanHeatFluxConvergence, "oceanHeatFluxConvergence");
    tessellation->register_output_field(airStressOceanU, "airStressOceanU");
    tessellation->register_output_field(airStressOceanV, "airStressOceanV");
    tessellation->register_output_field(testVar, "testVar");

    // init variables
    char paramNameTsmelt[7] = "Tsmelt";
    double snowSurfaceMeltTemp = column_get_constant(paramNameTsmelt);
    char paramNameTffresh[8] = "Tffresh";
    double KelvinToCelcius = column_get_constant(paramNameTffresh);
    int calcSurfaceTemperatureInt;
    char constNamecalc_Tsfc[10] = "calc_Tsfc";
    column_query_parameters_logical(constNamecalc_Tsfc, &calcSurfaceTemperatureInt);
    bool calcSurfaceTemperature = (bool) calcSurfaceTemperatureInt;

    for (int iParticleInit = 0 ; iParticleInit < tessellation->nParticlesInit ; iParticleInit++) {

      // zero out
      seaSurfaceTemperature[iParticleInit] = 0.0;
      seaSurfaceSalinity[iParticleInit] = 0.0;
      seaFreezingTemperature[iParticleInit] = 0.0;
      iceConcentration[iParticleInit] = 0.0;
      iceSurfaceTemperature[iParticleInit] = 0.0;
      oceanHeatFlux[iParticleInit] = 0.0;
      oceanShortwaveFlux[iParticleInit] = 0.0;
      freezingMeltingPotential[iParticleInit] = 0.0;
      latitude[iParticleInit] = 0.0;
      longitude[iParticleInit] = 0.0;
      airTemperature[iParticleInit] = 0.0;
      airSpecificHumidity[iParticleInit] = 0.0;
      shortwaveDown[iParticleInit] = 0.0;
      longwaveDown[iParticleInit] = 0.0;
      uAirVelocity[iParticleInit] = 0.0;
      vAirVelocity[iParticleInit] = 0.0;
      snowfallRate[iParticleInit] = 0.0;
      rainfallRate[iParticleInit] = 0.0;
      precipitationRate[iParticleInit] = 0.0;
      cloudFraction[iParticleInit] = 0.0;
      shortwaveVisibleDirectDown[iParticleInit] = 0.0;
      shortwaveVisibleDiffuseDown[iParticleInit] = 0.0;
      shortwaveIRDirectDown[iParticleInit] = 0.0;
      shortwaveIRDiffuseDown[iParticleInit] = 0.0;
      airDensity[iParticleInit] = 0.0;
      windSpeed[iParticleInit] = 0.0;
      airPotentialTemperature[iParticleInit] = 0.0;
      sensibleTransferCoefficient[iParticleInit] = 0.0;
      latentTransferCoefficient[iParticleInit] = 0.0;
      uAirStress[iParticleInit] = 0.0;
      vAirStress[iParticleInit] = 0.0;
      airLevelHeight[iParticleInit] = 0.0;
      atmosReferenceTemperature2mOcean[iParticleInit] = 0.0;
      atmosReferenceHumidity2mOcean[iParticleInit] = 0.0;
      airDragCoefficient[iParticleInit] = 0.0;
      airOceanDragCoefficientRatio[iParticleInit] = 0.0;
      albedoVisibleDirectOcean[iParticleInit] = 0.0;
      albedoIRDirectOcean[iParticleInit] = 0.0;
      albedoVisibleDiffuseOcean[iParticleInit] = 0.0;
      albedoIRDiffuseOcean[iParticleInit] = 0.0;
      longwaveUpOcean[iParticleInit] = 0.0;
      sensibleHeatFluxOcean[iParticleInit] = 0.0;
      latentHeatFluxOcean[iParticleInit] = 0.0;
      evaporativeWaterFluxOcean[iParticleInit] = 0.0;
      uOceanVelocity[iParticleInit] = 0.0;
      vOceanVelocity[iParticleInit] = 0.0;
      seaSurfaceTiltU[iParticleInit] = 0.0;
      seaSurfaceTiltV[iParticleInit] = 0.0;
      oceanMixedLayerDepth[iParticleInit] = 0.0;
      oceanHeatFluxConvergence[iParticleInit] = 0.0;
      airStressOceanU[iParticleInit] = 0.0;
      airStressOceanV[iParticleInit] = 0.0;
      testVar[iParticleInit] = 0.0;

      // non-zero values
      airDensity             [iParticleInit] = 1.3;
      airLevelHeight         [iParticleInit] = 10.0;
      seaSurfaceSalinity     [iParticleInit] = 34.0;
      airTemperature         [iParticleInit] = 253.0;
      airPotentialTemperature[iParticleInit] = 253.0;
      airSpecificHumidity    [iParticleInit] = 0.0006;
      longwaveDown           [iParticleInit] = 180.0;
      oceanMixedLayerDepth   [iParticleInit] = 20.0;

      iceConcentration[iParticleInit] = 1.0;

      // ice surface temperature
      if (calcSurfaceTemperature) {
	iceSurfaceTemperature[iParticleInit] = std::min(snowSurfaceMeltTemp, airTemperature[iParticleInit] - KelvinToCelcius); // deg C
      } else {
	iceSurfaceTemperature[iParticleInit] = seaFreezingTemperature[iParticleInit]; // default
      } // calcSurfaceTemperature

    } // iParticleInit

    // initialize coupler fields
    bool doRestart;
    configs->get({"ConfigGroup:simulationTiming","Config:doRestart"}, doRestart);
    int doRestartIn = (int) doRestart;

    for (int iParticleInit = 0 ; iParticleInit < tessellation->nParticlesInit ; iParticleInit++) {

      column_initialize_ocean_coupler_fields(&doRestartIn,
					     &seaSurfaceTemperature[iParticleInit],
					     &seaFreezingTemperature[iParticleInit],
					     &seaSurfaceSalinity[iParticleInit]);

    } // iParticleInit

  } // oceanType == "init"

} // Ocean::Ocean

// Initialize the ocean SST
void Ocean::init_sst(void) {

  if (oceanType == "init") {

    // init sea surface temperature
    for (int iParticleInit = 0 ; iParticleInit < tessellation->nParticlesInit ; iParticleInit++) {

      // interpolation
      double x = tessellation->xInit[iParticleInit];
      double y = tessellation->yInit[iParticleInit];

      int i[4], j[4];
      double weight[4];
      grid->interpolation_weights(x, y, i, j, weight);

      interpolate_field(&(forcing->seaSurfaceTemperature), seaSurfaceTemperature[iParticleInit], weight, i, j);

      seaSurfaceTemperature [iParticleInit] = std::max(seaSurfaceTemperature[iParticleInit], seaFreezingTemperature[iParticleInit]);

    } // iParticleInit

  } // oceanType

} // Ocean::init_sst

// Interpolate forcing for the ocean
void Ocean::interpolate_forcing(void) {

  if (oceanType == "init") {

    int i[4], j[4];
    double weight[4];

    // atmosphere
    for (int iParticleInit = 0 ; iParticleInit < tessellation->nParticlesInit ; iParticleInit++) {

      // interpolation
      double x = tessellation->xInit[iParticleInit];
      double y = tessellation->yInit[iParticleInit];

      grid->interpolation_weights(x, y, i, j, weight);

      // grid
      interpolate_field(&(grid->latitude),  latitude[iParticleInit],  weight, i, j);
      interpolate_field(&(grid->longitude), longitude[iParticleInit], weight, i, j);

      // atmospheric forcing
      interpolate_field(&(forcing->airTemperature),    airTemperature[iParticleInit],      weight, i, j);
      interpolate_field(&(forcing->specificHumidity),  airSpecificHumidity[iParticleInit], weight, i, j);
      interpolate_field(&(forcing->shortwaveFluxDown), shortwaveDown[iParticleInit],       weight, i, j);
      interpolate_field(&(forcing->longwaveFluxDown),  longwaveDown[iParticleInit],        weight, i, j);
      interpolate_field(&(forcing->xAtmWind),          uAirVelocity[iParticleInit],        weight, i, j);
      interpolate_field(&(forcing->yAtmWind),          vAirVelocity[iParticleInit],        weight, i, j);
      interpolate_field(&(forcing->snowfallRate),      snowfallRate[iParticleInit],        weight, i, j);
      interpolate_field(&(forcing->rainfallRate),      rainfallRate[iParticleInit],        weight, i, j);
      interpolate_field(&(forcing->precipitationRate), precipitationRate[iParticleInit],   weight, i, j);
      interpolate_field(&(forcing->cloudFraction),     cloudFraction[iParticleInit],       weight, i, j);

    } // iParticleInit

    // limit the air temperature
    bool useLimitTemperature;
    configs->get({"ConfigGroup:forcing","Config:useLimitTemperature"}, useLimitTemperature);
    if (useLimitTemperature) {

      for (int iParticleInit = 0 ; iParticleInit < tessellation->nParticlesInit ; iParticleInit++) {
	column_limit_temperature(&iceConcentration[iParticleInit],
				 &airTemperature[iParticleInit]);
      } // iParticle

    } // useLimitTemperature

    // limit the air specific humidity
    bool useLimitHumidity;
    configs->get({"ConfigGroup:forcing","Config:useLimitHumidity"}, useLimitHumidity);
    if (useLimitHumidity) {

      for (int iParticleInit = 0 ; iParticleInit < tessellation->nParticlesInit ; iParticleInit++) {
	column_limit_specific_humidity(&airTemperature[iParticleInit],
				       &airSpecificHumidity[iParticleInit]);
      } // iParticle

    } // useLimitHumidity

    // calculate the downwave radiation from cloudiness, longitude and time of day
    bool useGetShortwaveFromCloudiness;
    configs->get({"ConfigGroup:forcing","Config:useGetShortwaveFromCloudiness"}, useGetShortwaveFromCloudiness);
    if (useGetShortwaveFromCloudiness) {

      DEMSI::Time currentTime = simulationClock->get_time();

      int dayOfYear = currentTime.day_of_year();
      double secondsIntoDay = (double) currentTime.seconds_into_day();

      for (int iParticleInit = 0 ; iParticleInit < tessellation->nParticlesInit ; iParticleInit++) {
	column_shortwave_from_cloud_fraction(&shortwaveDown[iParticleInit],
					     &longitude[iParticleInit],
					     &latitude[iParticleInit],
					     &cloudFraction[iParticleInit],
					     &airSpecificHumidity[iParticleInit],
					     &secondsIntoDay,
					     &dayOfYear);
      } // iParticle

    } // useGetShortwaveFromCloudiness

    // split shortwave radiation between bands
    bool useSplitShortwave;
    configs->get({"ConfigGroup:forcing","Config:useSplitShortwave"}, useSplitShortwave);
    if (useSplitShortwave) {

      double fracShortwaveVisibleDirect  = 0.28; // fraction of incoming shortwave in visible direct band
      double fracShortwaveVisibleDiffuse = 0.24; // fraction of incoming shortwave in visible diffuse band
      double fracShortwaveIRDirectDown   = 0.31; // fraction of incoming shortwave in near IR direct band
      double fracShortwaveIRDiffuseDown  = 0.17; // fraction of incoming shortwave in near IR diffuse band

      for (int iParticleInit = 0 ; iParticleInit < tessellation->nParticlesInit ; iParticleInit++) {
	shortwaveVisibleDirectDown[iParticleInit]  = shortwaveDown[iParticleInit] * fracShortwaveVisibleDirect;
	shortwaveVisibleDiffuseDown[iParticleInit] = shortwaveDown[iParticleInit] * fracShortwaveVisibleDiffuse;
	shortwaveIRDirectDown[iParticleInit]       = shortwaveDown[iParticleInit] * fracShortwaveIRDirectDown;
	shortwaveIRDiffuseDown[iParticleInit]      = shortwaveDown[iParticleInit] * fracShortwaveIRDiffuseDown;
      } // iParticleInit

    } // useSplitShortwave

    // apply physical limits to forcing variables
    bool useApplyPhysicalLimits;
    configs->get({"ConfigGroup:forcing","Config:useApplyPhysicalLimits"}, useApplyPhysicalLimits);
    if (useApplyPhysicalLimits) {

      for (int iParticleInit = 0 ; iParticleInit < tessellation->nParticlesInit ; iParticleInit++) {
	cloudFraction[iParticleInit]       = std::max(std::min(cloudFraction[iParticleInit],1.0),0.0);
	shortwaveDown[iParticleInit]       = std::max(shortwaveDown[iParticleInit],0.0);
	rainfallRate[iParticleInit]        = std::max(rainfallRate[iParticleInit],0.0);
	airSpecificHumidity[iParticleInit] = std::max(airSpecificHumidity[iParticleInit],0.0);
      } // iParticle

    } // useApplyPhysicalLimits

    // calculate downwelling longwave radiation
    std::string longwaveType;
    configs->get({"ConfigGroup:forcing","Config:longwaveType"}, longwaveType);
    if (longwaveType == "rosati_miyakoda") {

      for (int iParticleInit = 0 ; iParticleInit < tessellation->nParticlesInit ; iParticleInit++) {
	column_longwave_rosati_miyakoda(&longwaveDown[iParticleInit],
					&cloudFraction[iParticleInit],
					&iceConcentration[iParticleInit],
					&iceSurfaceTemperature[iParticleInit],
					&seaSurfaceTemperature[iParticleInit],
					&airSpecificHumidity[iParticleInit],
					&airTemperature[iParticleInit]);
      } // iParticleInit

    } else if (longwaveType == "parkinson_and_washington") {

      for (int iParticleInit = 0 ; iParticleInit < tessellation->nParticlesInit ; iParticleInit++) {
	column_longwave_parkinson_and_washington(&longwaveDown[iParticleInit],
						 &airTemperature[iParticleInit],
						 &cloudFraction[iParticleInit]);
      } // iParticleInit

    } // longwaveType

    // scale precipitation by factor
    std::string precipitationUnits;
    configs->get({"ConfigGroup:forcing","Config:precipitationUnits"}, precipitationUnits);

    // 'mm_per_month', 'mm_per_day','mm_per_sec', or 'mks'
    double precipitationFactor = 1.0;
    if (precipitationUnits == "mm_per_month") {

      const int secondsPerYear = 365 * 24 * 3600;
      precipitationFactor = 12.0 / secondsPerYear;

    } else if (precipitationUnits == "mm_per_day") {

      const int secondsPerDay = 24 * 3600;
      precipitationFactor = 1.0 / secondsPerDay;

    }

    for (int iParticleInit = 0 ; iParticleInit < tessellation->nParticlesInit ; iParticleInit++) {
      precipitationRate[iParticleInit] *= precipitationFactor;
    } // iParticleInit

    // split precipitation between snow and rain from air temperature
    bool useSplitPrecipitation;
    configs->get({"ConfigGroup:forcing","Config:useSplitPrecipitation"}, useSplitPrecipitation);
    if (useSplitPrecipitation) {

      for (int iParticleInit = 0 ; iParticleInit < tessellation->nParticlesInit ; iParticleInit++) {
	column_split_precipitation(&airTemperature[iParticleInit],
				   &precipitationRate[iParticleInit],
				   &snowfallRate[iParticleInit],
				   &rainfallRate[iParticleInit]);
      } // iParticleInit

    } // useSplitPrecipitation

    // post processing of forcing
    char airSpecificHeatStr[7] = "cp_air";
    double airSpecificHeat = column_get_constant(airSpecificHeatStr);
    char latentHeatSublimationStr[5] = "Lsub";
    double latentHeatSublimation = column_get_constant(latentHeatSublimationStr);

    for (int iParticleInit = 0 ; iParticleInit < tessellation->nParticlesInit ; iParticleInit++) {

      // wind speed
      windSpeed[iParticleInit] = std::sqrt(std::pow(uAirVelocity[iParticleInit],2) + std::pow(vAirVelocity[iParticleInit],2));

      // air potential temperature
      airPotentialTemperature[iParticleInit] = airTemperature[iParticleInit];

      // transfer coefficients
      sensibleTransferCoefficient[iParticleInit] = 1.20e-3 * airSpecificHeat       * airDensity[iParticleInit] * windSpeed[iParticleInit];
      latentTransferCoefficient[iParticleInit]   = 1.50e-3 * latentHeatSublimation * airDensity[iParticleInit] * windSpeed[iParticleInit];

      // air stresses
      double airStressCoefficient = 0.0012 * airDensity[iParticleInit] * windSpeed[iParticleInit];

      uAirStress[iParticleInit] = uAirVelocity[iParticleInit] * airStressCoefficient;
      vAirStress[iParticleInit] = vAirVelocity[iParticleInit] * airStressCoefficient;

    } // iParticleInit

    // interpolation of ocean variables to particles
    for (int iParticleInit = 0 ; iParticleInit < tessellation->nParticlesInit ; iParticleInit++) {

      // interpolation
      double x = tessellation->xInit[iParticleInit];
      double y = tessellation->yInit[iParticleInit];

      grid->interpolation_weights(x, y, i, j, weight);

      // oceanic forcing
      interpolate_field(&(forcing->xOcnCurrents), uOceanVelocity[iParticleInit], weight, i, j);
      interpolate_field(&(forcing->yOcnCurrents), vOceanVelocity[iParticleInit], weight, i, j);

    } // iParticleInit

    // ocean
    for (int iParticleInit = 0 ; iParticleInit < tessellation->nParticlesInit ; iParticleInit++) {

      // interpolation
      double x = tessellation->xInit[iParticleInit];
      double y = tessellation->yInit[iParticleInit];

      grid->interpolation_weights(x, y, i, j, weight);

      interpolate_field(&(forcing->seaSurfaceSalinity),       seaSurfaceSalinity[iParticleInit],       weight, i, j);
      interpolate_field(&(forcing->seaSurfaceTiltU),          seaSurfaceTiltU[iParticleInit],          weight, i, j);
      interpolate_field(&(forcing->seaSurfaceTiltV),          seaSurfaceTiltV[iParticleInit],          weight, i, j);
      interpolate_field(&(forcing->oceanMixedLayerDepth),     oceanMixedLayerDepth[iParticleInit],     weight, i, j);
      interpolate_field(&(forcing->oceanHeatFluxConvergence), oceanHeatFluxConvergence[iParticleInit], weight, i, j);

      column_postprocess_ocean_forcing(&seaSurfaceSalinity[iParticleInit],
				       &oceanMixedLayerDepth[iParticleInit],
				       &seaFreezingTemperature[iParticleInit]);

    } // iParticleInit

  } // oceanType

} // Ocean::interpolate_forcing

// update the slab ocean state for one time step
void Ocean::update(void) {

  if (oceanType == "init") {

    double thermoTimeStep = timeStepInterval->seconds();

    for (int iParticleInit = 0 ; iParticleInit < tessellation->nParticlesInit ; iParticleInit++) {

      if (tessellation->typeInit[iParticleInit] == 0) {

	double freezingMeltingPotentialTmp = 0.0;

	column_ocean_mixed_layer(&thermoTimeStep,
				 &seaSurfaceTemperature[iParticleInit],
				 &airPotentialTemperature[iParticleInit],
				 &uAirVelocity[iParticleInit],
				 &vAirVelocity[iParticleInit],
				 &windSpeed[iParticleInit],
				 &airLevelHeight[iParticleInit],
				 &airSpecificHumidity[iParticleInit],
				 &airDensity[iParticleInit],
				 &airStressOceanU[iParticleInit],
				 &airStressOceanV[iParticleInit],
				 &atmosReferenceTemperature2mOcean[iParticleInit],
				 &atmosReferenceHumidity2mOcean[iParticleInit],
				 &airDragCoefficient[iParticleInit],
				 &airOceanDragCoefficientRatio[iParticleInit],
				 &albedoVisibleDirectOcean[iParticleInit],
				 &shortwaveVisibleDirectDown[iParticleInit],
				 &albedoIRDirectOcean[iParticleInit],
				 &shortwaveIRDirectDown[iParticleInit],
				 &albedoVisibleDiffuseOcean[iParticleInit],
				 &shortwaveVisibleDiffuseDown[iParticleInit],
				 &albedoIRDiffuseOcean[iParticleInit],
				 &shortwaveIRDiffuseDown[iParticleInit],
				 &longwaveUpOcean[iParticleInit],
				 &sensibleHeatFluxOcean[iParticleInit],
				 &sensibleTransferCoefficient[iParticleInit],
				 &latentHeatFluxOcean[iParticleInit],
				 &evaporativeWaterFluxOcean[iParticleInit],
				 &longwaveDown[iParticleInit],
				 &iceConcentration[iParticleInit],
				 &oceanHeatFlux[iParticleInit],
				 &oceanShortwaveFlux[iParticleInit],
				 &oceanMixedLayerDepth[iParticleInit],
				 &seaFreezingTemperature[iParticleInit],
				 &oceanHeatFluxConvergence[iParticleInit],
				 &freezingMeltingPotentialTmp);

	// accumulate freezing melting temperature
	freezingMeltingPotential[iParticleInit] += freezingMeltingPotentialTmp;

      } // typeInit

    } // iParticleInit

  } // oceanType == "init"

} // Ocean::update

// update the slab ocean state for one time step
void Ocean::update_remap(const int nParticlesNew, const std::vector<int> iParticleInitFromNew) {

  if (oceanType == "init") {

    double thermoTimeStep = timeStepInterval->seconds();

    for (int iParticleNew = 0 ; iParticleNew < nParticlesNew ; iParticleNew++) {

      int iParticleInit = iParticleInitFromNew[iParticleNew];

      if (tessellation->typeInit[iParticleInit] == 0) {

	double freezingMeltingPotentialTmp = 0.0;

	column_ocean_mixed_layer(&thermoTimeStep,
				 &seaSurfaceTemperature[iParticleInit],
				 &airPotentialTemperature[iParticleInit],
				 &uAirVelocity[iParticleInit],
				 &vAirVelocity[iParticleInit],
				 &windSpeed[iParticleInit],
				 &airLevelHeight[iParticleInit],
				 &airSpecificHumidity[iParticleInit],
				 &airDensity[iParticleInit],
				 &airStressOceanU[iParticleInit],
				 &airStressOceanV[iParticleInit],
				 &atmosReferenceTemperature2mOcean[iParticleInit],
				 &atmosReferenceHumidity2mOcean[iParticleInit],
				 &airDragCoefficient[iParticleInit],
				 &airOceanDragCoefficientRatio[iParticleInit],
				 &albedoVisibleDirectOcean[iParticleInit],
				 &shortwaveVisibleDirectDown[iParticleInit],
				 &albedoIRDirectOcean[iParticleInit],
				 &shortwaveIRDirectDown[iParticleInit],
				 &albedoVisibleDiffuseOcean[iParticleInit],
				 &shortwaveVisibleDiffuseDown[iParticleInit],
				 &albedoIRDiffuseOcean[iParticleInit],
				 &shortwaveIRDiffuseDown[iParticleInit],
				 &longwaveUpOcean[iParticleInit],
				 &sensibleHeatFluxOcean[iParticleInit],
				 &sensibleTransferCoefficient[iParticleInit],
				 &latentHeatFluxOcean[iParticleInit],
				 &evaporativeWaterFluxOcean[iParticleInit],
				 &longwaveDown[iParticleInit],
				 &iceConcentration[iParticleInit],
				 &oceanHeatFlux[iParticleInit],
				 &oceanShortwaveFlux[iParticleInit],
				 &oceanMixedLayerDepth[iParticleInit],
				 &seaFreezingTemperature[iParticleInit],
				 &oceanHeatFluxConvergence[iParticleInit],
				 &freezingMeltingPotentialTmp);

	// accumulate freezing melting temperature
	freezingMeltingPotential[iParticleInit] += freezingMeltingPotentialTmp;

      } // type

    } // iParticleInit

  } // oceanType == "init"

} // Ocean::update_remap

// get all new elements with frazil ice accumulated
std::vector<int> Ocean::get_new_ice_elements(void) {

  std::vector<int> newIceElements;

  if (oceanType == "init") {

    for (int iParticleInit = 0 ; iParticleInit < tessellation->nParticlesInit ; iParticleInit++) {

      if (freezingMeltingPotential[iParticleInit] > 0.0) {
	// new frazil ice formed
	newIceElements.push_back(iParticleInit);
      }

    } // iParticleInit

  } // oceanType

  return newIceElements;

} // Ocean::get_new_ice_elements

// couple ice to ocean
void Ocean::ice_to_ocean(const int nParticlesNew, const std::vector<int> iParticleInitFromNew) {

  if (oceanType == "init") {

    // reset ocean coupling variables
    for (int iParticleInit = 0 ; iParticleInit < tessellation->nParticlesInit ; iParticleInit++) {

      iceConcentration     [iParticleInit] = 0.0;
      iceSurfaceTemperature[iParticleInit] = 0.0;
      oceanHeatFlux        [iParticleInit] = 0.0;
      oceanShortwaveFlux   [iParticleInit] = 0.0;

    } // iParticleInit

    // set ocean coupling variables to ice values
    int nCategories = column->columnDimensions->size("nCategories");

    for (int iParticleNew = 0 ; iParticleNew < nParticlesNew ; iParticleNew++) {

      int iParticleInit = iParticleInitFromNew[iParticleNew];

      iceConcentration[iParticleInit] = 0.0;
      for (int iCategory = 0 ; iCategory < nCategories ; iCategory++) {
	iceConcentration[iParticleInit] += (*column->iceAreaCategory)(iParticleNew,iCategory);
      } // iCategory

      iceSurfaceTemperature[iParticleInit] = (*column->surfaceTemperature)(iParticleNew);
      oceanHeatFlux        [iParticleInit] = (*column->oceanHeatFlux)     (iParticleNew);
      oceanShortwaveFlux   [iParticleInit] = (*column->oceanShortwaveFlux)(iParticleNew);

    } // iParticle

    // reset ocean coupling fluxes
    column->reset_ocean_coupler_fluxes();

  } // oceanType

} // Ocean::ice_to_ocean

// couple the ocean to ice
void Ocean::ocean_to_ice(const int nParticlesNew, const std::vector<int> iParticleInitFromNew) {

  if (oceanType == "init") {

    for (int iParticleNew = 0 ; iParticleNew < nParticlesNew ; iParticleNew++) {

      int iParticleInit = iParticleInitFromNew[iParticleNew];

      // freezing melting potential
      (*column->freezingMeltingPotential)(iParticleNew) = freezingMeltingPotential[iParticleInit];

      // ocean state
      (*column->seaSurfaceTemperature) (iParticleNew) = seaSurfaceTemperature [iParticleInit];
      (*column->seaSurfaceSalinity)    (iParticleNew) = seaSurfaceSalinity    [iParticleInit];
      (*column->seaFreezingTemperature)(iParticleNew) = seaFreezingTemperature[iParticleInit];

      // other
      (*column->oceanMixedLayerDepth)(iParticleNew) = oceanMixedLayerDepth[iParticleInit];

    } // iParticle

    // reset fluxes
    for (int iParticleInit = 0 ; iParticleInit < tessellation->nParticlesInit ; iParticleInit++) {

      freezingMeltingPotential[iParticleInit] = 0.0;

    } // iParticleInit

  } // oceanType

} // Ocean::ocean_to_ice

} // namespace DEMSI
