  // tracers
  DEMSI::ColumnVariable<double>* iceAreaCategory = NULL;
  DEMSI::ColumnVariable<double>* iceVolumeCategory = NULL;
  DEMSI::ColumnVariable<double>* snowVolumeCategory = NULL;
  DEMSI::ColumnVariable<double>* surfaceTemperature = NULL;
  DEMSI::ColumnVariable<double>* iceEnthalpy = NULL;
  DEMSI::ColumnVariable<double>* iceSalinity = NULL;
  DEMSI::ColumnVariable<double>* snowEnthalpy = NULL;
  DEMSI::ColumnVariable<double>* iceAge = NULL;
  DEMSI::ColumnVariable<double>* firstYearIceArea = NULL;
  DEMSI::ColumnVariable<double>* levelIceArea = NULL;
  DEMSI::ColumnVariable<double>* levelIceVolume = NULL;
  DEMSI::ColumnVariable<double>* pondArea = NULL;
  DEMSI::ColumnVariable<double>* pondDepth = NULL;
  DEMSI::ColumnVariable<double>* pondLidThickness = NULL;
  DEMSI::ColumnVariable<double>* snowScatteringAerosol = NULL;
  DEMSI::ColumnVariable<double>* snowBodyAerosol = NULL;
  DEMSI::ColumnVariable<double>* iceScatteringAerosol = NULL;
  DEMSI::ColumnVariable<double>* iceBodyAerosol = NULL;

  // tracers aggregate
  DEMSI::ColumnVariable<double>* iceAreaCell = NULL;
  DEMSI::ColumnVariable<double>* iceVolumeCell = NULL;
  DEMSI::ColumnVariable<double>* snowVolumeCell = NULL;
  DEMSI::ColumnVariable<double>* surfaceTemperatureCell = NULL;
  DEMSI::ColumnVariable<double>* iceEnthalpyCell = NULL;
  DEMSI::ColumnVariable<double>* iceSalinityCell = NULL;
  DEMSI::ColumnVariable<double>* snowEnthalpyCell = NULL;
  DEMSI::ColumnVariable<double>* iceAgeCell = NULL;
  DEMSI::ColumnVariable<double>* firstYearIceAreaCell = NULL;
  DEMSI::ColumnVariable<double>* levelIceAreaCell = NULL;
  DEMSI::ColumnVariable<double>* levelIceVolumeCell = NULL;
  DEMSI::ColumnVariable<double>* pondAreaCell = NULL;
  DEMSI::ColumnVariable<double>* pondDepthCell = NULL;
  DEMSI::ColumnVariable<double>* pondLidThicknessCell = NULL;
  DEMSI::ColumnVariable<double>* snowScatteringAerosolCell = NULL;
  DEMSI::ColumnVariable<double>* snowBodyAerosolCell = NULL;
  DEMSI::ColumnVariable<double>* iceScatteringAerosolCell = NULL;
  DEMSI::ColumnVariable<double>* iceBodyAerosolCell = NULL;

  // icestate
  DEMSI::ColumnVariable<double>* iceAreaCellInitial = NULL;
  DEMSI::ColumnVariable<double>* iceAreaCategoryInitial = NULL;
  DEMSI::ColumnVariable<double>* iceVolumeCategoryInitial = NULL;
  DEMSI::ColumnVariable<double>* iceThicknessCategoryInitial = NULL;
  DEMSI::ColumnVariable<double>* snowVolumeCategoryInitial = NULL;
  DEMSI::ColumnVariable<double>* openWaterArea = NULL;
  DEMSI::ColumnVariable<double>* totalMassCell = NULL;

  // tracer tendencies
  DEMSI::ColumnVariable<double>* iceAreaCategoryTend = NULL;
  DEMSI::ColumnVariable<double>* iceVolumeCategoryTend = NULL;
  DEMSI::ColumnVariable<double>* snowVolumeCategoryTend = NULL;
  DEMSI::ColumnVariable<double>* surfaceTemperatureTend = NULL;
  DEMSI::ColumnVariable<double>* iceEnthalpyTend = NULL;
  DEMSI::ColumnVariable<double>* iceSalinityTend = NULL;
  DEMSI::ColumnVariable<double>* snowEnthalpyTend = NULL;

  // atmos coupling
  DEMSI::ColumnVariable<double>* airLevelHeight = NULL;
  DEMSI::ColumnVariable<double>* airPotentialTemperature = NULL;
  DEMSI::ColumnVariable<double>* airTemperature = NULL;
  DEMSI::ColumnVariable<double>* airSpecificHumidity = NULL;
  DEMSI::ColumnVariable<double>* airDensity = NULL;
  DEMSI::ColumnVariable<double>* shortwaveVisibleDirectDown = NULL;
  DEMSI::ColumnVariable<double>* shortwaveVisibleDiffuseDown = NULL;
  DEMSI::ColumnVariable<double>* shortwaveIRDirectDown = NULL;
  DEMSI::ColumnVariable<double>* shortwaveIRDiffuseDown = NULL;
  DEMSI::ColumnVariable<double>* longwaveDown = NULL;
  DEMSI::ColumnVariable<double>* rainfallRate = NULL;
  DEMSI::ColumnVariable<double>* snowfallRate = NULL;
  DEMSI::ColumnVariable<double>* uAirVelocity = NULL;
  DEMSI::ColumnVariable<double>* vAirVelocity = NULL;
  DEMSI::ColumnVariable<double>* atmosReferenceSpeed10m = NULL;
  DEMSI::ColumnVariable<double>* atmosReferenceTemperature2m = NULL;
  DEMSI::ColumnVariable<double>* atmosReferenceHumidity2m = NULL;

  // atmos forcing
  DEMSI::ColumnVariable<double>* windSpeed = NULL;
  DEMSI::ColumnVariable<double>* uAirStress = NULL;
  DEMSI::ColumnVariable<double>* vAirStress = NULL;
  DEMSI::ColumnVariable<double>* shortwaveDown = NULL;
  DEMSI::ColumnVariable<double>* cloudFraction = NULL;
  DEMSI::ColumnVariable<double>* sensibleTransferCoefficient = NULL;
  DEMSI::ColumnVariable<double>* latentTransferCoefficient = NULL;
  DEMSI::ColumnVariable<double>* precipitationRate = NULL;

  // alternative atmos forcing
  DEMSI::ColumnVariable<double>* latentHeatFluxCouple = NULL;
  DEMSI::ColumnVariable<double>* sensibleHeatFluxCouple = NULL;
  DEMSI::ColumnVariable<double>* surfaceHeatFluxCouple = NULL;
  DEMSI::ColumnVariable<double>* surfaceConductiveFluxCouple = NULL;
  DEMSI::ColumnVariable<double>* airStressForcingU = NULL;
  DEMSI::ColumnVariable<double>* airStressForcingV = NULL;

  // ocean coupling
  DEMSI::ColumnVariable<double>* seaSurfaceTemperature = NULL;
  DEMSI::ColumnVariable<double>* seaSurfaceSalinity = NULL;
  DEMSI::ColumnVariable<double>* seaFreezingTemperature = NULL;
  DEMSI::ColumnVariable<double>* freezingMeltingPotential = NULL;
  DEMSI::ColumnVariable<double>* uOceanVelocity = NULL;
  DEMSI::ColumnVariable<double>* vOceanVelocity = NULL;
  DEMSI::ColumnVariable<double>* seaSurfaceTiltU = NULL;
  DEMSI::ColumnVariable<double>* seaSurfaceTiltV = NULL;
  DEMSI::ColumnVariable<double>* oceanMixedLayerDepth = NULL;
  DEMSI::ColumnVariable<double>* oceanHeatFluxConvergence = NULL;
  DEMSI::ColumnVariable<int>* landIceMask = NULL;

  // ridging
  DEMSI::ColumnVariable<double>* ridgeConvergence = NULL;
  DEMSI::ColumnVariable<double>* ridgeShear = NULL;
  DEMSI::ColumnVariable<double>* areaLossRidge = NULL;
  DEMSI::ColumnVariable<double>* areaGainRidge = NULL;
  DEMSI::ColumnVariable<double>* iceVolumeRidged = NULL;
  DEMSI::ColumnVariable<double>* openingRateRidge = NULL;
  DEMSI::ColumnVariable<double>* ridgeParticipationFunction = NULL;
  DEMSI::ColumnVariable<double>* ratioRidgeThicknessToIce = NULL;
  DEMSI::ColumnVariable<double>* fractionNewRidgeArea = NULL;
  DEMSI::ColumnVariable<double>* fractionNewRidgeVolume = NULL;
  DEMSI::ColumnVariable<double>* areaLossRidgeCategory = NULL;
  DEMSI::ColumnVariable<double>* areaGainRidgeCategory = NULL;
  DEMSI::ColumnVariable<double>* iceVolumeRidgedCategory = NULL;
  DEMSI::ColumnVariable<double>* raftingIceArea = NULL;
  DEMSI::ColumnVariable<double>* raftingIceVolume = NULL;

  // melt growth rates
  DEMSI::ColumnVariable<double>* lateralIceMeltFraction = NULL;
  DEMSI::ColumnVariable<double>* surfaceIceMelt = NULL;
  DEMSI::ColumnVariable<double>* surfaceIceMeltCategory = NULL;
  DEMSI::ColumnVariable<double>* basalIceMelt = NULL;
  DEMSI::ColumnVariable<double>* basalIceMeltCategory = NULL;
  DEMSI::ColumnVariable<double>* lateralIceMelt = NULL;
  DEMSI::ColumnVariable<double>* snowMelt = NULL;
  DEMSI::ColumnVariable<double>* snowMeltCategory = NULL;
  DEMSI::ColumnVariable<double>* congelation = NULL;
  DEMSI::ColumnVariable<double>* congelationCategory = NULL;
  DEMSI::ColumnVariable<double>* snowiceFormation = NULL;
  DEMSI::ColumnVariable<double>* snowiceFormationCategory = NULL;
  DEMSI::ColumnVariable<double>* snowThicknessChange = NULL;
  DEMSI::ColumnVariable<double>* snowThicknessChangeCategory = NULL;
  DEMSI::ColumnVariable<double>* frazilFormation = NULL;
  DEMSI::ColumnVariable<double>* frazilGrowthDiagnostic = NULL;

  // atmos fluxes
  DEMSI::ColumnVariable<double>* surfaceHeatFlux = NULL;
  DEMSI::ColumnVariable<double>* surfaceHeatFluxCategory = NULL;
  DEMSI::ColumnVariable<double>* surfaceConductiveFlux = NULL;
  DEMSI::ColumnVariable<double>* surfaceConductiveFluxCategory = NULL;
  DEMSI::ColumnVariable<double>* longwaveUp = NULL;
  DEMSI::ColumnVariable<double>* sensibleHeatFlux = NULL;
  DEMSI::ColumnVariable<double>* sensibleHeatFluxCategory = NULL;
  DEMSI::ColumnVariable<double>* latentHeatFlux = NULL;
  DEMSI::ColumnVariable<double>* latentHeatFluxCategory = NULL;
  DEMSI::ColumnVariable<double>* evaporativeWaterFlux = NULL;

  // ocean fluxes
  DEMSI::ColumnVariable<double>* oceanFreshWaterFlux = NULL;
  DEMSI::ColumnVariable<double>* oceanSaltFlux = NULL;
  DEMSI::ColumnVariable<double>* oceanHeatFlux = NULL;
  DEMSI::ColumnVariable<double>* oceanShortwaveFlux = NULL;
  DEMSI::ColumnVariable<double>* oceanFreshWaterFluxArea = NULL;
  DEMSI::ColumnVariable<double>* oceanSaltFluxArea = NULL;
  DEMSI::ColumnVariable<double>* oceanHeatFluxArea = NULL;
  DEMSI::ColumnVariable<double>* oceanShortwaveFluxArea = NULL;
  DEMSI::ColumnVariable<double>* oceanHeatFluxIceBottom = NULL;

  // ocean atmosphere
  DEMSI::ColumnVariable<double>* airStressOceanU = NULL;
  DEMSI::ColumnVariable<double>* airStressOceanV = NULL;
  DEMSI::ColumnVariable<double>* atmosReferenceTemperature2mOcean = NULL;
  DEMSI::ColumnVariable<double>* atmosReferenceHumidity2mOcean = NULL;
  DEMSI::ColumnVariable<double>* albedoVisibleDirectOcean = NULL;
  DEMSI::ColumnVariable<double>* albedoVisibleDiffuseOcean = NULL;
  DEMSI::ColumnVariable<double>* albedoIRDirectOcean = NULL;
  DEMSI::ColumnVariable<double>* albedoIRDiffuseOcean = NULL;
  DEMSI::ColumnVariable<double>* longwaveUpOcean = NULL;
  DEMSI::ColumnVariable<double>* sensibleHeatFluxOcean = NULL;
  DEMSI::ColumnVariable<double>* latentHeatFluxOcean = NULL;
  DEMSI::ColumnVariable<double>* evaporativeWaterFluxOcean = NULL;

  // shortwave
  DEMSI::ColumnVariable<double>* solarZenithAngleCosine = NULL;
  DEMSI::ColumnVariable<double>* albedoVisibleDirectCategory = NULL;
  DEMSI::ColumnVariable<double>* albedoVisibleDiffuseCategory = NULL;
  DEMSI::ColumnVariable<double>* albedoIRDirectCategory = NULL;
  DEMSI::ColumnVariable<double>* albedoIRDiffuseCategory = NULL;
  DEMSI::ColumnVariable<double>* albedoVisibleDirectCell = NULL;
  DEMSI::ColumnVariable<double>* albedoVisibleDiffuseCell = NULL;
  DEMSI::ColumnVariable<double>* albedoIRDirectCell = NULL;
  DEMSI::ColumnVariable<double>* albedoIRDiffuseCell = NULL;
  DEMSI::ColumnVariable<double>* albedoVisibleDirectArea = NULL;
  DEMSI::ColumnVariable<double>* albedoVisibleDiffuseArea = NULL;
  DEMSI::ColumnVariable<double>* albedoIRDirectArea = NULL;
  DEMSI::ColumnVariable<double>* albedoIRDiffuseArea = NULL;
  DEMSI::ColumnVariable<double>* shortwaveScalingFactor = NULL;
  DEMSI::ColumnVariable<double>* surfaceShortwaveFlux = NULL;
  DEMSI::ColumnVariable<double>* interiorShortwaveFlux = NULL;
  DEMSI::ColumnVariable<double>* penetratingShortwaveFlux = NULL;
  DEMSI::ColumnVariable<double>* shortwaveLayerPenetration = NULL;
  DEMSI::ColumnVariable<double>* absorbedShortwaveFlux = NULL;
  DEMSI::ColumnVariable<double>* absorbedShortwaveSnowLayer = NULL;
  DEMSI::ColumnVariable<double>* absorbedShortwaveIceLayer = NULL;
  DEMSI::ColumnVariable<double>* bareIceAlbedoCategory = NULL;
  DEMSI::ColumnVariable<double>* snowAlbedoCategory = NULL;
  DEMSI::ColumnVariable<double>* pondAlbedoCategory = NULL;
  DEMSI::ColumnVariable<double>* bareIceAlbedoCell = NULL;
  DEMSI::ColumnVariable<double>* snowAlbedoCell = NULL;
  DEMSI::ColumnVariable<double>* pondAlbedoCell = NULL;
  DEMSI::ColumnVariable<double>* effectivePondAreaCategory = NULL;
  DEMSI::ColumnVariable<double>* effectivePondAreaCell = NULL;
  DEMSI::ColumnVariable<double>* snowFractionCategory = NULL;

  // drag
  DEMSI::ColumnVariable<double>* airOceanDragCoefficientRatio = NULL;
  DEMSI::ColumnVariable<double>* oceanDragCoefficient = NULL;
  DEMSI::ColumnVariable<double>* oceanDragCoefficientSkin = NULL;
  DEMSI::ColumnVariable<double>* oceanDragCoefficientFloe = NULL;
  DEMSI::ColumnVariable<double>* oceanDragCoefficientKeel = NULL;
  DEMSI::ColumnVariable<double>* airDragCoefficient = NULL;
  DEMSI::ColumnVariable<double>* airDragCoefficientSkin = NULL;
  DEMSI::ColumnVariable<double>* airDragCoefficientFloe = NULL;
  DEMSI::ColumnVariable<double>* airDragCoefficientPond = NULL;
  DEMSI::ColumnVariable<double>* airDragCoefficientRidge = NULL;
  DEMSI::ColumnVariable<double>* dragFreeboard = NULL;
  DEMSI::ColumnVariable<double>* dragIceSnowDraft = NULL;
  DEMSI::ColumnVariable<double>* dragRidgeHeight = NULL;
  DEMSI::ColumnVariable<double>* dragRidgeSeparation = NULL;
  DEMSI::ColumnVariable<double>* dragKeelDepth = NULL;
  DEMSI::ColumnVariable<double>* dragKeelSeparation = NULL;
  DEMSI::ColumnVariable<double>* dragFloeLength = NULL;
  DEMSI::ColumnVariable<double>* dragFloeSeparation = NULL;

  // ponds
  DEMSI::ColumnVariable<double>* pondFreshWaterFlux = NULL;
  DEMSI::ColumnVariable<double>* pondSnowDepthDifference = NULL;
  DEMSI::ColumnVariable<double>* pondLidMeltFluxFraction = NULL;

  // aerosols
  DEMSI::ColumnVariable<double>* atmosAerosolFlux = NULL;
  DEMSI::ColumnVariable<double>* oceanAerosolFlux = NULL;

  // initial
  DEMSI::ColumnVariable<double>* initialSalinityProfile = NULL;
  DEMSI::ColumnVariable<double>* initialMeltingTemperatureProfile = NULL;

  // dynamics
  DEMSI::ColumnVariable<double>* airStressCellU = NULL;
  DEMSI::ColumnVariable<double>* airStressCellV = NULL;
  DEMSI::ColumnVariable<double>* oceanStressCellU = NULL;
  DEMSI::ColumnVariable<double>* oceanStressCellV = NULL;

  // diagnostics
  DEMSI::ColumnVariable<double>* meltOnset = NULL;
  DEMSI::ColumnVariable<double>* freezeOnset = NULL;
  DEMSI::ColumnVariable<double>* freezingMeltingPotentialInitial = NULL;
  DEMSI::ColumnVariable<double>* iceAreaTendencyThermodynamics = NULL;
  DEMSI::ColumnVariable<double>* iceVolumeTendencyThermodynamics = NULL;
  DEMSI::ColumnVariable<double>* iceAgeTendencyThermodynamics = NULL;

  // location
  DEMSI::ColumnVariable<double>* latitude = NULL;
  DEMSI::ColumnVariable<double>* longitude = NULL;

  // initial particle polygon
  DEMSI::ColumnVariable<int>* nVerticesInit = NULL;
  DEMSI::ColumnVariable<double>* verticesInit = NULL;
  DEMSI::ColumnVariable<double>* areaInit = NULL;
  DEMSI::ColumnVariable<double>* areaInitRemap = NULL;

  // element areas
  DEMSI::ColumnVariable<double>* effectiveElementArea = NULL;
