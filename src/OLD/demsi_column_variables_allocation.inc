  // tracers
  surfaceTemperature = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"surfaceTemperature",
	{"nCategories"});
  columnVariables->add(surfaceTemperature);

  iceEnthalpy = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"iceEnthalpy",
	{"nCategories", "nIceLayers"});
  columnVariables->add(iceEnthalpy);

  iceSalinity = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"iceSalinity",
	{"nCategories", "nIceLayers"});
  columnVariables->add(iceSalinity);

  snowEnthalpy = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"snowEnthalpy",
	{"nCategories", "nSnowLayers"});
  columnVariables->add(snowEnthalpy);

  iceAge = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"iceAge",
	{"nCategories"},
	useIceAge);
  columnVariables->add(iceAge);

  firstYearIceArea = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"firstYearIceArea",
	{"nCategories"},
	useFirstYearIce);
  columnVariables->add(firstYearIceArea);

  levelIceArea = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"levelIceArea",
	{"nCategories"},
	useLevelIce);
  columnVariables->add(levelIceArea);

  levelIceVolume = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"levelIceVolume",
	{"nCategories"},
	useLevelIce);
  columnVariables->add(levelIceVolume);

  pondArea = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"pondArea",
	{"nCategories"},
	usePonds);
  columnVariables->add(pondArea);

  pondDepth = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"pondDepth",
	{"nCategories"},
	usePonds);
  columnVariables->add(pondDepth);

  pondLidThickness = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"pondLidThickness",
	{"nCategories"},
	useLidThickness);
  columnVariables->add(pondLidThickness);

  snowScatteringAerosol = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"snowScatteringAerosol",
	{"nCategories", "nAerosols"},
	useAerosols);
  columnVariables->add(snowScatteringAerosol);

  snowBodyAerosol = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"snowBodyAerosol",
	{"nCategories", "nAerosols"},
	useAerosols);
  columnVariables->add(snowBodyAerosol);

  iceScatteringAerosol = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"iceScatteringAerosol",
	{"nCategories", "nAerosols"},
	useAerosols);
  columnVariables->add(iceScatteringAerosol);

  iceBodyAerosol = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"iceBodyAerosol",
	{"nCategories", "nAerosols"},
	useAerosols);
  columnVariables->add(iceBodyAerosol);


  // tracers aggregate
  surfaceTemperatureCell = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"surfaceTemperatureCell",
	{});
  columnVariables->add(surfaceTemperatureCell);

  iceEnthalpyCell = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"iceEnthalpyCell",
	{"nIceLayers"});
  columnVariables->add(iceEnthalpyCell);

  iceSalinityCell = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"iceSalinityCell",
	{"nIceLayers"});
  columnVariables->add(iceSalinityCell);

  snowEnthalpyCell = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"snowEnthalpyCell",
	{"nSnowLayers"});
  columnVariables->add(snowEnthalpyCell);

  iceAgeCell = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"iceAgeCell",
	{},
	useIceAge);
  columnVariables->add(iceAgeCell);

  firstYearIceAreaCell = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"firstYearIceAreaCell",
	{},
	useFirstYearIce);
  columnVariables->add(firstYearIceAreaCell);

  levelIceAreaCell = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"levelIceAreaCell",
	{},
	useLevelIce);
  columnVariables->add(levelIceAreaCell);

  levelIceVolumeCell = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"levelIceVolumeCell",
	{},
	useLevelIce);
  columnVariables->add(levelIceVolumeCell);

  pondAreaCell = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"pondAreaCell",
	{},
	usePonds);
  columnVariables->add(pondAreaCell);

  pondDepthCell = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"pondDepthCell",
	{},
	usePonds);
  columnVariables->add(pondDepthCell);

  pondLidThicknessCell = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"pondLidThicknessCell",
	{},
	useLidThickness);
  columnVariables->add(pondLidThicknessCell);

  snowScatteringAerosolCell = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"snowScatteringAerosolCell",
	{"nAerosols"},
	useAerosols);
  columnVariables->add(snowScatteringAerosolCell);

  snowBodyAerosolCell = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"snowBodyAerosolCell",
	{"nAerosols"},
	useAerosols);
  columnVariables->add(snowBodyAerosolCell);

  iceScatteringAerosolCell = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"iceScatteringAerosolCell",
	{"nAerosols"},
	useAerosols);
  columnVariables->add(iceScatteringAerosolCell);

  iceBodyAerosolCell = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"iceBodyAerosolCell",
	{"nAerosols"},
	useAerosols);
  columnVariables->add(iceBodyAerosolCell);


  // icestate
  iceAreaCellInitial = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"iceAreaCellInitial",
	{});
  columnVariables->add(iceAreaCellInitial);

  iceAreaCategoryInitial = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"iceAreaCategoryInitial",
	{"nCategories"});
  columnVariables->add(iceAreaCategoryInitial);

  iceVolumeCategoryInitial = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"iceVolumeCategoryInitial",
	{"nCategories"});
  columnVariables->add(iceVolumeCategoryInitial);

  iceThicknessCategoryInitial = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"iceThicknessCategoryInitial",
	{"nCategories"});
  columnVariables->add(iceThicknessCategoryInitial);

  snowVolumeCategoryInitial = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"snowVolumeCategoryInitial",
	{"nCategories"});
  columnVariables->add(snowVolumeCategoryInitial);

  openWaterArea = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"openWaterArea",
	{});
  columnVariables->add(openWaterArea);

  totalMassCell = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"totalMassCell",
	{});
  columnVariables->add(totalMassCell);


  // tracer tendencies
  iceAreaCategoryTend = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"iceAreaCategoryTend",
	{"nCategories"});
  columnVariables->add(iceAreaCategoryTend);

  iceVolumeCategoryTend = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"iceVolumeCategoryTend",
	{"nCategories"});
  columnVariables->add(iceVolumeCategoryTend);

  snowVolumeCategoryTend = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"snowVolumeCategoryTend",
	{"nCategories"});
  columnVariables->add(snowVolumeCategoryTend);

  surfaceTemperatureTend = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"surfaceTemperatureTend",
	{"nCategories"});
  columnVariables->add(surfaceTemperatureTend);

  iceEnthalpyTend = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"iceEnthalpyTend",
	{"nCategories", "nIceLayers"});
  columnVariables->add(iceEnthalpyTend);

  iceSalinityTend = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"iceSalinityTend",
	{"nCategories", "nIceLayers"});
  columnVariables->add(iceSalinityTend);

  snowEnthalpyTend = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"snowEnthalpyTend",
	{"nCategories", "nSnowLayers"});
  columnVariables->add(snowEnthalpyTend);


  // atmos coupling
  airLevelHeight = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"airLevelHeight",
	{});
  columnVariables->add(airLevelHeight);

  airPotentialTemperature = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"airPotentialTemperature",
	{});
  columnVariables->add(airPotentialTemperature);

  airTemperature = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"airTemperature",
	{});
  columnVariables->add(airTemperature);

  airSpecificHumidity = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"airSpecificHumidity",
	{});
  columnVariables->add(airSpecificHumidity);

  airDensity = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"airDensity",
	{});
  columnVariables->add(airDensity);

  shortwaveVisibleDirectDown = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"shortwaveVisibleDirectDown",
	{});
  columnVariables->add(shortwaveVisibleDirectDown);

  shortwaveVisibleDiffuseDown = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"shortwaveVisibleDiffuseDown",
	{});
  columnVariables->add(shortwaveVisibleDiffuseDown);

  shortwaveIRDirectDown = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"shortwaveIRDirectDown",
	{});
  columnVariables->add(shortwaveIRDirectDown);

  shortwaveIRDiffuseDown = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"shortwaveIRDiffuseDown",
	{});
  columnVariables->add(shortwaveIRDiffuseDown);

  longwaveDown = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"longwaveDown",
	{});
  columnVariables->add(longwaveDown);

  rainfallRate = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"rainfallRate",
	{});
  columnVariables->add(rainfallRate);

  snowfallRate = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"snowfallRate",
	{});
  columnVariables->add(snowfallRate);

  uAirVelocity = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"uAirVelocity",
	{});
  columnVariables->add(uAirVelocity);

  vAirVelocity = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"vAirVelocity",
	{});
  columnVariables->add(vAirVelocity);

  atmosReferenceSpeed10m = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"atmosReferenceSpeed10m",
	{});
  columnVariables->add(atmosReferenceSpeed10m);

  atmosReferenceTemperature2m = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"atmosReferenceTemperature2m",
	{});
  columnVariables->add(atmosReferenceTemperature2m);

  atmosReferenceHumidity2m = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"atmosReferenceHumidity2m",
	{});
  columnVariables->add(atmosReferenceHumidity2m);


  // atmos forcing
  windSpeed = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"windSpeed",
	{});
  columnVariables->add(windSpeed);

  uAirStress = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"uAirStress",
	{});
  columnVariables->add(uAirStress);

  vAirStress = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"vAirStress",
	{});
  columnVariables->add(vAirStress);

  shortwaveDown = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"shortwaveDown",
	{});
  columnVariables->add(shortwaveDown);

  cloudFraction = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"cloudFraction",
	{});
  columnVariables->add(cloudFraction);

  sensibleTransferCoefficient = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"sensibleTransferCoefficient",
	{});
  columnVariables->add(sensibleTransferCoefficient);

  latentTransferCoefficient = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"latentTransferCoefficient",
	{});
  columnVariables->add(latentTransferCoefficient);


  // alternative atmos forcing
  latentHeatFluxCouple = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"latentHeatFluxCouple",
	{"nCategories"});
  columnVariables->add(latentHeatFluxCouple);

  sensibleHeatFluxCouple = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"sensibleHeatFluxCouple",
	{"nCategories"});
  columnVariables->add(sensibleHeatFluxCouple);

  surfaceHeatFluxCouple = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"surfaceHeatFluxCouple",
	{"nCategories"});
  columnVariables->add(surfaceHeatFluxCouple);

  surfaceConductiveFluxCouple = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"surfaceConductiveFluxCouple",
	{"nCategories"});
  columnVariables->add(surfaceConductiveFluxCouple);

  airStressForcingU = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"airStressForcingU",
	{});
  columnVariables->add(airStressForcingU);

  airStressForcingV = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"airStressForcingV",
	{});
  columnVariables->add(airStressForcingV);


  // ocean coupling
  seaSurfaceTemperature = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"seaSurfaceTemperature",
	{});
  columnVariables->add(seaSurfaceTemperature);

  seaSurfaceSalinity = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"seaSurfaceSalinity",
	{});
  columnVariables->add(seaSurfaceSalinity);

  seaFreezingTemperature = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"seaFreezingTemperature",
	{});
  columnVariables->add(seaFreezingTemperature);

  freezingMeltingPotential = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"freezingMeltingPotential",
	{});
  columnVariables->add(freezingMeltingPotential);

  uOceanVelocity = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"uOceanVelocity",
	{});
  columnVariables->add(uOceanVelocity);

  vOceanVelocity = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"vOceanVelocity",
	{});
  columnVariables->add(vOceanVelocity);

  seaSurfaceTiltU = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"seaSurfaceTiltU",
	{});
  columnVariables->add(seaSurfaceTiltU);

  seaSurfaceTiltV = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"seaSurfaceTiltV",
	{});
  columnVariables->add(seaSurfaceTiltV);

  oceanMixedLayerDepth = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"oceanMixedLayerDepth",
	{});
  columnVariables->add(oceanMixedLayerDepth);

  oceanHeatFluxConvergence = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"oceanHeatFluxConvergence",
	{});
  columnVariables->add(oceanHeatFluxConvergence);

  landIceMask = new DEMSI::ColumnVariable<int>(particles, columnDimensions,
	"landIceMask",
	{});
  columnVariables->add(landIceMask);


  // ridging
  ridgeConvergence = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"ridgeConvergence",
	{});
  columnVariables->add(ridgeConvergence);

  ridgeShear = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"ridgeShear",
	{});
  columnVariables->add(ridgeShear);

  areaLossRidge = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"areaLossRidge",
	{});
  columnVariables->add(areaLossRidge);

  areaGainRidge = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"areaGainRidge",
	{});
  columnVariables->add(areaGainRidge);

  iceVolumeRidged = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"iceVolumeRidged",
	{});
  columnVariables->add(iceVolumeRidged);

  openingRateRidge = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"openingRateRidge",
	{});
  columnVariables->add(openingRateRidge);

  ridgeParticipationFunction = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"ridgeParticipationFunction",
	{"nCategories"});
  columnVariables->add(ridgeParticipationFunction);

  ratioRidgeThicknessToIce = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"ratioRidgeThicknessToIce",
	{"nCategories"});
  columnVariables->add(ratioRidgeThicknessToIce);

  fractionNewRidgeArea = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"fractionNewRidgeArea",
	{"nCategories"});
  columnVariables->add(fractionNewRidgeArea);

  fractionNewRidgeVolume = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"fractionNewRidgeVolume",
	{"nCategories"});
  columnVariables->add(fractionNewRidgeVolume);

  areaLossRidgeCategory = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"areaLossRidgeCategory",
	{"nCategories"});
  columnVariables->add(areaLossRidgeCategory);

  areaGainRidgeCategory = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"areaGainRidgeCategory",
	{"nCategories"});
  columnVariables->add(areaGainRidgeCategory);

  iceVolumeRidgedCategory = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"iceVolumeRidgedCategory",
	{"nCategories"});
  columnVariables->add(iceVolumeRidgedCategory);

  raftingIceArea = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"raftingIceArea",
	{"nCategories"});
  columnVariables->add(raftingIceArea);

  raftingIceVolume = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"raftingIceVolume",
	{"nCategories"});
  columnVariables->add(raftingIceVolume);


  // melt growth rates
  lateralIceMeltFraction = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"lateralIceMeltFraction",
	{});
  columnVariables->add(lateralIceMeltFraction);

  surfaceIceMelt = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"surfaceIceMelt",
	{});
  columnVariables->add(surfaceIceMelt);

  surfaceIceMeltCategory = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"surfaceIceMeltCategory",
	{"nCategories"});
  columnVariables->add(surfaceIceMeltCategory);

  basalIceMelt = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"basalIceMelt",
	{});
  columnVariables->add(basalIceMelt);

  basalIceMeltCategory = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"basalIceMeltCategory",
	{"nCategories"});
  columnVariables->add(basalIceMeltCategory);

  lateralIceMelt = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"lateralIceMelt",
	{});
  columnVariables->add(lateralIceMelt);

  snowMelt = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"snowMelt",
	{});
  columnVariables->add(snowMelt);

  snowMeltCategory = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"snowMeltCategory",
	{"nCategories"});
  columnVariables->add(snowMeltCategory);

  congelation = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"congelation",
	{});
  columnVariables->add(congelation);

  congelationCategory = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"congelationCategory",
	{"nCategories"});
  columnVariables->add(congelationCategory);

  snowiceFormation = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"snowiceFormation",
	{});
  columnVariables->add(snowiceFormation);

  snowiceFormationCategory = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"snowiceFormationCategory",
	{"nCategories"});
  columnVariables->add(snowiceFormationCategory);

  snowThicknessChange = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"snowThicknessChange",
	{});
  columnVariables->add(snowThicknessChange);

  snowThicknessChangeCategory = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"snowThicknessChangeCategory",
	{"nCategories"});
  columnVariables->add(snowThicknessChangeCategory);

  frazilFormation = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"frazilFormation",
	{});
  columnVariables->add(frazilFormation);

  frazilGrowthDiagnostic = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"frazilGrowthDiagnostic",
	{});
  columnVariables->add(frazilGrowthDiagnostic);


  // atmos fluxes
  surfaceHeatFlux = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"surfaceHeatFlux",
	{});
  columnVariables->add(surfaceHeatFlux);

  surfaceHeatFluxCategory = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"surfaceHeatFluxCategory",
	{"nCategories"});
  columnVariables->add(surfaceHeatFluxCategory);

  surfaceConductiveFlux = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"surfaceConductiveFlux",
	{});
  columnVariables->add(surfaceConductiveFlux);

  surfaceConductiveFluxCategory = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"surfaceConductiveFluxCategory",
	{"nCategories"});
  columnVariables->add(surfaceConductiveFluxCategory);

  longwaveUp = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"longwaveUp",
	{});
  columnVariables->add(longwaveUp);

  sensibleHeatFlux = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"sensibleHeatFlux",
	{});
  columnVariables->add(sensibleHeatFlux);

  sensibleHeatFluxCategory = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"sensibleHeatFluxCategory",
	{"nCategories"});
  columnVariables->add(sensibleHeatFluxCategory);

  latentHeatFlux = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"latentHeatFlux",
	{});
  columnVariables->add(latentHeatFlux);

  latentHeatFluxCategory = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"latentHeatFluxCategory",
	{"nCategories"});
  columnVariables->add(latentHeatFluxCategory);

  evaporativeWaterFlux = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"evaporativeWaterFlux",
	{});
  columnVariables->add(evaporativeWaterFlux);


  // ocean fluxes
  oceanFreshWaterFlux = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"oceanFreshWaterFlux",
	{});
  columnVariables->add(oceanFreshWaterFlux);

  oceanSaltFlux = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"oceanSaltFlux",
	{});
  columnVariables->add(oceanSaltFlux);

  oceanHeatFlux = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"oceanHeatFlux",
	{});
  columnVariables->add(oceanHeatFlux);

  oceanShortwaveFlux = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"oceanShortwaveFlux",
	{});
  columnVariables->add(oceanShortwaveFlux);

  oceanFreshWaterFluxArea = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"oceanFreshWaterFluxArea",
	{});
  columnVariables->add(oceanFreshWaterFluxArea);

  oceanSaltFluxArea = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"oceanSaltFluxArea",
	{});
  columnVariables->add(oceanSaltFluxArea);

  oceanHeatFluxArea = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"oceanHeatFluxArea",
	{});
  columnVariables->add(oceanHeatFluxArea);

  oceanShortwaveFluxArea = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"oceanShortwaveFluxArea",
	{});
  columnVariables->add(oceanShortwaveFluxArea);

  oceanHeatFluxIceBottom = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"oceanHeatFluxIceBottom",
	{});
  columnVariables->add(oceanHeatFluxIceBottom);


  // ocean atmosphere
  airStressOceanU = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"airStressOceanU",
	{});
  columnVariables->add(airStressOceanU);

  airStressOceanV = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"airStressOceanV",
	{});
  columnVariables->add(airStressOceanV);

  atmosReferenceTemperature2mOcean = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"atmosReferenceTemperature2mOcean",
	{});
  columnVariables->add(atmosReferenceTemperature2mOcean);

  atmosReferenceHumidity2mOcean = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"atmosReferenceHumidity2mOcean",
	{});
  columnVariables->add(atmosReferenceHumidity2mOcean);

  albedoVisibleDirectOcean = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"albedoVisibleDirectOcean",
	{});
  columnVariables->add(albedoVisibleDirectOcean);

  albedoVisibleDiffuseOcean = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"albedoVisibleDiffuseOcean",
	{});
  columnVariables->add(albedoVisibleDiffuseOcean);

  albedoIRDirectOcean = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"albedoIRDirectOcean",
	{});
  columnVariables->add(albedoIRDirectOcean);

  albedoIRDiffuseOcean = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"albedoIRDiffuseOcean",
	{});
  columnVariables->add(albedoIRDiffuseOcean);

  longwaveUpOcean = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"longwaveUpOcean",
	{});
  columnVariables->add(longwaveUpOcean);

  sensibleHeatFluxOcean = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"sensibleHeatFluxOcean",
	{});
  columnVariables->add(sensibleHeatFluxOcean);

  latentHeatFluxOcean = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"latentHeatFluxOcean",
	{});
  columnVariables->add(latentHeatFluxOcean);

  evaporativeWaterFluxOcean = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"evaporativeWaterFluxOcean",
	{});
  columnVariables->add(evaporativeWaterFluxOcean);


  // shortwave
  solarZenithAngleCosine = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"solarZenithAngleCosine",
	{});
  columnVariables->add(solarZenithAngleCosine);

  albedoVisibleDirectCategory = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"albedoVisibleDirectCategory",
	{"nCategories"});
  columnVariables->add(albedoVisibleDirectCategory);

  albedoVisibleDiffuseCategory = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"albedoVisibleDiffuseCategory",
	{"nCategories"});
  columnVariables->add(albedoVisibleDiffuseCategory);

  albedoIRDirectCategory = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"albedoIRDirectCategory",
	{"nCategories"});
  columnVariables->add(albedoIRDirectCategory);

  albedoIRDiffuseCategory = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"albedoIRDiffuseCategory",
	{"nCategories"});
  columnVariables->add(albedoIRDiffuseCategory);

  albedoVisibleDirectCell = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"albedoVisibleDirectCell",
	{});
  columnVariables->add(albedoVisibleDirectCell);

  albedoVisibleDiffuseCell = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"albedoVisibleDiffuseCell",
	{});
  columnVariables->add(albedoVisibleDiffuseCell);

  albedoIRDirectCell = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"albedoIRDirectCell",
	{});
  columnVariables->add(albedoIRDirectCell);

  albedoIRDiffuseCell = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"albedoIRDiffuseCell",
	{});
  columnVariables->add(albedoIRDiffuseCell);

  albedoVisibleDirectArea = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"albedoVisibleDirectArea",
	{});
  columnVariables->add(albedoVisibleDirectArea);

  albedoVisibleDiffuseArea = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"albedoVisibleDiffuseArea",
	{});
  columnVariables->add(albedoVisibleDiffuseArea);

  albedoIRDirectArea = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"albedoIRDirectArea",
	{});
  columnVariables->add(albedoIRDirectArea);

  albedoIRDiffuseArea = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"albedoIRDiffuseArea",
	{});
  columnVariables->add(albedoIRDiffuseArea);

  shortwaveScalingFactor = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"shortwaveScalingFactor",
	{});
  columnVariables->add(shortwaveScalingFactor);

  surfaceShortwaveFlux = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"surfaceShortwaveFlux",
	{"nCategories"});
  columnVariables->add(surfaceShortwaveFlux);

  interiorShortwaveFlux = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"interiorShortwaveFlux",
	{"nCategories"});
  columnVariables->add(interiorShortwaveFlux);

  penetratingShortwaveFlux = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"penetratingShortwaveFlux",
	{"nCategories"});
  columnVariables->add(penetratingShortwaveFlux);

  shortwaveLayerPenetration = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"shortwaveLayerPenetration",
	{"nCategories", "nIceLayersP1"});
  columnVariables->add(shortwaveLayerPenetration);

  absorbedShortwaveFlux = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"absorbedShortwaveFlux",
	{});
  columnVariables->add(absorbedShortwaveFlux);

  absorbedShortwaveSnowLayer = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"absorbedShortwaveSnowLayer",
	{"nCategories", "nSnowLayers"});
  columnVariables->add(absorbedShortwaveSnowLayer);

  absorbedShortwaveIceLayer = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"absorbedShortwaveIceLayer",
	{"nCategories", "nIceLayers"});
  columnVariables->add(absorbedShortwaveIceLayer);

  bareIceAlbedoCategory = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"bareIceAlbedoCategory",
	{"nCategories"});
  columnVariables->add(bareIceAlbedoCategory);

  snowAlbedoCategory = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"snowAlbedoCategory",
	{"nCategories"});
  columnVariables->add(snowAlbedoCategory);

  pondAlbedoCategory = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"pondAlbedoCategory",
	{"nCategories"});
  columnVariables->add(pondAlbedoCategory);

  bareIceAlbedoCell = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"bareIceAlbedoCell",
	{});
  columnVariables->add(bareIceAlbedoCell);

  snowAlbedoCell = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"snowAlbedoCell",
	{});
  columnVariables->add(snowAlbedoCell);

  pondAlbedoCell = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"pondAlbedoCell",
	{});
  columnVariables->add(pondAlbedoCell);

  effectivePondAreaCategory = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"effectivePondAreaCategory",
	{"nCategories"});
  columnVariables->add(effectivePondAreaCategory);

  effectivePondAreaCell = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"effectivePondAreaCell",
	{});
  columnVariables->add(effectivePondAreaCell);

  snowFractionCategory = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"snowFractionCategory",
	{"nCategories"});
  columnVariables->add(snowFractionCategory);


  // drag
  airOceanDragCoefficientRatio = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"airOceanDragCoefficientRatio",
	{});
  columnVariables->add(airOceanDragCoefficientRatio);

  oceanDragCoefficient = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"oceanDragCoefficient",
	{});
  columnVariables->add(oceanDragCoefficient);

  oceanDragCoefficientSkin = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"oceanDragCoefficientSkin",
	{},
	useFormDrag);
  columnVariables->add(oceanDragCoefficientSkin);

  oceanDragCoefficientFloe = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"oceanDragCoefficientFloe",
	{},
	useFormDrag);
  columnVariables->add(oceanDragCoefficientFloe);

  oceanDragCoefficientKeel = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"oceanDragCoefficientKeel",
	{},
	useFormDrag);
  columnVariables->add(oceanDragCoefficientKeel);

  airDragCoefficient = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"airDragCoefficient",
	{});
  columnVariables->add(airDragCoefficient);

  airDragCoefficientSkin = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"airDragCoefficientSkin",
	{},
	useFormDrag);
  columnVariables->add(airDragCoefficientSkin);

  airDragCoefficientFloe = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"airDragCoefficientFloe",
	{},
	useFormDrag);
  columnVariables->add(airDragCoefficientFloe);

  airDragCoefficientPond = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"airDragCoefficientPond",
	{},
	useFormDrag);
  columnVariables->add(airDragCoefficientPond);

  airDragCoefficientRidge = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"airDragCoefficientRidge",
	{},
	useFormDrag);
  columnVariables->add(airDragCoefficientRidge);

  dragFreeboard = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"dragFreeboard",
	{},
	useFormDrag);
  columnVariables->add(dragFreeboard);

  dragIceSnowDraft = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"dragIceSnowDraft",
	{},
	useFormDrag);
  columnVariables->add(dragIceSnowDraft);

  dragRidgeHeight = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"dragRidgeHeight",
	{},
	useFormDrag);
  columnVariables->add(dragRidgeHeight);

  dragRidgeSeparation = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"dragRidgeSeparation",
	{},
	useFormDrag);
  columnVariables->add(dragRidgeSeparation);

  dragKeelDepth = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"dragKeelDepth",
	{},
	useFormDrag);
  columnVariables->add(dragKeelDepth);

  dragKeelSeparation = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"dragKeelSeparation",
	{},
	useFormDrag);
  columnVariables->add(dragKeelSeparation);

  dragFloeLength = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"dragFloeLength",
	{},
	useFormDrag);
  columnVariables->add(dragFloeLength);

  dragFloeSeparation = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"dragFloeSeparation",
	{},
	useFormDrag);
  columnVariables->add(dragFloeSeparation);


  // ponds
  pondFreshWaterFlux = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"pondFreshWaterFlux",
	{});
  columnVariables->add(pondFreshWaterFlux);

  pondSnowDepthDifference = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"pondSnowDepthDifference",
	{"nCategories"});
  columnVariables->add(pondSnowDepthDifference);

  pondLidMeltFluxFraction = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"pondLidMeltFluxFraction",
	{"nCategories"});
  columnVariables->add(pondLidMeltFluxFraction);


  // aerosols
  atmosAerosolFlux = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"atmosAerosolFlux",
	{"nAerosols"});
  columnVariables->add(atmosAerosolFlux);

  oceanAerosolFlux = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"oceanAerosolFlux",
	{"nAerosols"});
  columnVariables->add(oceanAerosolFlux);


  // initial
  initialSalinityProfile = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"initialSalinityProfile",
	{"nIceLayersP1"});
  columnVariables->add(initialSalinityProfile);

  initialMeltingTemperatureProfile = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"initialMeltingTemperatureProfile",
	{"nIceLayersP1"});
  columnVariables->add(initialMeltingTemperatureProfile);


  // dynamics
  airStressCellU = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"airStressCellU",
	{});
  columnVariables->add(airStressCellU);

  airStressCellV = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"airStressCellV",
	{});
  columnVariables->add(airStressCellV);

  oceanStressCellU = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"oceanStressCellU",
	{});
  columnVariables->add(oceanStressCellU);

  oceanStressCellV = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"oceanStressCellV",
	{});
  columnVariables->add(oceanStressCellV);


  // diagnostics
  meltOnset = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"meltOnset",
	{});
  columnVariables->add(meltOnset);

  freezeOnset = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"freezeOnset",
	{});
  columnVariables->add(freezeOnset);

  freezingMeltingPotentialInitial = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"freezingMeltingPotentialInitial",
	{});
  columnVariables->add(freezingMeltingPotentialInitial);

  iceAreaTendencyThermodynamics = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"iceAreaTendencyThermodynamics",
	{});
  columnVariables->add(iceAreaTendencyThermodynamics);

  iceVolumeTendencyThermodynamics = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"iceVolumeTendencyThermodynamics",
	{});
  columnVariables->add(iceVolumeTendencyThermodynamics);

  iceAgeTendencyThermodynamics = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"iceAgeTendencyThermodynamics",
	{});
  columnVariables->add(iceAgeTendencyThermodynamics);


  // location
  latitude = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"latitude",
	{});
  columnVariables->add(latitude);

  longitude = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"longitude",
	{});
  columnVariables->add(longitude);

  precipitationRate = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"precipitationRate",
	{});
  columnVariables->add(precipitationRate);
