module column_icepack

  implicit none

  private
  save

  ! time stepping
  public :: &
       column_prep_radiation, &
       column_step_radiation, &
       column_step_therm1, &
       column_step_therm2, &
       column_step_ridge, &
       column_coupling_prep, &
       column_ocean_mixed_layer

  ! init
  public :: &
       column_recompute_constants, &
       column_get_string_length, &
       column_init_thermo, &
       column_init_profiles, &
       column_init_itd, &
       column_init_category_areas, &
       column_init_particle_state, &
       column_init_orbit, &
       column_init_shortwave, &
       column_initialize_atmos_coupler_fields, &
       column_initialize_ocean_coupler_fields

  ! parameters
  public :: &
       column_configure, &
       column_init_parameters_real, &
       column_init_parameters_integer, &
       column_init_parameters_character, &
       column_init_parameters_logical, &
       column_query_parameters_real, &
       column_query_parameters_integer, &
       column_query_parameters_character, &
       column_query_parameters_logical

  ! tracers
  public :: &
       column_init_tracer_flags, &
       column_init_tracer_numbers, &
       column_init_tracer_indices, &
       column_aggregate, &
       column_print_tracer_object

  ! forcing
  public :: &
       column_limit_temperature, &
       column_limit_specific_humidity, &
       column_longwave_rosati_miyakoda, &
       column_longwave_parkinson_and_washington, &
       column_shortwave_from_cloud_fraction, &
       column_split_precipitation, &
       column_postprocess_ocean_forcing

  ! warnings
  public :: &
       column_warnings_reset, &
       column_warnings_number, &
       column_warnings_getone, &
       column_warnings_aborted

  ! misc
  public :: &
       column_initial_air_drag_coefficient, &
       column_get_constant, &
       column_liquidus_temperature

contains

  !-----------------------------------------------------------------------------
  ! time stepping
  !-----------------------------------------------------------------------------

  subroutine column_prep_radiation(&
       nCategories, &
       nIceLayers, &
       nIceLayersP1, &
       nSnowLayers, &
       iceAreaCell, &
       iceAreaCategory, &
       shortwaveVisibleDirectDown, &
       shortwaveVisibleDiffuseDown, &
       shortwaveIRDirectDown, &
       shortwaveIRDiffuseDown, &
       albedoVisibleDirectArea, &
       albedoVisibleDiffuseArea, &
       albedoIRDirectArea, &
       albedoIRDiffuseArea, &
       shortwaveScalingFactor, &
       surfaceShortwaveFlux, &
       interiorShortwaveFlux, &
       penetratingShortwaveFlux, &
       shortwaveLayerPenetration, &
       absorbedShortwaveSnowLayer, &
       absorbedShortwaveIceLayer) bind(C)

    use iso_c_binding, only: c_int, c_double

    use icepack_intfc, only: &
         icepack_prep_radiation

    integer(kind=c_int) :: &
         nCategories, &
         nIceLayers, &
         nIceLayersP1, &
         nSnowLayers

    real(kind=c_double) :: &
         iceAreaCell, &
         shortwaveVisibleDirectDown, &
         shortwaveVisibleDiffuseDown, &
         shortwaveIRDirectDown, &
         shortwaveIRDiffuseDown, &
         albedoVisibleDirectArea, &
         albedoVisibleDiffuseArea, &
         albedoIRDirectArea, &
         albedoIRDiffuseArea, &
         shortwaveScalingFactor

    real(kind=c_double), dimension(nCategories) :: &
         iceAreaCategory, &
         surfaceShortwaveFlux, &
         interiorShortwaveFlux, &
         penetratingShortwaveFlux

    real(kind=c_double), dimension(nIceLayersP1,nCategories) :: &
         shortwaveLayerPenetration

    real(kind=c_double), dimension(nSnowLayers,nCategories) :: &
         absorbedShortwaveSnowLayer

    real(kind=c_double), dimension(nIceLayers,nCategories) :: &
         absorbedShortwaveIceLayer

    call icepack_prep_radiation(&
         nCategories, &
         nIceLayers, &
         nSnowLayers, &
         iceAreaCell, &
         iceAreaCategory(:), &
         shortwaveVisibleDirectDown, &
         shortwaveVisibleDiffuseDown, &
         shortwaveIRDirectDown, &
         shortwaveIRDiffuseDown, &
         albedoVisibleDirectArea, &
         albedoVisibleDiffuseArea, &
         albedoIRDirectArea, &
         albedoIRDiffuseArea, &
         shortwaveScalingFactor, &
         surfaceShortwaveFlux(:), &
         interiorShortwaveFlux(:), &
         penetratingShortwaveFlux(:), &
         shortwaveLayerPenetration(:,:), &
         absorbedShortwaveSnowLayer(:,:), &
         absorbedShortwaveIceLayer(:,:))

  end subroutine column_prep_radiation

  !-----------------------------------------------------------------------------

  subroutine column_step_radiation(&
       nCategories, &
       nIceLayers, &
       nSnowLayers, &
       nAerosols, &
       nIceLayersP1, &
       timeStep, &
       calendarTypeIn, &
       daysInYear, &
       dayOfNextShortwaveCalculation, &
       dayOfYear, &
       secondsIntoDay, &
       latitude, &
       longitude, &
       tracerArrayCategory, &
       nTracers, &
       iceAreaCategory, &
       iceVolumeCategory, &
       snowVolumeCategory, &
       surfaceTemperature, &
       levelIceArea, &
       pondArea, &
       pondDepth, &
       pondLidThickness, &
       shortwaveVisibleDirectDown, &
       shortwaveVisibleDiffuseDown, &
       shortwaveIRDirectDown, &
       shortwaveIRDiffuseDown, &
       solarZenithAngleCosine, &
       snowfallRate, &
       albedoVisibleDirectCategory, &
       albedoVisibleDiffuseCategory, &
       albedoIRDirectCategory, &
       albedoIRDiffuseCategory, &
       surfaceShortwaveFlux, &
       interiorShortwaveFlux, &
       penetratingShortwaveFlux, &
       shortwaveLayerPenetration, &
       absorbedShortwaveSnowLayer, &
       absorbedShortwaveIceLayer, &
       bareIceAlbedoCategory, &
       snowAlbedoCategory, &
       pondAlbedoCategory, &
       effectivePondAreaCategory, &
       snowFractionCategory, &
       pondSnowDepthDifference, &
       pondLidMeltFluxFraction, &
       snowScatteringAerosol, &
       snowBodyAerosol, &
       iceScatteringAerosol, &
       iceBodyAerosol, &
       lInitializationIn, &
       nzAerosols, &
       nSpectralIntervals, &
       nModal1, &
       nModal2, &
       nShortwaveBio, &
       maxAerosolType, &
       indexShortwaveAerosol, &
       verticalShortwaveGrid, &
       verticalGrid, &
       brineFraction, &
       aerosolMassExtinctionCrossSection, &
       aerosolSingleScatterAlbedo, &
       aerosolAsymmetryParameter, &
       modalMassExtinctionCrossSection, &
       modalSingleScatterAlbedo, &
       modalAsymmetryParameter, &
       bioTracerShortwave, &
       modalBCabsorptionParameter) bind(C)

    use iso_c_binding, only: c_int, c_double, c_char

    use icepack_intfc, only: &
         icepack_step_radiation, &
         icepack_dbl_kind, &
         icepack_char_len_long

    use icepack_parameters, only: &
         pi

    integer(kind=c_int) :: &
         nCategories, &
         nIceLayers, &
         nSnowLayers, &
         nAerosols, &
         nIceLayersP1,&
         nTracers, &
         secondsIntoDay, &
         daysInYear, &
         lInitializationIn

    real(kind=c_double) :: &
         timeStep, &
         latitude, &
         longitude, &
         dayOfNextShortwaveCalculation, &
         dayOfYear

    character(kind=c_char, len=1), dimension(icepack_char_len_long), intent(in) :: &
         calendarTypeIn

    real(kind=c_double), dimension(nTracers,nCategories) :: &
         tracerArrayCategory

    real(kind=c_double), dimension(nCategories) :: &
         iceAreaCategory, &
         iceVolumeCategory, &
         snowVolumeCategory, &
         surfaceTemperature, &
         levelIceArea, &
         pondArea, &
         pondDepth, &
         pondLidThickness, &
         albedoVisibleDirectCategory, &
         albedoVisibleDiffuseCategory, &
         albedoIRDirectCategory, &
         albedoIRDiffuseCategory, &
         surfaceShortwaveFlux, &
         interiorShortwaveFlux, &
         penetratingShortwaveFlux, &
         bareIceAlbedoCategory, &
         snowAlbedoCategory, &
         pondAlbedoCategory, &
         effectivePondAreaCategory, &
         snowFractionCategory, &
         pondSnowDepthDifference, &
         pondLidMeltFluxFraction

    real(kind=c_double) :: &
         shortwaveVisibleDirectDown, &
         shortwaveVisibleDiffuseDown, &
         shortwaveIRDirectDown, &
         shortwaveIRDiffuseDown, &
         solarZenithAngleCosine, &
         snowfallRate

    real(kind=c_double), dimension(nIceLayersP1,nCategories) :: &
         shortwaveLayerPenetration

    real(kind=c_double), dimension(nIceLayers,nCategories) :: &
         absorbedShortwaveIceLayer

    real(kind=c_double), dimension(nSnowLayers,nCategories) :: &
         absorbedShortwaveSnowLayer

    real(kind=c_double), dimension(nAerosols,nCategories) :: &
         snowScatteringAerosol, &
         snowBodyAerosol, &
         iceScatteringAerosol, &
         iceBodyAerosol

    ! bgc aerosol dummy arrays
    integer(kind=c_int) :: &
         nzAerosols, &
         nSpectralIntervals, &
         nModal1, &
         nModal2, &
         nShortwaveBio, &
         maxAerosolType

    integer(kind=c_int), dimension(maxAerosolType) :: &
         indexShortwaveAerosol

    real(kind=c_double), dimension(nIceLayersP1) :: &
         verticalShortwaveGrid, &
         verticalGrid

    real(kind=c_double), dimension(nCategories) :: &
         brineFraction

    real(kind=c_double), dimension(nSpectralIntervals,maxAerosolType) :: &
         aerosolMassExtinctionCrossSection, &
         aerosolSingleScatterAlbedo, &
         aerosolAsymmetryParameter

    real(kind=c_double), dimension(nSpectralIntervals,nModal1) :: &
         modalMassExtinctionCrossSection, &
         modalSingleScatterAlbedo, &
         modalAsymmetryParameter

    real(kind=c_double), dimension(nShortwaveBio,nCategories) :: &
         bioTracerShortwave

    real(kind=c_double), dimension(nSpectralIntervals,nModal1,nModal2) :: &
         modalBCabsorptionParameter

    integer :: &
         indexChlorophyllShortwave = 0

    ! biogeochemistry
    logical :: &
         useZaerosols = .False., &
         useShortwaveBioabsorption = .False., &
         useModalAerosols = .False.

    integer :: &
         nAlgae = 3, &
         nBioLayers = 0, &
         nBioTracersShortwave = 0

    ! other
    real(kind=icepack_dbl_kind) :: &
         longitudeAdjust

    character(len=icepack_char_len_long) :: &
         calendarType

    logical :: &
         lInitialization

    integer :: &
         iCategory, &
         iAerosol

    ! aerosols array
    real(kind=icepack_dbl_kind), dimension(:,:), allocatable :: &
         aerosolsArray

    longitudeAdjust = longitude
    if (longitudeAdjust > pi) longitudeAdjust = longitudeAdjust - 2.0_icepack_dbl_kind * pi

    call convert_cstring_to_fstring(calendarTypeIn, calendarType)

    ! set aerosols array
    allocate(aerosolsArray(4*nAerosols,nCategories))
    do iCategory = 1, nCategories
       do iAerosol = 1, nAerosols

          aerosolsArray(1+4*(iAerosol-1), iCategory) = snowScatteringAerosol(iAerosol,iCategory)
          aerosolsArray(2+4*(iAerosol-1), iCategory) = snowBodyAerosol(iAerosol,iCategory)
          aerosolsArray(3+4*(iAerosol-1), iCategory) = iceScatteringAerosol(iAerosol,iCategory)
          aerosolsArray(4+4*(iAerosol-1), iCategory) = iceBodyAerosol(iAerosol,iCategory)

       enddo ! iAerosol
    enddo ! iCategory

    lInitialization = convert_clogical_to_flogical(lInitializationIn)

    call icepack_step_radiation(&
         timeStep, &
         nCategories, &
         nAlgae, &
         useZaerosols, &
         nBioLayers, &
         nTracers, &
         nBioTracersShortwave, &
         nIceLayers, &
         nSnowLayers, &
         nAerosols, &
         nzAerosols, &
         useShortwaveBioabsorption, &
         indexChlorophyllShortwave, &
         indexShortwaveAerosol, &
         verticalShortwaveGrid(:), &
         verticalGrid(:), &
         brineFraction(:), &
         iceAreaCategory(:), &
         iceVolumeCategory(:), &
         snowVolumeCategory(:), &
         surfaceTemperature(:), &
         levelIceArea(:), &
         pondArea(:), &
         pondDepth(:), &
         pondLidThickness(:), &
         aerosolsArray, &
         bioTracerShortwave(:,:), &
         tracerArrayCategory(:,:), &
         latitude, &
         longitudeAdjust, &
         calendarType, &
         daysInYear, &
         dayOfNextShortwaveCalculation, &
         dayOfYear, &
         secondsIntoDay, &
         aerosolMassExtinctionCrossSection(:,:), &
         aerosolSingleScatterAlbedo(:,:), &
         aerosolAsymmetryParameter(:,:), &
         modalMassExtinctionCrossSection(:,:), &
         modalSingleScatterAlbedo(:,:), &
         modalAsymmetryParameter(:,:), &
         modalBCabsorptionParameter(:,:,:), &
         useModalAerosols, &
         shortwaveVisibleDirectDown, &
         shortwaveVisibleDiffuseDown, &
         shortwaveIRDirectDown, &
         shortwaveIRDiffuseDown, &
         solarZenithAngleCosine, &
         snowfallRate, &
         albedoVisibleDirectCategory(:), &
         albedoVisibleDiffuseCategory(:), &
         albedoIRDirectCategory(:), &
         albedoIRDiffuseCategory(:), &
         surfaceShortwaveFlux(:), &
         interiorShortwaveFlux(:), &
         penetratingShortwaveFlux(:), &
         shortwaveLayerPenetration(:,:), &
         absorbedShortwaveSnowLayer(:,:), &
         absorbedShortwaveIceLayer(:,:), &
         bareIceAlbedoCategory(:), &
         snowAlbedoCategory(:), &
         pondAlbedoCategory(:), &
         effectivePondAreaCategory(:), &
         snowFractionCategory(:), &
         pondSnowDepthDifference(:), &
         pondLidMeltFluxFraction(:), &
         .false., &
         lInitialization)

    deallocate(aerosolsArray)

  end subroutine column_step_radiation

  !-----------------------------------------------------------------------------

  subroutine column_step_therm1(&
       nCategories, &
       nIceLayers, &
       nSnowLayers, &
       nAerosols, &
       timeStep, &
       useAerosolsIn, &
       latitude, &
       iceAreaCellInitial, &
       iceAreaCategoryInitial, &
       iceVolumeCategoryInitial, &
       snowVolumeCategoryInitial, &
       iceAreaCell, &
       iceAreaCategory, &
       iceVolumeCell, &
       iceVolumeCategory, &
       snowVolumeCell, &
       snowVolumeCategory, &
       uVelocity, &
       vVelocity, &
       surfaceTemperature, &
       snowEnthalpy, &
       iceEnthalpy, &
       iceSalinity, &
       levelIceArea, &
       levelIceVolume, &
       pondArea, &
       pondDepth, &
       pondLidThickness, &
       iceAge, &
       firstYearIceArea, &
       uAirVelocity, &
       vAirVelocity, &
       windSpeed, &
       airLevelHeight, &
       airSpecificHumidity, &
       airDensity, &
       airTemperature, &
       atmosReferenceTemperature2m, &
       atmosReferenceHumidity2m, &
       atmosReferenceSpeed10m, &
       airOceanDragCoefficientRatio, &
       oceanDragCoefficient, &
       oceanDragCoefficientSkin, &
       oceanDragCoefficientFloe, &
       oceanDragCoefficientKeel, &
       airDragCoefficient, &
       airDragCoefficientSkin, &
       airDragCoefficientFloe, &
       airDragCoefficientPond, &
       airDragCoefficientRidge, &
       dragFreeboard, &
       dragIceSnowDraft, &
       dragRidgeHeight, &
       dragRidgeSeparation, &
       dragKeelDepth, &
       dragKeelSeparation, &
       dragFloeLength, &
       dragFloeSeparation, &
       airStressForcingU, &
       airStressForcingV, &
       airStressCellU, &
       airStressCellV, &
       airPotentialTemperature, &
       seaSurfaceTemperature, &
       seaSurfaceSalinity, &
       seaFreezingTemperature, &
       oceanStressCellU, &
       oceanStressCellV, &
       oceanHeatFluxIceBottom, &
       freezingMeltingPotential, &
       lateralIceMeltFraction, &
       snowfallRate, &
       rainfallRate, &
       pondFreshWaterFlux, &
       surfaceHeatFlux, &
       surfaceHeatFluxCategory, &
       surfaceConductiveFlux, &
       surfaceConductiveFluxCategory, &
       surfaceShortwaveFlux, &
       interiorShortwaveFlux, &
       penetratingShortwaveFlux, &
       absorbedShortwaveFlux, &
       longwaveUp, &
       absorbedShortwaveSnowLayer, &
       absorbedShortwaveIceLayer, &
       longwaveDown, &
       sensibleHeatFlux, &
       sensibleHeatFluxCategory, &
       latentHeatFlux, &
       latentHeatFluxCategory, &
       evaporativeWaterFlux, &
       oceanFreshWaterFlux, &
       oceanSaltFlux, &
       oceanHeatFlux, &
       oceanShortwaveFlux, &
       latentHeatFluxCouple, &
       sensibleHeatFluxCouple, &
       surfaceHeatFluxCouple, &
       surfaceConductiveFluxCouple, &
       atmosAerosolFlux, &
       oceanAerosolFlux, &
       pondSnowDepthDifference, &
       pondLidMeltFluxFraction, &
       surfaceIceMelt, &
       surfaceIceMeltCategory, &
       basalIceMelt, &
       basalIceMeltCategory, &
       snowMelt, &
       snowMeltCategory, &
       congelation, &
       congelationCategory, &
       snowiceFormation, &
       snowiceFormationCategory, &
       snowThicknessChangeCategory, &
       meltOnset, &
       freezeOnset, &
       snowScatteringAerosol, &
       snowBodyAerosol, &
       iceScatteringAerosol, &
       iceBodyAerosol, &
       dayOfYear) bind(C)

    use iso_c_binding, only: c_int, c_double, c_char

    use icepack_intfc, only: &
         icepack_step_therm1, &
         icepack_dbl_kind, &
         icepack_char_len_long

    use icepack_parameters, only: &
         puny

    integer(c_int) :: &
         nCategories, &
         nIceLayers, &
         nSnowLayers, &
         nAerosols, &
         useAerosolsIn

    real(c_double) :: &
         timeStep, &
         latitude, &
         iceAreaCellInitial, &
         iceAreaCell, &
         iceVolumeCell, &
         snowVolumeCell, &
         uVelocity, &
         vVelocity, &
         uAirVelocity, &
         vAirVelocity, &
         windSpeed, &
         airLevelHeight, &
         airSpecificHumidity, &
         airDensity, &
         airTemperature, &
         atmosReferenceTemperature2m, &
         atmosReferenceHumidity2m, &
         atmosReferenceSpeed10m, &
         airOceanDragCoefficientRatio, &
         oceanDragCoefficient, &
         oceanDragCoefficientSkin, &
         oceanDragCoefficientFloe, &
         oceanDragCoefficientKeel, &
         airDragCoefficient, &
         airDragCoefficientSkin, &
         airDragCoefficientFloe, &
         airDragCoefficientPond, &
         airDragCoefficientRidge, &
         dragFreeboard, &
         dragIceSnowDraft, &
         dragRidgeHeight, &
         dragRidgeSeparation, &
         dragKeelDepth, &
         dragKeelSeparation, &
         dragFloeLength, &
         dragFloeSeparation, &
         airStressForcingU, &
         airStressForcingV, &
         airStressCellU, &
         airStressCellV, &
         airPotentialTemperature, &
         seaSurfaceTemperature, &
         seaSurfaceSalinity, &
         seaFreezingTemperature, &
         oceanStressCellU, &
         oceanStressCellV, &
         oceanHeatFluxIceBottom, &
         freezingMeltingPotential, &
         lateralIceMeltFraction, &
         snowfallRate, &
         rainfallRate, &
         pondFreshWaterFlux, &
         surfaceHeatFlux, &
         surfaceConductiveFlux, &
         absorbedShortwaveFlux, &
         longwaveUp, &
         longwaveDown, &
         sensibleHeatFlux, &
         latentHeatFlux, &
         evaporativeWaterFlux, &
         oceanFreshWaterFlux, &
         oceanSaltFlux, &
         oceanHeatFlux, &
         oceanShortwaveFlux, &
         surfaceIceMelt, &
         basalIceMelt, &
         snowMelt, &
         congelation, &
         snowiceFormation, &
         meltOnset, &
         freezeOnset, &
         dayOfYear

    real(c_double), dimension(nCategories) :: &
         iceAreaCategoryInitial, &
         iceVolumeCategoryInitial, &
         snowVolumeCategoryInitial, &
         iceAreaCategory, &
         iceVolumeCategory, &
         snowVolumeCategory, &
         surfaceTemperature, &
         levelIceArea, &
         levelIceVolume, &
         pondArea, &
         pondDepth, &
         pondLidThickness, &
         iceAge, &
         firstYearIceArea, &
         surfaceHeatFluxCategory, &
         surfaceConductiveFluxCategory, &
         surfaceShortwaveFlux, &
         interiorShortwaveFlux, &
         penetratingShortwaveFlux, &
         sensibleHeatFluxCategory, &
         latentHeatFluxCategory, &
         surfaceIceMeltCategory, &
         basalIceMeltCategory, &
         snowMeltCategory, &
         congelationCategory, &
         snowiceFormationCategory, &
         snowThicknessChangeCategory, &
         latentHeatFluxCouple, &
         sensibleHeatFluxCouple, &
         surfaceHeatFluxCouple, &
         surfaceConductiveFluxCouple, &
         pondSnowDepthDifference, &
         pondLidMeltFluxFraction

    real(c_double), dimension(nSnowLayers, nCategories) :: &
         snowEnthalpy, &
         absorbedShortwaveSnowLayer

    real(c_double), dimension(nIceLayers, nCategories) :: &
         iceEnthalpy, &
         iceSalinity, &
         absorbedShortwaveIceLayer

    real(c_double), dimension(nAerosols) :: &
         atmosAerosolFlux, &
         oceanAerosolFlux

    real(kind=c_double), dimension(nAerosols,nCategories) :: &
         snowScatteringAerosol, &
         snowBodyAerosol, &
         iceScatteringAerosol, &
         iceBodyAerosol

    ! local
    logical :: &
         northernHemisphereMask, &
         useAerosols

    real(kind=icepack_dbl_kind), dimension(:,:,:), allocatable :: &
         specificSnowAerosol, &
         specificIceAerosol

    integer :: &
         iCategory, &
         iAerosol

    ! convert flag from c_int to logical
    useAerosols = convert_clogical_to_flogical(useAerosolsIn)

    ! initial state values
    iceAreaCellInitial = iceAreaCell

    do iCategory = 1, nCategories

       iceAreaCategoryInitial(iCategory) = iceAreaCategory(iCategory)
       iceVolumeCategoryInitial(iCategory) = iceVolumeCategory(iCategory)
       snowVolumeCategoryInitial(iCategory) = snowVolumeCategory(iCategory)

    enddo ! iCategory

    ! aerosol
    if (useAerosols) then

       allocate(specificSnowAerosol(nAerosols, 2, nCategories))
       allocate(specificIceAerosol(nAerosols, 2, nCategories))

       do iCategory = 1, nCategories
          do iAerosol = 1, nAerosols

             specificSnowAerosol(iAerosol, 1, iCategory) = &
                  snowScatteringAerosol(iAerosol,iCategory) * snowVolumeCategoryInitial(iCategory)
             specificSnowAerosol(iAerosol, 2, iCategory) = &
                  snowBodyAerosol(iAerosol,iCategory)       * snowVolumeCategoryInitial(iCategory)

             specificIceAerosol(iAerosol, 1, iCategory) = &
                  iceScatteringAerosol(iAerosol,iCategory)   * iceVolumeCategoryInitial(iCategory)
             specificIceAerosol(iAerosol, 2, iCategory) = &
                  iceBodyAerosol(iAerosol,iCategory)         * iceVolumeCategoryInitial(iCategory)

          enddo ! iAerosol
       enddo ! iCategory

    else

       allocate(specificSnowAerosol(1, 1, 1))
       allocate(specificIceAerosol(1, 1, 1))
       specificSnowAerosol = 0.0_icepack_dbl_kind
       specificIceAerosol  = 0.0_icepack_dbl_kind

    end if

    ! hemisphere mask
    if (latitude > 0.0_c_double) then
       northernHemisphereMask = .true.
    else
       northernHemisphereMask = .false.
    endif

    call icepack_step_therm1(&
         timeStep, &
         nCategories, &
         nIceLayers, &
         nSnowLayers, &
         nAerosols, &
         iceAreaCategoryInitial(:), &
         iceVolumeCategoryInitial(:), &
         snowVolumeCategoryInitial(:), &
         iceAreaCell, &
         iceAreaCategory(:), &
         iceVolumeCell, &
         iceVolumeCategory(:), &
         snowVolumeCell, &
         snowVolumeCategory(:), &
         uVelocity, &
         vVelocity, &
         surfaceTemperature(:), &
         snowEnthalpy(:,:), &
         iceEnthalpy(:,:), &
         iceSalinity(:,:), &
         levelIceArea(:), &
         levelIceVolume(:), &
         pondArea(:), &
         pondDepth(:), &
         pondLidThickness(:), &
         iceAge(:), &
         firstYearIceArea(:), &
         specificSnowAerosol(:,:,:), &
         specificIceAerosol(:,:,:), &
         uAirVelocity, &
         vAirVelocity, &
         windSpeed, &
         airLevelHeight, &
         airSpecificHumidity, &
         airDensity, &
         airTemperature, &
         atmosReferenceTemperature2m, &
         atmosReferenceHumidity2m, &
         atmosReferenceSpeed10m, &
         airOceanDragCoefficientRatio, &
         oceanDragCoefficient, &
         oceanDragCoefficientSkin, &
         oceanDragCoefficientFloe, &
         oceanDragCoefficientKeel, &
         airDragCoefficient, &
         airDragCoefficientSkin, &
         airDragCoefficientFloe, &
         airDragCoefficientPond, &
         airDragCoefficientRidge, &
         dragFreeboard, &
         dragIceSnowDraft, &
         dragRidgeHeight, &
         dragRidgeSeparation, &
         dragKeelDepth, &
         dragKeelSeparation, &
         dragFloeLength, &
         dragFloeSeparation, &
         airStressForcingU, &
         airStressForcingV, &
         airStressCellU, &
         airStressCellV, &
         airPotentialTemperature, &
         seaSurfaceTemperature, &
         seaSurfaceSalinity, &
         seaFreezingTemperature, &
         oceanStressCellU, &
         oceanStressCellV, &
         oceanHeatFluxIceBottom, &
         freezingMeltingPotential, &
         lateralIceMeltFraction, &
         snowfallRate, &
         rainfallRate, &
         pondFreshWaterFlux, &
         surfaceHeatFlux, &
         surfaceHeatFluxCategory(:), &
         surfaceConductiveFlux, &
         surfaceConductiveFluxCategory(:), &
         surfaceShortwaveFlux(:), &
         interiorShortwaveFlux(:), &
         penetratingShortwaveFlux(:), &
         absorbedShortwaveFlux, &
         longwaveUp, &
         absorbedShortwaveSnowLayer(:,:), &
         absorbedShortwaveIceLayer(:,:), &
         longwaveDown, &
         sensibleHeatFlux, &
         sensibleHeatFluxCategory(:), &
         latentHeatFlux, &
         latentHeatFluxCategory(:), &
         evaporativeWaterFlux, &
         oceanFreshWaterFlux, &
         oceanSaltFlux, &
         oceanHeatFlux, &
         oceanShortwaveFlux, &
         latentHeatFluxCouple(:), &
         sensibleHeatFluxCouple(:), &
         surfaceHeatFluxCouple(:), &
         surfaceConductiveFluxCouple(:), &
         atmosAerosolFlux(:), &
         oceanAerosolFlux(:), &
         pondSnowDepthDifference(:), &
         pondLidMeltFluxFraction(:), &
         surfaceIceMelt, &
         surfaceIceMeltCategory(:), &
         basalIceMelt, &
         basalIceMeltCategory(:), &
         snowMelt, &
         snowMeltCategory(:), &
         congelation, &
         congelationCategory(:), &
         snowiceFormation, &
         snowiceFormationCategory(:), &
         snowThicknessChangeCategory(:), &
         northernHemisphereMask, &
         .not. northernHemisphereMask, &
         meltOnset, &
         freezeOnset, &
         dayOfYear)

    ! aerosol
    if (useAerosols) then

       do iCategory = 1, nCategories
          do iAerosol = 1, nAerosols

             if (snowVolumeCategory(iCategory) > puny) &
                  specificSnowAerosol(iAerosol, :, iCategory) = &
                  specificSnowAerosol(iAerosol, :, iCategory) / snowVolumeCategory(iCategory)

             if (iceVolumeCategory(iCategory) > puny) &
                  specificIceAerosol(iAerosol, :, iCategory)  = &
                  specificIceAerosol(iAerosol, :, iCategory)  / iceVolumeCategory(iCategory)

             snowScatteringAerosol(iAerosol,iCategory) = specificSnowAerosol(iAerosol, 1, iCategory)
             snowBodyAerosol(iAerosol,iCategory)       = specificSnowAerosol(iAerosol, 2, iCategory)

             iceScatteringAerosol(iAerosol,iCategory)  = specificIceAerosol(iAerosol, 1, iCategory)
             iceBodyAerosol(iAerosol,iCategory)        = specificIceAerosol(iAerosol, 2, iCategory)

          enddo ! iAerosol
       enddo ! iCategory

    endif

    ! aerosols
    deallocate(specificSnowAerosol)
    deallocate(specificIceAerosol)

  end subroutine column_step_therm1

  !-----------------------------------------------------------------------------

  subroutine column_step_therm2(&
       nCategories, &
       nCategoriesP1, &
       nTracers, &
       nBaseTracers, &
       nMaxAncestorTracers, &
       nAerosols, &
       nIceLayers, &
       nIceLayersP1, &
       nSnowLayers, &
       timeStep, &
       categoryThicknessLimits, &
       iceAreaCategory, &
       iceVolumeCategory, &
       snowVolumeCategory, &
       iceAreaCategoryInitial, &
       iceVolumeCategoryInitial, &
       tracerArrayCategory, &
       openWaterArea, &
       iceAreaCell, &
       parentIndex, &
       firstAncestorMask, &
       ancestorNumber, &
       ancestorIndices, &
       seaFreezingTemperature, &
       seaSurfaceSalinity, &
       initialSalinityProfile, &
       lateralIceMeltFraction, &
       lateralIceMelt, &
       freezingMeltingPotential, &
       frazilFormation, &
       rainfallRate, &
       pondFreshWaterFlux, &
       oceanFreshWaterFlux, &
       oceanSaltFlux, &
       oceanHeatFlux, &
       updateOceanFluxesIn, &
       oceanAerosolFlux, &
       frazilGrowthDiagnostic, &
       freezeOnset, &
       dayOfYear) bind(C)

    use iso_c_binding, only: c_int, c_double

    use icepack_intfc, only: &
         icepack_step_therm2, &
         icepack_dbl_kind

    integer(kind=c_int) :: &
         nCategories, &
         nCategoriesP1, &
         nTracers, &
         nBaseTracers, &
         nMaxAncestorTracers, &
         nAerosols, &
         nIceLayers, &
         nIceLayersP1, &
         nSnowLayers, &
         updateOceanFluxesIn

    real(kind=c_double) :: &
         timeStep, &
         openWaterArea, &
         iceAreaCell, &
         seaFreezingTemperature, &
         seaSurfaceSalinity, &
         lateralIceMeltFraction, &
         lateralIceMelt, &
         freezingMeltingPotential, &
         frazilFormation, &
         rainfallRate, &
         pondFreshWaterFlux, &
         oceanFreshWaterFlux, &
         oceanSaltFlux, &
         oceanHeatFlux, &
         frazilGrowthDiagnostic, &
         freezeOnset, &
         dayOfYear

    real(kind=c_double), dimension(nCategories) :: &
         iceAreaCategory, &
         iceVolumeCategory, &
         snowVolumeCategory, &
         iceAreaCategoryInitial, &
         iceVolumeCategoryInitial

    real(kind=c_double), dimension(nCategoriesP1) :: &
         categoryThicknessLimits

    real(kind=c_double), dimension(nAerosols) :: &
         oceanAerosolFlux

    real(kind=c_double), dimension(nTracers,nCategories) :: &
         tracerArrayCategory

    integer(kind=c_int), dimension(nTracers) :: &
         parentIndex, &
         ancestorNumber

    real(kind=c_double), dimension(nTracers, nBaseTracers) :: &
         firstAncestorMask

    integer(kind=c_int), dimension(nTracers, nMaxAncestorTracers) :: &
         ancestorIndices

    real(kind=c_double), dimension(nIceLayersP1) :: &
         initialSalinityProfile

    logical :: &
         updateOceanFluxes

    integer, parameter :: &
         nBioLayers = 0, &
         nBioLayersP2 = 0, &
         nBioLayersP1 = 0, &
         nZBGCTracers = 0, &
         nBioTracers = 0, &
         nBioTracersLayer = 0

    real(kind=icepack_dbl_kind), dimension(:), allocatable :: &
         biologyGrid, &
         interfaceBiologyGrid, &
         verticalGrid, &
         oceanBioFluxes, &
         oceanBioConcentrationsUsed

    logical, dimension(:), allocatable :: &
         newlyFormedIce

    real(kind=icepack_dbl_kind) :: &
         zSalinityFlux

    updateOceanFluxes = convert_clogical_to_flogical(updateOceanFluxesIn)

    allocate(biologyGrid(nBioLayersP2))
    allocate(interfaceBiologyGrid(nBioLayersP1))
    allocate(verticalGrid(0)) ! allocate(verticalGrid(nIceLayersP1))

    allocate(newlyFormedIce(nCategories))

    allocate(oceanBioFluxes(nZBGCTracers))
    allocate(oceanBioConcentrationsUsed(nBioTracers))

    call icepack_step_therm2(&
         timeStep, &
         nCategories, &
         nAerosols, &
         nBioTracersLayer, &
         nIceLayers, &
         nSnowLayers, &
         categoryThicknessLimits(:), &
         nBioLayers, &
         iceAreaCategory(:), &
         iceVolumeCategory(:), &
         snowVolumeCategory(:), &
         iceAreaCategoryInitial(:), &
         iceVolumeCategoryInitial(:), &
         tracerArrayCategory(:,:), &
         openWaterArea, &
         iceAreaCell, &
         parentIndex, & ! trcr_depend
         firstAncestorMask, & ! trcr_base
         ancestorNumber, & ! n_trcr_strata
         ancestorIndices, & ! nt_strata
         seaFreezingTemperature, &
         seaSurfaceSalinity, &
         initialSalinityProfile(:), &
         lateralIceMeltFraction, &
         lateralIceMelt, &
         freezingMeltingPotential, &
         frazilFormation, &
         rainfallRate, &
         pondFreshWaterFlux, &
         oceanFreshWaterFlux, &
         oceanSaltFlux, &
         oceanHeatFlux, &
         updateOceanFluxes, &
         biologyGrid(:), & ! bgrid, intent(in)
         verticalGrid(:), & ! cgrid, intent(in)
         interfaceBiologyGrid(:), & ! igrid, intent(in)
         oceanAerosolFlux(:), &
         newlyFormedIce, &
         zSalinityFlux, &
         oceanBioFluxes(:), &
         oceanBioConcentrationsUsed(:), &
         frazilGrowthDiagnostic, &
         freezeOnset, &
         dayOfYear)

    deallocate(biologyGrid)
    deallocate(interfaceBiologyGrid)
    deallocate(verticalGrid)

    deallocate(newlyFormedIce)

    deallocate(oceanBioFluxes)
    deallocate(oceanBioConcentrationsUsed)

  end subroutine column_step_therm2

  !-----------------------------------------------------------------------------

  subroutine column_step_ridge(&
       nCategories, &
       nCategoriesP1, &
       nTracers, &
       nBaseTracers, &
       nMaxAncestorTracers, &
       nAerosols, &
       nIceLayers, &
       nSnowLayers, &
       nDynamicsSubcycles, &
       timeStep, &
       categoryThicknessLimits, &
       ridgeConvergence, &
       ridgeShear, &
       iceAreaCategory, &
       iceVolumeCategory, &
       snowVolumeCategory, &
       tracerArrayCategory, &
       openWaterArea, &
       parentIndex, &
       firstAncestorMask, &
       ancestorNumber, &
       ancestorIndices, &
       areaLossRidge, &
       areaGainRidge, &
       iceVolumeRidged, &
       openingRateRidge, &
       pondFreshWaterFlux, &
       oceanFreshWaterFlux, &
       oceanHeatFlux, &
       oceanAerosolFlux, &
       iceAreaCell, &
       oceanSaltFlux) bind(C)

    use iso_c_binding, only: c_int, c_double

    use icepack_intfc, only: &
         icepack_step_ridge, &
         icepack_dbl_kind

    integer(kind=c_int) :: &
         nCategories, &
         nCategoriesP1, &
         nTracers, &
         nBaseTracers, &
         nMaxAncestorTracers, &
         nAerosols, &
         nIceLayers, &
         nSnowLayers, &
         nDynamicsSubcycles

    real(kind=c_double) :: &
         timeStep, &
         ridgeConvergence, &
         ridgeShear, &
         openWaterArea, &
         areaLossRidge, &
         areaGainRidge, &
         iceVolumeRidged, &
         openingRateRidge, &
         pondFreshWaterFlux, &
         oceanFreshWaterFlux, &
         oceanHeatFlux, &
         iceAreaCell, &
         oceanSaltFlux

    real(c_double), dimension(nCategories) :: &
         iceAreaCategory, &
         iceVolumeCategory, &
         snowVolumeCategory, &
         ridgeParticipationFunction, &
         ratioRidgeThicknessToIce, &
         fractionNewRidgeArea, &
         fractionNewRidgeVolume, &
         areaLossRidgeCategory, &
         areaGainRidgeCategory, &
         iceVolumeRidgedCategory, &
         raftingIceArea, &
         raftingIceVolume

    real(kind=c_double), dimension(nCategoriesP1) :: &
         categoryThicknessLimits

    real(kind=c_double), dimension(nTracers,nCategories) :: &
         tracerArrayCategory

    integer(kind=c_int), dimension(nTracers) :: &
         parentIndex, &
         ancestorNumber

    real(kind=c_double), dimension(nTracers, nBaseTracers) :: &
         firstAncestorMask

    integer(kind=c_int), dimension(nTracers, nMaxAncestorTracers) :: &
         ancestorIndices

    real(kind=c_double), dimension(nAerosols) :: &
         oceanAerosolFlux

    integer, parameter :: &
         nBioLayers = 0, &
         nZBGCTracers = 0

    logical, dimension(:), allocatable :: &
         newlyFormedIce

    real(kind=icepack_dbl_kind) :: &
         zSalinityFlux

    real(kind=icepack_dbl_kind), dimension(:), allocatable :: &
         oceanBioFluxes

    allocate(newlyFormedIce(nCategories))

    allocate(oceanBioFluxes(nZBGCTracers))

    call icepack_step_ridge(&
         timeStep, &
         nDynamicsSubcycles, &
         nIceLayers, &
         nSnowLayers, &
         nBioLayers, &
         nCategories, &
         categoryThicknessLimits, &
         ridgeConvergence, &
         ridgeShear, &
         iceAreaCategory(:), &
         tracerArrayCategory, &
         iceVolumeCategory(:), &
         snowVolumeCategory(:), &
         openWaterArea, &
         parentIndex, & ! trcr_depend
         firstAncestorMask, & ! trcr_base
         ancestorNumber, & ! n_trcr_strata
         ancestorIndices, & ! nt_strata
         areaLossRidge, &
         areaGainRidge, &
         iceVolumeRidged, &
         openingRateRidge, &
         pondFreshWaterFlux, &
         oceanFreshWaterFlux, &
         oceanHeatFlux, &
         nAerosols, &
         oceanAerosolFlux(:), &
         ridgeParticipationFunction(:), &
         ratioRidgeThicknessToIce(:), &
         fractionNewRidgeArea(:), &
         fractionNewRidgeVolume(:), &
         areaLossRidgeCategory(:), &
         areaGainRidgeCategory(:), &
         iceVolumeRidgedCategory(:), &
         raftingIceArea(:), &
         raftingIceVolume(:), &
         iceAreaCell, &
         oceanSaltFlux, &
         newlyFormedIce(:), &
         zSalinityFlux, &
         oceanBioFluxes(:))

    deallocate(newlyFormedIce)
    deallocate(oceanBioFluxes)

  end subroutine column_step_ridge

  !-----------------------------------------------------------------------------

  subroutine column_ocean_mixed_layer(&
       timeStep, &
       seaSurfaceTemperature, &
       airPotentialTemperature, &
       uAirVelocity, &
       vAirVelocity, &
       windSpeed, &
       airLevelHeight, &
       airSpecificHumidity, &
       airDensity, &
       airStressOceanU, &
       airStressOceanV, &
       atmosReferenceTemperature2mOcean, &
       atmosReferenceHumidity2mOcean, &
       airDragCoefficient, &
       airOceanDragCoefficientRatio, &
       albedoVisibleDirectOcean, &
       shortwaveVisibleDirectDown, &
       albedoIRDirectOcean, &
       shortwaveIRDirectDown, &
       albedoVisibleDiffuseOcean, &
       shortwaveVisibleDiffuseDown, &
       albedoIRDiffuseOcean, &
       shortwaveIRDiffuseDown, &
       longwaveUpOcean, &
       sensibleHeatFluxOcean, &
       sensibleTransferCoefficient, &
       latentHeatFluxOcean, &
       evaporativeWaterFluxOcean, &
       longwaveDown, &
       iceAreaCell, &
       oceanHeatFlux, &
       oceanShortwaveFlux, &
       oceanMixedLayerDepth, &
       seaFreezingTemperature, &
       oceanHeatFluxConvergence, &
       freezingMeltingPotential) bind(C)

    use iso_c_binding, only: c_double

    use icepack_intfc, only: &
         icepack_atm_boundary, &
         icepack_ocn_mixed_layer, &
         icepack_dbl_kind

    use icepack_parameters, only: &
         oceanAlbedo => albocn

    real(c_double) :: &
         timeStep, &
         seaSurfaceTemperature, &
         airPotentialTemperature, &
         uAirVelocity, &
         vAirVelocity, &
         windSpeed, &
         airLevelHeight, &
         airSpecificHumidity, &
         airDensity, &
         airStressOceanU, &
         airStressOceanV, &
         atmosReferenceTemperature2mOcean, &
         atmosReferenceHumidity2mOcean, &
         airDragCoefficient, &
         airOceanDragCoefficientRatio, &
         albedoVisibleDirectOcean, &
         shortwaveVisibleDirectDown, &
         albedoIRDirectOcean, &
         shortwaveIRDirectDown, &
         albedoVisibleDiffuseOcean, &
         shortwaveVisibleDiffuseDown, &
         albedoIRDiffuseOcean, &
         shortwaveIRDiffuseDown, &
         longwaveUpOcean, &
         sensibleHeatFluxOcean, &
         latentHeatFluxOcean, &
         evaporativeWaterFluxOcean, &
         longwaveDown, &
         iceAreaCell, &
         oceanHeatFlux, &
         oceanShortwaveFlux, &
         oceanMixedLayerDepth, &
         seaFreezingTemperature, &
         oceanHeatFluxConvergence, &
         freezingMeltingPotential

    real(kind=icepack_dbl_kind) :: &
         sensibleTransferCoefficient, &
         latentTransferCoefficient, &
         potentialTemperatureDifference, &
         specificHumidityDifference

    call icepack_atm_boundary(&
         'ocn', &
         seaSurfaceTemperature, &
         airPotentialTemperature, &
         uAirVelocity, &
         vAirVelocity, &
         windSpeed, &
         airLevelHeight, &
         airSpecificHumidity, &
         airDensity, &
         airStressOceanU, &
         airStressOceanV, &
         atmosReferenceTemperature2mOcean, &
         atmosReferenceHumidity2mOcean, &
         potentialTemperatureDifference, &
         specificHumidityDifference, &
         latentTransferCoefficient, &
         sensibleTransferCoefficient, &
         airDragCoefficient, &
         airOceanDragCoefficientRatio)

    albedoVisibleDirectOcean  = oceanAlbedo
    albedoIRDirectOcean       = oceanAlbedo
    albedoVisibleDiffuseOcean = oceanAlbedo
    albedoIRDiffuseOcean      = oceanAlbedo

    call icepack_ocn_mixed_layer(&
         albedoVisibleDirectOcean, &
         shortwaveVisibleDirectDown, &
         albedoIRDirectOcean, &
         shortwaveIRDirectDown, &
         albedoVisibleDiffuseOcean, &
         shortwaveVisibleDiffuseDown, &
         albedoIRDiffuseOcean, &
         shortwaveIRDiffuseDown, &
         seaSurfaceTemperature, &
         longwaveUpOcean, &
         sensibleHeatFluxOcean, &
         sensibleTransferCoefficient, &
         latentHeatFluxOcean, &
         latentTransferCoefficient, &
         evaporativeWaterFluxOcean, &
         longwaveDown, &
         potentialTemperatureDifference, &
         specificHumidityDifference, &
         iceAreaCell, &
         oceanHeatFlux, &
         oceanShortwaveFlux, &
         oceanMixedLayerDepth, &
         seaFreezingTemperature, &
         oceanHeatFluxConvergence, &
         freezingMeltingPotential, &
         timeStep)

  end subroutine column_ocean_mixed_layer

  !-----------------------------------------------------------------------------
  ! init
  !-----------------------------------------------------------------------------

  subroutine column_recompute_constants() bind(C)

    use icepack_intfc, only: &
         icepack_recompute_constants

    call icepack_recompute_constants

  end subroutine column_recompute_constants

  !-----------------------------------------------------------------------------

  function column_get_string_length() result(strLen) bind(C)

    use iso_c_binding, only: c_int

    use icepack_intfc, only: &
         icepack_char_len_long

    integer(c_int) :: strLen

    strLen = icepack_char_len_long

  end function column_get_string_length

  !-----------------------------------------------------------------------------

  subroutine column_init_thermo(nilyr, nilyrP1, sprofile) bind(C)

    use iso_c_binding, only: c_int, c_double

    use icepack_intfc, only: &
         icepack_init_thermo

    integer(c_int) :: &
         nilyr, &
         nilyrP1

    real(c_double), dimension(nilyrP1) :: &
         sprofile

    call icepack_init_thermo(nilyr, sprofile)

  end subroutine column_init_thermo

  !-----------------------------------------------------------------------------

  subroutine column_init_profiles(&
       nIceLayersP1, &
       initialSalinityProfileVertical, &
       initialSalinityProfile, &
       initialMeltingTemperatureProfile) bind(C)

    use iso_c_binding, only: c_int, c_double

    use icepack_intfc, only: &
         icepack_liquidus_temperature

    integer(c_int) :: &
         nIceLayersP1

    real(c_double), dimension(nIceLayersP1) :: &
         initialSalinityProfileVertical

    real(c_double), dimension(nIceLayersP1) :: &
         initialSalinityProfile, &
         initialMeltingTemperatureProfile

    integer :: &
         iIceLayer

    do iIceLayer = 1, nIceLayersP1

       ! these profiles are not used by mushy
       initialSalinityProfile(iIceLayer)           = initialSalinityProfileVertical(iIceLayer)
       initialMeltingTemperatureProfile(iIceLayer) = &
            icepack_liquidus_temperature(initialSalinityProfileVertical(iIceLayer))

    enddo ! iIceLayer

  end subroutine column_init_profiles

  !-----------------------------------------------------------------------------

  subroutine column_init_itd(ncat, hin_max) bind(C)

    use iso_c_binding, only: c_int, c_double

    use icepack_intfc, only: &
         icepack_init_itd

    integer(c_int), intent(in) :: &
         ncat ! number of categories

    real(c_double), intent(out), dimension(0:ncat) :: &
         hin_max  ! category limits (m)

    call icepack_init_itd(ncat, hin_max)

  end subroutine column_init_itd

  !-----------------------------------------------------------------------------

  subroutine column_init_category_areas(&
       nCategoriesP1, &
       nCategories, &
       categoryThicknessLimits, &
       initialCategoryIceThickness, &
       initialCategoryIceArea) bind(C)

    use iso_c_binding, only: c_int, c_double

    use icepack_intfc, only: &
         icepack_dbl_kind

    use icepack_parameters, only: &
         puny

    integer(c_int) :: &
         nCategoriesP1, &
         nCategories

    real(c_double), dimension(nCategoriesP1) :: &
         categoryThicknessLimits

    real(c_double), dimension(nCategories) :: &
         initialCategoryIceThickness, &
         initialCategoryIceArea

    integer :: &
         iCategory

    real(kind=icepack_dbl_kind) :: &
         areaCategorySum

    real(kind=icepack_dbl_kind), parameter :: &
         thicknessWithLargestArea = 3.0_icepack_dbl_kind ! initial ice thickness with greatest area

    !--------------------------------------------------------
    ! areas
    !--------------------------------------------------------

    ! initialize sum of areas in categories
    areaCategorySum = 0.0_icepack_dbl_kind

    do iCategory = 1, nCategories

       ! parabola, max at h=thicknessWithLargestArea, zero at h=0, 2*thicknessWithLargestArea
       initialCategoryIceArea(iCategory) = &
            max(0.0_icepack_dbl_kind, &
            (2.0_icepack_dbl_kind * thicknessWithLargestArea * initialCategoryIceThickness(iCategory) - &
             initialCategoryIceThickness(iCategory)**2))

       areaCategorySum = areaCategorySum + initialCategoryIceArea(iCategory)

    enddo ! iCategory

    ! normalize
    do iCategory = 1, nCategories

       initialCategoryIceArea(iCategory) = initialCategoryIceArea(iCategory) / &
            (areaCategorySum + puny / nCategories)

    enddo ! iCategory

  end subroutine column_init_category_areas

  !-----------------------------------------------------------------------------

  subroutine column_init_particle_state(&
       nCategories, &
       nIceLayers, &
       nSnowLayers, &
       surfaceTemperature, &
       airTemperature, &
       seaFreezingTemperature, &
       initialSalinityProfile, &
       initialMeltingTemperatureProfile, &
       iceEnthalpy, &
       iceSalinity, &
       snowEnthalpy) bind(C)

    use iso_c_binding, only: c_int, c_double

    use icepack_intfc, only: &
         icepack_init_trcr, &
         icepack_dbl_kind

    integer(c_int) :: &
         nCategories, &
         nIceLayers, &
         nSnowLayers

    real(c_double), dimension(nCategories) :: &
         surfaceTemperature

    real(c_double) :: &
         airTemperature, &
         seaFreezingTemperature

    real(c_double), dimension(nIceLayers) :: &
         initialSalinityProfile, &
         initialMeltingTemperatureProfile

    real(c_double), dimension(nIceLayers,nCategories) :: &
         iceEnthalpy, &
         iceSalinity

    real(c_double), dimension(nSnowLayers,nCategories) :: &
         snowEnthalpy

    integer :: &
         iCategory, &
         iIceLayer

    do iCategory = 1, nCategories

       call icepack_init_trcr(&
            airTemperature, &
            seaFreezingTemperature, &
            initialSalinityProfile(:), &
            initialMeltingTemperatureProfile(:), &
            surfaceTemperature(iCategory), &
            nIceLayers, &
            nSnowLayers, &
            iceEnthalpy(:,iCategory), &
            snowEnthalpy(:,iCategory))

       do iIceLayer = 1, nIceLayers
          iceSalinity(iIceLayer,iCategory) = initialSalinityProfile(iIceLayer)
       enddo ! iIceLayer

    enddo ! iCategory

  end subroutine column_init_particle_state

  !-----------------------------------------------------------------------------

  subroutine column_init_orbit() bind(C)

    use icepack_intfc, only: &
         icepack_init_orbit

    call icepack_init_orbit()

  end subroutine column_init_orbit

  !-----------------------------------------------------------------------------

  subroutine column_init_shortwave(&
       nCategories, &
       doRestartIn, &
       iceAreaCategory, &
       albedoVisibleDirectCategory, &
       albedoVisibleDiffuseCategory, &
       albedoIRDirectCategory, &
       albedoIRDiffuseCategory, &
       albedoVisibleDirectCell, &
       albedoVisibleDiffuseCell, &
       albedoIRDirectCell, &
       albedoIRDiffuseCell, &
       solarZenithAngleCosine, &
       bareIceAlbedoCategory, &
       snowAlbedoCategory, &
       pondAlbedoCategory, &
       bareIceAlbedoCell, &
       snowAlbedoCell, &
       pondAlbedoCell, &
       effectivePondAreaCategory, &
       effectivePondAreaCell, &
       albedoVisibleDirectArea, &
       albedoVisibleDiffuseArea, &
       albedoIRDirectArea, &
       albedoIRDiffuseArea, &
       shortwaveScalingFactor, &
       shortwaveVisibleDirectDown, &
       shortwaveVisibleDiffuseDown, &
       shortwaveIRDirectDown, &
       shortwaveIRDiffuseDown) bind(C)

    use iso_c_binding, only: c_int, c_double

    use icepack_intfc, only: &
         icepack_dbl_kind

    use icepack_parameters, only: &
         puny

    integer(c_int) :: &
         nCategories, &
         doRestartIn

    real(c_double), dimension(nCategories) :: &
         iceAreaCategory, &
         albedoVisibleDirectCategory, &
         albedoVisibleDiffuseCategory, &
         albedoIRDirectCategory, &
         albedoIRDiffuseCategory, &
         bareIceAlbedoCategory, &
         snowAlbedoCategory, &
         pondAlbedoCategory, &
         effectivePondAreaCategory

    real(c_double) :: &
         albedoVisibleDirectCell, &
         albedoVisibleDiffuseCell, &
         albedoIRDirectCell, &
         albedoIRDiffuseCell, &
         solarZenithAngleCosine, &
         bareIceAlbedoCell, &
         snowAlbedoCell, &
         pondAlbedoCell, &
         effectivePondAreaCell, &
         albedoVisibleDirectArea, &
         albedoVisibleDiffuseArea, &
         albedoIRDirectArea, &
         albedoIRDiffuseArea, &
         shortwaveScalingFactor, &
         shortwaveVisibleDirectDown, &
         shortwaveVisibleDiffuseDown, &
         shortwaveIRDirectDown, &
         shortwaveIRDiffuseDown

    integer :: &
         iCategory

    logical :: &
         doRestart

    doRestart = convert_clogical_to_flogical(doRestartIn)

    albedoVisibleDirectCell  = 0.0_icepack_dbl_kind
    albedoVisibleDiffuseCell = 0.0_icepack_dbl_kind
    albedoIRDirectCell       = 0.0_icepack_dbl_kind
    albedoIRDiffuseCell      = 0.0_icepack_dbl_kind

    do iCategory = 1, nCategories

       ! aggregate albedos
       if (iceAreaCategory(iCategory) > puny) then

          albedoVisibleDirectCell  = albedoVisibleDirectCell  + &
               albedoVisibleDirectCategory(iCategory)  * iceAreaCategory(iCategory)
          albedoVisibleDiffuseCell = albedoVisibleDiffuseCell + &
               albedoVisibleDiffuseCategory(iCategory) * iceAreaCategory(iCategory)
          albedoIRDirectCell       = albedoIRDirectCell       + &
               albedoIRDirectCategory(iCategory)       * iceAreaCategory(iCategory)
          albedoIRDiffuseCell      = albedoIRDiffuseCell      + &
               albedoIRDiffuseCategory(iCategory)      * iceAreaCategory(iCategory)

          if (solarZenithAngleCosine > puny) then ! sun above horizon

             bareIceAlbedoCell = bareIceAlbedoCell + &
                  bareIceAlbedoCategory(iCategory) * iceAreaCategory(iCategory)
             snowAlbedoCell    = snowAlbedoCell    + &
                  snowAlbedoCategory(iCategory)    * iceAreaCategory(iCategory)
             pondAlbedoCell    = pondAlbedoCell    + &
                  pondAlbedoCategory(iCategory)    * iceAreaCategory(iCategory)

          endif

          effectivePondAreaCell = effectivePondAreaCell + &
               effectivePondAreaCategory(iCategory) * iceAreaCategory(iCategory)

       endif

    enddo ! iCategory

    ! Store grid box mean albedos and fluxes before scaling by aice
    albedoVisibleDirectArea  = albedoVisibleDirectCell
    albedoVisibleDiffuseArea = albedoVisibleDiffuseCell
    albedoIRDirectArea       = albedoIRDirectCell
    albedoIRDiffuseArea      = albedoIRDiffuseCell

    ! Save net shortwave for scaling factor in scale_factor
    if (.not. doRestart) then
       shortwaveScalingFactor = &
            shortwaveVisibleDirectDown  * (1.0_icepack_dbl_kind - albedoVisibleDirectArea) + &
            shortwaveVisibleDiffuseDown * (1.0_icepack_dbl_kind - albedoVisibleDiffuseArea) + &
            shortwaveIRDirectDown       * (1.0_icepack_dbl_kind - albedoIRDirectArea) + &
            shortwaveIRDiffuseDown      * (1.0_icepack_dbl_kind - albedoIRDiffuseArea)
    endif

  end subroutine column_init_shortwave

  !-----------------------------------------------------------------------------

  subroutine column_initialize_atmos_coupler_fields(&
       useAerosolsIn, &
       windSpeed, &
       uAirVelocity, &
       vAirVelocity, &
       longwaveUp, &
       airDragCoefficient, &
       atmosAerosolFlux) bind(C)

    use iso_c_binding, only: c_int, c_double

    use icepack_intfc, only: &
         icepack_dbl_kind

    use icepack_parameters, only: &
         seaiceStefanBoltzmann => stefan_boltzmann, &
         seaiceFreshWaterFreezingPoint => Tffresh

    integer(c_int) :: &
         useAerosolsIn

    real(c_double) :: &
         windSpeed, &
         uAirVelocity, &
         vAirVelocity, &
         longwaveUp, &
         airDragCoefficient, &
         atmosAerosolFlux

    logical :: &
         useAerosols

    !-------------------------------------------------------------
    ! Physics fluxes received from atmosphere
    !-------------------------------------------------------------

    windSpeed = sqrt(uAirVelocity**2 + vAirVelocity**2)

    !-------------------------------------------------------------
    ! fluxes sent to atmosphere
    !-------------------------------------------------------------

    longwaveUp = -seaiceStefanBoltzmann * seaiceFreshWaterFreezingPoint**4

    airDragCoefficient = column_initial_air_drag_coefficient()

    !-------------------------------------------------------------
    ! Aerosol fluxes received from atmosphere
    !-------------------------------------------------------------

    useAerosols = convert_clogical_to_flogical(useAerosolsIn)
    if (useAerosols) then
       atmosAerosolFlux = 1.e-12_icepack_dbl_kind
    endif

  end subroutine column_initialize_atmos_coupler_fields

  !-----------------------------------------------------------------------------

  subroutine column_initialize_ocean_coupler_fields(&
       doRestartIn, &
       seaSurfaceTemperature, &
       seaFreezingTemperature, &
       seaSurfaceSalinity) bind(C)

    use iso_c_binding, only: c_int, c_double

    use icepack_intfc, only: &
         icepack_liquidus_temperature

    integer(c_int) :: &
         doRestartIn

    real(c_double) :: &
         seaSurfaceTemperature, &
         seaFreezingTemperature, &
         seaSurfaceSalinity

    logical :: &
         doRestart

    !-------------------------------------------------------------
    ! Physics fluxes received from ocean
    !-------------------------------------------------------------

    seaFreezingTemperature = icepack_liquidus_temperature(seaSurfaceSalinity)

    ! sea surface temperature is not initialized if we're restarting
    doRestart = convert_clogical_to_flogical(doRestartIn)
    if (.not. doRestart) then
       seaSurfaceTemperature = seaFreezingTemperature
    endif

  end subroutine column_initialize_ocean_coupler_fields

  !-----------------------------------------------------------------------------

  function column_initial_air_drag_coefficient() result(airDragCoefficient) bind(C)

    use iso_c_binding, only: c_double

    use icepack_intfc, only: &
         icepack_dbl_kind

    use icepack_parameters, only: &
         seaiceVonKarmanConstant => vonkar, &
         seaiceIceSurfaceRoughness => iceruf, &
         seaiceStabilityReferenceHeight => zref

    real(kind=c_double) :: airDragCoefficient

    ! atmo drag for RASM
    airDragCoefficient = (seaiceVonKarmanConstant/log(seaiceStabilityReferenceHeight/seaiceIceSurfaceRoughness)) &
                       * (seaiceVonKarmanConstant/log(seaiceStabilityReferenceHeight/seaiceIceSurfaceRoughness))

  end function column_initial_air_drag_coefficient

  !-----------------------------------------------------------------------------

  subroutine column_coupling_prep(&
       nCategories, &
       includePondFreshwaterFeedbackIn, &
       timeStep, &
       freezingMeltingPotentialInitial, &
       freezingMeltingPotential, &
       albedoVisibleDirectCell, &
       albedoVisibleDiffuseCell, &
       albedoIRDirectCell, &
       albedoIRDiffuseCell, &
       bareIceAlbedoCell, &
       snowAlbedoCell, &
       pondAlbedoCell, &
       effectivePondAreaCell, &
       solarZenithAngleCosine, &
       albedoVisibleDirectArea, &
       albedoVisibleDiffuseArea, &
       albedoIRDirectArea, &
       albedoIRDiffuseArea, &
       oceanFreshWaterFluxArea, &
       oceanSaltFluxArea, &
       oceanHeatFluxArea, &
       oceanShortwaveFluxArea, &
       oceanSaltFlux, &
       oceanHeatFlux, &
       oceanShortwaveFlux, &
       shortwaveScalingFactor, &
       pondFreshWaterFlux, &
       oceanFreshWaterFlux, &
       shortwaveVisibleDirectDown, &
       shortwaveVisibleDiffuseDown, &
       shortwaveIRDirectDown, &
       shortwaveIRDiffuseDown, &
       iceAreaCategory, &
       albedoVisibleDirectCategory, &
       albedoVisibleDiffuseCategory, &
       albedoIRDirectCategory, &
       albedoIRDiffuseCategory, &
       bareIceAlbedoCategory, &
       snowAlbedoCategory, &
       pondAlbedoCategory, &
       effectivePondAreaCategory) bind(C)

    use iso_c_binding, only: c_int, c_double

    use icepack_intfc, only: &
         icepack_dbl_kind

    use icepack_parameters, only: &
         puny, &
         densityFreshwater => rhofresh

    integer(c_int) :: &
         nCategories, &
         includePondFreshwaterFeedbackIn

    real(kind=c_double) :: &
         timeStep, &
         freezingMeltingPotentialInitial, &
         freezingMeltingPotential, &
         albedoVisibleDirectCell, &
         albedoVisibleDiffuseCell, &
         albedoIRDirectCell, &
         albedoIRDiffuseCell, &
         bareIceAlbedoCell, &
         snowAlbedoCell, &
         pondAlbedoCell, &
         effectivePondAreaCell, &
         solarZenithAngleCosine, &
         albedoVisibleDirectArea, &
         albedoVisibleDiffuseArea, &
         albedoIRDirectArea, &
         albedoIRDiffuseArea, &
         oceanFreshWaterFluxArea, &
         oceanSaltFluxArea, &
         oceanHeatFluxArea, &
         oceanShortwaveFluxArea, &
         oceanSaltFlux, &
         oceanHeatFlux, &
         oceanShortwaveFlux, &
         shortwaveScalingFactor, &
         pondFreshWaterFlux, &
         oceanFreshWaterFlux, &
         shortwaveVisibleDirectDown, &
         shortwaveVisibleDiffuseDown, &
         shortwaveIRDirectDown, &
         shortwaveIRDiffuseDown

    real(kind=c_double), dimension(nCategories) :: &
         iceAreaCategory, &
         albedoVisibleDirectCategory, &
         albedoVisibleDiffuseCategory, &
         albedoIRDirectCategory, &
         albedoIRDiffuseCategory, &
         bareIceAlbedoCategory, &
         snowAlbedoCategory, &
         pondAlbedoCategory, &
         effectivePondAreaCategory

    logical :: &
         includePondFreshwaterFeedback

    integer :: &
         iCategory

    ! convert flag from c_int to logical
    includePondFreshwaterFeedback = convert_clogical_to_flogical(includePondFreshwaterFeedbackIn)

    !-------------------------------------------------------------------
    ! store initial freezing melting potential
    !-------------------------------------------------------------------

    freezingMeltingPotentialInitial = freezingMeltingPotential

    !-------------------------------------------------------------------
    ! aggregate albedos
    !-------------------------------------------------------------------

    albedoVisibleDirectCell  = 0.0_icepack_dbl_kind
    albedoVisibleDiffuseCell = 0.0_icepack_dbl_kind
    albedoIRDirectCell       = 0.0_icepack_dbl_kind
    albedoIRDiffuseCell      = 0.0_icepack_dbl_kind

    bareIceAlbedoCell = 0.0_icepack_dbl_kind
    snowAlbedoCell    = 0.0_icepack_dbl_kind
    pondAlbedoCell    = 0.0_icepack_dbl_kind

    effectivePondAreaCell = 0.0_icepack_dbl_kind

    do iCategory = 1, nCategories

       albedoVisibleDirectCell = albedoVisibleDirectCell + &
            albedoVisibleDirectCategory(iCategory) * iceAreaCategory(iCategory)
       albedoVisibleDiffuseCell = albedoVisibleDiffuseCell + &
            albedoVisibleDiffuseCategory(iCategory) * iceAreaCategory(iCategory)
       albedoIRDirectCell = albedoIRDirectCell + &
            albedoIRDirectCategory(iCategory) * iceAreaCategory(iCategory)
       albedoIRDiffuseCell = albedoIRDiffuseCell + &
            albedoIRDiffuseCategory(iCategory) * iceAreaCategory(iCategory)

       ! sun above horizon
       if (solarZenithAngleCosine > puny) then

          bareIceAlbedoCell = bareIceAlbedoCell + &
               bareIceAlbedoCategory(iCategory) * iceAreaCategory(iCategory)
          snowAlbedoCell = snowAlbedoCell + &
               snowAlbedoCategory(iCategory) * iceAreaCategory(iCategory)
          pondAlbedoCell = pondAlbedoCell + &
               pondAlbedoCategory(iCategory) * iceAreaCategory(iCategory)

       endif

       effectivePondAreaCell = effectivePondAreaCell + &
            effectivePondAreaCategory(iCategory) * iceAreaCategory(iCategory)

    enddo ! iCategory

    !-------------------------------------------------------------------
    ! reduce oceanFreshWaterFlux by pondFreshWaterFlux for coupling
    !-------------------------------------------------------------------

    if (includePondFreshwaterFeedback) then
       pondFreshWaterFlux  = pondFreshWaterFlux * densityFreshwater / timeStep
       oceanFreshWaterFlux = oceanFreshWaterFlux - pondFreshWaterFlux
    endif

    !-------------------------------------------------------------------
    ! Store grid box mean albedos and fluxes before scaling by aice
    !-------------------------------------------------------------------

    albedoVisibleDirectArea  = albedoVisibleDirectCell
    albedoVisibleDiffuseArea = albedoVisibleDiffuseCell
    albedoIRDirectArea       = albedoIRDirectCell
    albedoIRDiffuseArea      = albedoIRDiffuseCell
    oceanFreshWaterFluxArea  = oceanFreshWaterFlux
    oceanSaltFluxArea        = oceanSaltFlux
    oceanHeatFluxArea        = oceanHeatFlux
    oceanShortwaveFluxArea   = oceanShortwaveFlux

    !-----------------------------------------------------------------
    ! Save net shortwave for scaling factor in shortwaveScalingFactor
    !-----------------------------------------------------------------

    shortwaveScalingFactor = &
         shortwaveVisibleDirectDown  * (1.0_icepack_dbl_kind - albedoVisibleDirectArea) + &
         shortwaveVisibleDiffuseDown * (1.0_icepack_dbl_kind - albedoVisibleDiffuseArea) + &
         shortwaveIRDirectDown       * (1.0_icepack_dbl_kind - albedoIRDirectArea) + &
         shortwaveIRDiffuseDown      * (1.0_icepack_dbl_kind - albedoIRDiffuseArea)

  end subroutine column_coupling_prep

  !-----------------------------------------------------------------------------
  ! parameters
  !-----------------------------------------------------------------------------

  subroutine column_configure() bind(C)

    use icepack_intfc, only: &
         icepack_configure

    call icepack_configure()

  end subroutine column_configure

  !-----------------------------------------------------------------------------

  subroutine column_init_parameters_real(paramNameIn, paramValue) bind(C)

    use iso_c_binding, only: c_char, c_double, c_null_char

    use icepack_intfc, only: &
         icepack_init_parameters, &
         icepack_char_len_long

    character(kind=c_char, len=1), dimension(icepack_char_len_long), intent(in) :: &
         paramNameIn

    real(kind=c_double), intent(in) :: &
         paramValue

    character(len=icepack_char_len_long) :: &
         paramName

    call convert_cstring_to_fstring(paramNameIn, paramName)

    if (trim(paramName) == "secday"           ) call icepack_init_parameters(secday_in            = paramValue)
    if (trim(paramName) == "puny"             ) call icepack_init_parameters(puny_in              = paramValue)
    if (trim(paramName) == "bignum"           ) call icepack_init_parameters(bignum_in            = paramValue)
    if (trim(paramName) == "pi"               ) call icepack_init_parameters(pi_in                = paramValue)
    if (trim(paramName) == "rhos"             ) call icepack_init_parameters(rhos_in              = paramValue)
    if (trim(paramName) == "rhoi"             ) call icepack_init_parameters(rhoi_in              = paramValue)
    if (trim(paramName) == "rhosi"            ) call icepack_init_parameters(rhosi_in             = paramValue)
    if (trim(paramName) == "rhow"             ) call icepack_init_parameters(rhow_in              = paramValue)
    if (trim(paramName) == "rhofresh"         ) call icepack_init_parameters(rhofresh_in          = paramValue)
    if (trim(paramName) == "cp_ice"           ) call icepack_init_parameters(cp_ice_in            = paramValue)
    if (trim(paramName) == "cp_ocn"           ) call icepack_init_parameters(cp_ocn_in            = paramValue)
    if (trim(paramName) == "depressT"         ) call icepack_init_parameters(depressT_in          = paramValue)
    if (trim(paramName) == "viscosity_dyn"    ) call icepack_init_parameters(viscosity_dyn_in     = paramValue)
    if (trim(paramName) == "Tocnfrz"          ) call icepack_init_parameters(Tocnfrz_in           = paramValue)
    if (trim(paramName) == "Tffresh"          ) call icepack_init_parameters(Tffresh_in           = paramValue)
    if (trim(paramName) == "Lsub"             ) call icepack_init_parameters(Lsub_in              = paramValue)
    if (trim(paramName) == "Lvap"             ) call icepack_init_parameters(Lvap_in              = paramValue)
    if (trim(paramName) == "Timelt"           ) call icepack_init_parameters(Timelt_in            = paramValue)
    if (trim(paramName) == "Tsmelt"           ) call icepack_init_parameters(Tsmelt_in            = paramValue)
    if (trim(paramName) == "ice_ref_salinity" ) call icepack_init_parameters(ice_ref_salinity_in  = paramValue)
    if (trim(paramName) == "kice"             ) call icepack_init_parameters(kice_in              = paramValue)
    if (trim(paramName) == "kseaice"          ) call icepack_init_parameters(kseaice_in           = paramValue)
    if (trim(paramName) == "ksno"             ) call icepack_init_parameters(ksno_in              = paramValue)
    if (trim(paramName) == "hs_min"           ) call icepack_init_parameters(hs_min_in            = paramValue)
    if (trim(paramName) == "snowpatch"        ) call icepack_init_parameters(snowpatch_in         = paramValue)
    if (trim(paramName) == "saltmax"          ) call icepack_init_parameters(saltmax_in           = paramValue)
    if (trim(paramName) == "phi_init"         ) call icepack_init_parameters(phi_init_in          = paramValue)
    if (trim(paramName) == "min_salin"        ) call icepack_init_parameters(min_salin_in         = paramValue)
    if (trim(paramName) == "salt_loss"        ) call icepack_init_parameters(salt_loss_in         = paramValue)
    if (trim(paramName) == "dSin0_frazil"     ) call icepack_init_parameters(dSin0_frazil_in      = paramValue)
    if (trim(paramName) == "dts_b"            ) call icepack_init_parameters(dts_b_in             = paramValue)
    if (trim(paramName) == "ustar_min"        ) call icepack_init_parameters(ustar_min_in         = paramValue)
    if (trim(paramName) == "a_rapid_mode"     ) call icepack_init_parameters(a_rapid_mode_in      = paramValue)
    if (trim(paramName) == "Rac_rapid_mode"   ) call icepack_init_parameters(Rac_rapid_mode_in    = paramValue)
    if (trim(paramName) == "aspect_rapid_mode") call icepack_init_parameters(aspect_rapid_mode_in = paramValue)
    if (trim(paramName) == "dSdt_slow_mode"   ) call icepack_init_parameters(dSdt_slow_mode_in    = paramValue)
    if (trim(paramName) == "phi_c_slow_mode"  ) call icepack_init_parameters(phi_c_slow_mode_in   = paramValue)
    if (trim(paramName) == "phi_i_mushy"      ) call icepack_init_parameters(phi_i_mushy_in       = paramValue)
    if (trim(paramName) == "emissivity"       ) call icepack_init_parameters(emissivity_in        = paramValue)
    if (trim(paramName) == "albocn"           ) call icepack_init_parameters(albocn_in            = paramValue)
    if (trim(paramName) == "vonkar"           ) call icepack_init_parameters(vonkar_in            = paramValue)
    if (trim(paramName) == "stefan_boltzmann" ) call icepack_init_parameters(stefan_boltzmann_in  = paramValue)
    if (trim(paramName) == "kappav"           ) call icepack_init_parameters(kappav_in            = paramValue)
    if (trim(paramName) == "hi_ssl"           ) call icepack_init_parameters(hi_ssl_in            = paramValue)
    if (trim(paramName) == "hs_ssl"           ) call icepack_init_parameters(hs_ssl_in            = paramValue)
    if (trim(paramName) == "awtvdr"           ) call icepack_init_parameters(awtvdr_in            = paramValue)
    if (trim(paramName) == "awtidr"           ) call icepack_init_parameters(awtidr_in            = paramValue)
    if (trim(paramName) == "awtvdf"           ) call icepack_init_parameters(awtvdf_in            = paramValue)
    if (trim(paramName) == "awtidf"           ) call icepack_init_parameters(awtidf_in            = paramValue)
    if (trim(paramName) == "albicev"          ) call icepack_init_parameters(albicev_in           = paramValue)
    if (trim(paramName) == "albicei"          ) call icepack_init_parameters(albicei_in           = paramValue)
    if (trim(paramName) == "albsnowv"         ) call icepack_init_parameters(albsnowv_in          = paramValue)
    if (trim(paramName) == "albsnowi"         ) call icepack_init_parameters(albsnowi_in          = paramValue)
    if (trim(paramName) == "ahmax"            ) call icepack_init_parameters(ahmax_in             = paramValue)
    if (trim(paramName) == "R_ice"            ) call icepack_init_parameters(R_ice_in             = paramValue)
    if (trim(paramName) == "R_pnd"            ) call icepack_init_parameters(R_pnd_in             = paramValue)
    if (trim(paramName) == "R_snw"            ) call icepack_init_parameters(R_snw_in             = paramValue)
    if (trim(paramName) == "dT_mlt"           ) call icepack_init_parameters(dT_mlt_in            = paramValue)
    if (trim(paramName) == "rsnw_mlt"         ) call icepack_init_parameters(rsnw_mlt_in          = paramValue)
    if (trim(paramName) == "kalg"             ) call icepack_init_parameters(kalg_in              = paramValue)
    if (trim(paramName) == "Cf"               ) call icepack_init_parameters(Cf_in                = paramValue)
    if (trim(paramName) == "Pstar"            ) call icepack_init_parameters(Pstar_in             = paramValue)
    if (trim(paramName) == "Cstar"            ) call icepack_init_parameters(Cstar_in             = paramValue)
    if (trim(paramName) == "dragio"           ) call icepack_init_parameters(dragio_in            = paramValue)
    if (trim(paramName) == "gravit"           ) call icepack_init_parameters(gravit_in            = paramValue)
    if (trim(paramName) == "iceruf"           ) call icepack_init_parameters(iceruf_in            = paramValue)
    if (trim(paramName) == "mu_rdg"           ) call icepack_init_parameters(mu_rdg_in            = paramValue)
    if (trim(paramName) == "cp_air"           ) call icepack_init_parameters(cp_air_in            = paramValue)
    if (trim(paramName) == "cp_wv"            ) call icepack_init_parameters(cp_wv_in             = paramValue)
    if (trim(paramName) == "zvir"             ) call icepack_init_parameters(zvir_in              = paramValue)
    if (trim(paramName) == "zref"             ) call icepack_init_parameters(zref_in              = paramValue)
    if (trim(paramName) == "qqqice"           ) call icepack_init_parameters(qqqice_in            = paramValue)
    if (trim(paramName) == "TTTice"           ) call icepack_init_parameters(TTTice_in            = paramValue)
    if (trim(paramName) == "qqqocn"           ) call icepack_init_parameters(qqqocn_in            = paramValue)
    if (trim(paramName) == "TTTocn"           ) call icepack_init_parameters(TTTocn_in            = paramValue)
    if (trim(paramName) == "grid_o"           ) call icepack_init_parameters(grid_o_in            = paramValue)
    if (trim(paramName) == "l_sk"             ) call icepack_init_parameters(l_sk_in              = paramValue)
    if (trim(paramName) == "initbio_frac"     ) call icepack_init_parameters(initbio_frac_in      = paramValue)
    if (trim(paramName) == "phi_snow"         ) call icepack_init_parameters(phi_snow_in          = paramValue)
    if (trim(paramName) == "grid_oS"          ) call icepack_init_parameters(grid_oS_in           = paramValue)
    if (trim(paramName) == "l_skS"            ) call icepack_init_parameters(l_skS_in             = paramValue)
    if (trim(paramName) == "fr_resp"          ) call icepack_init_parameters(fr_resp_in           = paramValue)
    if (trim(paramName) == "algal_vel"        ) call icepack_init_parameters(algal_vel_in         = paramValue)
    if (trim(paramName) == "R_dFe2dust"       ) call icepack_init_parameters(R_dFe2dust_in        = paramValue)
    if (trim(paramName) == "dustFe_sol"       ) call icepack_init_parameters(dustFe_sol_in        = paramValue)
    if (trim(paramName) == "T_max"            ) call icepack_init_parameters(T_max_in             = paramValue)
    if (trim(paramName) == "fsal"             ) call icepack_init_parameters(fsal_in              = paramValue)
    if (trim(paramName) == "op_dep_min"       ) call icepack_init_parameters(op_dep_min_in        = paramValue)
    if (trim(paramName) == "fr_graze_s"       ) call icepack_init_parameters(fr_graze_s_in        = paramValue)
    if (trim(paramName) == "fr_graze_e"       ) call icepack_init_parameters(fr_graze_e_in        = paramValue)
    if (trim(paramName) == "fr_mort2min"      ) call icepack_init_parameters(fr_mort2min_in       = paramValue)
    if (trim(paramName) == "fr_dFe"           ) call icepack_init_parameters(fr_dFe_in            = paramValue)
    if (trim(paramName) == "k_nitrif"         ) call icepack_init_parameters(k_nitrif_in          = paramValue)
    if (trim(paramName) == "t_iron_conv"      ) call icepack_init_parameters(t_iron_conv_in       = paramValue)
    if (trim(paramName) == "max_loss"         ) call icepack_init_parameters(max_loss_in          = paramValue)
    if (trim(paramName) == "max_dfe_doc1"     ) call icepack_init_parameters(max_dfe_doc1_in      = paramValue)
    if (trim(paramName) == "fr_resp_s"        ) call icepack_init_parameters(fr_resp_s_in         = paramValue)
    if (trim(paramName) == "y_sk_DMS"         ) call icepack_init_parameters(y_sk_DMS_in          = paramValue)
    if (trim(paramName) == "t_sk_conv"        ) call icepack_init_parameters(t_sk_conv_in         = paramValue)
    if (trim(paramName) == "t_sk_ox"          ) call icepack_init_parameters(t_sk_ox_in           = paramValue)
    if (trim(paramName) == "frazil_scav"      ) call icepack_init_parameters(frazil_scav_in       = paramValue)
    if (trim(paramName) == "sk_l"             ) call icepack_init_parameters(sk_l_in              = paramValue)
    if (trim(paramName) == "min_bgc"          ) call icepack_init_parameters(min_bgc_in           = paramValue)
    if (trim(paramName) == "hs0"              ) call icepack_init_parameters(hs0_in               = paramValue)
    if (trim(paramName) == "dpscale"          ) call icepack_init_parameters(dpscale_in           = paramValue)
    if (trim(paramName) == "rfracmin"         ) call icepack_init_parameters(rfracmin_in          = paramValue)
    if (trim(paramName) == "rfracmax"         ) call icepack_init_parameters(rfracmax_in          = paramValue)
    if (trim(paramName) == "pndaspect"        ) call icepack_init_parameters(pndaspect_in         = paramValue)
    if (trim(paramName) == "hs1"              ) call icepack_init_parameters(hs1_in               = paramValue)
    if (trim(paramName) == "hp1"              ) call icepack_init_parameters(hp1_in               = paramValue)

  end subroutine column_init_parameters_real

  !-----------------------------------------------------------------------------

  subroutine column_init_parameters_integer(paramNameIn, paramValue) bind(C)

    use iso_c_binding, only: c_char, c_int, c_null_char

    use icepack_intfc, only: &
         icepack_init_parameters, &
         icepack_char_len_long

    character(kind=c_char, len=1), dimension(icepack_char_len_long), intent(in) :: &
         paramNameIn

    integer(kind=c_int), intent(in) :: &
         paramValue

    character(len=icepack_char_len_long) :: &
         paramName

    call convert_cstring_to_fstring(paramNameIn, paramName)

    if (trim(paramName) == "ktherm"     ) call icepack_init_parameters(ktherm_in      = paramValue)
    if (trim(paramName) == "kstrength"  ) call icepack_init_parameters(kstrength_in   = paramValue)
    if (trim(paramName) == "krdg_partic") call icepack_init_parameters(krdg_partic_in = paramValue)
    if (trim(paramName) == "krdg_redist") call icepack_init_parameters(krdg_redist_in = paramValue)
    if (trim(paramName) == "natmiter"   ) call icepack_init_parameters(natmiter_in    = paramValue)
    if (trim(paramName) == "kitd"       ) call icepack_init_parameters(kitd_in        = paramValue)
    if (trim(paramName) == "kcatbound"  ) call icepack_init_parameters(kcatbound_in   = paramValue)

  end subroutine column_init_parameters_integer

  !-----------------------------------------------------------------------------

  subroutine column_init_parameters_character(paramNameIn, paramValue) bind(C)

    use iso_c_binding, only: c_char, c_double, c_null_char

    use icepack_intfc, only: &
         icepack_dbl_kind, &
         icepack_init_parameters, &
         icepack_char_len_long

    character(kind=c_char, len=1), dimension(icepack_char_len_long), intent(in) :: &
         paramNameIn

    character(kind=c_char, len=1), dimension(icepack_char_len_long), intent(in) :: &
         paramValue

    character(len=icepack_char_len_long) :: &
         paramName, &
         fparamValue

    call convert_cstring_to_fstring(paramValue, fparamValue)

    call convert_cstring_to_fstring(paramNameIn, paramName)

    if (trim(paramName) == "conduct"       ) call icepack_init_parameters(conduct_in        = fparamValue)
    if (trim(paramName) == "fbot_xfer_type") call icepack_init_parameters(fbot_xfer_type_in = fparamValue)
    if (trim(paramName) == "tfrz_option"   ) call icepack_init_parameters(tfrz_option_in    = fparamValue)
    if (trim(paramName) == "shortwave"     ) call icepack_init_parameters(shortwave_in      = fparamValue)
    if (trim(paramName) == "albedo_type"   ) call icepack_init_parameters(albedo_type_in    = fparamValue)
    if (trim(paramName) == "atmbndy"       ) call icepack_init_parameters(atmbndy_in        = fparamValue)
    if (trim(paramName) == "bgc_flux_type" ) call icepack_init_parameters(bgc_flux_type_in  = fparamValue)
    if (trim(paramName) == "frzpnd"        ) call icepack_init_parameters(frzpnd_in         = fparamValue)

  end subroutine column_init_parameters_character

  !-----------------------------------------------------------------------------

  subroutine column_init_parameters_logical(paramNameIn, paramValue) bind(C)

    use iso_c_binding, only: c_char, c_int, c_null_char

    use icepack_intfc, only: &
         icepack_dbl_kind, &
         icepack_init_parameters, &
         icepack_char_len_long

    character(kind=c_char, len=1), dimension(icepack_char_len_long), intent(in) :: &
         paramNameIn

    integer(kind=c_int), intent(in) :: &
         paramValue

    character(len=icepack_char_len_long) :: &
         paramName

    logical :: &
         fparamValue

    fparamValue = convert_clogical_to_flogical(paramValue)

    call convert_cstring_to_fstring(paramNameIn, paramName)

    if (trim(paramName) == "heat_capacity") call icepack_init_parameters(heat_capacity_in = fparamValue)
    if (trim(paramName) == "calc_Tsfc"    ) call icepack_init_parameters(calc_Tsfc_in     = fparamValue)
    if (trim(paramName) == "update_ocn_f" ) call icepack_init_parameters(update_ocn_f_in  = fparamValue)
    if (trim(paramName) == "calc_strair"  ) call icepack_init_parameters(calc_strair_in   = fparamValue)
    if (trim(paramName) == "formdrag"     ) call icepack_init_parameters(formdrag_in      = fparamValue)
    if (trim(paramName) == "highfreq"     ) call icepack_init_parameters(highfreq_in      = fparamValue)
    if (trim(paramName) == "z_tracers"    ) call icepack_init_parameters(z_tracers_in     = fparamValue)
    if (trim(paramName) == "scale_bgc"    ) call icepack_init_parameters(scale_bgc_in     = fparamValue)
    if (trim(paramName) == "solve_zbgc"   ) call icepack_init_parameters(solve_zbgc_in    = fparamValue)
    if (trim(paramName) == "dEdd_algae"   ) call icepack_init_parameters(dEdd_algae_in    = fparamValue)
    if (trim(paramName) == "modal_aero"   ) call icepack_init_parameters(modal_aero_in    = fparamValue)
    if (trim(paramName) == "skl_bgc"      ) call icepack_init_parameters(skl_bgc_in       = fparamValue)
    if (trim(paramName) == "solve_zsal"   ) call icepack_init_parameters(solve_zsal_in    = fparamValue)

  end subroutine column_init_parameters_logical

  !-----------------------------------------------------------------------------

  subroutine column_query_parameters_real(paramNameIn, paramValue) bind(C)

    use iso_c_binding, only: c_char, c_double, c_null_char

    use icepack_intfc, only: &
         icepack_query_parameters, &
         icepack_char_len_long

    character(kind=c_char, len=1), dimension(icepack_char_len_long), intent(in) :: &
         paramNameIn

    real(kind=c_double), intent(out) :: &
         paramValue

    character(len=icepack_char_len_long) :: &
         paramName

    call convert_cstring_to_fstring(paramNameIn, paramName)

    if (trim(paramName) == "c0"               ) call icepack_query_parameters(c0_out                = paramValue)
    if (trim(paramName) == "c1"               ) call icepack_query_parameters(c1_out                = paramValue)
    if (trim(paramName) == "c1p5"             ) call icepack_query_parameters(c1p5_out              = paramValue)
    if (trim(paramName) == "c2"               ) call icepack_query_parameters(c2_out                = paramValue)
    if (trim(paramName) == "c3"               ) call icepack_query_parameters(c3_out                = paramValue)
    if (trim(paramName) == "c4"               ) call icepack_query_parameters(c4_out                = paramValue)
    if (trim(paramName) == "c5"               ) call icepack_query_parameters(c5_out                = paramValue)
    if (trim(paramName) == "c6"               ) call icepack_query_parameters(c6_out                = paramValue)
    if (trim(paramName) == "c8"               ) call icepack_query_parameters(c8_out                = paramValue)
    if (trim(paramName) == "c10"              ) call icepack_query_parameters(c10_out               = paramValue)
    if (trim(paramName) == "c15"              ) call icepack_query_parameters(c15_out               = paramValue)
    if (trim(paramName) == "c16"              ) call icepack_query_parameters(c16_out               = paramValue)
    if (trim(paramName) == "c20"              ) call icepack_query_parameters(c20_out               = paramValue)
    if (trim(paramName) == "c25"              ) call icepack_query_parameters(c25_out               = paramValue)
    if (trim(paramName) == "c180"             ) call icepack_query_parameters(c180_out              = paramValue)
    if (trim(paramName) == "c100"             ) call icepack_query_parameters(c100_out              = paramValue)
    if (trim(paramName) == "c1000"            ) call icepack_query_parameters(c1000_out             = paramValue)
    if (trim(paramName) == "p001"             ) call icepack_query_parameters(p001_out              = paramValue)
    if (trim(paramName) == "p01"              ) call icepack_query_parameters(p01_out               = paramValue)
    if (trim(paramName) == "p1"               ) call icepack_query_parameters(p1_out                = paramValue)
    if (trim(paramName) == "p2"               ) call icepack_query_parameters(p2_out                = paramValue)
    if (trim(paramName) == "p4"               ) call icepack_query_parameters(p4_out                = paramValue)
    if (trim(paramName) == "p5"               ) call icepack_query_parameters(p5_out                = paramValue)
    if (trim(paramName) == "p6"               ) call icepack_query_parameters(p6_out                = paramValue)
    if (trim(paramName) == "p05"              ) call icepack_query_parameters(p05_out               = paramValue)
    if (trim(paramName) == "p15"              ) call icepack_query_parameters(p15_out               = paramValue)
    if (trim(paramName) == "p25"              ) call icepack_query_parameters(p25_out               = paramValue)
    if (trim(paramName) == "p75"              ) call icepack_query_parameters(p75_out               = paramValue)
    if (trim(paramName) == "p333"             ) call icepack_query_parameters(p333_out              = paramValue)
    if (trim(paramName) == "p666"             ) call icepack_query_parameters(p666_out              = paramValue)
    if (trim(paramName) == "spval_const"      ) call icepack_query_parameters(spval_const_out       = paramValue)
    if (trim(paramName) == "pih"              ) call icepack_query_parameters(pih_out               = paramValue)
    if (trim(paramName) == "piq"              ) call icepack_query_parameters(piq_out               = paramValue)
    if (trim(paramName) == "pi2"              ) call icepack_query_parameters(pi2_out               = paramValue)
    if (trim(paramName) == "secday"           ) call icepack_query_parameters(secday_out            = paramValue)
    if (trim(paramName) == "puny"             ) call icepack_query_parameters(puny_out              = paramValue)
    if (trim(paramName) == "bignum"           ) call icepack_query_parameters(bignum_out            = paramValue)
    if (trim(paramName) == "pi"               ) call icepack_query_parameters(pi_out                = paramValue)
    if (trim(paramName) == "rad_to_deg"       ) call icepack_query_parameters(rad_to_deg_out        = paramValue)
    if (trim(paramName) == "Lfresh"           ) call icepack_query_parameters(Lfresh_out            = paramValue)
    if (trim(paramName) == "cprho"            ) call icepack_query_parameters(cprho_out             = paramValue)
    if (trim(paramName) == "Cp"               ) call icepack_query_parameters(Cp_out                = paramValue)
    if (trim(paramName) == "rhos"             ) call icepack_query_parameters(rhos_out              = paramValue)
    if (trim(paramName) == "rhoi"             ) call icepack_query_parameters(rhoi_out              = paramValue)
    if (trim(paramName) == "rhosi"            ) call icepack_query_parameters(rhosi_out             = paramValue)
    if (trim(paramName) == "rhow"             ) call icepack_query_parameters(rhow_out              = paramValue)
    if (trim(paramName) == "rhofresh"         ) call icepack_query_parameters(rhofresh_out          = paramValue)
    if (trim(paramName) == "cp_ice"           ) call icepack_query_parameters(cp_ice_out            = paramValue)
    if (trim(paramName) == "cp_ocn"           ) call icepack_query_parameters(cp_ocn_out            = paramValue)
    if (trim(paramName) == "depressT"         ) call icepack_query_parameters(depressT_out          = paramValue)
    if (trim(paramName) == "viscosity_dyn"    ) call icepack_query_parameters(viscosity_dyn_out     = paramValue)
    if (trim(paramName) == "Tocnfrz"          ) call icepack_query_parameters(Tocnfrz_out           = paramValue)
    if (trim(paramName) == "Tffresh"          ) call icepack_query_parameters(Tffresh_out           = paramValue)
    if (trim(paramName) == "Lsub"             ) call icepack_query_parameters(Lsub_out              = paramValue)
    if (trim(paramName) == "Lvap"             ) call icepack_query_parameters(Lvap_out              = paramValue)
    if (trim(paramName) == "Timelt"           ) call icepack_query_parameters(Timelt_out            = paramValue)
    if (trim(paramName) == "Tsmelt"           ) call icepack_query_parameters(Tsmelt_out            = paramValue)
    if (trim(paramName) == "ice_ref_salinity" ) call icepack_query_parameters(ice_ref_salinity_out  = paramValue)
    if (trim(paramName) == "kice"             ) call icepack_query_parameters(kice_out              = paramValue)
    if (trim(paramName) == "kseaice"          ) call icepack_query_parameters(kseaice_out           = paramValue)
    if (trim(paramName) == "ksno"             ) call icepack_query_parameters(ksno_out              = paramValue)
    if (trim(paramName) == "hs_min"           ) call icepack_query_parameters(hs_min_out            = paramValue)
    if (trim(paramName) == "snowpatch"        ) call icepack_query_parameters(snowpatch_out         = paramValue)
    if (trim(paramName) == "saltmax"          ) call icepack_query_parameters(saltmax_out           = paramValue)
    if (trim(paramName) == "phi_init"         ) call icepack_query_parameters(phi_init_out          = paramValue)
    if (trim(paramName) == "min_salin"        ) call icepack_query_parameters(min_salin_out         = paramValue)
    if (trim(paramName) == "salt_loss"        ) call icepack_query_parameters(salt_loss_out         = paramValue)
    if (trim(paramName) == "dSin0_frazil"     ) call icepack_query_parameters(dSin0_frazil_out      = paramValue)
    if (trim(paramName) == "dts_b"            ) call icepack_query_parameters(dts_b_out             = paramValue)
    if (trim(paramName) == "ustar_min"        ) call icepack_query_parameters(ustar_min_out         = paramValue)
    if (trim(paramName) == "a_rapid_mode"     ) call icepack_query_parameters(a_rapid_mode_out      = paramValue)
    if (trim(paramName) == "Rac_rapid_mode"   ) call icepack_query_parameters(Rac_rapid_mode_out    = paramValue)
    if (trim(paramName) == "aspect_rapid_mode") call icepack_query_parameters(aspect_rapid_mode_out = paramValue)
    if (trim(paramName) == "dSdt_slow_mode"   ) call icepack_query_parameters(dSdt_slow_mode_out    = paramValue)
    if (trim(paramName) == "phi_c_slow_mode"  ) call icepack_query_parameters(phi_c_slow_mode_out   = paramValue)
    if (trim(paramName) == "phi_i_mushy"      ) call icepack_query_parameters(phi_i_mushy_out       = paramValue)
    if (trim(paramName) == "emissivity"       ) call icepack_query_parameters(emissivity_out        = paramValue)
    if (trim(paramName) == "albocn"           ) call icepack_query_parameters(albocn_out            = paramValue)
    if (trim(paramName) == "vonkar"           ) call icepack_query_parameters(vonkar_out            = paramValue)
    if (trim(paramName) == "stefan_boltzmann" ) call icepack_query_parameters(stefan_boltzmann_out  = paramValue)
    if (trim(paramName) == "kappav"           ) call icepack_query_parameters(kappav_out            = paramValue)
    if (trim(paramName) == "hi_ssl"           ) call icepack_query_parameters(hi_ssl_out            = paramValue)
    if (trim(paramName) == "hs_ssl"           ) call icepack_query_parameters(hs_ssl_out            = paramValue)
    if (trim(paramName) == "awtvdr"           ) call icepack_query_parameters(awtvdr_out            = paramValue)
    if (trim(paramName) == "awtidr"           ) call icepack_query_parameters(awtidr_out            = paramValue)
    if (trim(paramName) == "awtvdf"           ) call icepack_query_parameters(awtvdf_out            = paramValue)
    if (trim(paramName) == "awtidf"           ) call icepack_query_parameters(awtidf_out            = paramValue)
    if (trim(paramName) == "albicev"          ) call icepack_query_parameters(albicev_out           = paramValue)
    if (trim(paramName) == "albicei"          ) call icepack_query_parameters(albicei_out           = paramValue)
    if (trim(paramName) == "albsnowv"         ) call icepack_query_parameters(albsnowv_out          = paramValue)
    if (trim(paramName) == "albsnowi"         ) call icepack_query_parameters(albsnowi_out          = paramValue)
    if (trim(paramName) == "ahmax"            ) call icepack_query_parameters(ahmax_out             = paramValue)
    if (trim(paramName) == "R_ice"            ) call icepack_query_parameters(R_ice_out             = paramValue)
    if (trim(paramName) == "R_pnd"            ) call icepack_query_parameters(R_pnd_out             = paramValue)
    if (trim(paramName) == "R_snw"            ) call icepack_query_parameters(R_snw_out             = paramValue)
    if (trim(paramName) == "dT_mlt"           ) call icepack_query_parameters(dT_mlt_out            = paramValue)
    if (trim(paramName) == "rsnw_mlt"         ) call icepack_query_parameters(rsnw_mlt_out          = paramValue)
    if (trim(paramName) == "kalg"             ) call icepack_query_parameters(kalg_out              = paramValue)
    if (trim(paramName) == "Cf"               ) call icepack_query_parameters(Cf_out                = paramValue)
    if (trim(paramName) == "Pstar"            ) call icepack_query_parameters(Pstar_out             = paramValue)
    if (trim(paramName) == "Cstar"            ) call icepack_query_parameters(Cstar_out             = paramValue)
    if (trim(paramName) == "dragio"           ) call icepack_query_parameters(dragio_out            = paramValue)
    if (trim(paramName) == "gravit"           ) call icepack_query_parameters(gravit_out            = paramValue)
    if (trim(paramName) == "iceruf"           ) call icepack_query_parameters(iceruf_out            = paramValue)
    if (trim(paramName) == "mu_rdg"           ) call icepack_query_parameters(mu_rdg_out            = paramValue)
    if (trim(paramName) == "cp_air"           ) call icepack_query_parameters(cp_air_out            = paramValue)
    if (trim(paramName) == "cp_wv"            ) call icepack_query_parameters(cp_wv_out             = paramValue)
    if (trim(paramName) == "zvir"             ) call icepack_query_parameters(zvir_out              = paramValue)
    if (trim(paramName) == "zref"             ) call icepack_query_parameters(zref_out              = paramValue)
    if (trim(paramName) == "qqqice"           ) call icepack_query_parameters(qqqice_out            = paramValue)
    if (trim(paramName) == "TTTice"           ) call icepack_query_parameters(TTTice_out            = paramValue)
    if (trim(paramName) == "qqqocn"           ) call icepack_query_parameters(qqqocn_out            = paramValue)
    if (trim(paramName) == "TTTocn"           ) call icepack_query_parameters(TTTocn_out            = paramValue)
    if (trim(paramName) == "grid_o"           ) call icepack_query_parameters(grid_o_out            = paramValue)
    if (trim(paramName) == "l_sk"             ) call icepack_query_parameters(l_sk_out              = paramValue)
    if (trim(paramName) == "initbio_frac"     ) call icepack_query_parameters(initbio_frac_out      = paramValue)
    if (trim(paramName) == "phi_snow"         ) call icepack_query_parameters(phi_snow_out          = paramValue)
    if (trim(paramName) == "grid_oS"          ) call icepack_query_parameters(grid_oS_out           = paramValue)
    if (trim(paramName) == "l_skS"            ) call icepack_query_parameters(l_skS_out             = paramValue)
    if (trim(paramName) == "fr_resp"          ) call icepack_query_parameters(fr_resp_out           = paramValue)
    if (trim(paramName) == "algal_vel"        ) call icepack_query_parameters(algal_vel_out         = paramValue)
    if (trim(paramName) == "R_dFe2dust"       ) call icepack_query_parameters(R_dFe2dust_out        = paramValue)
    if (trim(paramName) == "dustFe_sol"       ) call icepack_query_parameters(dustFe_sol_out        = paramValue)
    if (trim(paramName) == "T_max"            ) call icepack_query_parameters(T_max_out             = paramValue)
    if (trim(paramName) == "fsal"             ) call icepack_query_parameters(fsal_out              = paramValue)
    if (trim(paramName) == "op_dep_min"       ) call icepack_query_parameters(op_dep_min_out        = paramValue)
    if (trim(paramName) == "fr_graze_s"       ) call icepack_query_parameters(fr_graze_s_out        = paramValue)
    if (trim(paramName) == "fr_graze_e"       ) call icepack_query_parameters(fr_graze_e_out        = paramValue)
    if (trim(paramName) == "fr_mort2min"      ) call icepack_query_parameters(fr_mort2min_out       = paramValue)
    if (trim(paramName) == "fr_dFe"           ) call icepack_query_parameters(fr_dFe_out            = paramValue)
    if (trim(paramName) == "k_nitrif"         ) call icepack_query_parameters(k_nitrif_out          = paramValue)
    if (trim(paramName) == "t_iron_conv"      ) call icepack_query_parameters(t_iron_conv_out       = paramValue)
    if (trim(paramName) == "max_loss"         ) call icepack_query_parameters(max_loss_out          = paramValue)
    if (trim(paramName) == "max_dfe_doc1"     ) call icepack_query_parameters(max_dfe_doc1_out      = paramValue)
    if (trim(paramName) == "fr_resp_s"        ) call icepack_query_parameters(fr_resp_s_out         = paramValue)
    if (trim(paramName) == "y_sk_DMS"         ) call icepack_query_parameters(y_sk_DMS_out          = paramValue)
    if (trim(paramName) == "t_sk_conv"        ) call icepack_query_parameters(t_sk_conv_out         = paramValue)
    if (trim(paramName) == "t_sk_ox"          ) call icepack_query_parameters(t_sk_ox_out           = paramValue)
    if (trim(paramName) == "frazil_scav"      ) call icepack_query_parameters(frazil_scav_out       = paramValue)
    if (trim(paramName) == "sk_l"             ) call icepack_query_parameters(sk_l_out              = paramValue)
    if (trim(paramName) == "min_bgc"          ) call icepack_query_parameters(min_bgc_out           = paramValue)
    if (trim(paramName) == "hs0"              ) call icepack_query_parameters(hs0_out               = paramValue)
    if (trim(paramName) == "dpscale"          ) call icepack_query_parameters(dpscale_out           = paramValue)
    if (trim(paramName) == "rfracmin"         ) call icepack_query_parameters(rfracmin_out          = paramValue)
    if (trim(paramName) == "rfracmax"         ) call icepack_query_parameters(rfracmax_out          = paramValue)
    if (trim(paramName) == "pndaspect"        ) call icepack_query_parameters(pndaspect_out         = paramValue)
    if (trim(paramName) == "hs1"              ) call icepack_query_parameters(hs1_out               = paramValue)
    if (trim(paramName) == "hp1"              ) call icepack_query_parameters(hp1_out               = paramValue)

  end subroutine column_query_parameters_real

  !-----------------------------------------------------------------------------

  subroutine column_query_parameters_integer(paramNameIn, paramValue) bind(C)

    use iso_c_binding, only: c_char, c_int, c_null_char

    use icepack_intfc, only: &
         icepack_query_parameters, &
         icepack_char_len_long

    character(kind=c_char, len=1), dimension(icepack_char_len_long), intent(in) :: &
         paramNameIn

    integer(kind=c_int), intent(out) :: &
         paramValue

    character(len=icepack_char_len_long) :: &
         paramName

    call convert_cstring_to_fstring(paramNameIn, paramName)

    if (trim(paramName) == "ktherm"     ) call icepack_query_parameters(ktherm_out      = paramValue)
    if (trim(paramName) == "kstrength"  ) call icepack_query_parameters(kstrength_out   = paramValue)
    if (trim(paramName) == "krdg_partic") call icepack_query_parameters(krdg_partic_out = paramValue)
    if (trim(paramName) == "krdg_redist") call icepack_query_parameters(krdg_redist_out = paramValue)
    if (trim(paramName) == "natmiter"   ) call icepack_query_parameters(natmiter_out    = paramValue)
    if (trim(paramName) == "kitd"       ) call icepack_query_parameters(kitd_out        = paramValue)
    if (trim(paramName) == "kcatbound"  ) call icepack_query_parameters(kcatbound_out   = paramValue)

  end subroutine column_query_parameters_integer

  !-----------------------------------------------------------------------------

  subroutine column_query_parameters_character(paramNameIn, paramValue) bind(C)

    use iso_c_binding, only: c_char, c_double, c_null_char

    use icepack_intfc, only: &
         icepack_dbl_kind, &
         icepack_query_parameters, &
         icepack_char_len_long

    character(kind=c_char, len=1), dimension(icepack_char_len_long), intent(in) :: &
         paramNameIn

    character(kind=c_char, len=1), dimension(icepack_char_len_long), intent(out) :: &
         paramValue

    character(len=icepack_char_len_long) :: &
         paramName, &
         fparamValue

    call convert_cstring_to_fstring(paramNameIn, paramName)

    if (trim(paramName) == "conduct"       ) call icepack_query_parameters(conduct_out        = fparamValue)
    if (trim(paramName) == "fbot_xfer_type") call icepack_query_parameters(fbot_xfer_type_out = fparamValue)
    if (trim(paramName) == "tfrz_option"   ) call icepack_query_parameters(tfrz_option_out    = fparamValue)
    if (trim(paramName) == "shortwave"     ) call icepack_query_parameters(shortwave_out      = fparamValue)
    if (trim(paramName) == "albedo_type"   ) call icepack_query_parameters(albedo_type_out    = fparamValue)
    if (trim(paramName) == "atmbndy"       ) call icepack_query_parameters(atmbndy_out        = fparamValue)
    if (trim(paramName) == "bgc_flux_type" ) call icepack_query_parameters(bgc_flux_type_out  = fparamValue)
    if (trim(paramName) == "frzpnd"        ) call icepack_query_parameters(frzpnd_out         = fparamValue)

    call convert_fstring_to_cstring(fparamValue, paramValue)

  end subroutine column_query_parameters_character

  !-----------------------------------------------------------------------------

  subroutine column_query_parameters_logical(paramNameIn, paramValue) bind(C)

    use iso_c_binding, only: c_char, c_int, c_null_char

    use icepack_intfc, only: &
         icepack_dbl_kind, &
         icepack_query_parameters, &
         icepack_char_len_long

    character(kind=c_char, len=1), dimension(icepack_char_len_long), intent(in) :: &
         paramNameIn

    integer(kind=c_int), intent(out) :: &
         paramValue

    character(len=icepack_char_len_long) :: &
         paramName

    logical :: &
         fparamValue

    call convert_cstring_to_fstring(paramNameIn, paramName)

    if (trim(paramName) == "heat_capacity") call icepack_query_parameters(heat_capacity_out = fparamValue)
    if (trim(paramName) == "calc_Tsfc"    ) call icepack_query_parameters(calc_Tsfc_out     = fparamValue)
    if (trim(paramName) == "update_ocn_f" ) call icepack_query_parameters(update_ocn_f_out  = fparamValue)
    if (trim(paramName) == "calc_strair"  ) call icepack_query_parameters(calc_strair_out   = fparamValue)
    if (trim(paramName) == "formdrag"     ) call icepack_query_parameters(formdrag_out      = fparamValue)
    if (trim(paramName) == "highfreq"     ) call icepack_query_parameters(highfreq_out      = fparamValue)
    if (trim(paramName) == "z_tracers"    ) call icepack_query_parameters(z_tracers_out     = fparamValue)
    if (trim(paramName) == "scale_bgc"    ) call icepack_query_parameters(scale_bgc_out     = fparamValue)
    if (trim(paramName) == "solve_zbgc"   ) call icepack_query_parameters(solve_zbgc_out    = fparamValue)
    if (trim(paramName) == "dEdd_algae"   ) call icepack_query_parameters(dEdd_algae_out    = fparamValue)
    if (trim(paramName) == "modal_aero"   ) call icepack_query_parameters(modal_aero_out    = fparamValue)
    if (trim(paramName) == "skl_bgc"      ) call icepack_query_parameters(skl_bgc_out       = fparamValue)
    if (trim(paramName) == "solve_zsal"   ) call icepack_query_parameters(solve_zsal_out    = fparamValue)

    paramValue = convert_flogical_to_clogical(fparamValue)

  end subroutine column_query_parameters_logical

  !-----------------------------------------------------------------------------

  function column_get_constant(paramNameIn) result(paramValue) bind(C)

    use iso_c_binding, only: c_char, c_double

    use icepack_intfc, only: &
         icepack_char_len_long

    use icepack_parameters, only: &
         dragio, &
         cp_air, &
         Lsub, &
         rhos, &
         rhoi, &
         Tsmelt, &
         Tffresh

    character(kind=c_char, len=1), dimension(icepack_char_len_long), intent(in) :: &
         paramNameIn

    real(c_double) :: &
         paramValue

    character(len=icepack_char_len_long) :: &
         paramName

    call convert_cstring_to_fstring(paramNameIn, paramName)

    if (paramName == "dragio") then
       paramValue = dragio
    else if (paramName == "cp_air") then
       paramValue = cp_air
    else if (paramName == "Lsub") then
       paramValue = Lsub
    else if (paramName == "rhos") then
       paramValue = rhos
    else if (paramName == "rhoi") then
       paramValue = rhoi
    else if (paramName == "Tsmelt") then
       paramValue = Tsmelt
    else if (paramName == "Tffresh") then
       paramValue = Tffresh
    endif

  end function column_get_constant

  !-----------------------------------------------------------------------------

  function column_liquidus_temperature(salinity) result(liquidusTemperature) bind(C)

    use iso_c_binding, only: c_double

    use icepack_intfc, only: &
         icepack_liquidus_temperature

    real(kind=c_double), intent(in) :: &
         salinity

    real(kind=c_double) :: &
         liquidusTemperature

    liquidusTemperature = icepack_liquidus_temperature(salinity)

  end function column_liquidus_temperature

  !-----------------------------------------------------------------------------
  ! tracers
  !-----------------------------------------------------------------------------

  subroutine column_init_tracer_flags(tracerNameIn, useTracerValue) bind(C)

    use iso_c_binding, only: c_char, c_int, c_null_char

    use icepack_intfc, only: &
         icepack_init_tracer_flags, &
         icepack_char_len_long

    character(kind=c_char, len=1), dimension(icepack_char_len_long), intent(in) :: &
         tracerNameIn

    integer(kind=c_int), intent(in) :: &
         useTracerValue

    character(len=icepack_char_len_long) :: &
         tracerName

    logical :: &
         fUseTracerValue

    fUseTracerValue = convert_clogical_to_flogical(useTracerValue)

    call convert_cstring_to_fstring(tracerNameIn, tracerName)

    if (trim(tracerName) == "useIceAge"        ) call icepack_init_tracer_flags(tr_iage_in      = fUseTracerValue)
    if (trim(tracerName) == "useFirstYearIce"  ) call icepack_init_tracer_flags(tr_FY_in        = fUseTracerValue)
    if (trim(tracerName) == "useLevelIce"      ) call icepack_init_tracer_flags(tr_lvl_in       = fUseTracerValue)
    if (trim(tracerName) == "useMeltPonds"     ) call icepack_init_tracer_flags(tr_pond_in      = fUseTracerValue)
    if (trim(tracerName) == "useCesmMeltponds" ) call icepack_init_tracer_flags(tr_pond_cesm_in = fUseTracerValue)
    if (trim(tracerName) == "useLevelMeltponds") call icepack_init_tracer_flags(tr_pond_lvl_in  = fUseTracerValue)
    if (trim(tracerName) == "useTopoMeltponds" ) call icepack_init_tracer_flags(tr_pond_topo_in = fUseTracerValue)
    if (trim(tracerName) == "useAerosols"      ) call icepack_init_tracer_flags(tr_aero_in      = fUseTracerValue)

  end subroutine column_init_tracer_flags

  !-----------------------------------------------------------------------------

  subroutine column_init_tracer_numbers(nTracers) bind(C)

    use iso_c_binding, only: c_int

    use icepack_intfc, only: &
         icepack_init_tracer_numbers

    integer(kind=c_int), intent(in) :: &
         nTracers

    call icepack_init_tracer_numbers(ntrcr_in = nTracers)

  end subroutine column_init_tracer_numbers

  !-----------------------------------------------------------------------------

  subroutine column_init_tracer_indices(&
       indexSurfaceTemperature, &
       indexIceEnthalpy, &
       indexSnowEnthalpy, &
       indexIceSalinity, &
       indexIceAge, &
       indexFirstYearIceArea, &
       indexLevelIceArea, &
       indexLevelIceVolume, &
       indexPondArea, &
       indexPondDepth, &
       indexPondLidThickness, &
       indexAerosols) bind(C)

    use iso_c_binding, only: c_int

    use icepack_intfc, only: &
         icepack_init_tracer_indices

    integer(kind=c_int), intent(in) :: &
         indexSurfaceTemperature, &
         indexIceEnthalpy, &
         indexSnowEnthalpy, &
         indexIceSalinity, &
         indexIceAge, &
         indexFirstYearIceArea, &
         indexLevelIceArea, &
         indexLevelIceVolume, &
         indexPondArea, &
         indexPondDepth, &
         indexPondLidThickness, &
         indexAerosols

    call icepack_init_tracer_indices(&
         nt_Tsfc_in = indexSurfaceTemperature, &
         nt_qice_in = indexIceEnthalpy, &
         nt_qsno_in = indexSnowEnthalpy, &
         nt_sice_in = indexIceSalinity, &
         nt_iage_in = indexIceAge, &
         nt_FY_in   = indexFirstYearIceArea, &
         nt_alvl_in = indexLevelIceArea, &
         nt_vlvl_in = indexLevelIceVolume, &
         nt_apnd_in = indexPondArea, &
         nt_hpnd_in = indexPondDepth, &
         nt_ipnd_in = indexPondLidThickness, &
         nt_aero_in = indexAerosols)

  end subroutine column_init_tracer_indices

  !-----------------------------------------------------------------------------

  subroutine column_aggregate(&
       nCategories, &
       nTracers, &
       nBaseTracers, &
       nMaxAncestorTracers, &
       parentIndex, & ! trcr_depend
       firstAncestorMask, & ! trcr_base
       ancestorNumber, & ! n_trcr_strata
       ancestorIndices, &
       tracerArrayCategory, & ! trcrn
       tracerArray, & ! trcr
       iceAreaCategory, &
       iceVolumeCategory, &
       snowVolumeCategory, &
       iceAreaCell, &
       iceVolumeCell, &
       snowVolumeCell, &
       openWaterArea) bind(C)

    use iso_c_binding, only: c_int, c_double

    use icepack_intfc, only: &
         icepack_aggregate

    integer(c_int) :: &
         nCategories, &
         nTracers, &
         nBaseTracers, &
         nMaxAncestorTracers

    real(kind=c_double), dimension(nTracers,nCategories) :: &
         tracerArrayCategory

    real(kind=c_double), dimension(nTracers) :: &
         tracerArray

    integer(c_int), dimension(nTracers) :: &
         parentIndex, &
         ancestorNumber

    real(kind=c_double), dimension(nTracers, nBaseTracers) :: &
         firstAncestorMask

    integer(c_int), dimension(nTracers, nMaxAncestorTracers) :: &
         ancestorIndices

    real(kind=c_double), dimension(nCategories) :: &
         iceAreaCategory, &
         iceVolumeCategory, &
         snowVolumeCategory

    real(kind=c_double) :: &
         iceAreaCell, &
         iceVolumeCell, &
         snowVolumeCell, &
         openWaterArea

    call icepack_aggregate(&
         nCategories, &
         iceAreaCategory(:), &
         tracerArrayCategory, & ! trcrn
         iceVolumeCategory(:), &
         snowVolumeCategory(:), &
         iceAreaCell, &
         tracerArray, & ! trcr
         iceVolumeCell, &
         snowVolumeCell, &
         openWaterArea, &
         nTracers, &
         parentIndex, & ! trcr_depend
         firstAncestorMask, & ! trcr_base
         ancestorNumber, & ! n_trcr_strata
         ancestorIndices) ! nt_strata

  end subroutine column_aggregate

  !-----------------------------------------------------------------------------

  subroutine column_print_tracer_object(&
       nTracers, &
       nBaseTracers, &
       nMaxAncestorTracers, &
       indexSurfaceTemperature, &
       indexIceEnthalpy, &
       indexSnowEnthalpy, &
       indexIceSalinity, &
       indexIceAge, &
       indexFirstYearIceArea, &
       indexLevelIceArea, &
       indexLevelIceVolume, &
       indexPondArea, &
       indexPondDepth, &
       indexPondLidThickness, &
       indexAerosols, &
       parentIndex, &
       firstAncestorMask, &
       ancestorIndices, &
       ancestorNumber) bind(C)

    use iso_c_binding, only: c_int, c_double

    integer(c_int), intent(in) :: &
         nTracers, &
         nBaseTracers, &
         nMaxAncestorTracers, &
         indexSurfaceTemperature, &
         indexIceEnthalpy, &
         indexSnowEnthalpy, &
         indexIceSalinity, &
         indexIceAge, &
         indexFirstYearIceArea, &
         indexLevelIceArea, &
         indexLevelIceVolume, &
         indexPondArea, &
         indexPondDepth, &
         indexPondLidThickness, &
         indexAerosols

    integer(c_int), dimension(nTracers), intent(in) :: &
         parentIndex, &
         ancestorNumber

    real(c_double), dimension(nTracers, nBaseTracers), intent(in) :: &
         firstAncestorMask

    integer(c_int), dimension(nTracers, nMaxAncestorTracers), intent(in) :: &
         ancestorIndices

    integer :: &
         i, j

    write(*,*) "----------------------------------------------------"
    write(*,*) "tracerObject:"
    write(*,*) ""

    write(*,*) "nTracers:            ", nTracers
    write(*,*) "nBaseTracers:        ", nBaseTracers
    write(*,*) "nMaxAncestorTracers: ", nMaxAncestorTracers
    write(*,*) ""

    write(*,*) "index_surfaceTemperature: ", indexSurfaceTemperature
    write(*,*) "index_iceEnthalpy:        ", indexIceEnthalpy
    write(*,*) "index_snowEnthalpy:       ", indexSnowEnthalpy
    write(*,*) "index_iceSalinity:        ", indexIceSalinity
    write(*,*) "index_iceAge:             ", indexIceAge
    write(*,*) "index_firstYearIceArea:   ", indexFirstYearIceArea
    write(*,*) "index_levelIceArea:       ", indexLevelIceArea
    write(*,*) "index_levelIceVolume:     ", indexLevelIceVolume
    write(*,*) "index_pondArea:           ", indexPondArea
    write(*,*) "index_pondDepth:          ", indexPondDepth
    write(*,*) "index_pondLidThickness:   ", indexPondLidThickness
    write(*,*) "index_aerosols:           ", indexAerosols
    write(*,*) ""

    write(*,*) "parentIndex: ", size(parentIndex,1)
    do i = 1, size(parentIndex,1)
       write(*,*) "  ", i, parentIndex(i)
    enddo ! i
    write(*,*) ""

    write(*,*) "firstAncestorMask: ", size(firstAncestorMask,1),size(firstAncestorMask,2)
    do i = 1, size(firstAncestorMask,1)
       do j = 1, size(firstAncestorMask,2)
          write(*,*) "  ", i, j, firstAncestorMask(i,j)
       enddo ! j
    enddo ! i
    write(*,*) ""

    write(*,*) "ancestorIndices: ", size(ancestorIndices,1),size(ancestorIndices,2)
    do i = 1, size(ancestorIndices,1)
       do j = 1, size(ancestorIndices,2)
          write(*,*) "  ", i, j, ancestorIndices(i,j)
       enddo ! j
    enddo ! i
    write(*,*) ""

    write(*,*) "ancestorNumber: ", size(ancestorNumber,1)
    do i = 1, size(ancestorNumber,1)
       write(*,*) "  ", i, ancestorNumber(i)
    enddo ! i
    write(*,*) ""

    write(*,*) "----------------------------------------------------"

  end subroutine column_print_tracer_object

  !-----------------------------------------------------------------------------
  ! forcing
  !-----------------------------------------------------------------------------

  subroutine column_limit_temperature(&
       iceAreaCell, &
       airTemperature) bind(C)

    use iso_c_binding, only: c_double

    use icepack_intfc, only: &
         icepack_dbl_kind

    use icepack_parameters, only: &
         freshWaterFreezingPoint => Tffresh

    real(c_double) :: &
         iceAreaCell, &
         airTemperature

    ! limit air temperature values where ice is present
    if (iceAreaCell > 0.1_icepack_dbl_kind) then
       airTemperature = min(airTemperature, freshWaterFreezingPoint + 0.1_icepack_dbl_kind)
    endif

  end subroutine column_limit_temperature

  !-----------------------------------------------------------------------------

  subroutine column_limit_specific_humidity(&
       airTemperature, &
       airSpecificHumidity) bind(C)

    use iso_c_binding, only: c_double

    use icepack_intfc, only: &
         icepack_dbl_kind

    use icepack_parameters, only: &
         freshWaterFreezingPoint => Tffresh, &
         puny

    real(c_double), intent(in) :: &
         airTemperature

    real(c_double), intent(inout) :: &
         airSpecificHumidity

    real(c_double) :: &
         airSpecificHumidityMax

    ! convert air temperature from Kelvin to Celcius
    airSpecificHumidityMax = airTemperature - freshWaterFreezingPoint

    airSpecificHumidityMax = 2.0_icepack_dbl_kind + &
         ((0.7859_icepack_dbl_kind + 0.03477_icepack_dbl_kind * airSpecificHumidityMax) / &
          (1.0_icepack_dbl_kind    + 0.00412_icepack_dbl_kind * airSpecificHumidityMax)) + &
          0.00422_icepack_dbl_kind * airSpecificHumidityMax

    ! vapor pressure
    airSpecificHumidityMax = 10.0_icepack_dbl_kind**airSpecificHumidityMax ! saturated
    airSpecificHumidityMax = max(airSpecificHumidityMax,puny) ! prevent division by zero

    ! specific humidity
    airSpecificHumidityMax = &
         (0.622_icepack_dbl_kind * airSpecificHumidityMax) / &
         (1.0e5_icepack_dbl_kind - 0.378_icepack_dbl_kind * airSpecificHumidityMax)

    ! limit the specific humidity
    airSpecificHumidity = min(airSpecificHumidity, airSpecificHumidityMax)

  end subroutine column_limit_specific_humidity

  !-----------------------------------------------------------------------------

  subroutine column_longwave_rosati_miyakoda(&
       longwaveDown, &
       cloudFraction, &
       iceAreaCell, &
       surfaceTemperature, &
       seaSurfaceTemperature, &
       airSpecificHumidity, &
       airTemperature) bind(C)

    ! Rosati, A. and K. Miyakoda (1988),
    ! A general-circulation model for upper ocean simulation,
    ! J. Physical Oceanography, 18, 1601-1626,
    ! doi:10.1175/1520-0485(1988)018<1601:AGCMFU>2.0.CO;2

    use iso_c_binding, only: c_double

    use icepack_intfc, only: &
         icepack_dbl_kind

    use icepack_parameters, only: &
         stefanBoltzmann => stefan_boltzmann, &
         freshWaterFreezingPoint => Tffresh, &
         iceSnowEmissivity => emissivity

    real(kind=c_double), intent(out) :: &
         longwaveDown

    real(kind=c_double), intent(in) :: &
         cloudFraction, &
         iceAreaCell, &
         surfaceTemperature, &
         seaSurfaceTemperature, &
         airSpecificHumidity, &
         airTemperature

    real(kind=c_double) :: &
         clearSkyFraction, &
         combinedSurfaceTemperature, &
         vapourPressureSqrt, &
         airPotentialTemperature, &
         airSeaTemperatureDifferenceTerm

    ! get a clear sky fraction
    clearSkyFraction = 1.0_icepack_dbl_kind - 0.8_icepack_dbl_kind * cloudFraction

    ! combined ice and ocean temperature in Kelvin
    combinedSurfaceTemperature = &
         surfaceTemperature    * iceAreaCell + &
         seaSurfaceTemperature * (1.0_icepack_dbl_kind - iceAreaCell) + &
         freshWaterFreezingPoint

    ! square root of the vapour pressure
    vapourPressureSqrt = sqrt((1000.0_icepack_dbl_kind * airSpecificHumidity) / &
                              (0.622_icepack_dbl_kind + 0.378_icepack_dbl_kind * airSpecificHumidity))

    ! potential temperature (CICE comment: get this from stability?)
    airPotentialTemperature = airTemperature

    ! unknown
    !airSeaTemperatureDifferenceTerm = airPotentialTemperature**3 * &
    !     (airPotentialTemperature * (0.39_icepack_dbl_kind - 0.05_icepack_dbl_kind * vapourPressureSqrt) * clearSkyFraction + &
    !     4.0_icepack_dbl_kind * (combinedSurfaceTemperature - airPotentialTemperature))

    ! Different version for bfb to CICE
    airSeaTemperatureDifferenceTerm = airPotentialTemperature*airPotentialTemperature*airPotentialTemperature * &
         (airPotentialTemperature * (0.39_icepack_dbl_kind - 0.05_icepack_dbl_kind * vapourPressureSqrt) * clearSkyFraction + &
         4.0_icepack_dbl_kind * (combinedSurfaceTemperature - airPotentialTemperature))

    ! final longwave calculation from stefan-boltzmann law
    longwaveDown = iceSnowEmissivity * stefanBoltzmann * &
         (combinedSurfaceTemperature**4 - airSeaTemperatureDifferenceTerm)

  end subroutine column_longwave_rosati_miyakoda

  !-----------------------------------------------------------------------------

  subroutine column_longwave_parkinson_and_washington(&
       longWavedown, &
       airTemperature, &
       cloudFraction) bind(C)

    use iso_c_binding, only: c_double

    use icepack_intfc, only: &
         icepack_dbl_kind

    use icepack_parameters, only: &
         stefanBoltzmann => stefan_boltzmann, &
         freshWaterFreezingPoint => Tffresh

    real(c_double), intent(out) :: &
         longWavedown

    real(c_double), intent(in) :: &
         airTemperature, &
         cloudFraction

    ! Longwave down
    ! Parkinson, C. L. and W. M. Washington (1979),
    ! Large-scale numerical-model of sea ice,
    ! JGR, 84, 311-337, doi:10.1029/JC084iC01p00311

    longWavedown = &
         stefanBoltzmann * airTemperature**4 * &
         (1.0_icepack_dbl_kind - 0.261_icepack_dbl_kind * exp(-7.77e-4_icepack_dbl_kind * (freshWaterFreezingPoint - airTemperature)**2)) * &
         (1.0_icepack_dbl_kind + 0.275_icepack_dbl_kind * cloudFraction)

  end subroutine column_longwave_parkinson_and_washington

  !-----------------------------------------------------------------------------

  subroutine column_shortwave_from_cloud_fraction(&
       shortwaveDown, &
       longitudeIn, &
       latitude, &
       cloudFraction, &
       airSpecificHumidity, &
       secondsToday, &
       dayOfYear) bind(C)

    use iso_c_binding, only: c_double, c_int

    use icepack_intfc, only: &
         icepack_dbl_kind

    use icepack_parameters, only: &
         pi

    real(c_double), intent(out) :: &
         shortwaveDown

    real(c_double), intent(in) :: &
         longitudeIn, &
         latitude, &
         cloudFraction, &
         airSpecificHumidity

    real(c_double), intent(in) :: &
         secondsToday

    integer(c_int), intent(in) :: &
         dayOfYear

    real(kind=c_double) :: &
         longitude, &
         solarTime, &
         hourAngle, &
         declination, &
         cosZ, &
         e, &
         d, &
         sw0, &
         degreesToRadians

    real(kind=c_double), parameter :: &
         secondsPerDay = 86400.0_c_double

    degreesToRadians = pi / 180.0_c_double

    ! longitude needs to be [-pi,pi] not [0,2pi]
    longitude = longitudeIn
    if (longitude > pi) longitude = longitude - 2.0_icepack_dbl_kind * pi

    solarTime = mod(real(secondsToday,kind=icepack_dbl_kind),secondsPerDay)/3600.0_icepack_dbl_kind + 12.0_icepack_dbl_kind*sin(0.5_icepack_dbl_kind*longitude)

    hourAngle = (12.0_icepack_dbl_kind - solarTime)*pi/12.0_icepack_dbl_kind

    ! solar declinatiom
    declination = 23.44_icepack_dbl_kind*cos((172.0_icepack_dbl_kind-real(dayOfYear,icepack_dbl_kind)) * 2.0_icepack_dbl_kind*pi/365.0_icepack_dbl_kind)*degreesToRadians

    ! solar zenith angle
    cosZ = sin(latitude)*sin(declination) + cos(latitude)*cos(declination)*cos(hourAngle)
    cosZ = max(cosZ,0.0_icepack_dbl_kind)

    e = 1.0e5_icepack_dbl_kind*airSpecificHumidity/(0.622_icepack_dbl_kind + 0.378_icepack_dbl_kind*airSpecificHumidity)

    d = (cosZ + 2.7_icepack_dbl_kind)*e*1.0e-5_icepack_dbl_kind+1.085_icepack_dbl_kind*cosZ+0.1_icepack_dbl_kind

    sw0 = 1353.0_icepack_dbl_kind*cosZ**2/d
    sw0 = max(sw0,0.0_icepack_dbl_kind)

    shortwaveDown = sw0 * (1.0_icepack_dbl_kind - 0.6_icepack_dbl_kind * cloudFraction**3)

  end subroutine column_shortwave_from_cloud_fraction

  !-----------------------------------------------------------------------------

  subroutine column_split_precipitation(&
       airTemperature, &
       precipitationRate, &
       snowfallRate, &
       rainfallRate) bind(C)

    use iso_c_binding, only: c_double

    use icepack_parameters, only: &
         freshWaterFreezingPoint => Tffresh

    real(c_double) :: &
         airTemperature, &
         precipitationRate, &
         snowfallRate, &
         rainfallRate

    if (airTemperature < freshWaterFreezingPoint) then

       snowfallRate = precipitationRate
       rainfallRate = 0.0_c_double

    else

       snowfallRate = 0.0_c_double
       rainfallRate = precipitationRate

    endif

  end subroutine column_split_precipitation

  !-----------------------------------------------------------------------------

  subroutine column_postprocess_ocean_forcing(&
       seaSurfaceSalinity, &
       oceanMixedLayerDepth, &
       seaFreezingTemperature) bind(C)

    use iso_c_binding, only: c_double

    use icepack_intfc, only: &
         icepack_sea_freezing_temperature

    real(kind=c_double) :: &
         seaSurfaceSalinity, &
         oceanMixedLayerDepth, &
         seaFreezingTemperature

    seaSurfaceSalinity   = max(seaSurfaceSalinity, 0.0_c_double)
    oceanMixedLayerDepth = max(oceanMixedLayerDepth, 0.0_c_double)

    ! sea freezing temperature
    seaFreezingTemperature = icepack_sea_freezing_temperature(seaSurfaceSalinity)

  end subroutine column_postprocess_ocean_forcing

  !-----------------------------------------------------------------------------
  ! warnings
  !-----------------------------------------------------------------------------

  subroutine column_warnings_reset() bind(C)

    use icepack_warnings, only: &
         icepack_warnings_resets

    call icepack_warnings_resets()

  end subroutine column_warnings_reset

  !-----------------------------------------------------------------------------

  function column_warnings_number() result(num) bind(C)

    use iso_c_binding, only: c_int

    use icepack_warnings, only: &
         icepack_warnings_number

    integer(c_int) :: num

    num = icepack_warnings_number()

  end function column_warnings_number

  !-----------------------------------------------------------------------------

  subroutine column_warnings_getone(iWarning, warningOut) bind(C)

    use iso_c_binding, only: c_int, c_char, c_null_char

    use icepack_intfc, only: &
         icepack_char_len_long

    use icepack_warnings, only: &
         icepack_warnings_getone

    integer(c_int), intent(in) :: iWarning
    character(len=1, kind=c_char), dimension(icepack_char_len_long), intent(out) :: warningOut

    character(len=icepack_char_len_long) :: warning

    integer :: &
         i

    warning = icepack_warnings_getone(iWarning+1) ! add one for c to fortran indices

    do i = 1, len_trim(warning)
       warningOut(i) = warning(i:i)
    enddo ! i
    warningOut(len_trim(warning)+1) = c_null_char

  end subroutine column_warnings_getone

  !-----------------------------------------------------------------------------

  function column_warnings_aborted() result(aborted) bind(C)

    use iso_c_binding, only: c_int

    use icepack_warnings, only: &
         icepack_warnings_aborted

    integer(c_int) :: aborted

    if (icepack_warnings_aborted()) then
       aborted = 1
    else
       aborted = 0
    endif

  end function column_warnings_aborted

  !-----------------------------------------------------------------------------
  ! type manipulations
  !-----------------------------------------------------------------------------

  subroutine convert_cstring_to_fstring(cstring, fstring)

    use iso_c_binding, only: &
         c_char, &
         c_null_char

    use icepack_intfc, only: &
         icepack_char_len_long

    character(kind=c_char, len=1), dimension(icepack_char_len_long), intent(in) :: &
         cstring

    character(len=icepack_char_len_long), intent(out) :: &
         fstring

    integer :: &
         i

    fstring = ""
    do i = 1, icepack_char_len_long
       if (cstring(i) == c_null_char) then
          exit
       else
          fstring(i:i) = cstring(i)
       endif
    enddo ! i

  end subroutine convert_cstring_to_fstring

  !-----------------------------------------------------------------------------

  subroutine convert_fstring_to_cstring(fstring, cstring)

    use iso_c_binding, only: &
         c_char, &
         c_null_char

    use icepack_intfc, only: &
         icepack_char_len_long

    character(len=icepack_char_len_long), intent(in) :: &
         fstring

    character(kind=c_char, len=1), dimension(icepack_char_len_long), intent(out) :: &
         cstring

    integer :: &
         i

    do i = 1, len_trim(fstring)
       cstring(i) = fstring(i:i)
    enddo ! i
    cstring(len_trim(fstring)+1) = c_null_char

  end subroutine convert_fstring_to_cstring

  !-----------------------------------------------------------------------------

  function convert_clogical_to_flogical(clogical) result(flogical)

    use iso_c_binding, only: c_int

    integer(c_int) :: &
         clogical

    logical :: &
         flogical

    flogical = .not. (clogical == 0)

  end function convert_clogical_to_flogical

  !-----------------------------------------------------------------------------

  function convert_flogical_to_clogical(flogical) result(clogical)

    use iso_c_binding, only: c_int

    logical :: &
         flogical

    integer(c_int) :: &
         clogical

    if (flogical) then
       clogical = 1
    else
       clogical = 0
    end if

  end function convert_flogical_to_clogical

  !-----------------------------------------------------------------------------

end module column_icepack
