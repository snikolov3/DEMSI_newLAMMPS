module column_ridging_dev

  implicit none

  private
  save

  public :: &
       ridging_set_variables, &
       ridging_ridge

contains

  !-----------------------------------------------------------------------------

  subroutine ridging_set_variables(&
       nCategories, &
       openWaterArea, &
       iceAreaCategory, &
       iceVolumeCategory, &
       ridgingIceThickness, &
       ridgingIceThicknessWeight, &
       netToGrossClosingRatio) bind(C)

    use iso_c_binding, only: c_int, c_double

    real(kind=c_double), parameter :: &
         participationFolding = 0.05_c_double, &
         puny = 1.0e-11_c_double

    integer(kind=c_int) :: &
         nCategories

    real(kind=c_double) :: &
         openWaterArea, &
         ridgingIceThickness, &
         ridgingIceThicknessWeight, &
         netToGrossClosingRatio

    real(kind=c_double), dimension(nCategories) :: &
         iceAreaCategory, &
         iceVolumeCategory

    real(kind=c_double), dimension(0:nCategories) :: &
         cumulativeITD

    real(kind=c_double), dimension(0:nCategories) :: &
         participationFunction

    real(kind=c_double) :: &
         iceThickness, &
         ridgeThicknessConst, &
         maxRaftThickness, &
         minRidgeThickness, &
         maxRidgeThickness, &
         meanRidgeThickness

    real(kind=c_double), dimension(nCategories) :: &
         thickeningRatio

    integer :: &
         iCategory

    ! calculate the cumulative ITD
    cumulativeITD = 0.0_c_double

    ! first entry is open water
    cumulativeITD(0) = openWaterArea

    do iCategory = 1, nCategories
       if (iceAreaCategory(iCategory) > puny) then
          cumulativeITD(iCategory) = cumulativeITD(iCategory-1) + iceAreaCategory(iCategory)
       else
          cumulativeITD(iCategory) = cumulativeITD(iCategory-1)
       endif
    enddo ! iCategory

    ! calculate the participation function
    ! b(h) = exp(-G(h)/astar)
    ! apartic(n) = [exp(-G(n-1)/astar - exp(-G(n)/astar] / [1-exp(-1/astar)].
    ! The expression for apartic is found by integrating b(h)g(h)
    ! between the category boundaries.
    participationFunction(0) = &
         (1.0_c_double - exp(-cumulativeITD(0) / participationFolding)) / &
         (1.0_c_double - exp(-1.0_c_double / participationFolding))

    do iCategory = 1, nCategories

       participationFunction(iCategory) = &
            (exp(-cumulativeITD(iCategory-1) / participationFolding) - exp(-cumulativeITD(iCategory) / participationFolding)) / &
            (1.0_c_double - exp(-1.0_c_double / participationFolding))

    enddo ! iCategory

    ! find weighted ice thickness
    ridgingIceThickness = 0.0_c_double
    ridgingIceThicknessWeight = 0.0_c_double

    do iCategory = 1, nCategories

       iceThickness = 0.0_c_double
       if (iceAreaCategory(iCategory) > puny) then
          iceThickness = iceVolumeCategory(iCategory) / iceAreaCategory(iCategory)
       endif

       ridgingIceThickness       = ridgingIceThickness       + participationFunction(iCategory) * iceThickness
       ridgingIceThicknessWeight = ridgingIceThicknessWeight + participationFunction(iCategory)

    enddo ! iCategory

    ! thickening ratio
    ridgeThicknessConst = 25.0_c_double
    maxRaftThickness = 1.0_c_double

    do iCategory = 1, nCategories

       if (iceAreaCategory(iCategory) > puny) then

          iceThickness = iceVolumeCategory(iCategory) / iceAreaCategory(iCategory)

          minRidgeThickness = min(2.0_c_double*iceThickness, iceThickness + maxRaftThickness)
          maxRidgeThickness = 2.0_c_double*sqrt(ridgeThicknessConst*iceThickness)
          maxRidgeThickness = max(minRidgeThickness,puny)

          meanRidgeThickness = 0.5 * (minRidgeThickness + maxRidgeThickness)

          thickeningRatio(iCategory) = meanRidgeThickness / iceThickness

       else

          thickeningRatio(iCategory) = 1.0_c_double ! avoid div by zero

       endif

    enddo ! iCategory

    ! netToGrossClosingRatio
    ! weighting in contact for this element
    netToGrossClosingRatio = participationFunction(0)
    do iCategory = 1, nCategories
       netToGrossClosingRatio = netToGrossClosingRatio + &
            participationFunction(iCategory) * (1.0_c_double - 1.0_c_double / thickeningRatio(iCategory))
    enddo ! iCategory

  end subroutine ridging_set_variables

  !-----------------------------------------------------------------------------

  subroutine ridging_ridge(&
       nCategories, &
       openWaterArea, &
       iceAreaCategory, &
       iceVolumeCategory, &
       categoryThicknessLimits, &
       netToGrossClosingRatio, &
       effectiveElementArea, &
       changeEffectiveElementArea) bind(C)

    use iso_c_binding, only: c_int, c_double

    integer(kind=c_int) :: &
         nCategories

    real(kind=c_double) :: &
         openWaterArea, &
         netToGrossClosingRatio, &
         effectiveElementArea, &
         changeEffectiveElementArea

    real(kind=c_double), dimension(0:nCategories) :: &
         categoryThicknessLimits

    real(kind=c_double), dimension(nCategories) :: &
         iceAreaCategory, &
         iceVolumeCategory

    real(kind=c_double), dimension(0:nCategories) :: &
         cumulativeITD

    real(kind=c_double), dimension(0:nCategories) :: &
         participationFunction

    real(kind=c_double), dimension(nCategories) :: &
         thickeningRatio

    real(kind=c_double) :: &
         elementAreaRatio, &
         netClosing, &
         opening

    real(kind=c_double), dimension(nCategories) :: &
         iceAreaCategoryInit, &
         iceVolumeCategoryInit

    real(kind=c_double) :: &
         iceThickness, &
         ridgeThicknessConst, &
         maxRaftThickness, &
         minRidgeThickness, &
         maxRidgeThickness, &
         meanRidgeThickness, &
         grossClosing, &
         iceAreaSrc, &
         iceVolumeSrc, &
         iceAreaDst, &
         iceVolumeDst, &
         ridgeDistWidth, &
         hLeft, &
         hRight, &
         expLeft, &
         expRight, &
         areaFractionDst, &
         volumeFractionDst, &
         advectiveDivergence, &
         areaSum, &
         iceAreaReduction, &
         openWaterReduction, &
         rescaleFactor, &
         weight1, &
         weight2, &
         areaDecrease

    integer(kind=c_int) :: &
         iCategory, &
         iCategorySrc, &
         iCategoryDst, &
         iteration

    integer(kind=c_int), parameter :: &
         nIterations = 100

    real(kind=c_double), parameter :: &
         participationFolding = 0.05_c_double, &
         puny = 1.0e-11_c_double, &
         ridgingFolding = 3.0_c_double

    ! calculate the cumulative ITD
    cumulativeITD = 0.0

    ! first entry is open water
    cumulativeITD(0) = openWaterArea

    do iCategory = 1, nCategories
       if (iceAreaCategory(iCategory) > puny) then
          cumulativeITD(iCategory) = cumulativeITD(iCategory-1) + iceAreaCategory(iCategory)
       else
          cumulativeITD(iCategory) = cumulativeITD(iCategory-1)
       endif
    enddo ! iCategory

    ! calculate the participation function
    ! b(h) = exp(-G(h)/astar)
    ! apartic(n) = [exp(-G(n-1)/astar - exp(-G(n)/astar] / [1-exp(-1/astar)].
    ! The expression for apartic is found by integrating b(h)g(h)
    ! between the category boundaries.
    participationFunction(0) = &
         (1.0_c_double - exp(-cumulativeITD(0) / participationFolding)) / &
         (1.0_c_double - exp(-1.0_c_double / participationFolding))

    do iCategory = 1, nCategories

       participationFunction(iCategory) = &
            (exp(-cumulativeITD(iCategory-1) / participationFolding) - exp(-cumulativeITD(iCategory) / participationFolding)) / &
            (1.0_c_double - exp(-1.0_c_double / participationFolding))

    enddo ! iCategory

    ! thickening ratio
    ridgeThicknessConst = 25.0_c_double
    maxRaftThickness = 1.0_c_double

    do iCategory = 1, nCategories

       if (iceAreaCategory(iCategory) > puny) then

          iceThickness = iceVolumeCategory(iCategory) / &
                         iceAreaCategory  (iCategory)

          minRidgeThickness = min(2.0_c_double*iceThickness, iceThickness + maxRaftThickness)
          maxRidgeThickness = 2.0_c_double*sqrt(ridgeThicknessConst*iceThickness)
          maxRidgeThickness = max(minRidgeThickness,puny)

          meanRidgeThickness = 0.5 * (minRidgeThickness + maxRidgeThickness)

          thickeningRatio(iCategory) = meanRidgeThickness / iceThickness

       else

          thickeningRatio(iCategory) = 1.0_c_double ! avoid div by zero

       endif

    enddo ! iCategory

    ! ratio of effective element areas before and after ridging
    elementAreaRatio = (effectiveElementArea + changeEffectiveElementArea) / effectiveElementArea

    ! change the effective element area due to ridging
    effectiveElementArea = effectiveElementArea + changeEffectiveElementArea

    ! modify area and thickness in ITD from element size change
    do iCategory = 1, nCategories

       iceAreaCategory  (iCategory) = iceAreaCategory  (iCategory) / elementAreaRatio
       iceVolumeCategory(iCategory) = iceVolumeCategory(iCategory) / elementAreaRatio

    enddo ! iCategory

    ! element net closing
    netClosing = (1.0_c_double / elementAreaRatio - 1.0_c_double)

    !-----------------------------------------------------------------
    ! Compute divu_adv, the divergence rate given by the transport/
    ! advection scheme, which may not be equal to divu as computed
    ! from the velocity field.
    !
    ! If divu_adv < 0, make sure the closing rate is large enough
    ! to give asum = 1.0 after ridging.
    !-----------------------------------------------------------------

    areaSum = openWaterArea
    do iCategory = 1, nCategories
       areaSum = areaSum + iceAreaCategory(iCategory)
    enddo ! iCategory

    advectiveDivergence = (1.0_c_double - areaSum)! / dt

    if (advectiveDivergence < 0.0) then
       netClosing = max(netClosing, -advectiveDivergence)
    endif

    !-----------------------------------------------------------------
    ! Compute the (non-negative) opening rate that will give
    ! asum = 1.0 after ridging.
    !-----------------------------------------------------------------

    opening = netClosing + advectiveDivergence

    ! iterate since might completely evacuate a category
    ridgingIteration: do iteration = 1, nIterations

       ! save initial values of iceAreaCategory, iceVolumeCategory
       do iCategory = 1, nCategories
          iceAreaCategoryInit  (iCategory) = iceAreaCategory  (iCategory)
          iceVolumeCategoryInit(iCategory) = iceVolumeCategory(iCategory)
       enddo ! iCategory

       ! gross closing from net closing
       grossClosing = netClosing / netToGrossClosingRatio

       !-----------------------------------------------------------------
       ! Reduce the closing rate if more than 100% of the open water
       ! would be removed.  Reduce the opening rate proportionately.
       !-----------------------------------------------------------------

       if (participationFunction(0) > 0.0_c_double) then

          openWaterReduction = participationFunction(0) * grossClosing! * dt

          if (openWaterReduction > openWaterArea) then

             rescaleFactor = openWaterArea / openWaterReduction

             grossClosing = grossClosing * rescaleFactor
             opening      = opening      * rescaleFactor

          endif
       endif

       !-----------------------------------------------------------------
       ! Reduce the closing rate if more than 100% of any ice category
       ! would be removed.  Reduce the opening rate proportionately.
       !-----------------------------------------------------------------

       do iCategory = 1, nCategories

          if (iceAreaCategory(iCategory) > puny .and. &
              participationFunction(iCategory) > 0.0_c_double) then

             iceAreaReduction = participationFunction(iCategory) * grossClosing + opening! * dt

             if (iceAreaReduction > iceAreaCategory(iCategory)) then

                rescaleFactor = iceAreaCategory(iCategory) / iceAreaReduction

                grossClosing = grossClosing * rescaleFactor
                opening      = opening      * rescaleFactor

             endif
          endif

       enddo ! iCategory

       ! change open water fraction
       openWaterArea = openWaterArea - participationFunction(0) * grossClosing! * timeStep

       ! loop over source categories
       do iCategorySrc = 1, nCategories

          ! only move ice from source if ice present, participation and closing is positive
          if (iceAreaCategoryInit(iCategorySrc) > puny .and. &
              participationFunction(iCategorySrc) > 0.0_c_double .and. &
              grossClosing > 0.0_c_double) then

             ! area and volume of source ice
             iceAreaSrc = participationFunction(iCategorySrc) * grossClosing! * timeStep
             iceVolumeSrc = iceVolumeCategoryInit(iCategorySrc) * (iceAreaSrc / iceAreaCategoryInit(iCategorySrc))

             ! Remove area and volume from source category
             iceAreaCategory  (iCategorySrc) = iceAreaCategory  (iCategorySrc) - iceAreaSrc
             iceVolumeCategory(iCategorySrc) = iceVolumeCategory(iCategorySrc) - iceVolumeSrc

             ! Area and volume of destination ridge
             iceAreaDst = iceAreaSrc / thickeningRatio(iCategorySrc)
             iceVolumeDst = iceVolumeSrc ! (=volume of source ridge from conservation)

             ! source ice thickness
             iceThickness = iceVolumeCategoryInit(iCategorySrc) / &
                            iceAreaCategoryInit  (iCategorySrc)

             ! minimum thickness of new ridge for this source thickness
             minRidgeThickness = min(2.0_c_double*iceThickness, iceThickness + maxRaftThickness)

             ! width of the thickness distribution for new ridges from this source thickness
             ridgeDistWidth = ridgingFolding * sqrt(iceThickness)

             ! loop over destination categories
             do iCategoryDst = 1, nCategories

                ! category integration limits
                hLeft  = max(categoryThicknessLimits(iCategoryDst-1),minRidgeThickness)
                hRight =     categoryThicknessLimits(iCategoryDst)

                if (iCategoryDst < nCategories) then
                   ! not last category

                   if (minRidgeThickness < hRight) then

                      expLeft  = exp(-(hLeft  - minRidgeThickness) / ridgeDistWidth)
                      expRight = exp(-(hRight - minRidgeThickness) / ridgeDistWidth)

                   else

                      expLeft  = 0.0_c_double
                      expRight = 0.0_c_double

                   endif

                else
                   ! last category

                   expLeft  = exp(-(hLeft - minRidgeThickness) / ridgeDistWidth)
                   expRight = 0.0_c_double

                endif

                areaFractionDst = expLeft - expRight
                volumeFractionDst = ((hLeft + ridgeDistWidth) * expLeft - (hRight + ridgeDistWidth) * expRight) / &
                                    (minRidgeThickness + ridgeDistWidth)

                ! Add area and volume to destination category
                iceAreaCategory  (iCategoryDst) = iceAreaCategory  (iCategoryDst) + areaFractionDst   * iceAreaDst
                iceVolumeCategory(iCategoryDst) = iceVolumeCategory(iCategoryDst) + volumeFractionDst * iceVolumeSrc

             enddo ! iCategoryDst

          endif

       enddo ! iCategorySrc

       ! total area after ridging (should == 1)
       areaSum = openWaterArea
       do iCategory = 1, nCategories
          areaSum = areaSum + iceAreaCategory(iCategory)
       enddo ! iCategory

       ! check if iterations complete
       if (abs(areaSum - 1.0_c_double) < puny) then
          exit ridgingIteration
       else

          ! ridging not complete
          advectiveDivergence = (1.0_c_double - areaSum)! / dt

          netClosing = max(0.0_c_double, -advectiveDivergence)
          opening    = max(0.0_c_double,  advectiveDivergence)

       endif

    enddo ridgingIteration

  end subroutine ridging_ridge

  !-----------------------------------------------------------------------------

end module column_ridging_dev
