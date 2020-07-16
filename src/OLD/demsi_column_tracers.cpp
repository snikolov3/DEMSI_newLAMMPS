#include "demsi_column_tracers.h"
#include "demsi_configs.h"
#include "demsi_column_variables.h"
#include <cmath>

// column tracer routines
extern "C" {void column_init_tracer_flags(char*, int*);}
extern "C" {void column_init_tracer_numbers(int*);}
extern "C" {void column_init_tracer_indices(int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*);}
extern "C" {void column_print_tracer_object(int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int[], double[], int[], int[]);}

namespace DEMSI {

// convert 2D indices to 1D index for fortran array ordering
int twoDimIndices(const int i, const int ni, const int j) {

  return j * ni + i;

}

// constructor for ColumnTracers object
ColumnTracers::ColumnTracers(DEMSI::ColumnDimensions* columnDimensionsIn, DEMSI::Configs* configsIn) {

  // set pointers to other objects
  configs = configsIn;
  columnDimensions = columnDimensionsIn;

  // set which tracers to use
  get_use_tracers();

  // get the number of tracers
  get_tracer_number();

  // allocate the tracer info arrays
  allocate_arrays();

  // init the child tracer indices
  init_child_indices();

  // init the parent tracer indices
  init_parent_indices();

  // set the first ancestor mask
  set_first_ancestor_mask();

  // set the ancester indices
  set_ancester_indices();

  // set the use tracer flags in icepack
  set_icepack_tracer_flags();

  // set the number of tracers in icepack
  set_icepack_tracer_numbers();

  // set the tracer indices in icepack
  set_icepack_tracer_indices();

}

// set which tracers to use
void ColumnTracers::get_use_tracers(void) {

  configs->get({"ConfigGroup:columnTracers","Config:useIceAge"},         useIceAge);
  configs->get({"ConfigGroup:columnTracers","Config:useFirstYearIce"},   useFirstYearIce);
  configs->get({"ConfigGroup:columnTracers","Config:useLevelIce"},       useLevelIce);
  configs->get({"ConfigGroup:columnTracers","Config:useCesmMeltponds"},  useCesmMeltponds);
  configs->get({"ConfigGroup:columnTracers","Config:useLevelMeltponds"}, useLevelMeltponds);
  configs->get({"ConfigGroup:columnTracers","Config:useTopoMeltponds"},  useTopoMeltponds);
  configs->get({"ConfigGroup:columnTracers","Config:useAerosols"},       useAerosols);

  useMeltPonds = useCesmMeltponds or useLevelMeltponds or useTopoMeltponds;

}

// set the total number of tracers
void ColumnTracers::get_tracer_number(void) {

  // layer numbers
  int nIceLayers  = columnDimensions->size("nIceLayers");
  int nSnowLayers = columnDimensions->size("nSnowLayers");
  int nAerosols   = columnDimensions->size("nAerosols");

  // surfaceTemperature
  nTracers = 1;

  // iceEnthalpy
  nTracers = nTracers + nIceLayers;

  // snowEnthalpy
  nTracers = nTracers + nSnowLayers;

  // ice Salinity
  nTracers = nTracers + nIceLayers;

  // iceAge
  if (useIceAge) {
    nTracers = nTracers + 1;
  }

  // firstYearIceArea
  if (useFirstYearIce) {
    nTracers = nTracers + 1;
  }

  // level ice tracers
  if (useLevelIce) {
    nTracers = nTracers + 2;
  }

  // pond tracers
  if (useMeltPonds) {
    nTracers = nTracers + 2;
  }

  // level or topo ponds
  if (useLevelMeltponds or
      useTopoMeltponds) {
    nTracers = nTracers + 1;
  }

  // aerosols
  if (useAerosols) {
    nTracers = nTracers + nAerosols*4;
  }

}

// allocate the tracer info arrays
void ColumnTracers::allocate_arrays(void) {

  int nCategories = columnDimensions->size("nCategories");

  // allocate the category tracer array
  //allocate(tracerObject % tracerArrayCategory(tracerObject % nTracers,nCategories))
  tracerArrayCategory = new double[nTracers*nCategories];

  // allocate the cell tracer array
  //allocate(tracerObject % tracerArray(tracerObject % nTracers))
  tracerArray = new double[nTracers];

  // allocate other arrays
  //allocate(tracerObject % parentIndex(tracerObject % nTracers))
  //allocate(tracerObject % firstAncestorMask(tracerObject % nTracers, tracerObject % nBaseTracers))
  //allocate(tracerObject % ancestorIndices(tracerObject % nTracers, tracerObject % nMaxAncestorTracers))
  //allocate(tracerObject % ancestorNumber(tracerObject % nTracers))
  parentIndex = new int[nTracers];
  firstAncestorMask = new double[nTracers*nBaseTracers];
  ancestorIndices = new int[nTracers*nMaxAncestorTracers];
  ancestorNumber = new int[nTracers];

}

void ColumnTracers::init_child_indices(void) {

  // layer numbers
  int nIceLayers  = columnDimensions->size("nIceLayers");
  int nSnowLayers = columnDimensions->size("nSnowLayers");
  int nAerosols   = columnDimensions->size("nAerosols");

  int nTracers;
  int indexMissingValue = 0;

  // ice/snow surface temperature
  indexSurfaceTemperature = 1;
  nTracers = 1;

  // ice enthalpy
  indexIceEnthalpy = nTracers + 1;
  nTracers = nTracers + nIceLayers;

  // snow enthalpy
  indexSnowEnthalpy = nTracers + 1;
  nTracers = nTracers + nSnowLayers;

  // ice salinity
  indexIceSalinity = nTracers + 1;
  nTracers = nTracers + nIceLayers;

  // ice age
  indexIceAge = indexMissingValue;
  if (useIceAge) {
    nTracers = nTracers + 1;
    indexIceAge = nTracers;
  }

  // first year ice
  indexFirstYearIceArea = indexMissingValue;
  if (useFirstYearIce) {
    nTracers = nTracers + 1;
    indexFirstYearIceArea = nTracers;
  }

  // level ice
  indexLevelIceArea   = indexMissingValue;
  indexLevelIceVolume = indexMissingValue;
  if (useLevelIce) {
    nTracers = nTracers + 1;
    indexLevelIceArea = nTracers;
    nTracers = nTracers + 1;
    indexLevelIceVolume = nTracers;
  }

  // ponds
  indexPondArea         = indexMissingValue;
  indexPondDepth        = indexMissingValue;
  indexPondLidThickness = indexMissingValue;

  if (useMeltPonds) {
    nTracers = nTracers + 1;
    indexPondArea = nTracers;
    nTracers = nTracers + 1;
    indexPondDepth = nTracers;
  }
  if (useLevelMeltponds) {
    nTracers = nTracers + 1;
    indexPondLidThickness = nTracers;
  }
  if (useTopoMeltponds) {
    nTracers = nTracers + 1;
    indexPondLidThickness = nTracers;
  }

  // aerosols
  indexAerosols = indexMissingValue;
  if (useAerosols) {
    indexAerosols = nTracers + 1;
  }

}

// init the parent tracer indices
void ColumnTracers::init_parent_indices(void) {

  // layer numbers
  int nIceLayers  = columnDimensions->size("nIceLayers");
  int nSnowLayers = columnDimensions->size("nSnowLayers");
  int nAerosols   = columnDimensions->size("nAerosols");

  // ice/snow surface temperature
  parentIndex[indexSurfaceTemperature - 1] = 0;

  // ice enthalpy and salinity
  for (int iIceLayer = 0 ; iIceLayer < nIceLayers ; iIceLayer++) {
    parentIndex[indexIceEnthalpy - 1 + iIceLayer] = 1;
    parentIndex[indexIceSalinity - 1 + iIceLayer] = 1;
  } // iIceLayer

  // snow enthalpy
  for (int iSnowLayer = 0 ; iSnowLayer < nSnowLayers ; iSnowLayer++) {
    parentIndex[indexSnowEnthalpy - 1] = 2;
  } // iSnowLayer

  // ice age
  if (useIceAge) {
    parentIndex[indexIceAge - 1] = 1;
  }

  // first year ice
  if (useFirstYearIce) {
    parentIndex[indexFirstYearIceArea - 1] = 0;
  }

  // level ice area
  if (useLevelIce) {
    parentIndex[indexLevelIceArea   - 1] = 0;
    parentIndex[indexLevelIceVolume - 1] = 1;
  }

  // cesm melt ponds
  if (useCesmMeltponds) {
    parentIndex[indexPondArea  - 1] = 0;
    parentIndex[indexPondDepth - 1] = 2 + indexPondArea;
  }

  // level ice ponds
  if (useLevelMeltponds) {
    parentIndex[indexPondArea         - 1] = 2 + indexLevelIceArea;
    parentIndex[indexPondDepth        - 1] = 2 + indexPondArea;
    parentIndex[indexPondLidThickness - 1] = 2 + indexPondArea;
  }

  // topo melt ponds
  if (useTopoMeltponds) {
    parentIndex[indexPondArea - 1]         = 0;
    parentIndex[indexPondDepth - 1]        = 2 + indexPondArea;
    parentIndex[indexPondLidThickness - 1] = 2 + indexPondArea;
  }

  // aerosols
  if (useAerosols) {
    for (int iAerosol = 0 ; iAerosol < nAerosols ; iAerosol++) {
      parentIndex[indexAerosols - 1 + iAerosol*4    ] = 2; // snow
      parentIndex[indexAerosols - 1 + iAerosol*4 + 1] = 2; // snow
      parentIndex[indexAerosols - 1 + iAerosol*4 + 2] = 1; // ice
      parentIndex[indexAerosols - 1 + iAerosol*4 + 3] = 1; // ice
    } // iAerosol
  }

}

// set the first ancestor mask
void ColumnTracers::set_first_ancestor_mask(void) {

  // mask for base quantity on which tracers are carried
  for (int i = 0 ; i < nTracers*nBaseTracers ; i++) {
    firstAncestorMask[i] = 0.0;
  }

  for (int iTracer = 0 ; iTracer < nTracers ; iTracer++) {

    if (parentIndex[iTracer] == 0) {

      // ice area
      firstAncestorMask[twoDimIndices(iTracer,nTracers,0)] = 1.0;

    } else if (parentIndex[iTracer] == 1) { // ice volume

      // ice volume
      firstAncestorMask[twoDimIndices(iTracer,nTracers,1)] = 1.0;

    } else if (parentIndex[iTracer] == 2) { // snow volume

      // snow volume
      firstAncestorMask[twoDimIndices(iTracer,nTracers,2)] = 1.0;

    } else {

      // default: ice area
      firstAncestorMask[twoDimIndices(iTracer,nTracers,0)] = 1.0;

    }

  } // iTracer

}

// set the ancester indices
void ColumnTracers::set_ancester_indices(void) {

  // initialize
  for (int i = 0 ; i < nTracers ; i++) {
    ancestorNumber[i] = 0;
  }
  for (int i = 0 ; i < nTracers*nMaxAncestorTracers ; i++) {
    ancestorIndices[i] = 0;
  }

  // cesm meltponds
  if (useCesmMeltponds) {

    // melt pond depth
    ancestorNumber [indexPondDepth - 1] = 1;
    ancestorIndices[twoDimIndices(indexPondDepth - 1, nTracers, 0)] = indexPondArea; // on melt pond area

  }

  // level melt ponds
  if (useLevelMeltponds) {

    // melt pond area
    ancestorNumber [indexPondArea - 1] = 1;
    ancestorIndices[twoDimIndices(indexPondArea - 1, nTracers, 0)] = indexLevelIceArea; // on level ice area

    // melt pond depth
    ancestorNumber [indexPondDepth - 1] = 2;
    ancestorIndices[twoDimIndices(indexPondDepth - 1, nTracers, 1)] = indexPondArea; // on melt pond area
    ancestorIndices[twoDimIndices(indexPondDepth - 1, nTracers, 0)] = indexLevelIceArea; // on level ice area

    // refrozen pond lid
    ancestorNumber [indexPondLidThickness - 1] = 2;
    ancestorIndices[twoDimIndices(indexPondLidThickness - 1, nTracers, 1)] = indexPondArea; // on melt pond area
    ancestorIndices[twoDimIndices(indexPondLidThickness - 1, nTracers, 0)] = indexLevelIceArea; // on level ice area

  }

  // topographic melt ponds
  if (useTopoMeltponds) {

    // melt pond depth
    ancestorNumber [indexPondDepth - 1] = 1;
    ancestorIndices[twoDimIndices(indexPondDepth - 1, nTracers, 0)] = indexPondArea; // on melt pond area

    // refrozen pond lid
    ancestorNumber [indexPondLidThickness - 1] = 1;
    ancestorIndices[twoDimIndices(indexPondLidThickness - 1, nTracers, 0)] = indexPondArea; // on melt pond area

  }

}

// print the tracer object
void ColumnTracers::print(void) {

  column_print_tracer_object(
	&nTracers,
	&nBaseTracers,
	&nMaxAncestorTracers,
	&indexSurfaceTemperature,
	&indexIceEnthalpy,
	&indexSnowEnthalpy,
	&indexIceSalinity,
	&indexIceAge,
	&indexFirstYearIceArea,
	&indexLevelIceArea,
	&indexLevelIceVolume,
	&indexPondArea,
	&indexPondDepth,
	&indexPondLidThickness,
	&indexAerosols,
	parentIndex,
	firstAncestorMask,
	ancestorIndices,
	ancestorNumber);

}

// set the use tracer flags in icepack
void ColumnTracers::set_icepack_tracer_flags(void) {

  int useIceAgeInt = (int) useIceAge;
  column_init_tracer_flags((char*) "useIceAge",         &useIceAgeInt);
  int useFirstYearIceInt = (int) useFirstYearIce;
  column_init_tracer_flags((char*) "useFirstYearIce",   &useFirstYearIceInt);
  int useLevelIceInt = (int) useLevelIce;
  column_init_tracer_flags((char*) "useLevelIce",       &useLevelIceInt);
  int useMeltPondsInt = (int) useMeltPonds;
  column_init_tracer_flags((char*) "useMeltPonds",      &useMeltPondsInt);
  int useCesmMeltpondsInt = (int) useCesmMeltponds;
  column_init_tracer_flags((char*) "useCesmMeltponds",  &useCesmMeltpondsInt);
  int useLevelMeltpondsInt = (int) useLevelMeltponds;
  column_init_tracer_flags((char*) "useLevelMeltponds", &useLevelMeltpondsInt);
  int useTopoMeltpondsInt = (int) useTopoMeltponds;
  column_init_tracer_flags((char*) "useTopoMeltponds",  &useTopoMeltpondsInt);
  int useAerosolsInt = (int) useAerosols;
  column_init_tracer_flags((char*) "useAerosols",       &useAerosolsInt);

}

// set the number of tracers in icepack
void ColumnTracers::set_icepack_tracer_numbers(void) {

  column_init_tracer_numbers(&nTracers);

}

// set the tracer indices in icepack
void ColumnTracers::set_icepack_tracer_indices(void) {

  column_init_tracer_indices(
       &indexSurfaceTemperature,
       &indexIceEnthalpy,
       &indexSnowEnthalpy,
       &indexIceSalinity,
       &indexIceAge,
       &indexFirstYearIceArea,
       &indexLevelIceArea,
       &indexLevelIceVolume,
       &indexPondArea,
       &indexPondDepth,
       &indexPondLidThickness,
       &indexAerosols);

}

// Set pointers to the column tracer variables
void ColumnTracers::set_tracer_variable_pointers(
	DEMSI::ColumnVariable<double>* iceAreaCategoryIn,
	DEMSI::ColumnVariable<double>* iceVolumeCategoryIn,
	DEMSI::ColumnVariable<double>* snowVolumeCategoryIn,
	DEMSI::ColumnVariable<double>* surfaceTemperatureIn,
	DEMSI::ColumnVariable<double>* iceEnthalpyIn,
	DEMSI::ColumnVariable<double>* iceSalinityIn,
	DEMSI::ColumnVariable<double>* snowEnthalpyIn,
	DEMSI::ColumnVariable<double>* iceAgeIn,
	DEMSI::ColumnVariable<double>* firstYearIceAreaIn,
	DEMSI::ColumnVariable<double>* levelIceAreaIn,
	DEMSI::ColumnVariable<double>* levelIceVolumeIn,
	DEMSI::ColumnVariable<double>* pondAreaIn,
	DEMSI::ColumnVariable<double>* pondDepthIn,
	DEMSI::ColumnVariable<double>* pondLidThicknessIn,
	DEMSI::ColumnVariable<double>* snowScatteringAerosolIn,
	DEMSI::ColumnVariable<double>* snowBodyAerosolIn,
	DEMSI::ColumnVariable<double>* iceScatteringAerosolIn,
	DEMSI::ColumnVariable<double>* iceBodyAerosolIn) {

  iceAreaCategory        = iceAreaCategoryIn;
  iceVolumeCategory      = iceVolumeCategoryIn;
  snowVolumeCategory     = snowVolumeCategoryIn;
  surfaceTemperature     = surfaceTemperatureIn;
  iceEnthalpy            = iceEnthalpyIn;
  iceSalinity            = iceSalinityIn;
  snowEnthalpy           = snowEnthalpyIn;
  iceAge                 = iceAgeIn;
  firstYearIceArea       = firstYearIceAreaIn;
  levelIceArea           = levelIceAreaIn;
  levelIceVolume         = levelIceVolumeIn;
  pondArea               = pondAreaIn;
  pondDepth              = pondDepthIn;
  pondLidThickness       = pondLidThicknessIn;
  snowScatteringAerosol  = snowScatteringAerosolIn;
  snowBodyAerosol        = snowBodyAerosolIn;
  iceScatteringAerosol   = iceScatteringAerosolIn;
  iceBodyAerosol         = iceBodyAerosolIn;

}

// Set pointers to the column aggregated tracer variables
void ColumnTracers::set_aggregated_tracer_variable_pointers(
	DEMSI::ColumnVariable<double>* iceAreaCellIn,
	DEMSI::ColumnVariable<double>* iceVolumeCellIn,
	DEMSI::ColumnVariable<double>* snowVolumeCellIn,
	DEMSI::ColumnVariable<double>* surfaceTemperatureCellIn,
	DEMSI::ColumnVariable<double>* iceEnthalpyCellIn,
	DEMSI::ColumnVariable<double>* iceSalinityCellIn,
	DEMSI::ColumnVariable<double>* snowEnthalpyCellIn,
	DEMSI::ColumnVariable<double>* iceAgeCellIn,
	DEMSI::ColumnVariable<double>* firstYearIceAreaCellIn,
	DEMSI::ColumnVariable<double>* levelIceAreaCellIn,
	DEMSI::ColumnVariable<double>* levelIceVolumeCellIn,
	DEMSI::ColumnVariable<double>* pondAreaCellIn,
	DEMSI::ColumnVariable<double>* pondDepthCellIn,
	DEMSI::ColumnVariable<double>* pondLidThicknessCellIn,
	DEMSI::ColumnVariable<double>* snowScatteringAerosolCellIn,
	DEMSI::ColumnVariable<double>* snowBodyAerosolCellIn,
	DEMSI::ColumnVariable<double>* iceScatteringAerosolCellIn,
	DEMSI::ColumnVariable<double>* iceBodyAerosolCellIn) {

  iceAreaCell               = iceAreaCellIn;
  iceVolumeCell             = iceVolumeCellIn;
  snowVolumeCell            = snowVolumeCellIn;
  surfaceTemperatureCell    = surfaceTemperatureCellIn;
  iceEnthalpyCell           = iceEnthalpyCellIn;
  iceSalinityCell           = iceSalinityCellIn;
  snowEnthalpyCell          = snowEnthalpyCellIn;
  iceAgeCell                = iceAgeCellIn;
  firstYearIceAreaCell      = firstYearIceAreaCellIn;
  levelIceAreaCell          = levelIceAreaCellIn;
  levelIceVolumeCell        = levelIceVolumeCellIn;
  pondAreaCell              = pondAreaCellIn;
  pondDepthCell             = pondDepthCellIn;
  pondLidThicknessCell      = pondLidThicknessCellIn;
  snowScatteringAerosolCell = snowScatteringAerosolCellIn;
  snowBodyAerosolCell       = snowBodyAerosolCellIn;
  iceScatteringAerosolCell  = iceScatteringAerosolCellIn;
  iceBodyAerosolCell        = iceBodyAerosolCellIn;

}

// Set the tracer category arrays for a particle
void ColumnTracers::set_icepack_tracer_array_category(const int iParticle) {

  // layer numbers
  int nCategories = columnDimensions->size("nCategories");
  int nIceLayers  = columnDimensions->size("nIceLayers");
  int nSnowLayers = columnDimensions->size("nSnowLayers");
  int nAerosols   = columnDimensions->size("nAerosols");

  int iTracer = 0;

  // surfaceTemperature
  for (int iCategory = 0 ; iCategory < nCategories ; iCategory++) {
    tracerArrayCategory[iTracer + iCategory*nTracers] = (*surfaceTemperature)(iParticle,iCategory);
  } // iCategory
  iTracer += 1;

  // iceEnthalpy
  for (int iCategory = 0 ; iCategory < nCategories ; iCategory++) {
    for (int iIceLayer = 0 ; iIceLayer < nIceLayers ; iIceLayer++) {
      tracerArrayCategory[iTracer + iIceLayer + iCategory*nTracers] = (*iceEnthalpy)(iParticle,iCategory,iIceLayer);
    } // iIceLayer
  } // iCategory
  iTracer += nIceLayers;

  // snowEnthalpy
  for (int iCategory = 0 ; iCategory < nCategories ; iCategory++) {
    for (int iSnowLayer = 0 ; iSnowLayer < nSnowLayers ; iSnowLayer++) {
      tracerArrayCategory[iTracer + iSnowLayer + iCategory*nTracers] = (*snowEnthalpy)(iParticle,iCategory,iSnowLayer);
    } // iSnowLayer
  } // iCategory
  iTracer += nSnowLayers;

  // ice Salinity
  for (int iCategory = 0 ; iCategory < nCategories ; iCategory++) {
    for (int iIceLayer = 0 ; iIceLayer < nIceLayers ; iIceLayer++) {
      tracerArrayCategory[iTracer + iIceLayer + iCategory*nTracers] = (*iceSalinity)(iParticle,iCategory,iIceLayer);
    } // iIceLayer
  } // iCategory
  iTracer += nIceLayers;

  // iceAge
  if (useIceAge) {
    for (int iCategory = 0 ; iCategory < nCategories ; iCategory++) {
      tracerArrayCategory[iTracer + iCategory*nTracers] = (*iceAge)(iParticle,iCategory);
    } // iCategory
    iTracer += 1;
  }

  // firstYearIceArea
  if (useFirstYearIce) {
    for (int iCategory = 0 ; iCategory < nCategories ; iCategory++) {
      tracerArrayCategory[iTracer + iCategory*nTracers] = (*firstYearIceArea)(iParticle,iCategory);
    } // iCategory
    iTracer += 1;
  }

  // level ice tracers
  if (useLevelIce) {
    for (int iCategory = 0 ; iCategory < nCategories ; iCategory++) {
      tracerArrayCategory[iTracer + iCategory*nTracers] = (*levelIceArea)(iParticle,iCategory);
    } // iCategory
    iTracer += 1;
    for (int iCategory = 0 ; iCategory < nCategories ; iCategory++) {
      tracerArrayCategory[iTracer + iCategory*nTracers] = (*levelIceVolume)(iParticle,iCategory);
    } // iCategory
    iTracer += 1;
  }

  // pond tracers
  if (useCesmMeltponds or
      useLevelMeltponds or
      useTopoMeltponds) {
    for (int iCategory = 0 ; iCategory < nCategories ; iCategory++) {
      tracerArrayCategory[iTracer + iCategory*nTracers] = (*pondArea)(iParticle,iCategory);
    } // iCategory
    iTracer += 1;
    for (int iCategory = 0 ; iCategory < nCategories ; iCategory++) {
      tracerArrayCategory[iTracer + iCategory*nTracers] = (*pondDepth)(iParticle,iCategory);
    } // iCategory
    iTracer += 1;
  }

  // level or topo ponds
  if (useLevelMeltponds or
      useTopoMeltponds) {
    for (int iCategory = 0 ; iCategory < nCategories ; iCategory++) {
      tracerArrayCategory[iTracer + iCategory*nTracers] = (*pondLidThickness)(iParticle,iCategory);
    } // iCategory
    iTracer += 1;
  }

  // aerosols
  if (useAerosols) {
    for (int iCategory = 0 ; iCategory < nCategories ; iCategory++) {
      for (int iAerosol = 0 ; iAerosol < nAerosols ; iAerosol++) {

	tracerArrayCategory[iTracer + 4*iAerosol   + iCategory*nTracers] = (*snowScatteringAerosol)(iParticle,iCategory,iAerosol);
	tracerArrayCategory[iTracer + 4*iAerosol+1 + iCategory*nTracers] = (*snowBodyAerosol      )(iParticle,iCategory,iAerosol);
	tracerArrayCategory[iTracer + 4*iAerosol+2 + iCategory*nTracers] = (*iceScatteringAerosol )(iParticle,iCategory,iAerosol);
	tracerArrayCategory[iTracer + 4*iAerosol+3 + iCategory*nTracers] = (*iceBodyAerosol       )(iParticle,iCategory,iAerosol);

      } // iAerosol
    } // iCategory
  }

}

// Get the tracer category arrays for a particle
void ColumnTracers::get_icepack_tracer_array_category(const int iParticle) {

    // layer numbers
  int nCategories = columnDimensions->size("nCategories");
  int nIceLayers  = columnDimensions->size("nIceLayers");
  int nSnowLayers = columnDimensions->size("nSnowLayers");
  int nAerosols   = columnDimensions->size("nAerosols");

  int iTracer = 0;

  // surfaceTemperature
  for (int iCategory = 0 ; iCategory < nCategories ; iCategory++) {
    (*surfaceTemperature)(iParticle,iCategory) = tracerArrayCategory[iTracer + iCategory*nTracers];
  } // iCategory
  iTracer += 1;

  // iceEnthalpy
  for (int iCategory = 0 ; iCategory < nCategories ; iCategory++) {
    for (int iIceLayer = 0 ; iIceLayer < nIceLayers ; iIceLayer++) {
      (*iceEnthalpy)(iParticle,iCategory,iIceLayer) = tracerArrayCategory[iTracer + iIceLayer + iCategory*nTracers];
    } // iIceLayer
  } // iCategory
  iTracer += nIceLayers;

  // snowEnthalpy
  for (int iCategory = 0 ; iCategory < nCategories ; iCategory++) {
    for (int iSnowLayer = 0 ; iSnowLayer < nSnowLayers ; iSnowLayer++) {
      (*snowEnthalpy)(iParticle,iCategory,iSnowLayer) = tracerArrayCategory[iTracer + iSnowLayer + iCategory*nTracers];
    } // iSnowLayer
  } // iCategory
  iTracer += nSnowLayers;

  // ice Salinity
  for (int iCategory = 0 ; iCategory < nCategories ; iCategory++) {
    for (int iIceLayer = 0 ; iIceLayer < nIceLayers ; iIceLayer++) {
      (*iceSalinity)(iParticle,iCategory,iIceLayer) = tracerArrayCategory[iTracer + iIceLayer + iCategory*nTracers];
    } // iIceLayer
  } // iCategory
  iTracer += nIceLayers;

  // iceAge
  if (useIceAge) {
    for (int iCategory = 0 ; iCategory < nCategories ; iCategory++) {
      (*iceAge)(iParticle,iCategory) = tracerArrayCategory[iTracer + iCategory*nTracers];
    } // iCategory
    iTracer += 1;
  }

  // firstYearIceArea
  if (useFirstYearIce) {
    for (int iCategory = 0 ; iCategory < nCategories ; iCategory++) {
      (*firstYearIceArea)(iParticle,iCategory) = tracerArrayCategory[iTracer + iCategory*nTracers];
    } // iCategory
    iTracer += 1;
  }

  // level ice tracers
  if (useLevelIce) {
    for (int iCategory = 0 ; iCategory < nCategories ; iCategory++) {
      (*levelIceArea)(iParticle,iCategory) = tracerArrayCategory[iTracer + iCategory*nTracers];
    } // iCategory
    iTracer += 1;
    for (int iCategory = 0 ; iCategory < nCategories ; iCategory++) {
      (*levelIceVolume)(iParticle,iCategory) = tracerArrayCategory[iTracer + iCategory*nTracers];
    } // iCategory
    iTracer += 1;
  }

  // pond tracers
  if (useCesmMeltponds or
      useLevelMeltponds or
      useTopoMeltponds) {
    for (int iCategory = 0 ; iCategory < nCategories ; iCategory++) {
      (*pondArea)(iParticle,iCategory) = tracerArrayCategory[iTracer + iCategory*nTracers];
    } // iCategory
    iTracer += 1;
    for (int iCategory = 0 ; iCategory < nCategories ; iCategory++) {
      (*pondDepth)(iParticle,iCategory) = tracerArrayCategory[iTracer + iCategory*nTracers];
    } // iCategory
    iTracer += 1;
  }

  // level or topo ponds
  if (useLevelMeltponds or
      useTopoMeltponds) {
    for (int iCategory = 0 ; iCategory < nCategories ; iCategory++) {
      (*pondLidThickness)(iParticle,iCategory) = tracerArrayCategory[iTracer + iCategory*nTracers];
    } // iCategory
    iTracer += 1;
  }

  // aerosols
  if (useAerosols) {
    for (int iCategory = 0 ; iCategory < nCategories ; iCategory++) {
      for (int iAerosol = 0 ; iAerosol < nAerosols ; iAerosol++) {

	(*snowScatteringAerosol)(iParticle,iCategory,iAerosol) = tracerArrayCategory[iTracer + 4*iAerosol   + iCategory*nTracers];
	(*snowBodyAerosol      )(iParticle,iCategory,iAerosol) = tracerArrayCategory[iTracer + 4*iAerosol+1 + iCategory*nTracers];
	(*iceScatteringAerosol )(iParticle,iCategory,iAerosol) = tracerArrayCategory[iTracer + 4*iAerosol+2 + iCategory*nTracers];
	(*iceBodyAerosol       )(iParticle,iCategory,iAerosol) = tracerArrayCategory[iTracer + 4*iAerosol+3 + iCategory*nTracers];

      } // iAerosol
    } // iCategory
  }

}

// Set the tracer particle arrays for a particle
void ColumnTracers::set_icepack_tracer_array(const int iParticle) {

  // layer numbers
  int nIceLayers  = columnDimensions->size("nIceLayers");
  int nSnowLayers = columnDimensions->size("nSnowLayers");
  int nAerosols   = columnDimensions->size("nAerosols");

  int iTracer = 0;

  // surfaceTemperature
  tracerArray[iTracer] = (*surfaceTemperatureCell)(iParticle);
  iTracer += 1;

  // iceEnthalpy
  for (int iIceLayer = 0 ; iIceLayer < nIceLayers ; iIceLayer++) {
    tracerArray[iTracer + iIceLayer] = (*iceEnthalpyCell)(iParticle,iIceLayer);
  } // iIceLayer
  iTracer += nIceLayers;

  // snowEnthalpy
  for (int iSnowLayer = 0 ; iSnowLayer < nSnowLayers ; iSnowLayer++) {
    tracerArray[iTracer + iSnowLayer] = (*snowEnthalpyCell)(iParticle,iSnowLayer);
  } // iSnowLayer
  iTracer += nSnowLayers;

  // ice Salinity
  for (int iIceLayer = 0 ; iIceLayer < nIceLayers ; iIceLayer++) {
    tracerArray[iTracer + iIceLayer] = (*iceSalinityCell)(iParticle,iIceLayer);
  } // iIceLayer
  iTracer += nIceLayers;

  // iceAge
  if (useIceAge) {
    tracerArray[iTracer] = (*iceAgeCell)(iParticle);
    iTracer += 1;
  }

  // firstYearIceArea
  if (useFirstYearIce) {
    tracerArray[iTracer] = (*firstYearIceAreaCell)(iParticle);
    iTracer += 1;
  }

  // level ice tracers
  if (useLevelIce) {
    tracerArray[iTracer] = (*levelIceAreaCell)(iParticle);
    iTracer += 1;
    tracerArray[iTracer] = (*levelIceVolumeCell)(iParticle);
    iTracer += 1;
  }

  // pond tracers
  if (useCesmMeltponds or
      useLevelMeltponds or
      useTopoMeltponds) {
    tracerArray[iTracer] = (*pondAreaCell)(iParticle);
    iTracer += 1;
    tracerArray[iTracer] = (*pondDepthCell)(iParticle);
    iTracer += 1;
  }

  // level or topo ponds
  if (useLevelMeltponds or
      useTopoMeltponds) {
    tracerArray[iTracer] = (*pondLidThicknessCell)(iParticle);
    iTracer += 1;
  }

  // aerosols
  if (useAerosols) {
    for (int iAerosol = 0 ; iAerosol < nAerosols ; iAerosol++) {

      tracerArray[iTracer + 4*iAerosol  ] = (*snowScatteringAerosolCell)(iParticle,iAerosol);
      tracerArray[iTracer + 4*iAerosol+1] = (*snowBodyAerosolCell      )(iParticle,iAerosol);
      tracerArray[iTracer + 4*iAerosol+2] = (*iceScatteringAerosolCell )(iParticle,iAerosol);
      tracerArray[iTracer + 4*iAerosol+3] = (*iceBodyAerosolCell       )(iParticle,iAerosol);

    } // iAerosol
  }

}

// Get the tracer particle arrays for a particle
void ColumnTracers::get_icepack_tracer_array(const int iParticle) {

  // layer numbers
  int nIceLayers  = columnDimensions->size("nIceLayers");
  int nSnowLayers = columnDimensions->size("nSnowLayers");
  int nAerosols   = columnDimensions->size("nAerosols");

  int iTracer = 0;

  // surfaceTemperature
  (*surfaceTemperatureCell)(iParticle) = tracerArray[iTracer];
  iTracer += 1;

  // iceEnthalpy
  for (int iIceLayer = 0 ; iIceLayer < nIceLayers ; iIceLayer++) {
    (*iceEnthalpyCell)(iParticle,iIceLayer) = tracerArray[iTracer + iIceLayer];
  } // iIceLayer
  iTracer += nIceLayers;

  // snowEnthalpy
  for (int iSnowLayer = 0 ; iSnowLayer < nSnowLayers ; iSnowLayer++) {
    (*snowEnthalpyCell)(iParticle,iSnowLayer) = tracerArray[iTracer + iSnowLayer];
  } // iSnowLayer
  iTracer += nSnowLayers;

  // ice Salinity
  for (int iIceLayer = 0 ; iIceLayer < nIceLayers ; iIceLayer++) {
    (*iceSalinityCell)(iParticle,iIceLayer) = tracerArray[iTracer + iIceLayer];
  } // iIceLayer
  iTracer += nIceLayers;

  // iceAge
  if (useIceAge) {
    (*iceAgeCell)(iParticle) = tracerArray[iTracer];
    iTracer += 1;
  }

  // firstYearIceArea
  if (useFirstYearIce) {
    (*firstYearIceAreaCell)(iParticle) = tracerArray[iTracer];
    iTracer += 1;
  }

  // level ice tracers
  if (useLevelIce) {
    (*levelIceAreaCell)(iParticle) = tracerArray[iTracer];
    iTracer += 1;
    (*levelIceVolumeCell)(iParticle) = tracerArray[iTracer];
    iTracer += 1;
  }

  // pond tracers
  if (useCesmMeltponds or
      useLevelMeltponds or
      useTopoMeltponds) {
    (*pondAreaCell)(iParticle) = tracerArray[iTracer];
    iTracer += 1;
    (*pondDepthCell)(iParticle) = tracerArray[iTracer];
    iTracer += 1;
  }

  // level or topo ponds
  if (useLevelMeltponds or
      useTopoMeltponds) {
    (*pondLidThicknessCell)(iParticle) = tracerArray[iTracer];
    iTracer += 1;
  }

  // aerosols
  if (useAerosols) {
    for (int iAerosol = 0 ; iAerosol < nAerosols ; iAerosol++) {

      (*snowScatteringAerosolCell)(iParticle,iAerosol) = tracerArray[iTracer + 4*iAerosol  ];
      (*snowBodyAerosolCell      )(iParticle,iAerosol) = tracerArray[iTracer + 4*iAerosol+1];
      (*iceScatteringAerosolCell )(iParticle,iAerosol) = tracerArray[iTracer + 4*iAerosol+2];
      (*iceBodyAerosolCell       )(iParticle,iAerosol) = tracerArray[iTracer + 4*iAerosol+3];

    } // iAerosol
  }

}

//------------------------------------------------------------------------------
// Tracer tree object describing tracer hierarchy
//------------------------------------------------------------------------------

// TracerTree constructor
TracerTree::TracerTree(DEMSI::ColumnVariable<double>* tracerIn) {
  tracer = tracerIn;
}

// TracerTree constructor
TracerTree::TracerTree(DEMSI::ColumnVariable<double>* tracerIn, DEMSI::ColumnVariable<double>* parentIn) {
  tracer = tracerIn;
  parent = parentIn;
}

// add a tracer node to the tree
void TracerTree::add(DEMSI::ColumnVariable<double>* tracer, DEMSI::ColumnVariable<double>* parent) {

  if (tracer != NULL and tracer->active()) {

    if (parent == this->tracer and tracer->active()) {
      TracerTree* newChild = new TracerTree(tracer, this->tracer);
      children.push_back(newChild);
    }

    for (int iChild = 0 ; iChild < this->children.size() ; iChild++) {
      this->children[iChild]->add(tracer, parent);
    }

  }

}

// print out the tracers in the tracer tree
void TracerTree::print(int level) {

  std::cout << std::endl << "Level: " << level << " " << (*tracer) << std::endl;

  for (int iChild = 0 ; iChild < this->children.size() ; iChild++) {
    this->children[iChild]->print(level+1);
  }

}

// get a vector of all the tracer variable pointers in the tree
std::vector<DEMSI::ColumnVariable<double>*> TracerTree::get_tracer_ptrs(void) {

  std::vector<DEMSI::ColumnVariable<double>*> ptrs;

  ptrs.push_back(this->tracer);

  for (int iChild = 0 ; iChild < this->children.size() ; iChild++) {
    std::vector<DEMSI::ColumnVariable<double>*> childPtrs;
    childPtrs = this->children[iChild]->get_tracer_ptrs();
    ptrs.insert(ptrs.end(), childPtrs.begin(), childPtrs.end());
  }

  return ptrs;

}

// convert tracers to their conserved versions
DEMSI::TracerTree* TracerTree::set_prev_pointer(DEMSI::TracerTree* prev) {

  this->prev = prev;

  DEMSI::TracerTree* prevPtr = this;
  for (int iChild = 0 ; iChild < this->children.size() ; iChild++) {
    prevPtr = this->children[iChild]->set_prev_pointer(prevPtr);
  }

  return prevPtr;

}

std::string TracerTree::name(void) {
  return this->tracer->name();
}

} // namespace DEMSI
