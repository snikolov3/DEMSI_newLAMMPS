#include <iostream>

#include "demsi_column.h"
#include "demsi_logging.h"
#include "demsi_configs.h"
#include "demsi_particles.h"
#include "demsi_partition.h"
#include "demsi_column_tracers.h"
#include "demsi_column_variables.h"
#include "demsi_time.h"
#include "demsi_grid.h"
#include "demsi_forcing.h"
#include "demsi_particles_io.h"

#include <Kokkos_Core.hpp>

#include <cmath>
#include <algorithm>
#include <iomanip>
#include <bitset>
#include <climits>
#include <fstream>

// icepack external routines
extern "C" {void column_recompute_constants(void);}
extern "C" {int column_get_string_length(void);}
extern "C" {void column_configure(void);}
extern "C" {void column_query_parameters_real(char*, double*);}
extern "C" {void column_query_parameters_integer(char*, int*);}
extern "C" {void column_query_parameters_character(char*, char*);}
extern "C" {void column_query_parameters_logical(char*, int*);}
extern "C" {void column_init_parameters_real(char*, double*);}
extern "C" {void column_init_parameters_integer(char*, int*);}
extern "C" {void column_init_parameters_character(char*, char*);}
extern "C" {void column_init_parameters_logical(char*, int*);}
extern "C" {void column_init_thermo(int*, int*, double[]);}
extern "C" {void column_init_profiles(int*, double*, double*, double*);}
extern "C" {void column_init_itd(int*, double[]);}
extern "C" {void column_init_category_areas(int*, int*, double*, double*, double*);}
extern "C" {void column_init_particle_state(int*, int*, int*, double*, double*, double*, double*, double*, double*, double*, double*);}
extern "C" {void column_init_orbit(void);}
extern "C" {void column_init_shortwave(int*, int*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*);}
extern "C" {void column_initialize_atmos_coupler_fields(int*, double*, double*, double*, double*, double*, double*);}
extern "C" {void column_initialize_ocean_coupler_fields(int*, double*, double*, double*);}
extern "C" {void column_aggregate(int*, int*, int*, int*, int*, double*, int*, int*, double*, double*, double*, double*, double*, double*, double*, double*, double*);}
extern "C" {void column_prep_radiation(int*, int*, int*, int*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*);}
extern "C" {void column_step_radiation(int*, int*, int*, int*, int*, double*, char*, int*, double*, double*, int*, double*, double*, double*, int*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, int*, int*, int*, int*, int*, int*, int*, int*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*);}
extern "C" {void column_step_therm1(int*, int*, int*, int*, double*, int*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*);}
extern "C" {void column_step_therm2(int*, int*, int*, int*, int*, int*, int*, int*, int*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, int*, double*, int*, int*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, int*, double*, double*, double*, double*);}
extern "C" {void column_step_ridge(int*, int*, int*, int*, int*, int*, int*, int*, int*, double*, double*, double*, double*, double*, double*, double*, double*, double*, int*, double*, int*, int*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*);}
extern "C" {void column_coupling_prep(int*, int*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*);}
extern "C" {void column_ocean_mixed_layer(double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*);}
extern "C" {double column_initial_air_drag_coefficient(void);}
extern "C" {double column_get_constant(char*);}
extern "C" {void column_limit_temperature(double*, double*);}
extern "C" {void column_limit_specific_humidity(double*, double*);}
extern "C" {void column_longwave_rosati_miyakoda(double*, double*, double*, double*, double*, double*, double*);}
extern "C" {void column_longwave_parkinson_and_washington(double*, double*, double*);}
extern "C" {void column_shortwave_from_cloud_fraction(double*, double*, double*, double*, double*, double*, int*);}
extern "C" {void column_split_precipitation(double*, double*, double*, double*);}
extern "C" {void column_postprocess_ocean_forcing(double*, double*, double*);}

// icepack warnings functions
extern "C" {void column_warnings_reset(void);}
extern "C" {int column_warnings_number(void);}
extern "C" {void column_warnings_getone(int*, char*);}
extern "C" {int column_warnings_aborted(void);}

namespace DEMSI {

//------------------------------------------------------------------------------
// Utils
//------------------------------------------------------------------------------

void convert_string_to_array(const std::string str, char* charArray) {

  size_t len = str.copy(charArray, str.length());
  charArray[len] = '\0';

}

// interpolate a forcing field for one particle
void interpolate_field(const Kokkos::View<double**>* fieldEulerian, DEMSI::ColumnVariable<double>* fieldParticle, const double weight[], const int i[], const int j[], const int iParticle) {

  if (fieldEulerian->use_count() > 0) {
    (*fieldParticle)(iParticle) =
      (*fieldEulerian)(i[0],j[0]) * weight[0] +
      (*fieldEulerian)(i[1],j[1]) * weight[1] +
      (*fieldEulerian)(i[2],j[2]) * weight[2] +
      (*fieldEulerian)(i[3],j[3]) * weight[3];
  }

}

template<typename T>
void print_binary_representation(const T& a) {
  const char* beg = reinterpret_cast<const char*>(&a);
  const char* end = beg + sizeof(a);
  while(beg != end)
    std::cout << std::bitset<CHAR_BIT>(*beg++) << ' ';
  std::cout << '\n';
}

//------------------------------------------------------------------------------
// Init
//------------------------------------------------------------------------------

// constructor for Column class
Column::Column(DEMSI::Partition* partitionIn, DEMSI::Log *logIn, DEMSI::Configs* configsIn, DEMSI::LammpsInstance* lammpsInstanceIn, DEMSI::Grid* gridIn, DEMSI::Particles *particlesIn, DEMSI::Tessellation* tessellationIn, DEMSI::Forcing* forcingIn, DEMSI::Clock* simulationClockIn, DEMSI::TimeInterval* timeStepIntervalIn) {

  // set pointers
  partition = partitionIn;
  log = logIn;
  configs = configsIn;
  lammpsInstance = lammpsInstanceIn;
  grid = gridIn;
  particles = particlesIn;
  tessellation = tessellationIn;
  forcing = forcingIn;
  simulationClock = simulationClockIn;
  timeStepInterval = timeStepIntervalIn;

  // compute some constants
  column_recompute_constants();

  // Are we using the column physics?
  configs->get({"ConfigGroup:useSections","Config:useColumn"}, useColumn);

  if (useColumn) {

    configs->get({"ConfigGroup:columnUseSections","Config:useColumnShortwave"}, useColumnShortwave);
    configs->get({"ConfigGroup:columnUseSections","Config:useColumnVerticalThermodynamics"}, useColumnVerticalThermodynamics);
    configs->get({"ConfigGroup:columnUseSections","Config:useColumnItdThermodynamics"}, useColumnItdThermodynamics);
    configs->get({"ConfigGroup:columnUseSections","Config:useColumnRidging"}, useColumnRidging);

    // set column configs
    (*log)(DEMSI::LOG::DEBUG) << "Column: set configs" << std::endl;
    column_configs();

  } // useColumn

}

// init the column object
void Column::init(void) {

  // what is the ocean type?
  oceanType = "embedded";
  if (configs->exists({"ConfigGroup:ocean","Config:oceanType"})) {
    configs->get({"ConfigGroup:ocean","Config:oceanType"}, oceanType);
  } // oceanType option exists

  // set up empty column variables object
  columnVariables = new DEMSI::ColumnVariables(log, partition, particles);

  // init DEMSI column dimensions
  (*log)(DEMSI::LOG::DEBUG) << "Column: Init column dimensions" << std::endl;
  init_demsi_column_dimensions();

  // allocate variables that are always used
  (*log)(DEMSI::LOG::DEBUG) << "Column: Allocate column mass variables" << std::endl;
  allocate_demsi_mass_variables();

  // init column variable processor transfer test
  (*log)(DEMSI::LOG::DEBUG) << "Column: Init column transfer test" << std::endl;
  init_processor_transfer();

  // init element area
  (*log)(DEMSI::LOG::DEBUG) << "Column: Initial element area" << std::endl;
  effective_element_area();

  // Are we using the column physics?
  if (useColumn) {

    // column tracers object
    (*log)(DEMSI::LOG::DEBUG) << "Column: Init tracer object" << std::endl;
    columnTracers = new ColumnTracers(columnDimensions, configs);

    // init DEMSI column variables
    (*log)(DEMSI::LOG::DEBUG) << "Column: Init column variables" << std::endl;
    init_demsi_column_variables();

    // get the column string length
    columnStrLen = column_get_string_length();
    (*log)(DEMSI::LOG::DEBUG) << "Column: columnStrLen: " << columnStrLen << std::endl;

    // configure icepack
    (*log)(DEMSI::LOG::DEBUG) << "Column: configure icepack" << std::endl;
    column_configure();
    check_warnings("column_configure");

    // initialize atmospheric coupler fields
    (*log)(DEMSI::LOG::DEBUG) << "Column: initialize atmos coupler fields" << std::endl;
    initialize_atmos_coupler_fields();

    // initialize oceanic coupler fields
    (*log)(DEMSI::LOG::DEBUG) << "Column: initialize ocean coupler fields" << std::endl;
    initialize_ocean_coupler_fields();

    // init vertical thermodynamics
    (*log)(DEMSI::LOG::DEBUG) << "Column: init vertical thermodynamics" << std::endl;
    init_thermo_vertical();

    // init ice thickness distribution
    (*log)(DEMSI::LOG::DEBUG) << "Column: init ice thickness distribution" << std::endl;
    init_ice_thickness_distribution();

  } else {

    // init ice thickness distribution
    (*log)(DEMSI::LOG::DEBUG) << "Column: init ice thickness distribution no column" << std::endl;
    init_ice_thickness_distribution_no_column();

  } // useColumn

  // set up mass related variables
  (*log)(DEMSI::LOG::DEBUG) << "Column: Init column mass variables" << std::endl;
  init_demsi_mass_variables();

  // Are we using the column physics?
  if (useColumn) {

    // initialize thermodynamic state by particle
    (*log)(DEMSI::LOG::DEBUG) << "Column: init ice state" << std::endl;
    init_column_state();

    // reset atmospheric coupler fluxes
    (*log)(DEMSI::LOG::DEBUG) << "Column: reset atmospheric coupler fluxes" << std::endl;
    reset_atmospheric_coupler_fluxes();

    // reset oceanic coupler fluxes
    (*log)(DEMSI::LOG::DEBUG) << "Column: reset oceanic coupler fluxes" << std::endl;
    reset_ocean_coupler_fluxes();

  } // useColumn

  // create tracer tree
  create_tracer_tree();

}

// init column physics after remapping
void Column::init_remap(void) {

  if (useColumn) {

    ratioRidgeThicknessToIce->set(1.0);
    airDensity->set(1.3);
    airLevelHeight->set(10.0);

  } // useColumn

} // Column::init_remap

// initialize column dimensions
void Column::init_demsi_column_dimensions(void) {

  columnDimensions = new DEMSI::ColumnDimensions(configs, log, useColumn);

}

// allocate mass variables that are always used
void Column::allocate_demsi_mass_variables(void) {

  // effective areas
  effectiveElementArea = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"effectiveElementArea",
	{});
  columnVariables->add(effectiveElementArea);

  // category variables
  iceAreaCategory = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"iceAreaCategory",
	{"nCategories"});
  columnVariables->add(iceAreaCategory);

  iceVolumeCategory = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"iceVolumeCategory",
	{"nCategories"});
  columnVariables->add(iceVolumeCategory);

  snowVolumeCategory = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"snowVolumeCategory",
	{"nCategories"});
  columnVariables->add(snowVolumeCategory);

  // aggregate variables
  iceAreaCell = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"iceAreaCell",
	{});
  columnVariables->add(iceAreaCell);

  iceVolumeCell = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"iceVolumeCell",
	{});
  columnVariables->add(iceVolumeCell);

  snowVolumeCell = new DEMSI::ColumnVariable<double>(particles, columnDimensions,
	"snowVolumeCell",
	{});
  columnVariables->add(snowVolumeCell);

}

// Initialize area and volume variables from input file
void Column::init_demsi_mass_variables_file(void) {

  // from particle input file
  std::string particleInputFilename;
  configs->get({"ConfigGroup:particleInput","Config:particleInputFile"}, particleInputFilename);
  DEMSI::ParticlesColumnRead* fileIn = new DEMSI::ParticlesColumnRead(partition, log, particles, particleInputFilename);

  log->check(fileIn->has_variable("iceFraction") and fileIn->has_variable("iceThickness"),
	     "file option for massInitType requires iceFraction and iceThickness in particle input file.");

  fileIn->get_variable("iceFraction", iceAreaCell);
  fileIn->get_variable("iceThickness", iceVolumeCell);
  if (fileIn->has_variable("snowThickness")) {
    fileIn->get_variable("snowThickness", snowVolumeCell);
  } else {
    snowVolumeCell->set(0.0);
  }

  delete fileIn;

  // check input values
  for (int iParticle = 0 ; iParticle < *(particles->nParticles) ; iParticle++) {
    log->check((*iceAreaCell)(iParticle) > 0.0, "Input ice area is equal to zero");
    log->check((*iceVolumeCell)(iParticle) > 0.0, "Input ice area is equal to zero");
  } // iParticle

  // set category values
  int nCategories = columnDimensions->size("nCategories");
  for (int iParticle = 0 ; iParticle < *(particles->nParticles) ; iParticle++) {

    for (int iCategory = 0 ; iCategory < nCategories ; iCategory++) {

      if ((*iceVolumeCell)(iParticle) >= categoryThicknessLimits[iCategory] and
	  (*iceVolumeCell)(iParticle) <  categoryThicknessLimits[iCategory+1]) {

	(*iceAreaCategory)   (iParticle,iCategory) = (*iceAreaCell)(iParticle);
	(*iceVolumeCategory) (iParticle,iCategory) = (*iceAreaCell)(iParticle) * (*iceVolumeCell) (iParticle);
	(*snowVolumeCategory)(iParticle,iCategory) = (*iceAreaCell)(iParticle) * (*snowVolumeCell)(iParticle);

      } else {

	(*iceAreaCategory)   (iParticle,iCategory) = 0.0;
	(*iceVolumeCategory) (iParticle,iCategory) = 0.0;
	(*snowVolumeCategory)(iParticle,iCategory) = 0.0;

      }

    } // iCategory

  } // iParticle

} // Column::init_demsi_mass_variables_file

// Initialize area and volume variables from constant config values
void Column::init_demsi_mass_variables_const(void) {

  // from config file
  double iceFractionInitial;
  configs->get({"ConfigGroup:initialization","Config:iceFraction"}, iceFractionInitial);
  iceAreaCell->set(iceFractionInitial);

  double iceThicknessInitial;
  configs->get({"ConfigGroup:initialization","Config:iceThickness"}, iceThicknessInitial);
  iceVolumeCell->set(iceThicknessInitial);

  double snowThicknessInitial = 0.0;
  if (configs->exists({"ConfigGroup:initialization","Config:snowThickess"})) {
    configs->get({"ConfigGroup:initialization","Config:snowThickess"}, snowThicknessInitial);
  }
  snowVolumeCell->set(snowThicknessInitial);

  // check input values
  for (int iParticle = 0 ; iParticle < *(particles->nParticles) ; iParticle++) {
    log->check((*iceAreaCell)(iParticle) > 0.0, "Input ice area is equal to zero");
    log->check((*iceVolumeCell)(iParticle) > 0.0, "Input ice area is equal to zero");
  } // iParticle

  // set category values
  int nCategories = columnDimensions->size("nCategories");
  for (int iParticle = 0 ; iParticle < *(particles->nParticles) ; iParticle++) {

    if (particles->type(iParticle) == 1) {

      for (int iCategory = 0 ; iCategory < nCategories ; iCategory++) {

	if ((*iceVolumeCell)(iParticle) >= categoryThicknessLimits[iCategory] and
	    (*iceVolumeCell)(iParticle) <  categoryThicknessLimits[iCategory+1]) {

	  (*iceAreaCategory)   (iParticle,iCategory) = (*iceAreaCell)(iParticle);
	  (*iceVolumeCategory) (iParticle,iCategory) = (*iceAreaCell)(iParticle) * (*iceVolumeCell) (iParticle);
	  (*snowVolumeCategory)(iParticle,iCategory) = (*iceAreaCell)(iParticle) * (*snowVolumeCell)(iParticle);

	} else {

	  (*iceAreaCategory)   (iParticle,iCategory) = 0.0;
	  (*iceVolumeCategory) (iParticle,iCategory) = 0.0;
	  (*snowVolumeCategory)(iParticle,iCategory) = 0.0;

	}

      } // iCategory

    } else if (particles->type(iParticle) == 2) {

      // coastal element
      (*iceAreaCategory)  (iParticle,nCategories-1) =   1.0;
      (*iceVolumeCategory)(iParticle,nCategories-1) = 100.0;

    } // particle type

  } // iParticle

} // Column::init_demsi_mass_variables_const

// Initialize area and volume variables from config thickness distribution
void Column::init_demsi_mass_variables_thickness_dist(void) {

  log->check(useColumn, "Column physics needed for thickness_dist mass inout option");

  int nCategories = columnDimensions->size("nCategories");

  // input areas
  std::vector<double> iceAreaCategoryIn = configs->get_double_array({"ConfigGroup:initialization","Array:iceAreaCategoryInit"}, "iceAreaCategory");
  log->check(iceAreaCategoryIn.size() == nCategories, "Wrong number of input ice areas (!= nCategories) from configs.");

  // thickness
  if (configs->exists({"ConfigGroup:initialization","Array:iceThicknessCategoryInit"})) {

    std::vector<double> iceThicknessCategoryIn = configs->get_double_array({"ConfigGroup:initialization","Array:iceThicknessCategoryInit"}, "iceThicknessCategory");
    log->check(iceThicknessCategoryIn.size() == nCategories, "Wrong number of input ice thicknesses (!= nCategories) from configs.");

    // check thicknesses in bounds
    for (int iParticle = 0 ; iParticle < *(particles->nParticles) ; iParticle++) {
      for (int iCategory = 0 ; iCategory < nCategories-1 ; iCategory++) {
	(*iceVolumeCategory)(iParticle,iCategory) = iceThicknessCategoryIn[iCategory];
      } // iCategory
    } // iParticle

  } else {

    for (int iParticle = 0 ; iParticle < *(particles->nParticles) ; iParticle++) {
      for (int iCategory = 0 ; iCategory < nCategories-1 ; iCategory++) {
	(*iceVolumeCategory)(iParticle,iCategory) = 0.5 * (categoryThicknessLimits[iCategory] + categoryThicknessLimits[iCategory+1]);
      } // iCategory
      (*iceVolumeCategory)(iParticle,nCategories-1) = categoryThicknessLimits[nCategories-1] + 1.0;
    } // iParticle

  }

  // config defines thickness categories
  // area
  for (int iParticle = 0 ; iParticle < *(particles->nParticles) ; iParticle++) {
    for (int iCategory = 0 ; iCategory < nCategories ; iCategory++) {
      (*iceAreaCategory)(iParticle,iCategory) = iceAreaCategoryIn[iCategory];
    } // iCategory
  } // iParticle

  // volume
  for (int iParticle = 0 ; iParticle < *(particles->nParticles) ; iParticle++) {
    for (int iCategory = 0 ; iCategory < nCategories ; iCategory++) {
      (*iceVolumeCategory)(iParticle,iCategory) *= (*iceAreaCategory)(iParticle,iCategory);
    } // iCategory
  } // iParticle

  // coastal elements
  for (int iParticle = 0 ; iParticle < *(particles->nParticles) ; iParticle++) {
    if (particles->type(iParticle) == 2) {
      for (int iCategory = 0 ; iCategory < nCategories-1 ; iCategory++) {
	(*iceAreaCategory)   (iParticle,iCategory) = 0.0;
	(*iceVolumeCategory) (iParticle,iCategory) = 0.0;
	(*snowVolumeCategory)(iParticle,iCategory) = 0.0;
      } // iCategory
      (*iceAreaCategory)   (iParticle,nCategories-1) =   1.0;
      (*iceVolumeCategory) (iParticle,nCategories-1) = 100.0;
      (*snowVolumeCategory)(iParticle,nCategories-1) =   0.0;
    } // coastal element
  } // iParticle

} // Column::init_demsi_mass_variables_thickness_dist

// Initialize area and volume variables from cice method
void Column::init_demsi_mass_variables_cice(void) {

  log->check(useColumn, "Column physics needed for cice mass inout option");

  int nCategories = columnDimensions->size("nCategories");

  // thickness
  for (int iParticle = 0 ; iParticle < *(particles->nParticles) ; iParticle++) {
    for (int iCategory = 0 ; iCategory < nCategories-1 ; iCategory++) {
      (*iceVolumeCategory)(iParticle,iCategory) = 0.5 * (categoryThicknessLimits[iCategory] + categoryThicknessLimits[iCategory+1]);
    } // iCategory
    (*iceVolumeCategory)(iParticle,nCategories-1) = categoryThicknessLimits[nCategories-1] + 1.0;
  } // iParticle

  // areas
  double thicknessWithLargestArea = 3.0;
  for (int iParticle = 0 ; iParticle < *(particles->nParticles) ; iParticle++) {

    // initialize sum of areas in categories
    double areaCategorySum = 0.0;
    for (int iCategory = 0 ; iCategory < nCategories-1 ; iCategory++) {

      // parabola, max at h=thicknessWithLargestArea, zero at h=0, 2*thicknessWithLargestArea
      (*iceAreaCategory)(iParticle,iCategory) =
	std::max(0.0, (2.0 * thicknessWithLargestArea * (*iceVolumeCategory)(iParticle,iCategory) - std::pow((*iceVolumeCategory)(iParticle,iCategory),2)));

      areaCategorySum += (*iceAreaCategory)(iParticle,iCategory);

    } // iCategory

    // normalize
    double puny = 1.0e-11;
    for (int iCategory = 0 ; iCategory < nCategories-1 ; iCategory++) {
      (*iceAreaCategory)(iParticle,iCategory) /= (areaCategorySum + puny / double(nCategories));
    } // iCategory

  } // iParticle

  // volume
  for (int iParticle = 0 ; iParticle < *(particles->nParticles) ; iParticle++) {
    for (int iCategory = 0 ; iCategory < nCategories-1 ; iCategory++) {
      (*iceVolumeCategory)(iParticle,iCategory) *= (*iceAreaCategory)(iParticle,iCategory);
    } // iCategory
  } // iParticle

  // snow
  double initialCategorySnowThickness = 0.2;
  for (int iParticle = 0 ; iParticle < *(particles->nParticles) ; iParticle++) {
    for (int iCategory = 0 ; iCategory < nCategories-1 ; iCategory++) {
      (*snowVolumeCategory)(iParticle,iCategory) = std::min((*iceAreaCategory)(iParticle,iCategory) * initialCategorySnowThickness, 0.2 * (*iceVolumeCategory)(iParticle,iCategory));
    } // iCategory
  } // iParticle

  // coastal elements
  for (int iParticle = 0 ; iParticle < *(particles->nParticles) ; iParticle++) {
    if (particles->type(iParticle) == 2) {
      for (int iCategory = 0 ; iCategory < nCategories-1 ; iCategory++) {
	(*iceAreaCategory)   (iParticle,iCategory) = 0.0;
	(*iceVolumeCategory) (iParticle,iCategory) = 0.0;
	(*snowVolumeCategory)(iParticle,iCategory) = 0.0;
      } // iCategory
      (*iceAreaCategory)   (iParticle,nCategories-1) =   1.0;
      (*iceVolumeCategory) (iParticle,nCategories-1) = 100.0;
      (*snowVolumeCategory)(iParticle,nCategories-1) =   0.0;
    } // coastal element
  } // iParticle

} // Column::init_demsi_mass_variables_cice

// initialize mass variables that are always used
void Column::init_demsi_mass_variables(void) {

  std::string massInitType; // 'file', 'const', 'thickness_dist', or 'cice'
  configs->get({"ConfigGroup:initialization","Config:massInitType"}, massInitType);

  if (massInitType == "file") {
    init_demsi_mass_variables_file();

  } else if (massInitType == "const") {
    init_demsi_mass_variables_const();

  } else if (massInitType == "thickness_dist") {
    init_demsi_mass_variables_thickness_dist();

  } else if (massInitType == "cice") {
    init_demsi_mass_variables_cice();

  } else {
    log->abort("Unknown config value: ConfigGroup:initialization, Config:massInitType: ", massInitType);

  }

} // Column::init_demsi_mass_variables

// init element area
void Column::effective_element_area(void) {

  bool useRemapping;
  configs->get({"ConfigGroup:remapping","Config:useRemapping"}, useRemapping, true);
  if ((useColumn and useColumnRidging) or useRemapping) {

    // initialize effective areas
    for (int iParticle = 0 ; iParticle < *(particles->nParticles) ; iParticle++) {

      int iParticleInitClosest;
      double minDistance = 1e30;

      for (int iParticleInit = 0 ; iParticleInit < tessellation->nParticlesInit ; iParticleInit++) {

	double dx = particles->x(iParticle,0) - tessellation->xInit[iParticleInit];
	double dy = particles->x(iParticle,1) - tessellation->yInit[iParticleInit];

	double distance = std::sqrt(dx*dx + dy*dy);

	if (distance < minDistance) {
	  minDistance = distance;
	  iParticleInitClosest = iParticleInit;
	}

      } // iParticleInit
      int iParticleInit = iParticleInitClosest;

      log->check(minDistance < 10.0, "Ridging: Did not find new particle amongst init", std::to_string(minDistance));

      (*effectiveElementArea)(iParticle) = tessellation->initPolygons[iParticleInit].area();

    } // iParticle

  } else {

    auto PI = LAMMPS_NS::MathConst::MY_PI;

    for (int iParticle = 0 ; iParticle < *(particles->nParticles) ; iParticle++) {
      (*effectiveElementArea)(iParticle) = PI * std::pow(particles->radius[iParticle],2);
    } // iParticle

  } // effective element type

} // Column::effective_element_area

// initialize DEMSI column variables
void Column::init_demsi_column_variables(void) {

  // flags for whether variables will be used
  bool useIceAge;
  bool useFirstYearIce;
  bool useLevelIce;
  bool usePonds;
  bool useLidThickness;
  bool useAerosols;
  bool useFormDrag;

  bool useCesmMeltponds;
  bool useLevelMeltponds;
  bool useTopoMeltponds;

  configs->get({"ConfigGroup:columnTracers","Config:useIceAge"},         useIceAge);
  configs->get({"ConfigGroup:columnTracers","Config:useFirstYearIce"},   useFirstYearIce);
  configs->get({"ConfigGroup:columnTracers","Config:useLevelIce"},       useLevelIce);
  configs->get({"ConfigGroup:columnTracers","Config:useCesmMeltponds"},  useCesmMeltponds);
  configs->get({"ConfigGroup:columnTracers","Config:useLevelMeltponds"}, useLevelMeltponds);
  configs->get({"ConfigGroup:columnTracers","Config:useTopoMeltponds"},  useTopoMeltponds);
  configs->get({"ConfigGroup:columnTracers","Config:useAerosols"},       useAerosols);

  usePonds = (useCesmMeltponds or useLevelMeltponds or useTopoMeltponds);
  useLidThickness = (useLevelMeltponds or useTopoMeltponds);

  configs->get({"ConfigGroup:atmosphere","Config:useFormDrag"}, useFormDrag);

  // allocate variables
#include "demsi_column_variables_allocation.inc"

  // set array pointers for inactive variables
  columnVariables->set_arrays_inactive_variables();

  // initialize to non-zero constant values
  levelIceArea->set(1.0);
  levelIceVolume->set(1.0);
  ratioRidgeThicknessToIce->set(1.0);
  airDensity->set(1.3);
  airLevelHeight->set(10.0);
  seaSurfaceSalinity->set(34.0);
  airTemperature->set(253.0);
  airPotentialTemperature->set(253.0);
  airSpecificHumidity->set(0.0006);
  longwaveDown->set(180.0);
  oceanMixedLayerDepth->set(20.0);

  // set tracers pointer in tracer object
  columnTracers->set_tracer_variable_pointers(
	iceAreaCategory,
	iceVolumeCategory,
	snowVolumeCategory,
	surfaceTemperature,
	iceEnthalpy,
	iceSalinity,
	snowEnthalpy,
	iceAge,
	firstYearIceArea,
	levelIceArea,
	levelIceVolume,
	pondArea,
	pondDepth,
	pondLidThickness,
	snowScatteringAerosol,
	snowBodyAerosol,
	iceScatteringAerosol,
	iceBodyAerosol);

  // set aggregated tracers pointer in aggregated tracer object
  columnTracers->set_aggregated_tracer_variable_pointers(
	iceAreaCell,
	iceVolumeCell,
	snowVolumeCell,
	surfaceTemperatureCell,
	iceEnthalpyCell,
	iceSalinityCell,
	snowEnthalpyCell,
	iceAgeCell,
	firstYearIceAreaCell,
	levelIceAreaCell,
	levelIceVolumeCell,
	pondAreaCell,
	pondDepthCell,
	pondLidThicknessCell,
	snowScatteringAerosolCell,
	snowBodyAerosolCell,
	iceScatteringAerosolCell,
	iceBodyAerosolCell);

}

// init the shortwave calculation
void Column::init_shortwave(void) {

  if (useColumn) {

    bool doRestart;
    configs->get({"ConfigGroup:simulationTiming","Config:doRestart"}, doRestart);

    if (useColumnShortwave and not doRestart) {

      // init orbital data
      std::string shortwaveType;
      configs->get({"ConfigGroup:shortwave","Config:shortwaveType"}, shortwaveType);

      if (shortwaveType == "dEdd") {
	column_init_orbit();
	check_warnings("column_init_orbit");
      }

      // run radiation
      run_step_radiation(true);

      // aggregate albedos etc
      int nCategories  = columnDimensions->size("nCategories");
      int doRestartIn = (int) doRestart;

      for (int iParticle = 0 ; iParticle < *(particles->nParticles) ; iParticle++) {

	if (particles->type(iParticle) == 1) {

	  column_init_shortwave(&nCategories,
				&doRestartIn,
				iceAreaCategory->get(iParticle),
				albedoVisibleDirectCategory->get(iParticle),
				albedoVisibleDiffuseCategory->get(iParticle),
				albedoIRDirectCategory->get(iParticle),
				albedoIRDiffuseCategory->get(iParticle),
				albedoVisibleDirectCell->get(iParticle),
				albedoVisibleDiffuseCell->get(iParticle),
				albedoIRDirectCell->get(iParticle),
				albedoIRDiffuseCell->get(iParticle),
				solarZenithAngleCosine->get(iParticle),
				bareIceAlbedoCategory->get(iParticle),
				snowAlbedoCategory->get(iParticle),
				pondAlbedoCategory->get(iParticle),
				bareIceAlbedoCell->get(iParticle),
				snowAlbedoCell->get(iParticle),
				pondAlbedoCell->get(iParticle),
				effectivePondAreaCategory->get(iParticle),
				effectivePondAreaCell->get(iParticle),
				albedoVisibleDirectArea->get(iParticle),
				albedoVisibleDiffuseArea->get(iParticle),
				albedoIRDirectArea->get(iParticle),
				albedoIRDiffuseArea->get(iParticle),
				shortwaveScalingFactor->get(iParticle),
				shortwaveVisibleDirectDown->get(iParticle),
				shortwaveVisibleDiffuseDown->get(iParticle),
				shortwaveIRDirectDown->get(iParticle),
				shortwaveIRDiffuseDown->get(iParticle));

	}

      } // iParticle

    } // useColumnShortwave and not doRestart

  } // useColumn

}

// initialize the sst
void Column::init_sst(void) {

  if (useColumn and oceanType == "embedded") {

    bool doRestart;
    configs->get({"ConfigGroup:simulationTiming","Config:doRestart"}, doRestart);
    if (not doRestart) {

      int i[4], j[4];
      double weight[4];

      // interpolation to particles
      for (int iParticle = 0 ; iParticle < *(particles->nParticles) ; iParticle++) {

	if (particles->type(iParticle) == 1) {

	  double x = particles->x(iParticle,0);
	  double y = particles->x(iParticle,1);

	  grid->interpolation_weights(x, y, i, j, weight);

	  // atmospheric forcing
	  interpolate_field(&(forcing->seaSurfaceTemperature), seaSurfaceTemperature, weight, i, j, iParticle);

	  (*seaSurfaceTemperature)(iParticle) = std::max((*seaSurfaceTemperature)(iParticle), (*seaFreezingTemperature)(iParticle));

	}

      } // iParticle

    } // not doRestart

  } // useColumn

}

// create the tracer tree
void Column::create_tracer_tree(void) {

  // NB need to correct for other melt pond options!

  // ice area as root parent of all tracers
  tracerTree = new DEMSI::TracerTree(iceAreaCategory);

  // children of iceAreaCategory
  tracerTree->add(iceVolumeCategory, iceAreaCategory);
  tracerTree->add(snowVolumeCategory, iceAreaCategory);
  tracerTree->add(surfaceTemperature, iceAreaCategory);
  tracerTree->add(firstYearIceArea, iceAreaCategory);
  tracerTree->add(levelIceArea, iceAreaCategory);

  // children of iceVolumeCategory
  tracerTree->add(iceEnthalpy, iceVolumeCategory);
  tracerTree->add(iceSalinity, iceVolumeCategory);
  tracerTree->add(iceAge, iceVolumeCategory);
  tracerTree->add(levelIceVolume, iceVolumeCategory);
  tracerTree->add(iceScatteringAerosol, iceVolumeCategory);
  tracerTree->add(iceBodyAerosol, iceVolumeCategory);

  // children of snowVolumeCategory
  tracerTree->add(snowEnthalpy, snowVolumeCategory);
  tracerTree->add(snowScatteringAerosol, snowVolumeCategory);
  tracerTree->add(snowBodyAerosol, snowVolumeCategory);

  // children of level ice area
  tracerTree->add(pondArea, levelIceArea);

  // children of pond area
  tracerTree->add(pondDepth, pondArea);
  tracerTree->add(pondLidThickness, pondArea);

  // get vector of tracer pointers
  columnTracerPtrs = tracerTree->get_tracer_ptrs();

  // set prev pointers
  DEMSI::TracerTree* tracerPtr = tracerTree->set_prev_pointer(NULL);

  // set next pointers
  tracerPtr->next = NULL;
  while (tracerPtr != NULL) {
    DEMSI::TracerTree* tracerPrev = tracerPtr->prev;
    if (tracerPrev != NULL) {
      tracerPrev->next = tracerPtr;
    }
    tracerPtr = tracerPtr->prev;
  }

#ifdef DEBUG
  std::cout << "Forward tracers: " << std::endl;
  DEMSI::TracerTree* tracer = tracerTree;
  while (tracer != NULL) {
    std::cout << tracer << " " << tracer->name() << std::endl;
    tracer = tracer->next;
  }
  exit(0);
#endif

}

//------------------------------------------------------------------------------
// Time stepping
//------------------------------------------------------------------------------

// column time stepping before dynamics
void Column::pre_dynamics(void) {

  if (useColumn) {

    // Reinitialize diagnostic variables
    (*log)(DEMSI::LOG::DEBUG) << "Column: reinitialize diagnostics" << std::endl << std::flush;
    reinitialize_diagnostics();

    // Prepare shortwave radiation calculation
    (*log)(DEMSI::LOG::DEBUG) << "Column: prepare radiation calculation" << std::endl << std::flush;
    run_prep_radiation();

    // Vertical thermodynamics
    (*log)(DEMSI::LOG::DEBUG) << "Column: vertical thermodynamics" << std::endl << std::flush;
    run_step_therm1();

    // ITD thermodynamics
    (*log)(DEMSI::LOG::DEBUG) << "Column: ITD thermodynamics" << std::endl << std::flush;
    run_step_therm2();

    // Aggregate tracers
    (*log)(DEMSI::LOG::DEBUG) << "Column: pre dynamics aggregate tracers" << std::endl << std::flush;
    aggregate();

  } // useColumn

}

// column time stepping during dynamics
void Column::dynamics(void) {

  if (useColumn) {

    // Run ridging
    (*log)(DEMSI::LOG::DEBUG) << "Column: ridging" << std::endl << std::flush;
    run_step_ridge();

    // Aggregate tracers
    (*log)(DEMSI::LOG::DEBUG) << "Column: dynamics aggregate tracers" << std::endl << std::flush;
    aggregate();


  } // useColumn

}

// column time stepping after dynamics
void Column::post_dynamics(void) {

  if (useColumn) {

    // Calculate the ocean SST for an embedded ocean
    if (oceanType == "embedded") {
      ocean_mixed_layer();
    }

    // Calculate shortwave radiation
    (*log)(DEMSI::LOG::DEBUG) << "Column: shortwave radiation" << std::endl << std::flush;
    run_step_radiation(false);

    // Prepare coupling
    (*log)(DEMSI::LOG::DEBUG) << "Column: prepare coupling" << std::endl << std::flush;
    coupling_prep();

  } // useColumn

}

// Interpolate column forcing to particles and post-process
void Column::interpolate_forcing(void) {

  if (useColumn) {

    int i[4], j[4];
    double weight[4];

    // interpolation to particles
    for (int iParticle = 0 ; iParticle < *(particles->nParticles) ; iParticle++) {

      if (particles->type(iParticle) == 1) {

	double x = particles->x(iParticle,0);
	double y = particles->x(iParticle,1);

	grid->interpolation_weights(x, y, i, j, weight);

	// grid
	interpolate_field(&(grid->latitude),  latitude,  weight, i, j, iParticle);
	interpolate_field(&(grid->longitude), longitude, weight, i, j, iParticle);

	// atmospheric forcing
	interpolate_field(&(forcing->airTemperature),    airTemperature,      weight, i, j, iParticle);
	interpolate_field(&(forcing->specificHumidity),  airSpecificHumidity, weight, i, j, iParticle);
	interpolate_field(&(forcing->shortwaveFluxDown), shortwaveDown,       weight, i, j, iParticle);
	interpolate_field(&(forcing->longwaveFluxDown),  longwaveDown,        weight, i, j, iParticle);
	interpolate_field(&(forcing->xAtmWind),          uAirVelocity,        weight, i, j, iParticle);
	interpolate_field(&(forcing->yAtmWind),          vAirVelocity,        weight, i, j, iParticle);
	interpolate_field(&(forcing->snowfallRate),      snowfallRate,        weight, i, j, iParticle);
	interpolate_field(&(forcing->rainfallRate),      rainfallRate,        weight, i, j, iParticle);
	interpolate_field(&(forcing->precipitationRate), precipitationRate,   weight, i, j, iParticle);
	interpolate_field(&(forcing->cloudFraction),     cloudFraction,       weight, i, j, iParticle);

      }

    } // iParticle

    // limit the air temperature
    limit_temperature();

    // limit the air specific humidity
    limit_specific_humidity();

    // calculate the downwave radiation from cloudiness, longitude and time of day
    shortwave_from_cloud_fraction();

    // split shortwave radiation between bands
    split_shortwave();

    // apply physical limits to forcing variables
    apply_physical_limits();

    // calculate downwelling longwave radiation
    calculate_longwave_down();

    // scale precipitation by factor
    scale_precipitation();

    // split precipitation between snow and rain from air temperature
    split_precipitation();

    // post processing of forcing
    char airSpecificHeatStr[7] = "cp_air";
    double airSpecificHeat = column_get_constant(airSpecificHeatStr);
    char latentHeatSublimationStr[5] = "Lsub";
    double latentHeatSublimation = column_get_constant(latentHeatSublimationStr);

    for (int iParticle = 0 ; iParticle < *(particles->nParticles) ; iParticle++) {

      if (particles->type(iParticle) == 1) {

	// wind speed
	(*windSpeed)(iParticle) = std::sqrt(std::pow((*uAirVelocity)(iParticle),2) + std::pow((*vAirVelocity)(iParticle),2));

	// air potential temperature
	(*airPotentialTemperature)(iParticle) = (*airTemperature)(iParticle);

	// transfer coefficients
	(*sensibleTransferCoefficient)(iParticle) = 1.20e-3 * airSpecificHeat       * (*airDensity)(iParticle) * (*windSpeed)(iParticle);
	(*latentTransferCoefficient)(iParticle)   = 1.50e-3 * latentHeatSublimation * (*airDensity)(iParticle) * (*windSpeed)(iParticle);

	// air stresses
	double airStressCoefficient = 0.0012 * (*airDensity)(iParticle) * (*windSpeed)(iParticle);

	(*uAirStress)(iParticle) = (*uAirVelocity)(iParticle) * airStressCoefficient;
	(*vAirStress)(iParticle) = (*vAirVelocity)(iParticle) * airStressCoefficient;

      }

    } // iParticle

    // interpolation of ocean variables to particles
    for (int iParticle = 0 ; iParticle < *(particles->nParticles) ; iParticle++) {

      if (particles->type(iParticle) == 1) {

	double x = particles->x(iParticle,0);
	double y = particles->x(iParticle,1);

	grid->interpolation_weights(x, y, i, j, weight);

	// oceanic forcing
	interpolate_field(&(forcing->xOcnCurrents), uOceanVelocity, weight, i, j, iParticle);
	interpolate_field(&(forcing->yOcnCurrents), vOceanVelocity, weight, i, j, iParticle);

      }

    } // iParticle

    // interpolation of embedded ocean variables to particles
    for (int iParticle = 0 ; iParticle < *(particles->nParticles) ; iParticle++) {

      if (particles->type(iParticle) == 1) {

	double x = particles->x(iParticle,0);
	double y = particles->x(iParticle,1);

	grid->interpolation_weights(x, y, i, j, weight);

	interpolate_field(&(forcing->seaSurfaceSalinity),       seaSurfaceSalinity,       weight, i, j, iParticle);
	interpolate_field(&(forcing->seaSurfaceTiltU),          seaSurfaceTiltU,          weight, i, j, iParticle);
	interpolate_field(&(forcing->seaSurfaceTiltV),          seaSurfaceTiltV,          weight, i, j, iParticle);
	interpolate_field(&(forcing->oceanMixedLayerDepth),     oceanMixedLayerDepth,     weight, i, j, iParticle);
	interpolate_field(&(forcing->oceanHeatFluxConvergence), oceanHeatFluxConvergence, weight, i, j, iParticle);

	column_postprocess_ocean_forcing(seaSurfaceSalinity->get(iParticle),
					 oceanMixedLayerDepth->get(iParticle),
					 seaFreezingTemperature->get(iParticle));

      }

    } // iParticle

  } // useColumn

}

// Time stepping of the icepack radiation preparation subroutine
void Column::run_prep_radiation(void) {

  bool calcSurfaceTemperature;
  configs->get({"ConfigGroup:atmosphere","Config:calcSurfaceTemperature"}, calcSurfaceTemperature);
  if (useColumnShortwave and calcSurfaceTemperature) {

    int nCategories  = columnDimensions->size("nCategories");
    int nIceLayers   = columnDimensions->size("nIceLayers");
    int nIceLayersP1 = columnDimensions->size("nIceLayersP1");
    int nSnowLayers  = columnDimensions->size("nSnowLayers");

    for (int iParticle = 0 ; iParticle < *(particles->nParticles) ; iParticle++) {

      if (particles->type(iParticle) == 1) {

	column_prep_radiation(&nCategories,
			      &nIceLayers,
			      &nIceLayersP1,
			      &nSnowLayers,
			      iceAreaCell->get(iParticle),
			      iceAreaCategory->get(iParticle),
			      shortwaveVisibleDirectDown->get(iParticle),
			      shortwaveVisibleDiffuseDown->get(iParticle),
			      shortwaveIRDirectDown->get(iParticle),
			      shortwaveIRDiffuseDown->get(iParticle),
			      albedoVisibleDirectArea->get(iParticle),
			      albedoVisibleDiffuseArea->get(iParticle),
			      albedoIRDirectArea->get(iParticle),
			      albedoIRDiffuseArea->get(iParticle),
			      shortwaveScalingFactor->get(iParticle),
			      surfaceShortwaveFlux->get(iParticle),
			      interiorShortwaveFlux->get(iParticle),
			      penetratingShortwaveFlux->get(iParticle),
			      shortwaveLayerPenetration->get(iParticle),
			      absorbedShortwaveSnowLayer->get(iParticle),
			      absorbedShortwaveIceLayer->get(iParticle));

      }

    } // iParticle

  } // useColumnShortwave

}

// Time stepping of the icepack radiation subroutine
void Column::run_step_radiation(const bool lInitialization) {

  if (useColumnShortwave) {

    int nCategories  = columnDimensions->size("nCategories");
    int nIceLayers   = columnDimensions->size("nIceLayers");
    int nIceLayersP1 = columnDimensions->size("nIceLayersP1");
    int nSnowLayers  = columnDimensions->size("nSnowLayers");
    int nAerosols    = columnDimensions->size("nAerosols");

    int lInitializationIn = (int) lInitialization;

    DEMSI::Time currentTime = simulationClock->get_time();

    char calendarTypeIn[columnStrLen];
    int daysInYear = currentTime.days_in_year();
    double dayOfNextShortwaveCalculation = 0.0;
    double dayOfYear = (double) currentTime.day_of_year();
    int secondsIntoDay = currentTime.seconds_into_day();

    double thermoTimeStep = timeStepInterval->seconds();

    // aerosol/bgc dummy arrays
    int nzAerosols = 3;
    int nSpectralIntervals = 3;
    int nModal1 = 10;
    int nModal2 = 9;

    int nShortwaveBio = (nIceLayers+nSnowLayers+2)*(1+nzAerosols);

    int maxBCType = 2;
    int maxDustType = 4;
    int maxAerosolType = maxBCType+maxDustType;

    int indexShortwaveAerosol[maxAerosolType];
    for (int i = 0 ; i < maxAerosolType ; i++) indexShortwaveAerosol[i] = 1;
    double verticalShortwaveGrid[nIceLayersP1];
    for (int i = 0 ; i < nIceLayersP1 ; i++) verticalShortwaveGrid[i] = 0.0;
    double verticalGrid[nIceLayersP1];
    for (int i = 0 ; i < nIceLayersP1 ; i++) verticalGrid[i] = 0.0;
    double brineFraction[nCategories];
    for (int i = 0 ; i < nCategories ; i++) brineFraction[i] = 0.0;
    double aerosolMassExtinctionCrossSection[nSpectralIntervals*maxAerosolType];
    for (int i = 0 ; i < nSpectralIntervals*maxAerosolType ; i++) aerosolMassExtinctionCrossSection[i] = 0.0;
    double aerosolSingleScatterAlbedo[nSpectralIntervals*maxAerosolType];
    for (int i = 0 ; i < nSpectralIntervals*maxAerosolType ; i++) aerosolSingleScatterAlbedo[i] = 0.0;
    double aerosolAsymmetryParameter[nSpectralIntervals*maxAerosolType];
    for (int i = 0 ; i < nSpectralIntervals*maxAerosolType ; i++) aerosolAsymmetryParameter[i] = 0.0;
    double modalMassExtinctionCrossSection[nSpectralIntervals*nModal1];
    for (int i = 0 ; i < nSpectralIntervals*nModal1 ; i++) modalMassExtinctionCrossSection[i] = 0.0;
    double modalSingleScatterAlbedo[nSpectralIntervals*nModal1];
    for (int i = 0 ; i < nSpectralIntervals*nModal1 ; i++) modalSingleScatterAlbedo[i] = 0.0;
    double modalAsymmetryParameter[nSpectralIntervals*nModal1];
    for (int i = 0 ; i < nSpectralIntervals*nModal1 ; i++) modalAsymmetryParameter[i] = 0.0;
    double bioTracerShortwave[nShortwaveBio*nCategories];
    for (int i = 0 ; i < nShortwaveBio*nCategories ; i++) bioTracerShortwave[i] = 0.0;
    double modalBCabsorptionParameter[nSpectralIntervals*nModal1*nModal2];
    for (int i = 0 ; i < nSpectralIntervals*nModal1*nModal2 ; i++) modalBCabsorptionParameter[i] = 0.0;

    for (int iParticle = 0 ; iParticle < *(particles->nParticles) ; iParticle++) {

      if (particles->type(iParticle) == 1) {

	columnTracers->set_icepack_tracer_array_category(iParticle);

	column_step_radiation(&nCategories,
			      &nIceLayers,
			      &nSnowLayers,
			      &nAerosols,
			      &nIceLayersP1,
			      &thermoTimeStep,
			      calendarTypeIn,
			      &daysInYear,
			      &dayOfNextShortwaveCalculation,
			      &dayOfYear,
			      &secondsIntoDay,
			      latitude->get(iParticle),
			      longitude->get(iParticle),
			      columnTracers->tracerArrayCategory,
			      &columnTracers->nTracers,
			      iceAreaCategory->get(iParticle),
			      iceVolumeCategory->get(iParticle),
			      snowVolumeCategory->get(iParticle),
			      surfaceTemperature->get(iParticle),
			      levelIceArea->get(iParticle),
			      pondArea->get(iParticle),
			      pondDepth->get(iParticle),
			      pondLidThickness->get(iParticle),
			      shortwaveVisibleDirectDown->get(iParticle),
			      shortwaveVisibleDiffuseDown->get(iParticle),
			      shortwaveIRDirectDown->get(iParticle),
			      shortwaveIRDiffuseDown->get(iParticle),
			      solarZenithAngleCosine->get(iParticle),
			      snowfallRate->get(iParticle),
			      albedoVisibleDirectCategory->get(iParticle),
			      albedoVisibleDiffuseCategory->get(iParticle),
			      albedoIRDirectCategory->get(iParticle),
			      albedoIRDiffuseCategory->get(iParticle),
			      surfaceShortwaveFlux->get(iParticle),
			      interiorShortwaveFlux->get(iParticle),
			      penetratingShortwaveFlux->get(iParticle),
			      shortwaveLayerPenetration->get(iParticle),
			      absorbedShortwaveSnowLayer->get(iParticle),
			      absorbedShortwaveIceLayer->get(iParticle),
			      bareIceAlbedoCategory->get(iParticle),
			      snowAlbedoCategory->get(iParticle),
			      pondAlbedoCategory->get(iParticle),
			      effectivePondAreaCategory->get(iParticle),
			      snowFractionCategory->get(iParticle),
			      pondSnowDepthDifference->get(iParticle),
			      pondLidMeltFluxFraction->get(iParticle),
			      snowScatteringAerosol->get(iParticle),
			      snowBodyAerosol->get(iParticle),
			      iceScatteringAerosol->get(iParticle),
			      iceBodyAerosol->get(iParticle),
			      &lInitializationIn,
			      &nzAerosols,
			      &nSpectralIntervals,
			      &nModal1,
			      &nModal2,
			      &nShortwaveBio,
			      &maxAerosolType,
			      indexShortwaveAerosol,
			      verticalShortwaveGrid,
			      verticalGrid,
			      brineFraction,
			      aerosolMassExtinctionCrossSection,
			      aerosolSingleScatterAlbedo,
			      aerosolAsymmetryParameter,
			      modalMassExtinctionCrossSection,
			      modalSingleScatterAlbedo,
			      modalAsymmetryParameter,
			      bioTracerShortwave,
			      modalBCabsorptionParameter);

	check_warnings("column_step_radiation", iParticle);

	columnTracers->get_icepack_tracer_array_category(iParticle);

      }

      } // iParticle

    } // useColumnShortwave

}

// Time stepping of the icepack vertical thermodyanmic subroutine
void Column::run_step_therm1(void) {

  if (useColumnVerticalThermodynamics) {

    int nCategories = columnDimensions->size("nCategories");
    int nIceLayers  = columnDimensions->size("nIceLayers");
    int nSnowLayers = columnDimensions->size("nSnowLayers");
    int nAerosols   = columnDimensions->size("nAerosols");

    double thermoTimeStep = timeStepInterval->seconds();

    DEMSI::Time currentTime = simulationClock->get_time();

    double dayOfYear = (double) currentTime.day_of_year();

    int useAerosolsIn = (int) columnTracers->useAerosols;

    for (int iParticle = 0 ; iParticle < *(particles->nParticles) ; iParticle++) {

      if (particles->type(iParticle) == 1) {

	double uVelocity = particles->v(iParticle,0);
	double vVelocity = particles->v(iParticle,1);

	column_step_therm1(&nCategories,
			   &nIceLayers,
			   &nSnowLayers,
			   &nAerosols,
			   &thermoTimeStep,
			   &useAerosolsIn,
			   latitude->get(iParticle),
			   iceAreaCellInitial->get(iParticle),
			   iceAreaCategoryInitial->get(iParticle),
			   iceVolumeCategoryInitial->get(iParticle),
			   snowVolumeCategoryInitial->get(iParticle),
			   iceAreaCell->get(iParticle),
			   iceAreaCategory->get(iParticle),
			   iceVolumeCell->get(iParticle),
			   iceVolumeCategory->get(iParticle),
			   snowVolumeCell->get(iParticle),
			   snowVolumeCategory->get(iParticle),
			   &uVelocity,
			   &vVelocity,
			   surfaceTemperature->get(iParticle),
			   snowEnthalpy->get(iParticle),
			   iceEnthalpy->get(iParticle),
			   iceSalinity->get(iParticle),
			   levelIceArea->get(iParticle),
			   levelIceVolume->get(iParticle),
			   pondArea->get(iParticle),
			   pondDepth->get(iParticle),
			   pondLidThickness->get(iParticle),
			   iceAge->get(iParticle),
			   firstYearIceArea->get(iParticle),
			   uAirVelocity->get(iParticle),
			   vAirVelocity->get(iParticle),
			   windSpeed->get(iParticle),
			   airLevelHeight->get(iParticle),
			   airSpecificHumidity->get(iParticle),
			   airDensity->get(iParticle),
			   airTemperature->get(iParticle),
			   atmosReferenceTemperature2m->get(iParticle),
			   atmosReferenceHumidity2m->get(iParticle),
			   atmosReferenceSpeed10m->get(iParticle),
			   airOceanDragCoefficientRatio->get(iParticle),
			   oceanDragCoefficient->get(iParticle),
			   oceanDragCoefficientSkin->get(iParticle),
			   oceanDragCoefficientFloe->get(iParticle),
			   oceanDragCoefficientKeel->get(iParticle),
			   airDragCoefficient->get(iParticle),
			   airDragCoefficientSkin->get(iParticle),
			   airDragCoefficientFloe->get(iParticle),
			   airDragCoefficientPond->get(iParticle),
			   airDragCoefficientRidge->get(iParticle),
			   dragFreeboard->get(iParticle),
			   dragIceSnowDraft->get(iParticle),
			   dragRidgeHeight->get(iParticle),
			   dragRidgeSeparation->get(iParticle),
			   dragKeelDepth->get(iParticle),
			   dragKeelSeparation->get(iParticle),
			   dragFloeLength->get(iParticle),
			   dragFloeSeparation->get(iParticle),
			   airStressForcingU->get(iParticle),
			   airStressForcingV->get(iParticle),
			   airStressCellU->get(iParticle),
			   airStressCellV->get(iParticle),
			   airPotentialTemperature->get(iParticle),
			   seaSurfaceTemperature->get(iParticle),
			   seaSurfaceSalinity->get(iParticle),
			   seaFreezingTemperature->get(iParticle),
			   oceanStressCellU->get(iParticle),
			   oceanStressCellV->get(iParticle),
			   oceanHeatFluxIceBottom->get(iParticle),
			   freezingMeltingPotential->get(iParticle),
			   lateralIceMeltFraction->get(iParticle),
			   snowfallRate->get(iParticle),
			   rainfallRate->get(iParticle),
			   pondFreshWaterFlux->get(iParticle),
			   surfaceHeatFlux->get(iParticle),
			   surfaceHeatFluxCategory->get(iParticle),
			   surfaceConductiveFlux->get(iParticle),
			   surfaceConductiveFluxCategory->get(iParticle),
			   surfaceShortwaveFlux->get(iParticle),
			   interiorShortwaveFlux->get(iParticle),
			   penetratingShortwaveFlux->get(iParticle),
			   absorbedShortwaveFlux->get(iParticle),
			   longwaveUp->get(iParticle),
			   absorbedShortwaveSnowLayer->get(iParticle),
			   absorbedShortwaveIceLayer->get(iParticle),
			   longwaveDown->get(iParticle),
			   sensibleHeatFlux->get(iParticle),
			   sensibleHeatFluxCategory->get(iParticle),
			   latentHeatFlux->get(iParticle),
			   latentHeatFluxCategory->get(iParticle),
			   evaporativeWaterFlux->get(iParticle),
			   oceanFreshWaterFlux->get(iParticle),
			   oceanSaltFlux->get(iParticle),
			   oceanHeatFlux->get(iParticle),
			   oceanShortwaveFlux->get(iParticle),
			   latentHeatFluxCouple->get(iParticle),
			   sensibleHeatFluxCouple->get(iParticle),
			   surfaceHeatFluxCouple->get(iParticle),
			   surfaceConductiveFluxCouple->get(iParticle),
			   atmosAerosolFlux->get(iParticle),
			   oceanAerosolFlux->get(iParticle),
			   pondSnowDepthDifference->get(iParticle),
			   pondLidMeltFluxFraction->get(iParticle),
			   surfaceIceMelt->get(iParticle),
			   surfaceIceMeltCategory->get(iParticle),
			   basalIceMelt->get(iParticle),
			   basalIceMeltCategory->get(iParticle),
			   snowMelt->get(iParticle),
			   snowMeltCategory->get(iParticle),
			   congelation->get(iParticle),
			   congelationCategory->get(iParticle),
			   snowiceFormation->get(iParticle),
			   snowiceFormationCategory->get(iParticle),
			   snowThicknessChangeCategory->get(iParticle),
			   meltOnset->get(iParticle),
			   freezeOnset->get(iParticle),
			   snowScatteringAerosol->get(iParticle),
			   snowBodyAerosol->get(iParticle),
			   iceScatteringAerosol->get(iParticle),
			   iceBodyAerosol->get(iParticle),
			   &dayOfYear);

	check_warnings("column_step_therm1", iParticle);

      }

    } // iParticle

  } // useColumnVerticalThermodynamics

}

// Time stepping of the icepack itd thermodyanmic subroutine
void Column::run_step_therm2(void) {

  if (useColumnItdThermodynamics) {

    int nCategories   = columnDimensions->size("nCategories");
    int nCategoriesP1 = columnDimensions->size("nCategoriesP1");
    int nIceLayers    = columnDimensions->size("nIceLayers");
    int nIceLayersP1  = columnDimensions->size("nIceLayersP1");
    int nSnowLayers   = columnDimensions->size("nSnowLayers");
    int nAerosols     = columnDimensions->size("nAerosols");

    double thermoTimeStep = timeStepInterval->seconds();

    DEMSI::Time currentTime = simulationClock->get_time();

    double dayOfYear = (double) currentTime.day_of_year();

    bool updateOceanFluxes;
    configs->get({"ConfigGroup:ocean","Config:updateOceanFluxes"}, updateOceanFluxes);
    int updateOceanFluxesIn = (int) updateOceanFluxes;

    for (int iParticle = 0 ; iParticle < *(particles->nParticles) ; iParticle++) {

      if (particles->type(iParticle) == 1) {

	columnTracers->set_icepack_tracer_array_category(iParticle);

	column_step_therm2(&nCategories,
			   &nCategoriesP1,
			   &columnTracers->nTracers,
			   &columnTracers->nBaseTracers,
			   &columnTracers->nMaxAncestorTracers,
			   &nAerosols,
			   &nIceLayers,
			   &nIceLayersP1,
			   &nSnowLayers,
			   &thermoTimeStep,
			   categoryThicknessLimits,
			   iceAreaCategory->get(iParticle),
			   iceVolumeCategory->get(iParticle),
			   snowVolumeCategory->get(iParticle),
			   iceAreaCategoryInitial->get(iParticle),
			   iceVolumeCategoryInitial->get(iParticle),
			   columnTracers->tracerArrayCategory,
			   openWaterArea->get(iParticle),
			   iceAreaCell->get(iParticle),
			   columnTracers->parentIndex,
			   columnTracers->firstAncestorMask,
			   columnTracers->ancestorNumber,
			   columnTracers->ancestorIndices,
			   seaFreezingTemperature->get(iParticle),
			   seaSurfaceSalinity->get(iParticle),
			   initialSalinityProfile->get(iParticle),
			   lateralIceMeltFraction->get(iParticle),
			   lateralIceMelt->get(iParticle),
			   freezingMeltingPotential->get(iParticle),
			   frazilFormation->get(iParticle),
			   rainfallRate->get(iParticle),
			   pondFreshWaterFlux->get(iParticle),
			   oceanFreshWaterFlux->get(iParticle),
			   oceanSaltFlux->get(iParticle),
			   oceanHeatFlux->get(iParticle),
			   &updateOceanFluxesIn,
			   oceanAerosolFlux->get(iParticle),
			   frazilGrowthDiagnostic->get(iParticle),
			   freezeOnset->get(iParticle),
			   &dayOfYear);

	check_warnings("column_step_therm2", iParticle);

	columnTracers->get_icepack_tracer_array_category(iParticle);

      }

    } // iParticle

  } // useColumnVerticalThermodynamics

}

// Time stepping of the icepack ridging subroutine
void Column::run_step_ridge(void) {

  if (useColumnRidging) {

    int nCategories   = columnDimensions->size("nCategories");
    int nCategoriesP1 = columnDimensions->size("nCategoriesP1");
    int nIceLayers    = columnDimensions->size("nIceLayers");
    int nSnowLayers   = columnDimensions->size("nSnowLayers");
    int nAerosols     = columnDimensions->size("nAerosols");

    double thermoTimeStep = timeStepInterval->seconds();

    int nDynamicsSubcycles;
    configs->get({"ConfigGroup:simulationTiming","Config:nDynamicsSubcycles"}, nDynamicsSubcycles);

    double dynamicsTimeStep = thermoTimeStep / double(nDynamicsSubcycles);

    auto changeEffectiveElementArea_particles = particles->changeEffectiveElementArea;

    for (int iParticle = 0 ; iParticle < *(particles->nParticles) ; iParticle++) {

      if (particles->type(iParticle) == 1) {

	double changeEffectiveElementArea = changeEffectiveElementArea_particles(iParticle);

	// ratio of effective element areas before and after ridging
	double elementAreaRatio = ((*effectiveElementArea)(iParticle) + changeEffectiveElementArea) / (*effectiveElementArea)(iParticle);

	// icepack convergence
	(*ridgeConvergence)(iParticle) = (1.0 / elementAreaRatio - 1.0) / dynamicsTimeStep;
	(*ridgeShear)(iParticle) = 0.0;

	// modify area and thickness in ITD from element size change
	for (int iCategory = 0 ; iCategory < nCategories ; iCategory++) {
	  (*iceAreaCategory  )(iParticle,iCategory) = (*iceAreaCategory  )(iParticle,iCategory) / elementAreaRatio;
	  (*iceVolumeCategory)(iParticle,iCategory) = (*iceVolumeCategory)(iParticle,iCategory) / elementAreaRatio;
	} // iCategory

	// change the effective element area due to ridging
	(*effectiveElementArea)(iParticle) = (*effectiveElementArea)(iParticle) + changeEffectiveElementArea;

	columnTracers->set_icepack_tracer_array_category(iParticle);

	column_step_ridge(&nCategories,
			  &nCategoriesP1,
			  &columnTracers->nTracers,
			  &columnTracers->nBaseTracers,
			  &columnTracers->nMaxAncestorTracers,
			  &nAerosols,
			  &nIceLayers,
			  &nSnowLayers,
			  &nDynamicsSubcycles,
			  &dynamicsTimeStep,
			  categoryThicknessLimits,
			  ridgeConvergence->get(iParticle),
			  ridgeShear->get(iParticle),
			  iceAreaCategory->get(iParticle),
			  iceVolumeCategory->get(iParticle),
			  snowVolumeCategory->get(iParticle),
			  columnTracers->tracerArrayCategory,
			  openWaterArea->get(iParticle),
			  columnTracers->parentIndex,
			  columnTracers->firstAncestorMask,
			  columnTracers->ancestorNumber,
			  columnTracers->ancestorIndices,
			  areaLossRidge->get(iParticle),
			  areaGainRidge->get(iParticle),
			  iceVolumeRidged->get(iParticle),
			  openingRateRidge->get(iParticle),
			  pondFreshWaterFlux->get(iParticle),
			  oceanFreshWaterFlux->get(iParticle),
			  oceanHeatFlux->get(iParticle),
			  oceanAerosolFlux->get(iParticle),
			  iceAreaCell->get(iParticle),
			  oceanSaltFlux->get(iParticle));

	check_warnings("column_step_ridge", iParticle);

	columnTracers->get_icepack_tracer_array_category(iParticle);

      }

    } // iParticle

  } // useColumnVerticalThermodynamics

}

// coupling prep
void Column::coupling_prep(void) {

  int nCategories = columnDimensions->size("nCategories");

  bool includePondFreshwaterFeedback = false;
  int includePondFreshwaterFeedbackIn = (int) includePondFreshwaterFeedback;

  double thermoTimeStep = timeStepInterval->seconds();

  for (int iParticle = 0 ; iParticle < *(particles->nParticles) ; iParticle++) {

    if (particles->type(iParticle) == 1) {

      column_coupling_prep(&nCategories,
			   &includePondFreshwaterFeedbackIn,
			   &thermoTimeStep,
			   freezingMeltingPotentialInitial->get(iParticle),
			   freezingMeltingPotential->get(iParticle),
			   albedoVisibleDirectCell->get(iParticle),
			   albedoVisibleDiffuseCell->get(iParticle),
			   albedoIRDirectCell->get(iParticle),
			   albedoIRDiffuseCell->get(iParticle),
			   bareIceAlbedoCell->get(iParticle),
			   snowAlbedoCell->get(iParticle),
			   pondAlbedoCell->get(iParticle),
			   effectivePondAreaCell->get(iParticle),
			   solarZenithAngleCosine->get(iParticle),
			   albedoVisibleDirectArea->get(iParticle),
			   albedoVisibleDiffuseArea->get(iParticle),
			   albedoIRDirectArea->get(iParticle),
			   albedoIRDiffuseArea->get(iParticle),
			   oceanFreshWaterFluxArea->get(iParticle),
			   oceanSaltFluxArea->get(iParticle),
			   oceanHeatFluxArea->get(iParticle),
			   oceanShortwaveFluxArea->get(iParticle),
			   oceanSaltFlux->get(iParticle),
			   oceanHeatFlux->get(iParticle),
			   oceanShortwaveFlux->get(iParticle),
			   shortwaveScalingFactor->get(iParticle),
			   pondFreshWaterFlux->get(iParticle),
			   oceanFreshWaterFlux->get(iParticle),
			   shortwaveVisibleDirectDown->get(iParticle),
			   shortwaveVisibleDiffuseDown->get(iParticle),
			   shortwaveIRDirectDown->get(iParticle),
			   shortwaveIRDiffuseDown->get(iParticle),
			   iceAreaCategory->get(iParticle),
			   albedoVisibleDirectCategory->get(iParticle),
			   albedoVisibleDiffuseCategory->get(iParticle),
			   albedoIRDirectCategory->get(iParticle),
			   albedoIRDiffuseCategory->get(iParticle),
			   bareIceAlbedoCategory->get(iParticle),
			   snowAlbedoCategory->get(iParticle),
			   pondAlbedoCategory->get(iParticle),
			   effectivePondAreaCategory->get(iParticle));

    }

  } // iParticle

}

// Ocean mixed layer
void Column::ocean_mixed_layer(void) {

  double thermoTimeStep = timeStepInterval->seconds();

  for (int iParticle = 0 ; iParticle < *(particles->nParticles) ; iParticle++) {

    if (particles->type(iParticle) == 1) {

      double airStressOceanU = 0.0;
      double airStressOceanV = 0.0;

      column_ocean_mixed_layer(&thermoTimeStep,
			       seaSurfaceTemperature->get(iParticle),
			       airPotentialTemperature->get(iParticle),
			       uAirVelocity->get(iParticle),
			       vAirVelocity->get(iParticle),
			       windSpeed->get(iParticle),
			       airLevelHeight->get(iParticle),
			       airSpecificHumidity->get(iParticle),
			       airDensity->get(iParticle),
			       &airStressOceanU,
			       &airStressOceanV,
			       atmosReferenceTemperature2mOcean->get(iParticle),
			       atmosReferenceHumidity2mOcean->get(iParticle),
			       airDragCoefficient->get(iParticle),
			       airOceanDragCoefficientRatio->get(iParticle),
			       albedoVisibleDirectOcean->get(iParticle),
			       shortwaveVisibleDirectDown->get(iParticle),
			       albedoIRDirectOcean->get(iParticle),
			       shortwaveIRDirectDown->get(iParticle),
			       albedoVisibleDiffuseOcean->get(iParticle),
			       shortwaveVisibleDiffuseDown->get(iParticle),
			       albedoIRDiffuseOcean->get(iParticle),
			       shortwaveIRDiffuseDown->get(iParticle),
			       longwaveUpOcean->get(iParticle),
			       sensibleHeatFluxOcean->get(iParticle),
			       sensibleTransferCoefficient->get(iParticle),
			       latentHeatFluxOcean->get(iParticle),
			       evaporativeWaterFluxOcean->get(iParticle),
			       longwaveDown->get(iParticle),
			       iceAreaCell->get(iParticle),
			       oceanHeatFlux->get(iParticle),
			       oceanShortwaveFlux->get(iParticle),
			       oceanMixedLayerDepth->get(iParticle),
			       seaFreezingTemperature->get(iParticle),
			       oceanHeatFluxConvergence->get(iParticle),
			       freezingMeltingPotential->get(iParticle));

    }

  } // iParticle

}

// reinitialize fluxes at end of timestep
void Column::finalize_timestep(void) {

  if (useColumn) {

    // reinitialize atmospheric fluxes
    reset_atmospheric_coupler_fluxes();

    // reinitialize oceanic fluxes
    if (oceanType == "embedded") {
      reset_ocean_coupler_fluxes();
    }

  } // useColumn

}

//------------------------------------------------------------------------------
// Column parameters
//------------------------------------------------------------------------------

// set column configs in Icepack
void Column::column_configs(void) {

  set_column_parameter("ocean", "minFrictionVelocity");
  set_column_parameter("shortwave", "visibleIceAlbedo");
  set_column_parameter("shortwave", "infraredIceAlbedo");
  set_column_parameter("shortwave", "visibleSnowAlbedo");
  set_column_parameter("shortwave", "infraredSnowAlbedo");
  set_column_parameter("atmosphere", "boundaryLayerIterationNumber");
  set_column_parameter("shortwave", "variableAlbedoThicknessLimit");
  set_column_parameter("shortwave", "shortwaveType");
  set_column_parameter("shortwave", "albedoType");
  set_column_parameter("shortwave", "iceShortwaveTuningParameter");
  set_column_parameter("shortwave", "pondShortwaveTuningParameter");
  set_column_parameter("shortwave", "snowShortwaveTuningParameter");
  set_column_parameter("shortwave", "tempChangeSnowGrainRadiusChange");
  set_column_parameter("shortwave", "maxMeltingSnowGrainRadius");
  set_column_parameter("ridging", "iceStrengthFormulation");
  set_column_parameter("ridging", "ridgingParticipationFunction");
  set_column_parameter("ridging", "ridgingRedistributionFunction");
  set_column_parameter("ridging", "ridigingEfoldingScale");
  set_column_parameter("atmosphere", "atmosBoundaryMethod");
  set_column_parameter("atmosphere", "calcSurfaceStresses");
  set_column_parameter("atmosphere", "useFormDrag");
  set_column_parameter("atmosphere", "useHighFrequencyCoupling");
  set_column_parameter("itd", "itdConversionType");
  set_column_parameter("itd", "categoryBoundsType");
  set_column_parameter("meltponds", "snowToIceTransitionDepth");
  set_column_parameter("meltponds", "pondFlushingTimescale");
  set_column_parameter("meltponds", "pondRefreezingType");
  set_column_parameter("meltponds", "minMeltwaterRetainedFraction");
  set_column_parameter("meltponds", "maxMeltwaterRetainedFraction");
  set_column_parameter("meltponds", "pondDepthToFractionRatio");
  set_column_parameter("meltponds", "snowOnPondIceTaperingParameter");
  set_column_parameter("meltponds", "criticalPondIceThickness");
  set_column_parameter("thermodynamics", "thermodynamicsType");
  set_column_parameter("atmosphere", "calcSurfaceTemperature");
  set_column_parameter("thermodynamics", "heatConductivityType");
  set_column_parameter("thermodynamics", "rapidModeChannelRadius");
  set_column_parameter("thermodynamics", "rapidModelCriticalRa");
  set_column_parameter("thermodynamics", "rapidModeAspectRatio");
  set_column_parameter("thermodynamics", "slowModeDrainageStrength");
  set_column_parameter("thermodynamics", "slowModeCriticalPorosity");
  set_column_parameter("thermodynamics", "congelationIcePorosity");
  set_column_parameter("ocean", "seaFreezingTemperatureType");
  set_column_parameter("shortwave", "algaeAbsorptionCoefficient");
  set_column_parameter("ocean", "oceanHeatTransferType");

}

// replace DEMSI string config values for Icepack integer ones, return true if config needs conversion
bool replace_icepack_int_values(const std::string configGroupName, const std::string configName, const std::string configValue, int &configValueIcepack) {

  configValueIcepack = -999;

  // ktherm
  if (configGroupName == "thermodynamics" and configName == "thermodynamicsType") {

    if (configValue == "zero layer") {
      configValueIcepack = 0;
    } else if (configValue == "BL99") {
      configValueIcepack = 1;
    } else if (configValue == "mushy") {
      configValueIcepack = 2;
    }
    return true;

  // kitd
  } else if (configGroupName == "itd" and configName == "itdConversionType") {

    if (configValue == "delta function") {
      configValueIcepack = 0;
    } else if ("linear remap") {
      configValueIcepack = 1;
    }
    return true;

  // kcatbound
  } else if (configGroupName == "itd" and configName == "categoryBoundsType") {

    if (configValue == "single category") {
      configValueIcepack = -1;
    } else if (configValue == "original") {
      configValueIcepack = 0;
    } else if (configValue == "new") {
      configValueIcepack = 1;
    } else if (configValue == "WMO") {
      configValueIcepack = 2;
    } else if (configValue == "asymptotic") {
      configValueIcepack = 3;
    }
    return true;

  // kstrength
  } else if (configGroupName == "ridging" and configName == "iceStrengthFormulation") {

    if (configValue == "Hibler79") {
      configValueIcepack = 0;
    } else if (configValue == "Rothrock75") {
      configValueIcepack = 1;
    }
    return true;

  // krdg_partic
  } else if (configGroupName == "ridging" and configName == "ridgingParticipationFunction") {

    if (configValue == "Thorndike75") {
      configValueIcepack = 0;
    } else if (configValue == "exponential") {
      configValueIcepack = 1;
    }
    return true;

  // krdg_redist
  } else if (configGroupName == "ridging" and configName == "ridgingRedistributionFunction") {

    if (configValue == "Hibler80") {
      configValueIcepack = 0;
    } else if (configValue == "exponential") {
      configValueIcepack = 1;
    }
    return true;

  }

  return false;

}

// set an individual column parameter in icepack
void Column::set_column_parameter(const std::string configGroupName, const std::string configName) {

  std::string icePackConfigName;
  configs->get_attribute({"ConfigGroup:"+configGroupName,"Config:"+configName}, "icepack_name", icePackConfigName);
  char cIcePackConfigName[1024];
  convert_string_to_array(icePackConfigName, cIcePackConfigName);

  DEMSI::CONFIGTYPE configType = configs->get_config_type({"ConfigGroup:"+configGroupName,"Config:"+configName});
  if (configType == DEMSI::CONFIGTYPE::INTEGER) {

    int paramValue;
    configs->get({"ConfigGroup:"+configGroupName,"Config:"+configName}, paramValue);
    column_init_parameters_integer(cIcePackConfigName, &paramValue);

  } else if (configType == DEMSI::CONFIGTYPE::DOUBLE) {

    double paramValue;
    configs->get({"ConfigGroup:"+configGroupName,"Config:"+configName}, paramValue);
    column_init_parameters_real(cIcePackConfigName, &paramValue);

  } else if (configType == DEMSI::CONFIGTYPE::STRING) {

    std::string paramValue;
    configs->get({"ConfigGroup:"+configGroupName,"Config:"+configName}, paramValue);

    // check if we need to convert from a DEMSI string config to Icepack integer one
    int paramValueIcepack;
    if (replace_icepack_int_values(configGroupName, configName, paramValue, paramValueIcepack)) {
      column_init_parameters_integer(cIcePackConfigName, &paramValueIcepack);
    } else {
      char cParamValue[1024];
      convert_string_to_array(paramValue, cParamValue);
      column_init_parameters_character(cIcePackConfigName, cParamValue);
    }

  } else if (configType == DEMSI::CONFIGTYPE::BOOLEAN) {

    bool paramValue;
    configs->get({"ConfigGroup:"+configGroupName,"Config:"+configName}, paramValue);
    int paramValueInt = (int) paramValue;
    column_init_parameters_logical(cIcePackConfigName, &paramValueInt);

  }

}

// Initialize the vertical thermo profiles: init_thermo_vertical
void Column::init_thermo_vertical(void) {

  int nIceLayersP1 = columnDimensions->size("nIceLayersP1");
  double salinityProfileTemplate[nIceLayersP1]; // vertical salinity profile

  int nIceLayers = columnDimensions->size("nIceLayers");
  column_init_thermo(&nIceLayers, &nIceLayersP1, &salinityProfileTemplate[0]);
  check_warnings("column_init_thermo");

  for (int iParticle = 0 ; iParticle < *(particles->nParticles) ; iParticle++) {

    if (particles->type(iParticle) == 1) {

      column_init_profiles(&nIceLayersP1,
			   salinityProfileTemplate,
			   initialSalinityProfile->get(iParticle),
			   initialMeltingTemperatureProfile->get(iParticle));

    }

  } // iParticle

}

void Column::init_ice_thickness_distribution(void) {

  int nCategories   = columnDimensions->size("nCategories");
  int nCategoriesP1 = columnDimensions->size("nCategoriesP1");

  categoryThicknessLimits = new double[nCategoriesP1];

  if (configs->exists({"ConfigGroup:initialization","Array:categoryThicknessLimits"})) {

    // config defines thickness categories
    std::vector<double> categoryThicknessLimitsIn = configs->get_double_array({"ConfigGroup:initialization","Array:categoryThicknessLimits"}, "categoryThicknessLimit");
    log->check(categoryThicknessLimitsIn.size() == nCategories-1, "Wrong number of input category thickness limits from configs.");

    categoryThicknessLimits[0] = 0.0;
    for (int iCategory = 1 ; iCategory < nCategories ; iCategory++) {
      categoryThicknessLimits[iCategory] = categoryThicknessLimitsIn[iCategory-1];
    } // iCategory
    categoryThicknessLimits[nCategories] = 1.0e8;

  } else {

    // Icepack set up standard thickness categories
    column_init_itd(&nCategories, categoryThicknessLimits);
    check_warnings("init_ice_thickness_distribution");

    // set upper bound to large number
    categoryThicknessLimits[nCategories] = 999.9;

  }

  for (int iCategory = 0 ; iCategory < nCategoriesP1 ; iCategory++) {
    (*log)(DEMSI::LOG::DEBUG) << iCategory << " " << categoryThicknessLimits[iCategory] << std::endl;
  } // iCategory

}

void Column::init_ice_thickness_distribution_no_column(void) {

  int nCategories = columnDimensions->size("nCategories");
  log->check(nCategories == 1, "nCategories must be 1 for no column physics.");

  categoryThicknessLimits = new double[2];

  categoryThicknessLimits[0] = 0.0;
  categoryThicknessLimits[1] = 999.9;

} // Column::init_ice_thickness_distribution_no_column

void Column::init_column_state(void) {

  // init state variables
  int nCategories   = columnDimensions->size("nCategories");
  int nIceLayers  = columnDimensions->size("nIceLayers");
  int nSnowLayers = columnDimensions->size("nSnowLayers");

  for (int iParticle = 0 ; iParticle < *(particles->nParticles) ; iParticle++) {

    if (particles->type(iParticle) == 1) {

      column_init_particle_state(&nCategories,
				 &nIceLayers,
				 &nSnowLayers,
				 surfaceTemperature->get(iParticle),
				 airTemperature->get(iParticle),
				 seaFreezingTemperature->get(iParticle),
				 initialSalinityProfile->get(iParticle),
				 initialMeltingTemperatureProfile->get(iParticle),
				 iceEnthalpy->get(iParticle),
				 iceSalinity->get(iParticle),
				 snowEnthalpy->get(iParticle));

    } // particle type

  } // iParticle

  // aggregate tracers
  aggregate();

}

void Column::aggregate(void) {

  int nCategories = columnDimensions->size("nCategories");

  for (int iParticle = 0 ; iParticle < *(particles->nParticles) ; iParticle++) {

    if (particles->type(iParticle) == 1) {

      // set the category tracer array
      columnTracers->set_icepack_tracer_array_category(iParticle);

      column_aggregate(&nCategories,
		       &columnTracers->nTracers,
		       &columnTracers->nBaseTracers,
		       &columnTracers->nMaxAncestorTracers,
		       columnTracers->parentIndex, // trcr_depend
		       columnTracers->firstAncestorMask, // trcr_base
		       columnTracers->ancestorNumber, // n_trcr_strata
		       columnTracers->ancestorIndices,
		       columnTracers->tracerArrayCategory, // trcrn
		       columnTracers->tracerArray, // trcr
		       iceAreaCategory->get(iParticle),
		       iceVolumeCategory->get(iParticle),
		       snowVolumeCategory->get(iParticle),
		       iceAreaCell->get(iParticle),
		       iceVolumeCell->get(iParticle),
		       snowVolumeCell->get(iParticle),
		       openWaterArea->get(iParticle));

      // set the cell tracer array
      columnTracers->get_icepack_tracer_array(iParticle);

    }

  } // iParticle

}

// initialize the atmospheric coupler fields
void Column::initialize_atmos_coupler_fields(void) {

  int useAerosolsIn = (int) columnTracers->useAerosols;

  for (int iParticle = 0 ; iParticle < *(particles->nParticles) ; iParticle++) {

    if (particles->type(iParticle) == 1) {

      column_initialize_atmos_coupler_fields(&useAerosolsIn,
					     windSpeed->get(iParticle),
					     uAirVelocity->get(iParticle),
					     vAirVelocity->get(iParticle),
					     longwaveUp->get(iParticle),
					     airDragCoefficient->get(iParticle),
					     atmosAerosolFlux->get(iParticle));

    }

  } // iParticle

}

// initialize the oceanic coupler fields
void Column::initialize_ocean_coupler_fields(void) {

  bool doRestart;
  configs->get({"ConfigGroup:simulationTiming","Config:doRestart"}, doRestart);
  int doRestartIn = (int) doRestart;

  for (int iParticle = 0 ; iParticle < *(particles->nParticles) ; iParticle++) {

    if (particles->type(iParticle) == 1) {

      column_initialize_ocean_coupler_fields(&doRestartIn,
					     seaSurfaceTemperature->get(iParticle),
					     seaFreezingTemperature->get(iParticle),
					     seaSurfaceSalinity->get(iParticle));

    }

  } // iParticle

}

// reset atmospheric coupler fluxes
void Column::reset_atmospheric_coupler_fluxes(void) {

  sensibleHeatFlux->set(0.0);
  latentHeatFlux->set(0.0);
  longwaveUp->set(0.0);
  evaporativeWaterFlux->set(0.0);

  airStressCellU->set(0.0);
  airStressCellV->set(0.0);

  absorbedShortwaveFlux->set(0.0);

  atmosReferenceSpeed10m->set(0.0);
  atmosReferenceTemperature2m->set(0.0);
  atmosReferenceHumidity2m->set(0.0);

}

// reset oceanic coupler fluxes
void Column::reset_ocean_coupler_fluxes(void) {

  oceanFreshWaterFlux->set(0.0);
  oceanSaltFlux->set(0.0);
  oceanHeatFlux->set(0.0);
  oceanShortwaveFlux->set(0.0);

}

// reinitialize diagnostic variables
void Column::reinitialize_diagnostics(void) {

  // forcing
  airDensity->set(1.3);
  airLevelHeight->set(10.0);

  // atmospheric fluxes
  surfaceHeatFlux->set(0.0);
  surfaceConductiveFlux->set(0.0);
  surfaceHeatFluxCategory->set(0.0);
  surfaceConductiveFluxCategory->set(0.0);
  latentHeatFluxCategory->set(0.0);
  sensibleHeatFluxCategory->set(0.0);

  // melt growth rates
  congelation->set(0.0);
  frazilFormation->set(0.0);
  snowiceFormation->set(0.0);
  snowThicknessChange->set(0.0);
  surfaceIceMelt->set(0.0);
  snowMelt->set(0.0);
  basalIceMelt->set(0.0);
  lateralIceMelt->set(0.0);

  // thermodynamic tendencies
  iceAreaTendencyThermodynamics->set(*iceAreaCell);
  iceVolumeTendencyThermodynamics->set(*iceVolumeCell);
  if (columnTracers->useIceAge) {
    iceAgeTendencyThermodynamics->set(*iceAgeCell);
  } else {
    iceAgeTendencyThermodynamics->set(0.0);
  }

  // ponds
  pondFreshWaterFlux->set(0.0);

  // shortwave
  bareIceAlbedoCell->set(0.0);
  snowAlbedoCell->set(0.0);
  pondAlbedoCell->set(0.0);

  // form drag
  char paramName[7] = "dragio";
  double iceOceanDragCoefficient = column_get_constant(paramName);
  oceanDragCoefficient->set(iceOceanDragCoefficient);
  airDragCoefficient->set(column_initial_air_drag_coefficient());

  bool useFormDrag;
  configs->get({"ConfigGroup:atmosphere","Config:useFormDrag"}, useFormDrag);
  if (useFormDrag) {

    airOceanDragCoefficientRatio->set(0.0);
    oceanDragCoefficientSkin->set(0.0);
    oceanDragCoefficientFloe->set(0.0);
    oceanDragCoefficientKeel->set(0.0);
    airDragCoefficientSkin->set(0.0);
    airDragCoefficientFloe->set(0.0);
    airDragCoefficientPond->set(0.0);
    airDragCoefficientRidge->set(0.0);
    dragFreeboard->set(0.0);
    dragIceSnowDraft->set(0.0);
    dragRidgeHeight->set(0.0);
    dragRidgeSeparation->set(0.0);
    dragKeelDepth->set(0.0);
    dragKeelSeparation->set(0.0);
    dragFloeLength->set(0.0);
    dragFloeSeparation->set(0.0);

  } // useFormDrag

}

//------------------------------------------------------------------------------
// Forcing
//------------------------------------------------------------------------------

// limit air temperature over ice
void Column::limit_temperature(void) {

  bool useLimitTemperature;
  configs->get({"ConfigGroup:forcing","Config:useLimitTemperature"}, useLimitTemperature);
  if (useLimitTemperature) {

    for (int iParticle = 0 ; iParticle < *(particles->nParticles) ; iParticle++) {

      if (particles->type(iParticle) == 1) {

	column_limit_temperature(iceAreaCell->get(iParticle),
				 airTemperature->get(iParticle));

      }

    } // iParticle

  } // useLimitTemperature

}

// limit specific humidity
void Column::limit_specific_humidity(void) {

  bool useLimitHumidity;
  configs->get({"ConfigGroup:forcing","Config:useLimitHumidity"}, useLimitHumidity);
  if (useLimitHumidity) {

    for (int iParticle = 0 ; iParticle < *(particles->nParticles) ; iParticle++) {

      if (particles->type(iParticle) == 1) {

	column_limit_specific_humidity(airTemperature->get(iParticle),
				       airSpecificHumidity->get(iParticle));

      }

    } // iParticle

  } // useLimitHumidity

}

// get shortwave radiation from cloudiness and location/time
void Column::shortwave_from_cloud_fraction(void) {

  bool useGetShortwaveFromCloudiness;
  configs->get({"ConfigGroup:forcing","Config:useGetShortwaveFromCloudiness"}, useGetShortwaveFromCloudiness);
  if (useGetShortwaveFromCloudiness) {

    DEMSI::Time currentTime = simulationClock->get_time();

    int dayOfYear = currentTime.day_of_year();
    double secondsIntoDay = (double) currentTime.seconds_into_day();

    for (int iParticle = 0 ; iParticle < *(particles->nParticles) ; iParticle++) {

      if (particles->type(iParticle) == 1) {

	column_shortwave_from_cloud_fraction(shortwaveDown->get(iParticle),
					     longitude->get(iParticle),
					     latitude->get(iParticle),
					     cloudFraction->get(iParticle),
					     airSpecificHumidity->get(iParticle),
					     &secondsIntoDay,
					     &dayOfYear);

      }

    } // iParticle

  } // useGetShortwaveFromCloudiness

}

// split shortwave between bands
void Column::split_shortwave(void) {

  bool useSplitShortwave;
  configs->get({"ConfigGroup:forcing","Config:useSplitShortwave"}, useSplitShortwave);
  if (useSplitShortwave) {

    double fracShortwaveVisibleDirect  = 0.28; // fraction of incoming shortwave in visible direct band
    double fracShortwaveVisibleDiffuse = 0.24; // fraction of incoming shortwave in visible diffuse band
    double fracShortwaveIRDirectDown   = 0.31; // fraction of incoming shortwave in near IR direct band
    double fracShortwaveIRDiffuseDown  = 0.17; // fraction of incoming shortwave in near IR diffuse band

    for (int iParticle = 0 ; iParticle < *(particles->nParticles) ; iParticle++) {

      if (particles->type(iParticle) == 1) {

	(*shortwaveVisibleDirectDown)(iParticle)  = (*shortwaveDown)(iParticle) * fracShortwaveVisibleDirect;
	(*shortwaveVisibleDiffuseDown)(iParticle) = (*shortwaveDown)(iParticle) * fracShortwaveVisibleDiffuse;
	(*shortwaveIRDirectDown)(iParticle)       = (*shortwaveDown)(iParticle) * fracShortwaveIRDirectDown;
	(*shortwaveIRDiffuseDown)(iParticle)      = (*shortwaveDown)(iParticle) * fracShortwaveIRDiffuseDown;

      }

    } // iParticle

  } // useSplitShortwave

}

// limit forcing variables to physically realistic values
void Column::apply_physical_limits(void) {

  bool useApplyPhysicalLimits;
  configs->get({"ConfigGroup:forcing","Config:useApplyPhysicalLimits"}, useApplyPhysicalLimits);
  if (useApplyPhysicalLimits) {

    for (int iParticle = 0 ; iParticle < *(particles->nParticles) ; iParticle++) {

      if (particles->type(iParticle) == 1) {

	if (cloudFraction       != NULL) (*cloudFraction)(iParticle)       = std::max(std::min((*cloudFraction)(iParticle),1.0),0.0);
	if (shortwaveDown       != NULL) (*shortwaveDown)(iParticle)       = std::max((*shortwaveDown)(iParticle),0.0);
	if (rainfallRate        != NULL) (*rainfallRate)(iParticle)        = std::max((*rainfallRate)(iParticle),0.0);
	if (airSpecificHumidity != NULL) (*airSpecificHumidity)(iParticle) = std::max((*airSpecificHumidity)(iParticle),0.0);

      }

    } // iParticle

  } // useApplyPhysicalLimits

}

// calculate the downwelling longwave radiation
void Column::calculate_longwave_down(void) {

  std::string longwaveType;
  configs->get({"ConfigGroup:forcing","Config:longwaveType"}, longwaveType);
  if (longwaveType == "rosati_miyakoda") {
    longwave_rosati_miyakoda();
  } else if (longwaveType == "parkinson_and_washington") {
    longwave_parkinson_and_washington();
  }

}

// calculate the downwelling longwave radiation from Rosati and Miyakoda
void Column::longwave_rosati_miyakoda(void) {

  for (int iParticle = 0 ; iParticle < *(particles->nParticles) ; iParticle++) {

    if (particles->type(iParticle) == 1) {

      column_longwave_rosati_miyakoda(longwaveDown->get(iParticle),
				      cloudFraction->get(iParticle),
				      iceAreaCell->get(iParticle),
				      surfaceTemperature->get(iParticle),
				      seaSurfaceTemperature->get(iParticle),
				      airSpecificHumidity->get(iParticle),
				      airTemperature->get(iParticle));

    }

  } // iParticle

}

// calculate the downwelling longwave radiation from Parkinson and Washington
void Column::longwave_parkinson_and_washington(void) {

  for (int iParticle = 0 ; iParticle < *(particles->nParticles) ; iParticle++) {

    if (particles->type(iParticle) == 1) {

      column_longwave_parkinson_and_washington(longwaveDown->get(iParticle),
					       airTemperature->get(iParticle),
					       cloudFraction->get(iParticle));

    }

  } // iParticle

}

// scale precipitation according to units
void Column::scale_precipitation(void) {

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

  (*precipitationRate) *= precipitationFactor;

}

// split precipitation between snow and rain by air temperature
void Column::split_precipitation(void) {

  bool useSplitPrecipitation;
  configs->get({"ConfigGroup:forcing","Config:useSplitPrecipitation"}, useSplitPrecipitation);
  if (useSplitPrecipitation) {

    for (int iParticle = 0 ; iParticle < *(particles->nParticles) ; iParticle++) {

      if (particles->type(iParticle) == 1) {

	column_split_precipitation(airTemperature->get(iParticle),
				   precipitationRate->get(iParticle),
				   snowfallRate->get(iParticle),
				   rainfallRate->get(iParticle));

      }

    } // iParticle

  } // useSplitPrecipitation

}

//------------------------------------------------------------------------------
// Connection with dynamics
//------------------------------------------------------------------------------

// update particle mass
void Column::update_mass(void) {

  int nCategories = columnDimensions->size("nCategories");

  char snowDensityStr[5] = "rhos";
  double snowDensity = column_get_constant(snowDensityStr);
  char iceDensityStr[5] = "rhoi";
  double iceDensity = 900.0;//column_get_constant(iceDensityStr);

  auto particles_mass = particles->mass;
  auto particles_ice_area = particles->ice_area;

  // temporary copy until we get Icepack variables on device
  Kokkos::View<double**> iceVolumeCategory_device("ice volume category device", *(particles->nParticles), nCategories);
  Kokkos::View<double**> snowVolumeCategory_device("snow volume category device", *(particles->nParticles), nCategories);
  Kokkos::View<double*> effectiveElementArea_device("New effective element area device", *(particles->nParticles));
  Kokkos::View<double*> iceAreaCell_device("Ice area cell device", *(particles->nParticles));
  Kokkos::View<double**>::HostMirror iceVolumeCategory_host  = Kokkos::create_mirror(iceVolumeCategory_device);
  Kokkos::View<double**>::HostMirror snowVolumeCategory_host = Kokkos::create_mirror(snowVolumeCategory_device);
  Kokkos::View<double*>::HostMirror effectiveElementArea_host = Kokkos::create_mirror(effectiveElementArea_device);
  Kokkos::View<double*>::HostMirror iceAreaCell_host = Kokkos::create_mirror(iceAreaCell_device);

  Kokkos::parallel_for("Column:ice and snow volume category copy", Kokkos::RangePolicy<host_execution_space>
      (0,*(particles->nParticles)), KOKKOS_LAMBDA (const int iParticle) {

    for (int iCategory = 0 ; iCategory < nCategories ; iCategory++) {
      iceVolumeCategory_host(iParticle,iCategory)  = (*iceVolumeCategory)(iParticle,iCategory);
      snowVolumeCategory_host(iParticle,iCategory) = (*snowVolumeCategory)(iParticle,iCategory);
    } // iCategory

    effectiveElementArea_host(iParticle) = (*effectiveElementArea)(iParticle);
    iceAreaCell_host(iParticle) = (*iceAreaCell)(iParticle);

  }); // iParticle

  Kokkos::deep_copy(iceVolumeCategory_device, iceVolumeCategory_host);
  Kokkos::deep_copy(snowVolumeCategory_device, snowVolumeCategory_host);
  Kokkos::deep_copy(effectiveElementArea_device, effectiveElementArea_host);
  Kokkos::deep_copy(iceAreaCell_device, iceAreaCell_host);

  double minimumIceThickness;
  configs->get({"ConfigGroup:dynamics","Config:minimumIceThickness"}, minimumIceThickness, true);
  double minimumIceConcentration;
  configs->get({"ConfigGroup:dynamics","Config:minimumIceConcentration"}, minimumIceConcentration, true);

  Kokkos::parallel_for("Column::update_mass", Kokkos::RangePolicy<device_execution_space>
      (0,*(particles->nParticles)), KOKKOS_LAMBDA (const int iParticle) {

    double minimumConcentration = std::max(iceAreaCell_device(iParticle),minimumIceConcentration);
    particles_ice_area[iParticle] = effectiveElementArea_device(iParticle) * minimumConcentration;

    double iceVolume = 0.0;
    double snowVolume = 0.0;
    for (int iCategory = 0 ; iCategory < nCategories ; iCategory++) {
      iceVolume  += iceVolumeCategory_device(iParticle,iCategory);
      snowVolume += snowVolumeCategory_device(iParticle,iCategory);
    } // iCategory

    iceVolume = std::max(iceVolume,minimumIceThickness);

    particles_mass[iParticle] = (iceVolume * iceDensity + snowVolume * snowDensity) *
      effectiveElementArea_device(iParticle);

  }); // iParticle

}

// Compute mean/minimum ice thickness, needed by contact model
void Column::compute_mean_min_thickness(void) {

  int nCategories = columnDimensions->size("nCategories");

  // temporary copy until we get Icepack variables on device
  Kokkos::View<double**> iceVolumeCategory_device("ice volume category device", *(particles->nParticles), nCategories);
  Kokkos::View<double**> iceAreaCategory_device("ice area category device", *(particles->nParticles), nCategories);
  Kokkos::View<double**>::HostMirror iceVolumeCategory_host  = Kokkos::create_mirror(iceVolumeCategory_device);
  Kokkos::View<double**>::HostMirror iceAreaCategory_host  = Kokkos::create_mirror(iceAreaCategory_device);

  Kokkos::parallel_for("Column:ice area and volume category copy for compute_mean_min_thickness", Kokkos::RangePolicy<host_execution_space>
      (0,*(particles->nParticles)), KOKKOS_LAMBDA (const int iParticle) {

    for (int iCategory = 0 ; iCategory < nCategories ; iCategory++) {
      iceVolumeCategory_host(iParticle,iCategory)  = (*iceVolumeCategory)(iParticle,iCategory);
      iceAreaCategory_host(iParticle,iCategory)    = (*iceAreaCategory)(iParticle,iCategory);
    } // iCategory

  }); // iParticle

  Kokkos::deep_copy(iceVolumeCategory_device, iceVolumeCategory_host);
  Kokkos::deep_copy(iceAreaCategory_device, iceAreaCategory_host);

  double minimumIceThickness;
  configs->get({"ConfigGroup:dynamics","Config:minimumIceThickness"}, minimumIceThickness, true);

  // this ensure KOKKOS_LAMBDA capture of these variables
  auto particles_minThickness = particles->minThickness;
  auto particles_meanThickness = particles->meanThickness;

  Kokkos::parallel_for("Column::compute_mean_min_thickness",
      Kokkos::RangePolicy<device_execution_space>
      (0,*(particles->nParticles)), KOKKOS_LAMBDA (const int iParticle) {

    //For now, assume lowest category is the minimum thickness, and
    // that categories are sorted in ascending order of thickness
    particles_minThickness(iParticle) = iceVolumeCategory_device(iParticle,0)
            / iceAreaCategory_device(iParticle,0);

    //Fractions will probably add up to one, but just in case:
    double meanThickness = 0;
    for (int iCategory = 0; iCategory < nCategories; iCategory++){
      meanThickness += iceVolumeCategory_device(iParticle, iCategory);
    }
    meanThickness = std::max(meanThickness,minimumIceThickness);
    particles_meanThickness(iParticle) = meanThickness;

  }); // iParticle

}

// convert ice/snow volumes to thicknesses
void Column::to_thickness(void) {

  int nCategories = columnDimensions->size("nCategories");

  for (int iParticle = 0 ; iParticle < *(particles->nParticles) ; iParticle++) {
    for (int iCategory = 0; iCategory < nCategories; iCategory++){
      if ((*iceAreaCategory)(iParticle, iCategory) > 0.0) {
	(*iceVolumeCategory) (iParticle, iCategory) /= (*iceAreaCategory)(iParticle, iCategory);
	(*snowVolumeCategory)(iParticle, iCategory) /= (*iceAreaCategory)(iParticle, iCategory);
      }
    } // iCategory
  } // iParticle

}

// convert ice/snow thicknesses to volumes
void Column::from_thickness(void) {

  int nCategories = columnDimensions->size("nCategories");

  for (int iParticle = 0 ; iParticle < *(particles->nParticles) ; iParticle++) {
    for (int iCategory = 0; iCategory < nCategories; iCategory++){
      (*iceVolumeCategory) (iParticle, iCategory) *= (*iceAreaCategory)(iParticle, iCategory);
      (*snowVolumeCategory)(iParticle, iCategory) *= (*iceAreaCategory)(iParticle, iCategory);
    } // iCategory
  } // iParticle

}

// delete particles that have too little sea ice
void Column::delete_particles(void) {

  auto particles_type = particles->type;

  Kokkos::View<double*> iceAreaCell_device("Ice area cell device", *(particles->nParticles));
  Kokkos::View<double*>::HostMirror iceAreaCell_host = Kokkos::create_mirror(iceAreaCell_device);

  Kokkos::parallel_for("Column:ice area copy", Kokkos::RangePolicy<host_execution_space>
      (0,*(particles->nParticles)), KOKKOS_LAMBDA (const int iParticle) {

    iceAreaCell_host(iParticle) = (*iceAreaCell)(iParticle);

  }); // iParticle

  Kokkos::deep_copy(iceAreaCell_device, iceAreaCell_host);

  double iceConcentrationRemovalLimit;
  configs->get({"ConfigGroup:","Config:iceConcentrationRemovalLimit"}, iceConcentrationRemovalLimit, true);

  Kokkos::parallel_for("Column::delete_particles set type", Kokkos::RangePolicy<device_execution_space>
      (0,*(particles->nParticles)), KOKKOS_LAMBDA (const int iParticle) {

      if (iceAreaCell_device(iParticle) < iceConcentrationRemovalLimit and
	  particles_type[iParticle] == 1) {
	particles_type[iParticle] = 3;
      } // particles to delete

  }); // iParticle

  // create group of atoms to delete
  lammpsInstance->one("group particlesToDelete type 3");

  // delete this group
  lammpsInstance->one("delete_atoms group particlesToDelete compress no bond yes");

}

//------------------------------------------------------------------------------
// Processor transfer
//------------------------------------------------------------------------------

// init the transfer across processors for column variables
void Column::init_processor_transfer(void) {

  bool testColumnProcTransfer;
  configs->get({"ConfigGroup:testing","Config:testColumnProcTransfer"}, testColumnProcTransfer);
  if (testColumnProcTransfer) {

    columnVariablesTestTransfer = new DEMSI::ColumnVariable<int>(particles, columnDimensions,
          "columnVariablesTestTransfer", {}, true);
    columnVariables->add(columnVariablesTestTransfer);

    for (int iParticle=0 ; iParticle < *(particles->nParticles) ; iParticle++) {
      (*columnVariablesTestTransfer)(iParticle) = (int) particles->globalID[iParticle];
    }

  } // testColumnProcTransfer

}

// Initialize transfer column data between processors after remapping.
void Column::reset_processor_transfer_remap(const int nParticles) {

  columnVariables->init_processor_transfer_remap(nParticles);

}

// reset the transfer across processors for column variables
void Column::reset_processor_transfer(void) {

  columnVariables->init_processor_transfer();

  bool testColumnProcTransfer;
  configs->get({"ConfigGroup:testing","Config:testColumnProcTransfer"}, testColumnProcTransfer);
  if (testColumnProcTransfer) {

    columnVariablesTestTransfer->resize(*(particles->nParticles));

    for (int iParticle=0 ; iParticle < *(particles->nParticles) ; iParticle++) {
      (*columnVariablesTestTransfer)(iParticle) = (int) particles->globalID[iParticle];
    } // iParticle

  } // testColumnProcTransfer

}

// perform the transfer across processors for column variables
void Column::processor_transfer(void) {

  columnVariables->processor_transfer();

  bool testColumnProcTransfer;
  configs->get({"ConfigGroup:testing","Config:testColumnProcTransfer"}, testColumnProcTransfer);
  if (testColumnProcTransfer) {

      bool fail = false;

      for (int iParticle=0 ; iParticle < *(particles->nParticles) ; iParticle++) {

	// check that the testParticleData contains the globalID of the particles
	if ((int) particles->globalID[iParticle] != (*columnVariablesTestTransfer)(iParticle)) {
	  (*log)(DEMSI::LOG::ERROR) << "Column communication test failed: " << iParticle << " " << particles->globalID[iParticle] << " " << (*columnVariablesTestTransfer)(iParticle) << std::endl;
	  fail = true;
	}

      }

      MPI_Barrier(partition->comm());

      if (fail) {
	log->abort("Column communication test failed!");
      }

  } // testColumnProcTransfer

}

//------------------------------------------------------------------------------
// Warnings
//------------------------------------------------------------------------------

// print icepack warning messages
bool Column::check_warnings(const std::string message) {

  int nWarnings = column_warnings_number();
  if (nWarnings > 0) (*log)(DEMSI::LOG::WARNING) << "Column warning: " << message << ": nWarnings: " << nWarnings << std::endl;
  for (int i = 0 ; i < nWarnings ; i++) {
    char warning[columnStrLen];
    column_warnings_getone(&i, warning);
    std::string warningStr(warning);
    (*log)(DEMSI::LOG::WARNING) << "Column warning: " << message << ": " << warningStr << std::endl;
  }
  column_warnings_reset();
  log->check(not column_warnings_aborted(), "Column: Icepack abort");
  if (nWarnings > 0) {
    return true;
  } else {
    return false;
  }

}

// print icepack warning messages from within particle loop
bool Column::check_warnings(const std::string message, const int iParticle) {

  int nWarnings = column_warnings_number();
  if (nWarnings > 0) (*log)(DEMSI::LOG::WARNING) << "Column warning: " << message << ": for particle: " << iParticle << ": nWarnings: " << nWarnings << std::endl;
  for (int i = 0 ; i < nWarnings ; i++) {
    char warning[columnStrLen];
    column_warnings_getone(&i, warning);
    std::string warningStr(warning);
    (*log)(DEMSI::LOG::WARNING) << "Column warning: " << message << ": for particle: " << iParticle << ": " << warningStr << " " << particles->type(iParticle) << std::endl;
  }
  column_warnings_reset();
  log->check(not column_warnings_aborted(), "Column: Icepack abort for particle: ", std::to_string(iParticle));
  if (nWarnings > 0) {
    return true;
  } else {
    return false;
  }

}

} // namespace DEMSI
