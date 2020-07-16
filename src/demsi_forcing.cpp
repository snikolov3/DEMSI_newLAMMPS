#include "demsi_forcing.h"
#include "demsi_logging.h"
#include "demsi_file_utils.h"
#include "demsi_time.h"

#include <vector>
#include <list>
#include <netcdf.h>

namespace DEMSI {

//------------------------------------------------------------------------------
// FieldIO class
//------------------------------------------------------------------------------

// FieldIO class constructor.
FieldIO::FieldIO(const std::string filenameTemplateIn, DEMSI::Log* logIn) {

  filenameTemplate = filenameTemplateIn;
  filename = "";
  log = logIn;

}

// FieldIO class destructor.
FieldIO::~FieldIO() {

  for (int i = 0; i < times.size(); ++i) {
    delete times[i];
  }
  times.clear();

}

// read in a single time of an input field.
void FieldIO::read_field(const std::string fieldName, const Time readTime, Grid* grid, Kokkos::View<double**>* field) {

  std::string filenameNew = expand_filename_template(filenameTemplate, readTime);

  int err;

  if (filenameNew != filename) {

    if (filename != "") {

      err = nc_close(nc_id);
      log->check(err == NC_NOERR, "Couldn't close file: ", filename);

    }

    filename = filenameNew;

    err = nc_open(filename.c_str(), 0, &nc_id);
    log->check(err == NC_NOERR, "Couldn't open file: ", filename);

    // get whether climatology
    int climatologyTmp;
    err = nc_get_att_int(nc_id, NC_GLOBAL, "climatology", &climatologyTmp);
    log->check(err == NC_NOERR, "Climatology attribute not present for: ", filename);
    climatology = (bool) climatologyTmp;

    // read times
    int dimidTime;
    err = nc_inq_dimid(nc_id, "Time", &dimidTime);
    log->check(err == NC_NOERR, "Can't get dimid for Time: ", filename);

    size_t nTimesTmp;
    err = nc_inq_dimlen(nc_id, dimidTime, &nTimesTmp);
    log->check(err == NC_NOERR, "Can't get Time dimension: ", filename);
    int nTimes = (int) nTimesTmp;

    int varidTime;
    err = nc_inq_varid(nc_id, "Time", &varidTime);
    log->check(err == NC_NOERR, "Can't get varid for Time: ", filename);

    size_t startTime[2] = {0, 0};
    size_t countTime[2] = {(size_t) nTimes, 64};

    char timesIn[nTimes][64];
    err = nc_get_vara_text(nc_id, varidTime, startTime, countTime, &timesIn[0][0]);
    log->check(err == NC_NOERR, "Can't read in times: ", filename);

    for (int i = 0; i < times.size(); ++i) {
      delete times[i];
    }
    times.clear();

    for (int i=0 ; i<nTimes ; i++) {
      Time *timeIn = new Time(readTime.get_calendar(), timesIn[i], log);
      times.push_back(timeIn);
    }

  }

  // find date time index
  int timeIndex = -1;
  if (not climatology) {
    for (unsigned i=0; i < times.size(); i++) {
      if (readTime == *(times[i])) {
	timeIndex = i;
	break;
      }
    }
  } else {
    for (unsigned i=0; i < times.size(); i++) {
      if (readTime.climatology() == *(times[i])) {
	timeIndex = i;
	break;
      }
    }
  }
  log->check(timeIndex >= 0, "Can't find time in file: ", filename, " for time: ", readTime.get());

  // read in data
  grid->read_field_varying(field, nc_id, fieldName, timeIndex);

}

//------------------------------------------------------------------------------
// Forcing field class
//------------------------------------------------------------------------------

// Forcing field class constructor.
ForcingField::ForcingField(std::string fieldNameIn, Kokkos::View<double**>* fieldOutIn, FieldIO* fieldIOIn, Grid* gridIn) {

  fieldName = fieldNameIn;

  grid = gridIn;

  fieldInPtr1 = new Kokkos::View<double **>("fieldInPtr1"+fieldName,gridIn->nx(),grid->ny());
  fieldInPtr2 = new Kokkos::View<double **>("fieldInPtr2"+fieldName,gridIn->nx(),grid->ny());

  fieldOut = fieldOutIn;

  fieldIO = fieldIOIn;

}

// Forcing field class destructor.
ForcingField::~ForcingField() {

  delete fieldInPtr1;
  delete fieldInPtr2;

}

// over load << operator.
std::ostream & operator<<(std::ostream & os, const ForcingField & forcingField)
{
  os << "ForcingField: {" << forcingField.get_name() << "}";
  return os;
}

// Return the forcing field name.
std::string ForcingField::get_name(void) const {
  return fieldName;
}

// Swap the input data in the two forcing time slots.
void ForcingField::shift_data(void) {

  Kokkos::View<double**>* tmp = fieldInPtr1;
  fieldInPtr1 = fieldInPtr2;
  fieldInPtr2 = tmp;

}

// Load new data into one of the input data slots.
void ForcingField::load_new_data(const Time readTime, const int slot) {

  if (slot == 1) {
    fieldIO->read_field(fieldName, readTime, grid, fieldInPtr1);
  } else if (slot == 2) {
    fieldIO->read_field(fieldName, readTime, grid, fieldInPtr2);
  }

}

// Interpolate the input data from the two time slots to the output forcing field.
void ForcingField::interpolate(double interpolants[2]) {

  // copy interpolants to device
  Kokkos::View<double[2], Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> > interpolants_host(&interpolants[0]);
  Kokkos::View<double[2]> interpolants_device("interpolants");
  Kokkos::deep_copy(interpolants_device, interpolants_host);

  typedef typename Kokkos::Experimental::MDRangePolicy< Kokkos::Experimental::Rank<2> > MDPolicyType_2D;
  MDPolicyType_2D mdpolicy_2d( {{0,0}}, {{ grid->nx(), grid->ny() }} );

  // prepare variables used in loop for capture on device
  auto fieldOut = *(this->fieldOut);
  auto fieldIn1 = *(this->fieldInPtr1);
  auto fieldIn2 = *(this->fieldInPtr2);

  Kokkos::parallel_for( "ForcingField::interpolate", mdpolicy_2d,
      KOKKOS_LAMBDA( const int i, const int j ) {
          fieldOut(i,j) = fieldIn1(i,j) * interpolants_device[0] + fieldIn2(i,j) * interpolants_device[1];
  });
  Kokkos::fence();

}

// Return the field name
std::string ForcingField::get_name(void) {
  return fieldName;
}

// Return forcing field
Kokkos::View<double**>* ForcingField::get_field(void) {
  return fieldOut;
}

//------------------------------------------------------------------------------
// Forcing group class
//------------------------------------------------------------------------------

// Forcing group class constructor
ForcingGroup::ForcingGroup(std::string nameIn, Grid* gridIn, std::string filenameTemplateIn, DEMSI::Clock* simulationClockIn, const DEMSI::Time forcingReferenceTimeIn, const DEMSI::TimeInterval forcingIntervalIn, const bool timeVaryingForcingIn, DEMSI::Log* logIn) {

  name = nameIn;

  grid = gridIn;
  log = logIn;

  filenameTemplate = filenameTemplateIn;

  simulationClock = simulationClockIn;
  forcingReferenceTime = forcingReferenceTimeIn;
  forcingInterval = forcingIntervalIn;

  timeVaryingForcing = timeVaryingForcingIn;

  // time varying forcing
  if (timeVaryingForcing) {

    // set the initial forcing times
    Time currentTime = simulationClock->get_time();

    // initialize the forcing data time
    forcingTimes[0] = forcingReferenceTime;
    forcingTimes[1] = forcingReferenceTime + forcingInterval;

    if (currentTime >= forcingTimes[0] and
	currentTime <  forcingTimes[1]) {
      // current forcing clock time is in the current interval

    } else if (currentTime < forcingTimes[0]) {
      // current forcing clock time earlier than the current interval

      while (true) {

	forcingTimes[1] = forcingTimes[0];
	forcingTimes[0] = forcingTimes[1] - forcingInterval;

	if (currentTime >= forcingTimes[0] and
	    currentTime <  forcingTimes[1]) {
	  break;
	}
      }

    } else if (currentTime >= forcingTimes[1]) {
      // current forcing clock time later than the current interval

      while (true) {

	forcingTimes[0] = forcingTimes[1];
	forcingTimes[1] = forcingTimes[0] + forcingInterval;

	if (currentTime >= forcingTimes[0] and
	    currentTime <  forcingTimes[1]) {
	  break;
	}
      }
    }

    // set the forcing group alarm
    forcingAlarm.set(simulationClock, forcingTimes[1], forcingInterval);

  } else { // timeVaryingForcing

    forcingTimes[0] = forcingReferenceTime;
    forcingTimes[1] = forcingReferenceTime;

  } // timeVaryingForcing

  // set up filename template
  fieldIO = new FieldIO(filenameTemplate, log);

}

// Forcing group class destructor
ForcingGroup::~ForcingGroup() {

  delete fieldIO;

  for (int i = 0; i < forcingFields.size(); ++i) {
    delete forcingFields[i];
  }
  forcingFields.clear();

}

// Initialize the input forcing data for both time slots.
void ForcingGroup::init_field_data(void) {

  // load initial data
  std::vector<ForcingField*>::iterator forcingFieldItr;
  for ( forcingFieldItr = forcingFields.begin(); forcingFieldItr != forcingFields.end(); ++forcingFieldItr ) {
    (*forcingFieldItr)->load_new_data(forcingTimes[0],1);
    (*forcingFieldItr)->load_new_data(forcingTimes[1],2);
  }

}

// over load << operator
std::ostream & operator<<(std::ostream & os, const ForcingGroup & forcingGroup)
{
  os << "ForcingGroup: {" << forcingGroup.get_name() << std::endl;
  for ( int i=0 ; i<forcingGroup.get_field_number() ; ++i) {
    os << forcingGroup.get_field(i) << std::endl;
  }
  os << "}";
  return os;
}

// Add a forcing field to the forcing group.
void ForcingGroup::add_forcing_field(std::string fieldName, Kokkos::View<double**>* fieldOut) {

  // check name isnt already in use
  std::vector<ForcingField*>::iterator forcingFieldItr;
  for ( forcingFieldItr = forcingFields.begin(); forcingFieldItr != forcingFields.end(); ++forcingFieldItr ) {
    log->check((*forcingFieldItr)->get_name() != fieldName, "Forcing field name already in use: ", fieldName);
  }

  ForcingField *forcingField = new ForcingField(fieldName, fieldOut, fieldIO, grid);

  forcingFields.push_back(forcingField);

  // add field pointer to fields to output
  forcingWriteFields.push_back(fieldOut);

}

// Return the forcing group name.
std::string ForcingGroup::get_name(void) const {
  return name;
}

// Return the number of forcing fields in the forcing group.
int ForcingGroup::get_field_number(void) const {
  return forcingFields.size();
}

// Return a forcing field from a forcing group given its index.
DEMSI::ForcingField ForcingGroup::get_field(int fieldIndex) const {
  return (*forcingFields[fieldIndex]);
}

// Calculate the linear interpolants given the current time and times of the input data slots
void get_interpolants(const Time forcingTimes[2], const Time currentTime, double (&interpolants)[2]) {

  long long int diff0 = forcingTimes[1] - forcingTimes[0];
  long long int diff1 = forcingTimes[1] - currentTime;

  interpolants[0] = (double) diff1 / (double) diff0;
  interpolants[1] = 1.0 - interpolants[0];

}

// Set fixed interpolants for time constant forcing
void get_interpolants_fixed(double (&interpolants)[2]) {

  interpolants[0] = 1.0;
  interpolants[1] = 0.0;

}

// Update the output forcing fields for the current simulation time for a forcing group.
void ForcingGroup::update(void) {

  // check if the alarm is ringing for loading new data
  if (timeVaryingForcing and forcingAlarm.is_ringing()) {

    // update forcing times
    forcingTimes[0] = forcingTimes[1];
    forcingTimes[1] = forcingTimes[0] + forcingInterval;

    std::vector<ForcingField*>::iterator forcingFieldItr;
    for ( forcingFieldItr = forcingFields.begin(); forcingFieldItr != forcingFields.end(); ++forcingFieldItr ) {
      (*forcingFieldItr)->shift_data();
      (*forcingFieldItr)->load_new_data(forcingTimes[1],2);
    }

    // reset forcing alarm
    forcingAlarm.reset();

  }

  // get the linear interpolants
  double interpolants[2];

  if (timeVaryingForcing) {
    get_interpolants(forcingTimes, simulationClock->get_time(), interpolants);
  } else {
    get_interpolants_fixed(interpolants);
  }

  // iterate over forcing fields and perform interpolation
  std::vector<ForcingField*>::iterator forcingFieldItr;
  for ( forcingFieldItr = forcingFields.begin(); forcingFieldItr != forcingFields.end(); ++forcingFieldItr ) {
    (*forcingFieldItr)->interpolate(interpolants);
  }

}

//------------------------------------------------------------------------------
// Forcing instance class
//------------------------------------------------------------------------------

// Forcing Instance class constructor
ForcingInstance::ForcingInstance(Grid* gridIn, DEMSI::Log* logIn) {

  grid = gridIn;
  log = logIn;

}

// Forcing Instance class destructor
ForcingInstance::~ForcingInstance() {

  for (int i = 0; i < forcingGroups.size(); ++i) {
    delete forcingGroups[i];
  }
  forcingGroups.clear();

}

// over load << operator
std::ostream & operator<<(std::ostream & os, const ForcingInstance & forcingInstance)
{
  os << "ForcingInstance: {" << std::endl;
  for ( int i=0 ; i<forcingInstance.get_group_number() ; ++i) {
    os << forcingInstance.get_group(i) << std::endl;
  }
  os << "}";
  return os;
}

// Add a forcing group to this forcing instance object.
  void ForcingInstance::add_forcing_group(std::string name, std::string filenameTemplate, DEMSI::Clock* simulationClock, const DEMSI::Time forcingTimeReference, const DEMSI::TimeInterval forcingInterval, const bool timeVaryingForcing) {

  // check name isnt already in use
  std::vector<ForcingGroup*>::iterator forcingGroupItr;
  for ( forcingGroupItr = forcingGroups.begin(); forcingGroupItr != forcingGroups.end(); ++forcingGroupItr ) {
    log->check((*forcingGroupItr)->get_name() != name, "Forcing group name already in use: ", name);
  }

  ForcingGroup *forcingGroup = new ForcingGroup(name, grid, filenameTemplate, simulationClock, forcingTimeReference, forcingInterval, timeVaryingForcing, log);

  forcingGroups.push_back(forcingGroup);

}

// Add a forcing field to a forcing group of this forcing instance object.
void ForcingInstance::add_forcing_field(std::string groupName, std::string fieldName, Kokkos::View<double**>* fieldOut) {

  std::vector<ForcingGroup*>::iterator forcingGroupItr;
  for ( forcingGroupItr = forcingGroups.begin(); forcingGroupItr != forcingGroups.end(); ++forcingGroupItr ) {
    if ((*forcingGroupItr)->get_name() == groupName) {
      (*forcingGroupItr)->add_forcing_field(fieldName, fieldOut);
    }
  }

}

// Init the field data for all forcing fields in all forcing groups for this forcing instance.
void ForcingInstance::init_field_data(void) {

  std::vector<ForcingGroup*>::iterator forcingGroupItr;
  for ( forcingGroupItr = forcingGroups.begin(); forcingGroupItr != forcingGroups.end(); ++forcingGroupItr ) {
    (*forcingGroupItr)->init_field_data();
  }

}

// Return the number of forcing groups in this forcing instance.
int ForcingInstance::get_group_number(void) const {
  return forcingGroups.size();
}

// Return a forcing group from this forcing instance given its index.
DEMSI::ForcingGroup ForcingInstance::get_group(int groupIndex) const {
  return (*forcingGroups[groupIndex]);
}

// Update all output forcing fields in this forcing instance.
void ForcingInstance::update(void) {

  std::vector<ForcingGroup*>::iterator forcingGroupItr;
  for ( forcingGroupItr = forcingGroups.begin(); forcingGroupItr != forcingGroups.end(); ++forcingGroupItr ) {
    (*forcingGroupItr)->update();
  }

}

//------------------------------------------------------------------------------
// Forcing class
//------------------------------------------------------------------------------

// Forcing class constructor.
Forcing::Forcing(DEMSI::Log* logIn, DEMSI::Configs *configs, DEMSI::Grid* gridIn, Calendar* calendar, Clock* clock) {

  log = logIn;
  grid = gridIn;

  // use forcing
  configs->get({"ConfigGroup:useSections","Config:useForcing"}, useForcing);
  if (useForcing) {

    // register atmospheric forcing fields
    forcingFields["xAtmWind"]          = &xAtmWind;
    forcingFields["yAtmWind"]          = &yAtmWind;
    forcingFields["longwaveFluxDown"]  = &longwaveFluxDown;
    forcingFields["shortwaveFluxDown"] = &shortwaveFluxDown;
    forcingFields["airTemperature"]    = &airTemperature;
    forcingFields["specificHumidity"]  = &specificHumidity;
    forcingFields["precipitationRate"] = &precipitationRate;
    forcingFields["cloudFraction"]     = &cloudFraction;

    // register oceanic forcing fields
    forcingFields["xOcnCurrents"]             = &xOcnCurrents;
    forcingFields["yOcnCurrents"]             = &yOcnCurrents;
    forcingFields["seaSurfaceTemperature"]    = &seaSurfaceTemperature;
    forcingFields["seaSurfaceSalinity"]       = &seaSurfaceSalinity;
    forcingFields["seaSurfaceTiltU"]          = &seaSurfaceTiltU;
    forcingFields["seaSurfaceTiltV"]          = &seaSurfaceTiltV;
    forcingFields["oceanMixedLayerDepth"]     = &oceanMixedLayerDepth;
    forcingFields["oceanHeatFluxConvergence"] = &oceanHeatFluxConvergence;

    // register forcing variables for possible gridded write
    grid->register_field_for_write(&xAtmWind);
    grid->register_field_for_write(&yAtmWind);
    grid->register_field_for_write(&longwaveFluxDown);
    grid->register_field_for_write(&shortwaveFluxDown);
    grid->register_field_for_write(&airTemperature);
    grid->register_field_for_write(&specificHumidity);
    grid->register_field_for_write(&precipitationRate);
    grid->register_field_for_write(&cloudFraction);
    grid->register_field_for_write(&xOcnCurrents);
    grid->register_field_for_write(&yOcnCurrents);
    grid->register_field_for_write(&seaSurfaceTemperature);
    grid->register_field_for_write(&seaSurfaceSalinity);
    grid->register_field_for_write(&seaSurfaceTiltU);
    grid->register_field_for_write(&seaSurfaceTiltV);
    grid->register_field_for_write(&oceanMixedLayerDepth);
    grid->register_field_for_write(&oceanHeatFluxConvergence);

    // create forcing instance
    forcingInstance = new DEMSI::ForcingInstance(grid, log);

    // get forcing group names
    std::vector<std::string> forcingGroupNames = configs->get_array({"ConfigGroup:forcing","Array:ForcingGroups"}, "ForcingGroup", "name");

    // iterate over forcing groups
    for (int iForcingGroup=0 ; iForcingGroup<forcingGroupNames.size() ; iForcingGroup++) {

      // forcing group configs
      std::string forcingFilenameTemplate;
      bool timeVaryingForcing;
      std::string forcingReferenceTimeStr;
      std::string forcingIntervalStr;

      // forcing group configs location
      std::list<std::string> configLocationRoot, configLocation;
      configLocationRoot.push_back("ConfigGroup:forcing");
      configLocationRoot.push_back("Array:ForcingGroups");
      configLocationRoot.push_back("ForcingGroup:" + forcingGroupNames[iForcingGroup]);

      // get forcing group configs
      configLocation = configLocationRoot;
      configLocation.push_back("Config:forcingFilenameTemplate");
      configs->get(configLocation, forcingFilenameTemplate);

      if (forcingFilenameTemplate != "none") {

	configLocation = configLocationRoot;
	configLocation.push_back("Config:timeVaryingForcing");
	configs->get(configLocation, timeVaryingForcing);

	configLocation = configLocationRoot;
	configLocation.push_back("Config:forcingReferenceTime");
	configs->get(configLocation, forcingReferenceTimeStr);

	configLocation = configLocationRoot;
	configLocation.push_back("Config:forcingInterval");
	configs->get(configLocation, forcingIntervalStr);

	DEMSI::Time forcingReferenceTime(calendar, forcingReferenceTimeStr, log);
	DEMSI::TimeInterval forcingInterval(forcingIntervalStr, log);

	// add forcing group to forcing instance
	forcingInstance->add_forcing_group(forcingGroupNames[iForcingGroup], forcingFilenameTemplate, clock, forcingReferenceTime, forcingInterval, timeVaryingForcing);

      }

      // get forcing group fields
      std::list<std::string> forcingFieldsLocation = configLocationRoot;
      forcingFieldsLocation.push_back("Array:ForcingFields");
      std::vector<std::string> forcingFieldNames = configs->get_array(forcingFieldsLocation, "ForcingField");

      // iterate over forcing group fields
      for (int iForcingField = 0 ; iForcingField < forcingFieldNames.size() ; iForcingField++) {

	Kokkos::View<double**>* forcingFieldPtr;
	try {
	  forcingFieldPtr = forcingFields.at(forcingFieldNames[iForcingField]);
	} catch (const std::out_of_range& e) {
	  log->abort("Input forcing field: ", forcingFieldNames[iForcingField], " not found in DEMSI");
	}
	if (forcingFilenameTemplate != "none") {
	  forcingInstance->add_forcing_field(forcingGroupNames[iForcingGroup], forcingFieldNames[iForcingField], forcingFieldPtr);
	}

	(*forcingFieldPtr) = Kokkos::View<double **>(forcingFieldNames[iForcingField], grid->nx(), grid->ny());

      }

    }

    forcingInstance->init_field_data();

  } // useForcing

}

// Forcing class destructor.
Forcing::~Forcing(void) {

  if (useForcing) delete forcingInstance;

}

// Update the forcing fields.
void Forcing::update(void) {

  if (useForcing) {

    // update input fields
    forcingInstance->update();

  } // useForcing

}

} // namespace DEMSI
