#ifndef DEMSI_FORCING_H_
#define DEMSI_FORCING_H_

#include <vector>
#include <map> 

#include "demsi_logging.h"
#include "demsi_time.h"
#include "demsi_configs.h"
#include "demsi_grid.h"

/*! \file   demsi_forcing.h
    \brief  Header file for the DEMSI::Forcing classes
*/

namespace DEMSI {

//------------------------------------------------------------------------------
/*! \class FieldIO
    \brief Class for reading in two dimensional forcing fields.

    This class has the ability to read in time slices of 2D fields defined on
    the demsi grid. The files are in netcdf format with a 'Time', 'nx', 'ny' and
    'strLen' dimensions. Files also contain a 'Time' variable of dimension
    (Time,strLen) with times in format 'YYYY-MM-DD_hh:mm:ss', and other fields
    of dimension (Time,nx,ny). This class reads in individual times from the
    defined inoput fields and allows the netcdf to remain open until a new file
    must be opened. The class operates on filename templates, such as
    'filename_$Y.nc' where the tokens are expanded with the desired read time.
*/
//------------------------------------------------------------------------------
class FieldIO{
public:

  /*! \brief FieldIO constructor
      \param filenameTemplateIn Forcing filename template to read data with.
      \param logIn Pointer to log object.
   */
  FieldIO(const std::string filenameTemplateIn, DEMSI::Log* logIn);

  /*! \brief Destructor.*/
  ~FieldIO();

  /*! \brief Read a single time of a forcing field.
      \param varname The input field name in the netcdf file.
      \param readTime The input time to read data from.
      \param grid The grid the input data is defined on.
      \param field The field pointer to write input data to.
  */
  void read_field(const std::string varname, const Time readTime, Grid* grid, Kokkos::View<double**>* field);

private:

  std::string filenameTemplate; /*! forcing filename template. */
  std::string filename; /*! current forcing filename. */

  int nc_id; /*! netcdf fileID. */

  bool climatology; /*! whether file is a climatology. */

  std::vector<DEMSI::Time*> times; /*! times of input data in current file. */

  DEMSI::Log* log; /*! Pointer to log object. */

};

//------------------------------------------------------------------------------
/*! \class ForcingField
    \brief Class that updates a single forcing field.

    This class updates a single forcing fields. It contains two input field
    arrays for the two input forcing time slots. It performs linear
    interpolation of these two fields to generate the output forcing field data
    and uses interpolants provides by the owning forcing group to do it.
*/
//------------------------------------------------------------------------------
class ForcingField{
public:

  /*! \brief ForcingField constructor
      \param fieldNameIn Name of the field in the input netcdf file.
      \param fieldOutIn Pointer to the output forcing field to be written to.
      \param fieldIOIn Pointer to the FieldIO object for this forcing group.
      \param gridIn Pointer to the Grid object for this forcing.
   */
  ForcingField(std::string fieldNameIn, Kokkos::View<double**>* fieldOutIn, FieldIO* fieldIOIn, Grid* gridIn);

  /*! \brief Destructor.*/
  ~ForcingField();

  /*! \brief Write the forcing field object to a stream.
      \param os Stream to write to.
      \param forcingField Forcing field object to write to stream.
  */
  friend std::ostream & operator<<(std::ostream & os, const ForcingField & forcingField);

  /*! \brief Returns the forcing field name.
      \return The field name in the input netcdf file of this forcing field.
  */
  std::string get_name(void) const;

  /*! \brief Swap the input data stored in the two input forcing time slots.*/
  void shift_data(void);

  /*! \brief Loads new data into a desired time slot.
      \param readTime Time to read input data from.
      \param slot The time slot to read input data into (1,2).
  */
  void load_new_data(const Time readTime, const int slot);

  /*! \brief Perform linear interpolation on the input fields and write to the output forcing field.
      \param interpolants The linear interpolants to use to interpolate.
  */
  void interpolate(double interpolants[2]);

  /*! \brief Return the field name
      \return Name of the field
  */
  std::string get_name(void);

  /*! \brief Return a pointer to the output field object. 
      \return Pointer to the output field object.
   */
  Kokkos::View<double**>* get_field(void);

private:

  std::string fieldName; /*! Name of field in inout netcdf file. */

  Kokkos::View<double**>* fieldOut; /*! Pointer to output forcing field. */

  Kokkos::View<double**>* fieldInPtr1; /*! Input field at time slot 1. */
  Kokkos::View<double**>* fieldInPtr2; /*! Input field at time slot 2. */

  FieldIO* fieldIO; /*! Pointer to FieldIO object for this forcing field. */

  Grid* grid; /*! Pointer to grid object for this forcing input. */

};

//------------------------------------------------------------------------------
/*! \class ForcingGroup
    \brief Class representing a single forcing input file template with multiple
    fields.

    This class represents a single forcing file template, e.g. 'forcing.$Y.nc',
    which may contain multiple fields and multiple forcing times. Forcing times
    must be separated by a constant time interval, as defined in the
    TimeInterval class, so that forcing times in the files are given by
    forcingTime = forcingReferenceTime + n*forcingInterval, where n is an
    integer. Filename template are expanded by the desired read time, so that
    'forcing.$Y.nc' would be expanded to 'forcing.2000.nc' is the desired read
    time was '2000-01-10_00:00:00'. Time token are '$Y' for year, '$M' for
    months, '$D' for days, '$h' for hours, '$m' for minutes, and  '$s' for
    seconds.
*/
//------------------------------------------------------------------------------
class ForcingGroup{
public:

  /*! \brief ForcingGroup constructor
      \param nameIn Identifying name of the forcing group.
      \param grid Pointer to the Grid object for this forcing group.
      \param filenameTemplateIn Input filename template for this forcing group.
      \param simulationClock Pointer to the simulation clock.
      \param forcingReferenceTimeIn Reference time for the forcing input data times.
      \param forcingIntervalIn Time interval between forcing input times.
      \param timeVaryingForcingIn If true forcing is time varying.
      \param logIn Pointer to log object.
   */
  ForcingGroup(std::string nameIn, Grid* grid, std::string filenameTemplateIn, DEMSI::Clock* simulationClock, const DEMSI::Time forcingReferenceTimeIn, const DEMSI::TimeInterval forcingIntervalIn, const bool timeVaryingForcingIn, DEMSI::Log* logIn);

  /*! \brief Destructor.*/
  ~ForcingGroup();

  /*! \brief Write the forcing group object to a stream.
      \param os Stream to write to.
      \param forcingGroup Forcing group object to write to stream.
  */
  friend std::ostream & operator<<(std::ostream & os, const ForcingGroup & forcingGroup);

  /*! \brief Add a forcing field to this forcing group.
      \param fieldName The name in the input netcdf file of the forcing field.
      \param fieldOut Pointer to the output forcing field.
  */
  void add_forcing_field(std::string fieldName, Kokkos::View<double**>* fieldOut);

  /*! \brief Read in the input data fields for all fields for the initial forcing times. */
  void init_field_data(void);

  /*! \brief Return the name of the forcing group.
      \return The name of the forcing group.
  */
  std::string get_name(void) const;

  /*! \brief Return the number of forcing fields in this forcing group.
      \return The number of forcing fields in this forcing group.
  */
  int get_field_number(void) const;

  /*! \brief Return a forcing field of the forcing group given an index.
      \param fieldIndex The index of the forcing field to return.
      \return The desired forcing field.
  */
  DEMSI::ForcingField get_field(int fieldIndex) const;

  /*! \brief Update all the output forcing fields for this forcing group. */
  void update(void);

private:

  std::string name; /*! Forcing group identifying name. */

  std::vector<ForcingField*> forcingFields; /*! Vector of forcing fields in this forcing group. */

  Clock* simulationClock; /*! Pointer to simulation clock. */
  Time forcingReferenceTime; /*! Forcing reference time for this group. */
  TimeInterval forcingInterval; /*! Interval between forcing input data. */

  Time forcingTimes[2]; /*! Current forcing data times before and after the current simulation time. */
  Alarm forcingAlarm; /*! Alarm signals when to read next input data time slice. */

  std::string filenameTemplate; /*! Filename template of this forcing group. */
  FieldIO* fieldIO; /*! FieldIO object for this forcing group. */

  bool timeVaryingForcing; /*! If true forcing is time varying. */

  Grid* grid; /*! Pointer to grid object for this forcing group. */

  Log* log; /*! Pointer to log. */

  std::vector<Kokkos::View<double**>*> forcingWriteFields; /*! Vector of field pointers to write out. */

};

//------------------------------------------------------------------------------
/*! \class ForcingInstance
    \brief Class representing all the specified forcing groups and their
    specifications.

    This class contains a vector of forcing groups, and member functions to
    create these. This object is created by the Forcing object using the desired
    forcing parameters specified in the Configs object.
*/
//------------------------------------------------------------------------------
class ForcingInstance{
public:

  /*! \brief ForcingInstance constructor
      \param grid Pointer to the Grid object for this forcing.
      \param logIn Pointer to log object.
   */
  ForcingInstance(Grid* grid, DEMSI::Log* logIn);

  /*! \brief Destructor.*/
  ~ForcingInstance();

  /*! \brief Write the forcing instance object to a stream.
      \param os Stream to write to.
      \param forcingInstance Forcing instance object to write to stream.
  */
  friend std::ostream & operator<<(std::ostream & os, const ForcingInstance & forcingInstance);

  /*! \brief Add a forcing group to the forcing instance.
      \param name The identifying name of the new forcing group.
      \param filenameTemplate The forcing filename template for the new forcing group.
      \param simulationClock Pointer to the simulation clock.
      \param forcingTimeReference Reference time for the forcing input data times for the new forcing group.
      \param forcingInterval Time interval between forcing input times for the new forcing group.
      \param timeVaryingForcing If true forcing is time varying.
  */
  void add_forcing_group(std::string name, std::string filenameTemplate, DEMSI::Clock* simulationClock, const DEMSI::Time forcingTimeReference, const DEMSI::TimeInterval forcingInterval, const bool timeVaryingForcing);

  /*! \brief Add a forcing field to a forcing group of the forcing instance.
      \param groupName The forcing group to add to.
      \param fieldName The netcdf input file fieldname to add.
      \param fieldOut Pointer to the output forcing field.
  */
  void add_forcing_field(std::string groupName, std::string fieldName, Kokkos::View<double**>* fieldOut);

  /*! \brief Read in the input data fields for all fields for the initial forcing times. */
  void init_field_data(void);

  /*! \brief Returns the number of forcing groups in this forcing instance.
      \return The number of forcing groups in this forcing instance.
   */
  int get_group_number(void) const;

  /*! \brief Return a forcing group of the forcing instance based on an index.
      \param groupIndex The index of the forcing group to return.
      \return The desired forcing group.
  */
  DEMSI::ForcingGroup get_group(int groupIndex) const;

  /*! \brief Update all forcing fields of the forcing instance based on the current time. */
  void update(void);

private:

  Grid* grid; /*! Pointer to grid object for this forcing instance. */

  Log* log; /*! Pointer to log object. */

  std::vector<DEMSI::ForcingGroup*> forcingGroups; /*! Forcing groups belonging to the forcing instance. */

};

//------------------------------------------------------------------------------
/*! \class Forcing
    \brief Class representing input forcing for a given grid.

    This class contains the input forcing fields as public members. These can
    be updated to the current simulation time using the update method. The
    configuration of the forcing input is contained in the Configs object.
*/
//------------------------------------------------------------------------------
class Forcing{
public:

  /*! \brief Forcing constructor
      \param logIn Pointer to log object.
      \param configs Simulation configs object.
      \param gridIn Grid object for this forcing object.
      \param calendar Pointer to the Calendar object for this forcing object.
      \param clock Pointer to the simulation clock.
   */
  Forcing(DEMSI::Log* logIn, DEMSI::Configs *configs, DEMSI::Grid* gridIn, Calendar* calendar, Clock* clock);

  /*! \brief Destructor.*/
  ~Forcing(void);

  /*! \brief Update the input forcing fields for the current simulation time. */
  void update(void);

  // forcing atmospheric fields
  /*! Input atmospheric wind speed in the x direction (m/s).*/
  Kokkos::View<double**> xAtmWind;
  /*! Input atmospheric wind speed in the y direction (m/s).*/
  Kokkos::View<double**> yAtmWind;
  /*! Input downwelling long wave radiation (W/m^2).*/
  Kokkos::View<double**> longwaveFluxDown;
  /*! Input downwelling  shortwave radiation (W/m^2).*/
  Kokkos::View<double**> shortwaveFluxDown;
  /*! Input air temperature (K).*/
  Kokkos::View<double**> airTemperature;
  /*! Input specific humidity (kg/kg).*/
  Kokkos::View<double**> specificHumidity;
  /*! Input precipitation rate (m/s).*/
  Kokkos::View<double**> precipitationRate;
  /*! Input rainfall rate (m/s).*/
  Kokkos::View<double**> rainfallRate;
  /*! Input snowfallRate rate (m/s).*/
  Kokkos::View<double**> snowfallRate;
  /*! Input cloud fraction (-)*/
  Kokkos::View<double**> cloudFraction;

  // forcing oceanic fields
  /*! Input oceanic current speed in the x direction (m/s).*/
  Kokkos::View<double**> xOcnCurrents;
  /*! Input oceanic current speed in the y direction (m/s).*/
  Kokkos::View<double**> yOcnCurrents;
  /*! Input sea surface temperature (C) */
  Kokkos::View<double**> seaSurfaceTemperature;
  /*! Input sea surface salinity (PSU) */
  Kokkos::View<double**> seaSurfaceSalinity;
  /*! Input sea surface tilt in x direction */
  Kokkos::View<double**> seaSurfaceTiltU;
  /*! Input sea surface tilt in y direction */
  Kokkos::View<double**> seaSurfaceTiltV;
  /*! Input ocean mixed layer depth */
  Kokkos::View<double**> oceanMixedLayerDepth;
  /*! Input ocean heat flux convergence */
  Kokkos::View<double**> oceanHeatFluxConvergence;

private:

  /*! Pointer to log object. */
  DEMSI::Log* log;

  /*! Pointer to grid object for this forcing object. */
  DEMSI::Grid* grid;

  /*! Forcing instance object for this forcing object. */
  DEMSI::ForcingInstance *forcingInstance;

  /*! map of forcing field names and pointers to the variable */
  std::map<std::string, Kokkos::View<double**>*> forcingFields;

  /*! if true use forcing system, otherwise dont */
  bool useForcing;

};

} // namespace DEMSI

#endif /* FORCING_H_ */
