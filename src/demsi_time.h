#ifndef DEMSI_TIME_H_
#define DEMSI_TIME_H_

#include <string>
#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
#include <algorithm>
#include <vector>

#include "demsi_logging.h"

/*! \file   demsi_time.h
    \brief  Header file for the DEMSI::Time classes
*/

namespace DEMSI {

//------------------------------------------------------------------------------
/*! \class Calendar
    \brief Class defining the calendar used by model times.

    The class specifies the calendar used by the other time classes. Calendar
    types defined are 'c360_DAY', 'cGREGORIAN_NO_LEAP', and 'cGREGORIAN'.
    - 'c360_DAY': 360 day years each with twelve thirty days months.
    - 'cNOLEAP': 365 day years with normal month lengths and no leap
    years.
    - 'cPROLEPTIC_GREGORIAN': Standard Gregorian calendar with leap years. This
    calendar is proleptic and continues with leap years before the start of the
    Gregorian calendar.
*/
//------------------------------------------------------------------------------
class Calendar
{
 private:

  enum {c360_DAY, cNOLEAP, cPROLEPTIC_GREGORIAN} calendarType;

  // cumulated days in year before month
  int daysInYearPrev360Days  [13] = {0, 30, 60, 90,120,150,180,210,240,270,300,330,360};
  int daysInYearPrevGregorian[13] = {0, 31, 59, 90,120,151,181,212,243,273,304,334,365};

  int daysInYearPrevGregorianLeap[13] = {0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

  /*! Pointer to the log object */
  Log* log;

 public:

  /*! \brief Calendar constructor
      \param calendarType The calendar type to create.
      \param log Pointer to the logging object.
      This can be "360_day" (for c360_DAY), "noleap" or "365_day" (for cNOLEAP),
      or "proleptic_gregorian" (for cPROLEPTIC_GREGORIAN).
   */
  Calendar(const std::string calendarType, DEMSI::Log* log);

  /*! \brief Default destructor.*/
  ~Calendar();

  /*! \brief Return the type of calendar
      \return The calendar type as a string.
  */
  std::string type(void) const;

  /*! \brief Determine if a given year in the calendar is a leap year.
      \param year The year to test if it is a leap year.
   */
  bool is_leap_year(const int year) const;

  /*! \brief Return the number of days in a year before a given month.
      \param year Given year to find the number of days of.
      \param month The month for which the number of days in the year before it
      is given.
   */
  int cumul_prev_days_in_month(const int year, const int month) const;

  /*! \brief Return the number of days in a given year for a calendar.
      \param year The year to determine the number of days of.
   */
  int days_in_year(const int year) const;

  /*! \brief Write the calendar object to a stream.
      \param os Stream to write to.
      \param calendar Calendar object to write to stream.

      The calendar is written to a stream in the format:
      "Calendar: {std::string calendarType}".
  */
  friend std::ostream & operator<<(std::ostream & os, const Calendar & calendar);

};

//------------------------------------------------------------------------------
/*! \class TimeInterval
    \brief Class defining the difference in time between two times.

    This class defines the difference in time between two times. This is stored
    as a unit of interval ('YEARS', 'MONTHS', 'DAYS', 'HOURS', 'MINUTES',
    'SECONDS', or 'MILLISECONDS') and a number of that interval.
*/
//------------------------------------------------------------------------------
class TimeInterval
{
 public:

  /*! \enum interval_type
      \brief Allowed values of the interval type.
   */
  enum interval_type {YEARS, MONTHS, DAYS, HOURS, MINUTES, SECONDS, MILLISECONDS};

 private:

  TimeInterval::interval_type intervalType;

  int intervalNumber;

  DEMSI::Log* log;

 public:

  /*! \brief Default constructor.*/
  TimeInterval();

  /*! \brief TimeInterval constructor
      \param intervalNumberIn The number of interval units.
      \param intervalTypeIn The interval type as an interval_type.
      \param logIn Pointer to the logging object.
   */
  TimeInterval(const int intervalNumberIn, const interval_type intervalTypeIn, DEMSI::Log* logIn);

  /*! \brief TimeInterval constructor
      \param intervalNumberIn The number of interval units.
      \param intervalTypeStr The interval type as a string. Allowed values:
      - 'YEARS'
      - 'MONTHS'
      - 'DAYS'
      - 'HOURS'
      - 'MINUTES'
      - 'SECONDS'
      - 'MILLISECONDS'
      \param logIn Pointer to the logging object.
   */
  TimeInterval(const int intervalNumberIn, const std::string intervalTypeStr, DEMSI::Log* logIn);

  /*! \brief TimeInterval constructor
      \param intervalStr The number of interval units and interval type in the
      format 'XXXX:NNNN' where XXXX is the interval type, and NNNN is the integer
      number of those units.
      \param logIn Pointer to the logging object.
   */
  TimeInterval(const std::string intervalStr, DEMSI::Log* logIn);

  /*! \brief Default destructor.*/
  ~TimeInterval();

  /*! \brief Write the time interval object to a stream.
      \param os Stream to write to.
      \param timeInterval Time interval object to write to stream.

      The time interval is written to a stream in the format:
      "TimeInterval: {IntervalNumber IntervalType}".
  */
  friend std::ostream & operator<<(std::ostream & os, const TimeInterval & timeInterval);

  /*! \brief Multiply the interval number by a factor but keep the same units.
      \param number Multiplication factor of the interval number.
  */
  TimeInterval operator*(const int &number) const;

  /*! \brief Return the interval type of the time interval.*/
  TimeInterval::interval_type get_interval_type() const;

  /*! \brief Return the number of units of the interval type.*/
  int get_interval_number() const;

  /*! \brief Return the interval type as a string.*/
  std::string get_interval_type_str() const;

  /*! \brief Get the interval size in seconds
      \return The interval size in second
   */
  double seconds() const;

};

//------------------------------------------------------------------------------
/*! \class Time
    \brief Class defining a single date and time.

    This class defines a basic date and time type. Internally the time is stored
    as an integer number of milliseconds since '0000-01-01_00:00:00'. This class
    can handle years in the range ?? to ??.
*/
//------------------------------------------------------------------------------
class Time
{
 private:

  // milliseconds since year 0
  long long int timeCounter;

  // pointer to the calendar object used by this time
  const Calendar* calendar;

  // Pointer to the logging object.
  DEMSI::Log* log;

 public:

  /*! \brief Default constructor.*/
  Time();

  /*! \brief Time constructor
      \param calendar The calendar to use with this time.
      \param logIn Pointer to the logging object.
   */
  Time(const Calendar* calendar, DEMSI::Log* logIn);

  /*! \brief Time constructor
      \param calendar The calendar to use with this time.
      \param timeStr The initial time in format 'YYYY-MM-DD_hh:mm:ss' or
      'YYYY-MM-DD_hh:mm:ss_xxx' where xxx is the number of milliseconds.
      \param logIn Pointer to the logging object.
   */
  Time(const Calendar* calendar, const std::string timeStr, DEMSI::Log* logIn);

  /*! \brief Time constructor
      \param calendar The calendar to use with this time.
      \param years The initial number of years.
      \param months The initial number of months.
      \param days The initial number of days of the month.
      \param hours The initial number of hours.
      \param minutes The initial number of minutes.
      \param seconds The initial number of seconds.
      \param milliseconds The initial number of milliseconds.
      \param logIn Pointer to the logging object.
   */
  Time(const Calendar* calendar, const int years, const int months, const int days, const int hours, const int minutes, const int seconds, const int milliseconds, DEMSI::Log* logIn);

  /*! \brief Default destructor.*/
  ~Time();

  /*! \brief Write the time object to a stream.
      \param os Stream to write to.
      \param time Time object to write to stream.

      The time object is written to a stream in the format:
      "Time: {YYYY-MM-DD_hh:mm:ss}" or "Time: {YYYY-MM-DD_hh:mm:ss_xxx}" if
      milliseconds are present.
  */
  friend std::ostream & operator<<(std::ostream & os, const Time & time);

  /*! \brief Calculate the difference in milliseconds between two times.
      \param time The time object to subtract from the calling time object.
  */
  long long int operator-(const Time & time) const;

  /*! \brief Add an interval to the time.
      \param timeInterval The time interval to add to the time.
   */
  void add_interval(const TimeInterval timeInterval);

  /*! \brief Subtract an interval from the time.
      \param timeInterval The time interval to subtract from the time.
   */
  void subtract_interval(const TimeInterval timeInterval);

  /*! \brief Add a time interval to the calling time object.
      \param timeInterval The time interval to add to the time object.
  */
  Time operator+(const TimeInterval & timeInterval) const;

  /*! \brief Subtract a time interval from the calling time object.
      \param timeInterval The time interval to subtract from the time object.
  */
  Time operator-(const TimeInterval & timeInterval) const;

  /*! \brief Determine if the time object is same as another time.
      \param time The other time to test against.
  */
  bool operator==(const Time & time) const;

  /*! \brief Determine if the time object is before another time.
      \param time The other time to test against.
  */
  bool operator<(const Time & time) const;

  /*! \brief Determine if the time object is after another time.
      \param time The other time to test against.
  */
  bool operator>(const Time & time) const;

  /*! \brief Determine if the time object is before or equal to another time.
      \param time The other time to test against.
  */
  bool operator<=(const Time & time) const;

  /*! \brief Determine if the time object is after or equal to another time.
      \param time The other time to test against.
  */
  bool operator>=(const Time & time) const;

  /*! \brief Return the current time.

    Return the current time as a string in either 'YYYY-MM-DD_hh:mm:ss' format
    if milliseconds is zero, or 'YYYY-MM-DD_hh:mm:ss_xxx' otherwise where xxx
    is the number of milliseconds.
   */
  std::string get() const;

  /*! \brief Return the current time as a char[64].

    Return the current time as a char[64] in either 'YYYY-MM-DD_hh:mm:ss' format
    if milliseconds is zero, or 'YYYY-MM-DD_hh:mm:ss_xxx' otherwise where xxx
    is the number of milliseconds.
   */
  char* get_char64() const;

  /*! \brief Set the current time.

    Set the current time from a string in either 'YYYY-MM-DD_hh:mm:ss'
    format if milliseconds is zero, or 'YYYY-MM-DD_hh:mm:ss_xxx' otherwise,
    where xxx is the number of milliseconds.
   */
  void set(std::string timeStr);

  /*! \brief Return the current time value in cf format.

    Return the current time value in cf format,
    e.g. "days since 01-07-015 00:00:00"
   */
  double get_cf_value(const Time referenceTime, const TimeInterval::interval_type intervalType) const;

  /*! \brief Return the current time in cf format.

    Return the current time units in cf format,
    e.g. "days since 01-07-015 00:00:00"
   */
  std::string get_cf_units(const Time referenceTime, const TimeInterval::interval_type intervalType) const;

  /*! \brief Set the current time in cf format.

    Set the current time from a string in cf format,
    e.g. "days since 01-07-015 00:00:00"
   */
  void set_cf(double timeDouble, std::string units);

  /*! \brief Return the individual time components of the time. */
  void get_components(int &years, int &months, int &days, int &hours, int &minutes, int &seconds, int &milliseconds) const;

  /*! \brief Return a pointer to the calendar used by the time. */
  const DEMSI::Calendar* get_calendar(void) const;

  /*! \brief Return the time counter value of the time. */
  long long int get_time_counter(void) const;

  /*! \brief Return the time converted into a climatological time by setting the year to 0. */
  Time climatology(void) const;

  /*! \brief Return the current day of year
      \return The current day of year
  */
  int day_of_year(void);

  /*! \brief Get the number of days in the current year
      \return The number of days in the current year
  */
  int days_in_year(void);

  /*! \brief Get the number of seconds into the current day
      \return The number of seconds into the current day
  */
  int seconds_into_day(void);

};

//------------------------------------------------------------------------------
/*! \class Clock
    \brief Class defining a simulation clock.

    This class defines a clock type for the simulation to keep track of time.
    The class has a current time and a set interval between states of the clock
    time and a method to advance the clock by the clock interval.
*/
//------------------------------------------------------------------------------
class Clock
{
 private:

  // current time of the clock
  Time clockTime;

  // time interval that will be added to the clock time when clock time is advanced
  TimeInterval clockInterval;

 public:

  /*! \brief Clock constructor
      \param startTime The initial time of the clock.
      \param clockIntervalIn The time interval to advance the clock by.
   */
  Clock(const Time startTime, const TimeInterval clockIntervalIn);

  /*! \brief Default destructor.*/
  ~Clock();

  /*! \brief Write the clock object to a stream.
      \param os Stream to write to.
      \param clock Clock object to write to stream.

      The clock object is written to a stream in the format:
      "Clock: {clockTime clockInterval}".
  */
  friend std::ostream & operator<<(std::ostream & os, const Clock & clock);

  /*! \brief Advance the clock time by the clock interval.*/
  void advance();

  /*! \brief Return the current clock time.*/
  Time get_time();

};

//------------------------------------------------------------------------------
/*! \class Alarm
    \brief Class defining a simulation alarm.

    This class defines an alarm type that allows the simulation to determine if
    fixed intervals of time have passed. This class is useful for allowing the
    simulation to read data in and out on set time schedules.
*/
//------------------------------------------------------------------------------
class Alarm
{
 private:

  // pointer to the clock that the alarm is tied to
  Clock* alarmClock;

  // time the alarm will next ring
  Time alarmTime;

  // time interval that will be added to the alarm time when alarm reset
  TimeInterval alarmInterval;

 public:

  /*! \brief Default constructor.*/
  Alarm();

  /*! \brief Alarm constructor
      \param clock The simulation clock to associate with the alarm.
      \param startTime The initial internal alarm time that is compared
      against the clock time to determine if it is ringing.
      \param alarmIntervalIn The time interval to increase alarm time by
      when the alarm is reset. The alarm time is increased by multiple units
      of the alarm interval until the alarm time is greater than the
      alarm clock time.
   */
  Alarm(Clock* clock, Time startTime, TimeInterval alarmIntervalIn);

  /*! \brief Default destructor.*/
  ~Alarm();

  /*! \brief Write the alarm object to a stream.
      \param os Stream to write to.
      \param alarm Alarm object to write to stream.

      The alarm object is written to a stream in the format:
      "Alarm: {alarmTime {alarmClock} alarmInterval}".
  */
  friend std::ostream & operator<<(std::ostream & os, const Alarm & alarm);

  /*! \brief Set the alarm private members.
      \param clock The simulation clock to associate with the alarm.
      \param startTime The initial internal alarm time that is compared
      against the clock time to determine if it is ringing.
      \param alarmIntervalIn The time interval to increase alarm time by
      when the alarm is reset. The alarm time is increased by multiple units
      of the alarm interval until the alarm time is greater than the
      alarm clock time.
  */
  void set(Clock* clock, const Time startTime, const TimeInterval alarmIntervalIn);

  /*! \brief Reset the alarm.

    Increase the internal alarm time by multiples of the alarm interval until
    the alarm time is greater than the alarm clock time.
   */
  void reset();

  /*! \brief Check to see if the alarm is ringing.

    Determine if the clock time is greater or equal to the alarm time and if so
    return true.
   */
  bool is_ringing() const; // check if the alarm is ringing

  /*! \brief Return the alarm time as a Time object.*/
  Time get_alarm_time() const;

  /*! \brief Return the clock associated with the alarm.*/
  Clock get_alarm_clock() const;

  /*! \brief Return the alarm interval associated with the alarm.*/
  TimeInterval get_alarm_interval() const;

};

} // namespace DEMSI

#endif
