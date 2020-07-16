#include "demsi_time.h"

namespace DEMSI {

//------------------------------------------------------------------------------
// time based constants
//------------------------------------------------------------------------------

const long long int millisecondsInSecond = 1000LL;
const long long int millisecondsInMinute = 60000LL;
const long long int millisecondsInHour   = 3600000LL;
const long long int millisecondsInDay    = 86400000LL;

//------------------------------------------------------------------------------
// utils
//------------------------------------------------------------------------------

// convert a string to an int, ensuring all string is consumed
int str_to_int(std::string intStr, DEMSI::Log* log)
{
  std::size_t size;
  int intOut = std::stoi(intStr, &size);
  log->check(size == intStr.length(), "Error: Time str_to_int failed to convert string to int: ", intStr);
  return intOut;
}

// split a string into tokens with another string as the delimiter
std::vector<std::string> split_string(const std::string strInput, const std::string delimiter)
{
  std::vector<std::string> list;
  std::string s = strInput;
  size_t pos = 0;
  std::string token;
  while ((pos = s.find(delimiter)) != std::string::npos) {
    token = s.substr(0, pos);
    list.push_back(token);
    s.erase(0, pos + delimiter.length());
  }
  list.push_back(s);
  return list;
}

// divide two long long ints and return a double
double division_long_long_int(const long long int numerator, const long long int denominator)
{
  long long int remainder      = std::abs(numerator);
  long long int denominatorAbs = std::abs(denominator);

  long long int result = 0;
  while (remainder >= denominatorAbs) {
    result = result + 1LL;
    remainder = remainder - denominatorAbs;
  }
  double resultDouble = (double) result + (double) remainder / (double) denominatorAbs;

  if (numerator < 0LL xor denominator < 0LL) {
    return -resultDouble;
  } else {
    return resultDouble;
  }
}

//------------------------------------------------------------------------------
// Calendar class
//------------------------------------------------------------------------------

// constructor
Calendar::Calendar(const std::string calendarTypeStr, DEMSI::Log* logIn)
{
  log = logIn;

  if (calendarTypeStr == "360_day") {
    calendarType = Calendar::c360_DAY;
  } else if (calendarTypeStr == "noleap" or calendarTypeStr == "365_day") {
    calendarType = Calendar::cNOLEAP;
  } else if (calendarTypeStr == "proleptic_gregorian") {
    calendarType = Calendar::cPROLEPTIC_GREGORIAN;
  } else {
    log->abort("Calendar: Unknown calendar input type: ", calendarTypeStr);
  }
}

// destructor
Calendar::~Calendar() {};

// Return the type of calendar
std::string Calendar::type(void) const
{
  if (calendarType == Calendar::c360_DAY) {
    return (std::string) "360_day";
  } else if (calendarType == Calendar::cNOLEAP) {
    return (std::string) "noleap/365_day";
  } else if (calendarType == Calendar::cPROLEPTIC_GREGORIAN) {
    return (std::string) "proleptic_gregorian";
  } else {
    log->abort("Calendar: Unknown calendar type: ", std::to_string(calendarType));
  }
}

bool Calendar::is_leap_year(const int year) const
{
  if (calendarType == Calendar::cPROLEPTIC_GREGORIAN) {
    if (year % 4 == 0) {
      // year is divisible by 4
      if (year % 100 == 0) {
	// year is divisible by 100
	if (year % 400 == 0) {
	  // year is divisible by 400
	  return true; // leap year
	} else {
	  return false; // common year
	}
      } else {
	return true; // leap year
      }
    } else {
      return false; // common year
    }
  } else {
    return false; // always common year
  }
}

int Calendar::cumul_prev_days_in_month(const int year, const int month) const
{

  if (calendarType == Calendar::c360_DAY) {
    return daysInYearPrev360Days[month];
  } else if (calendarType == Calendar::cNOLEAP) {
    return daysInYearPrevGregorian[month];
  } else {
    // proleptic_gregorian
    return daysInYearPrevGregorian[month] + (int) this->is_leap_year(year) * daysInYearPrevGregorianLeap[month];
  }
}

int Calendar::days_in_year(const int year) const
{
  if (calendarType == Calendar::c360_DAY) {
    return 360;
  } else if (calendarType == Calendar::cNOLEAP) {
    return 365;
  } else {
    if (this->is_leap_year(year)) {
      return 366;
    } else {
      return 365;
    }
  }
}

// over load << operator
std::ostream & operator<<(std::ostream & os, const Calendar & calendar)
{
  os << "Calendar: {" << calendar.type() << "}";
  return os;
}

//------------------------------------------------------------------------------
// time utility functions
//------------------------------------------------------------------------------

// get the months and days of month from a days of year
void get_months_and_days_from_days_of_year(const Calendar calendar, const int daysOfYear, const int years, int &months, int &days)
{
  // month
  for (int iMonth = 0 ; iMonth < 12 ; iMonth++) {
    if (daysOfYear >= calendar.cumul_prev_days_in_month(years,iMonth) and
	daysOfYear <  calendar.cumul_prev_days_in_month(years,iMonth+1)) {
      months = iMonth;
      break;
    }
  }

  // day
  days = daysOfYear - calendar.cumul_prev_days_in_month(years,months);
}

// get the time components from an input string
void get_time_components_from_string(const std::string timeStr, int &years, int &months, int &days, int &hours, int &minutes, int &seconds, int &milliseconds, DEMSI::Log* log)
{
  // accepts strings of formats "YYYY-MM-DD_hh:mm:ss" or "YYYY-MM-DD_hh:mm:ss_xxxx"
  // check input string length
  log->check(timeStr.length() == 19 or timeStr.length() == 24, "Error: timeStr wrong length in get_time_components_from_string");

  std::string yearsStr   = timeStr.substr(0, 4);
  std::string monthsStr  = timeStr.substr(5, 2);
  std::string daysStr    = timeStr.substr(8, 2);
  std::string hoursStr   = timeStr.substr(11, 2);
  std::string minutesStr = timeStr.substr(14, 2);
  std::string secondsStr = timeStr.substr(17, 2);

  years   = str_to_int(yearsStr,   log);
  months  = str_to_int(monthsStr,  log);
  days    = str_to_int(daysStr,    log);
  hours   = str_to_int(hoursStr,   log);
  minutes = str_to_int(minutesStr, log);
  seconds = str_to_int(secondsStr, log);

  if (timeStr.length() == 24) {
    std::string millisecondsStr = timeStr.substr(20, 4);
    milliseconds = str_to_int(millisecondsStr, log);
  } else {
    milliseconds = 0;
  }

  // store everything from index 0
  months = months - 1;
  days   = days   - 1;
}

// from time components get the time counter
long long int get_time_from_components(const Calendar calendar, const int years, const int months, const int days, const int hours, const int minutes, const int seconds, const int milliseconds)
{
  long long int integerTime = 0;

  // seconds in preceding years since reference time
  for (int iYear = 0 ; iYear < years ; iYear++) {
    integerTime = integerTime + (long long) calendar.days_in_year(iYear) * millisecondsInDay;
  }

  // seconds in preceding months
  integerTime = integerTime + (long long) calendar.cumul_prev_days_in_month(years,months) * millisecondsInDay;

  // remaining time
  integerTime = integerTime +
    (long long) days         * millisecondsInDay    +
    (long long) hours        * millisecondsInHour   +
    (long long) minutes      * millisecondsInMinute +
    (long long) seconds      * millisecondsInSecond +
    (long long) milliseconds;

  return integerTime;
}

// convert the time counter to time components
void get_components_from_time(const Calendar calendar, const long long int integerTime, int &years, int &months, int &days, int &hours, int &minutes, int &seconds, int &milliseconds)
{
  long long int integerTimeRemainder = integerTime;

  // year
  years = 0;
  if (integerTimeRemainder > 0) {
    while (integerTimeRemainder >= (long long) calendar.days_in_year(years) * millisecondsInDay)
      {
	integerTimeRemainder = integerTimeRemainder - (long long) calendar.days_in_year(years) * millisecondsInDay;
	years = years + 1;
      }
  } else {
    while (integerTimeRemainder < 0) {
      years = years - 1;
      integerTimeRemainder = integerTimeRemainder + (long long) calendar.days_in_year(years) * millisecondsInDay;
    }
  }

  // days
  int daysTmp = int(integerTimeRemainder / millisecondsInDay);
  integerTimeRemainder = integerTimeRemainder - (long long) daysTmp * millisecondsInDay;

  // hours
  hours = int(integerTimeRemainder / millisecondsInHour);
  integerTimeRemainder = integerTimeRemainder - (long long) hours   * millisecondsInHour;

  // minutes
  minutes = int(integerTimeRemainder / millisecondsInMinute);
  integerTimeRemainder = integerTimeRemainder - (long long) minutes * millisecondsInMinute;

  // seconds
  seconds = int(integerTimeRemainder / millisecondsInSecond);
  integerTimeRemainder = integerTimeRemainder - (long long) seconds * millisecondsInSecond;

  // milliseconds
  milliseconds = int(integerTimeRemainder);

  // month and day
  get_months_and_days_from_days_of_year(calendar, daysTmp, years, months, days);
}

//------------------------------------------------------------------------------
// Time interval class
//------------------------------------------------------------------------------

TimeInterval::interval_type get_interval_type_from_str(const std::string intervalTypeStr, DEMSI::Log* log)
{
  if (intervalTypeStr == "YEARS") {
    return TimeInterval::interval_type::YEARS;
  } else if (intervalTypeStr == "MONTHS") {
    return TimeInterval::interval_type::MONTHS;
  } else if (intervalTypeStr == "DAYS") {
    return TimeInterval::interval_type::DAYS;
  } else if (intervalTypeStr == "HOURS") {
    return TimeInterval::interval_type::HOURS;
  } else if (intervalTypeStr == "MINUTES") {
    return TimeInterval::interval_type::MINUTES;
  } else if (intervalTypeStr == "SECONDS") {
    return TimeInterval::interval_type::SECONDS;
  } else if (intervalTypeStr == "MILLISECONDS") {
    return TimeInterval::interval_type::MILLISECONDS;
  } else {
    log->abort("Error: Unknown interval type in get_interval_type_from_str: ", intervalTypeStr);
  }
}

// constructor
TimeInterval::TimeInterval() {};

TimeInterval::TimeInterval(const int intervalNumberIn, const interval_type intervalTypeIn, DEMSI::Log* logIn)
{
  log = logIn;
  intervalType = intervalTypeIn;
  intervalNumber = intervalNumberIn;
}

TimeInterval::TimeInterval(const int intervalNumberIn, const std::string intervalTypeStr, DEMSI::Log* logIn)
{
  log = logIn;
  intervalType = get_interval_type_from_str(intervalTypeStr, log);
  intervalNumber = intervalNumberIn;
}

TimeInterval::TimeInterval(const std::string intervalStr, DEMSI::Log* logIn)
{
  log = logIn;

  // read in an interval in string form e.g. "DAYS:20"

  // read in string to string stream
  std::istringstream intervalStrStream(intervalStr);
  log->check(!intervalStrStream.fail(), "Error: problem reading string in TimeInterval string constructor: ", intervalStr);

  // get the interval type from the string stream
  std::string intervalTypeStr;
  std::getline(intervalStrStream, intervalTypeStr, ':');
  intervalType = get_interval_type_from_str(intervalTypeStr, log);

  // get the interval number from the string stream
  std::string intervalNumberStr;
  std::getline(intervalStrStream, intervalNumberStr);
  intervalNumber = str_to_int(intervalNumberStr, log);
}

TimeInterval::~TimeInterval() {};

// over load << operator
std::ostream & operator<<(std::ostream & os, const TimeInterval & timeInterval)
{
  os << "TimeInterval: {" << timeInterval.get_interval_number() << " " << timeInterval.get_interval_type_str() << "}";
  return os;
}

// operators
TimeInterval TimeInterval::operator*(const int &number) const
{
  TimeInterval timeInterval = *this;
  timeInterval.intervalNumber = timeInterval.intervalNumber * number;
  return timeInterval;
}

// member functions
TimeInterval::interval_type TimeInterval::get_interval_type() const
{
  return intervalType;
}

// Return the number of units of the interval type
int TimeInterval::get_interval_number() const
{
  return intervalNumber;
}

// Return the interval type as a string
std::string TimeInterval::get_interval_type_str() const
{
  if (intervalType == interval_type::YEARS) {
    return "YEARS";
  } else if (intervalType == interval_type::MONTHS) {
    return "MONTHS";
  } else if (intervalType == interval_type::DAYS) {
    return "DAYS";
  } else if (intervalType == interval_type::HOURS) {
    return "HOURS";
  } else if (intervalType == interval_type::MINUTES) {
    return "MINUTES";
  } else if (intervalType == interval_type::SECONDS) {
    return "SECONDS";
  } else if (intervalType == interval_type::MILLISECONDS) {
    return "MILLISECONDS";
  } else {
    log->abort("Error: unknown interval type in get_interval_type_str");
  }
}

// Get the interval size in seconds
double TimeInterval::seconds() const {
  if (intervalType == interval_type::DAYS) {
    return (double) (intervalNumber * (millisecondsInDay / millisecondsInDay));
  } else if (intervalType == interval_type::HOURS) {
    return (double) (intervalNumber * (millisecondsInHour / millisecondsInDay));
  } else if (intervalType == interval_type::MINUTES) {
    return (double) (intervalNumber * (millisecondsInMinute / millisecondsInDay));
  } else if (intervalType == interval_type::SECONDS) {
    return (double) intervalNumber;
  } else {
    log->abort("Error: conversion to double seconds not supported for ", this->get_interval_type_str(), "interval type");
  }
}

//------------------------------------------------------------------------------
// Time class
//------------------------------------------------------------------------------

// constructors
Time::Time() {};

Time::Time(const Calendar* calendarIn, DEMSI::Log* logIn)
{
  calendar = calendarIn;
  log = logIn;
}

Time::Time(const Calendar* calendarIn, const std::string timeStr, DEMSI::Log* logIn)
{
  calendar = calendarIn;
  log = logIn;

  int years, months, days, hours, minutes, seconds, milliseconds;
  get_time_components_from_string(timeStr, years, months, days, hours, minutes, seconds, milliseconds, log);

  timeCounter = get_time_from_components(*calendar, years, months, days, hours, minutes, seconds, milliseconds);
}

  Time::Time(const Calendar* calendarIn, const int years, const int months, const int days, const int hours, const int minutes, const int seconds, const int milliseconds, DEMSI::Log* logIn)
{
  calendar = calendarIn;
  log = logIn;

  timeCounter = get_time_from_components(*calendar, years, months, days, hours, minutes, seconds, milliseconds);
}

// destructor
Time::~Time() {};

// over load << operator
std::ostream & operator<<(std::ostream & os, const Time & time)
{
  os << "Time: {" << time.get() << "}";
  return os;
}

// arithmetic

// subtraction
long long int Time::operator-(const Time & time) const
{
  long long int timeCounterDiff = timeCounter - time.get_time_counter();

  return timeCounterDiff;
}

// addition
Time Time::operator+(const TimeInterval & timeInterval) const
{
  Time time = *this;
  time.add_interval(timeInterval);
  return time;
}

// subtraction
Time Time::operator-(const TimeInterval & timeInterval) const
{
  Time time = *this;
  time.subtract_interval(timeInterval);
  return time;
}

// comparison
bool Time::operator==(const Time & time) const
{
  return timeCounter == time.timeCounter;
}
bool Time::operator<(const Time & time) const
{
  return timeCounter < time.timeCounter;
}
bool Time::operator>(const Time & time) const
{
  return timeCounter > time.timeCounter;
}
bool Time::operator<=(const Time & time) const
{
  return timeCounter <= time.timeCounter;
}
bool Time::operator>=(const Time & time) const
{
  return timeCounter >= time.timeCounter;
}

// add interval
void Time::add_interval(const TimeInterval timeInterval)
{
  if (timeInterval.get_interval_type() == TimeInterval::DAYS) {
    timeCounter += ((long long) timeInterval.get_interval_number() *
		    millisecondsInDay);
  } else if (timeInterval.get_interval_type() == TimeInterval::HOURS) {
    timeCounter += ((long long) timeInterval.get_interval_number() *
		    millisecondsInHour);
  } else if (timeInterval.get_interval_type() == TimeInterval::MINUTES) {
    timeCounter += ((long long) timeInterval.get_interval_number() *
		    millisecondsInMinute);
  } else if (timeInterval.get_interval_type() == TimeInterval::SECONDS) {
    timeCounter += ((long long) timeInterval.get_interval_number() *
		    millisecondsInSecond);
  } else if (timeInterval.get_interval_type() == TimeInterval::MILLISECONDS) {
    timeCounter +=  (long long) timeInterval.get_interval_number();
  } else if (timeInterval.get_interval_type() == TimeInterval::MONTHS) {

    int years, months, days, hours, minutes, seconds, milliseconds;
    get_components_from_time(*calendar, timeCounter, years, months, days, hours, minutes, seconds, milliseconds);

    months = months + timeInterval.get_interval_number();

    while(months > 11) {
      months -= 12;
      years += 1;
    }
    while(months < 0) {
      months += 12;
      years -= 1;
    }

    timeCounter = get_time_from_components(*calendar, years, months, days, hours, minutes, seconds, milliseconds);
  } else if (timeInterval.get_interval_type() == TimeInterval::YEARS) {
    int years, months, days, hours, minutes, seconds, milliseconds;
    get_components_from_time(*calendar, timeCounter, years, months, days, hours, minutes, seconds, milliseconds);
    years += timeInterval.get_interval_number();
    timeCounter = get_time_from_components(*calendar, years, months, days, hours, minutes, seconds, milliseconds);
  }
}

// subtract interval
void Time::subtract_interval(const TimeInterval timeInterval)
{
  if (timeInterval.get_interval_type() == TimeInterval::DAYS) {
    timeCounter -= ((long long) timeInterval.get_interval_number() *
		    millisecondsInDay);
  } else if (timeInterval.get_interval_type() == TimeInterval::HOURS) {
    timeCounter -= ((long long) timeInterval.get_interval_number() *
		    millisecondsInHour);
  } else if (timeInterval.get_interval_type() == TimeInterval::MINUTES) {
    timeCounter -= ((long long) timeInterval.get_interval_number() *
		    millisecondsInMinute);
  } else if (timeInterval.get_interval_type() == TimeInterval::SECONDS) {
    timeCounter -= ((long long) timeInterval.get_interval_number() *
		    millisecondsInSecond);
  } else if (timeInterval.get_interval_type() == TimeInterval::MILLISECONDS) {
    timeCounter -=  (long long) timeInterval.get_interval_number();
  } else if (timeInterval.get_interval_type() == TimeInterval::MONTHS) {

    int years, months, days, hours, minutes, seconds, milliseconds;
    get_components_from_time(*calendar, timeCounter, years, months, days, hours, minutes, seconds, milliseconds);

    months = months - timeInterval.get_interval_number();

    while(months > 11) {
      months -= 12;
      years += 1;
    }
    while(months < 0) {
      months += 12;
      years -= 1;
    }

    timeCounter = get_time_from_components(*calendar, years, months, days, hours, minutes, seconds, milliseconds);
  } else if (timeInterval.get_interval_type() == TimeInterval::YEARS) {
    int years, months, days, hours, minutes, seconds, milliseconds;
    get_components_from_time(*calendar, timeCounter, years, months, days, hours, minutes, seconds, milliseconds);
    years -= timeInterval.get_interval_number();
    timeCounter = get_time_from_components(*calendar, years, months, days, hours, minutes, seconds, milliseconds);
  }
}

// return time in string format
std::string Time::get() const
{
  int years, months, days, hours, minutes, seconds, milliseconds;

  get_components_from_time(*calendar, timeCounter, years, months, days, hours, minutes, seconds, milliseconds);

  months = months + 1;
  days   = days   + 1;

  char timeStr[25];
  if (milliseconds == 0) {
    sprintf(timeStr, "%04i-%02i-%02i_%02i:%02i:%02i", years, months, days, hours, minutes, seconds);
  } else {
    sprintf(timeStr, "%04i-%02i-%02i_%02i:%02i:%02i_%04i", years, months, days, hours, minutes, seconds, milliseconds);
  }

  return (std::string) timeStr;
}

// return a string in char[64] format
char* Time::get_char64() const {

  char* strOut = new char[64];

  std::string timeStr = this->get();
  for (int i=0 ; i < timeStr.length() ; i++) {
    strOut[i] = timeStr[i];
  }
  strOut[timeStr.length()] = '\0';

  return strOut;

}

// set the time
void Time::set(std::string timeStr)
{
  int years, months, days, hours, minutes, seconds, milliseconds;
  get_time_components_from_string(timeStr, years, months, days, hours, minutes, seconds, milliseconds, log);

  timeCounter = get_time_from_components(*calendar, years, months, days, hours, minutes, seconds, milliseconds);
}

// return time in cf string format
double Time::get_cf_value(const Time referenceTime, const TimeInterval::interval_type intervalType) const
{
  long long int timeDifference = this->timeCounter - referenceTime.timeCounter;

  if (intervalType == TimeInterval::interval_type::DAYS) {
    return division_long_long_int(timeDifference, millisecondsInDay);
  } else if (intervalType == TimeInterval::interval_type::HOURS) {
    return division_long_long_int(timeDifference, millisecondsInHour);
  } else if (intervalType == TimeInterval::interval_type::MINUTES) {
    return division_long_long_int(timeDifference, millisecondsInMinute);
  } else if (intervalType == TimeInterval::interval_type::SECONDS) {
    return division_long_long_int(timeDifference, millisecondsInDay);
  } else if (intervalType == TimeInterval::interval_type::MILLISECONDS) {
    return division_long_long_int(timeDifference, 1LL);
  } else {
    log->abort("Error: Incompatable interval type with cf format.");
  }
}

// return time in cf string format
std::string Time::get_cf_units(const Time referenceTime, const TimeInterval::interval_type intervalType) const
{
  std::string intervalTypeStr;
  if (intervalType == TimeInterval::interval_type::DAYS) {
    intervalTypeStr = "days";
  } else if (intervalType == TimeInterval::interval_type::HOURS) {
    intervalTypeStr = "hours";
  } else if (intervalType == TimeInterval::interval_type::MINUTES) {
    intervalTypeStr = "minutes";
  } else if (intervalType == TimeInterval::interval_type::SECONDS) {
    intervalTypeStr = "seconds";
  } else if (intervalType == TimeInterval::interval_type::MILLISECONDS) {
    intervalTypeStr = "milliseconds";
  } else {
    log->abort("Error: Incompatable interval type with cf format.");
  }

  std::string timeStr = intervalTypeStr + " since " + referenceTime.get();

  return timeStr;
}

// set the time from a string in cf format
void Time::set_cf(double timeDouble, std::string units)
{
  // parse the units string
  std::vector<std::string> unitsSplit = split_string(units, " since ");
  log->check(unitsSplit.size() == 2, "Error: cf units cant be split: ", units);
  std::string unitsType        = unitsSplit[0];
  std::string referenceTimeStr = unitsSplit[1];
  std::transform(unitsType.begin(), unitsType.end(), unitsType.begin(), ::toupper);

  // set the time to the reference time
  int years, months, days, hours, minutes, seconds, milliseconds;
  get_time_components_from_string(referenceTimeStr, years, months, days, hours, minutes, seconds, milliseconds, log);

  timeCounter = get_time_from_components(*calendar, years, months, days, hours, minutes, seconds, milliseconds);

  long long int conversionFactor;
  if (unitsType == "DAYS") {
    conversionFactor = millisecondsInDay;
  } else if (unitsType == "HOURS") {
    conversionFactor = millisecondsInHour;
  } else if (unitsType == "MINUTES") {
    conversionFactor = millisecondsInMinute;
  } else if (unitsType == "SECONDS") {
    conversionFactor = millisecondsInSecond;
  } else if (unitsType == "MILLISECONDS") {
    conversionFactor = 1LL;
  } else {
    log->abort("Error: cannot use units type in Time.set_cf");
  }

  // set milliseconds from whole part of input double
  double timeDoubleRound;
  double remainder = std::modf(timeDouble, &timeDoubleRound);
  timeCounter = timeCounter + conversionFactor * (long long int) timeDoubleRound;

  // get remainder
  double remainderWhole;
  remainder = remainder * (double) conversionFactor;
  if (std::modf(remainder, &remainderWhole) != 0.0) {
    log->abort("Error: Not a round number of milliseconds in cf time.");
  }
  timeCounter = timeCounter + (long long int) remainder;

}

// Return the individual time components of the time.
void Time::get_components(int &years, int &months, int &days, int &hours, int &minutes, int &seconds, int &milliseconds) const {
  get_components_from_time(*calendar, timeCounter, years, months, days, hours, minutes, seconds, milliseconds);
}

// Return a pointer to the calendar used by the time.
const DEMSI::Calendar* Time::get_calendar(void) const {
  return calendar;
}

// Return the time counter value of the time.
long long int Time::get_time_counter(void) const {
  return timeCounter;
}

// Return the time converted into a climatological time by setting the year to 0.
Time Time::climatology(void) const {

  int years, months, days, hours, minutes, seconds, milliseconds;
  this->get_components(years, months, days, hours, minutes, seconds, milliseconds);
  years = 0; // climatological year is zero

  Time climatologicalTime(this->get_calendar(), years, months, days, hours, minutes, seconds, milliseconds, log);
  return climatologicalTime;

}

// Get the current day of year
int Time::day_of_year(void) {
  int years, months, days, hours, minutes, seconds, milliseconds;
  this->get_components(years, months, days, hours, minutes, seconds, milliseconds);
  return days + 1 + calendar->cumul_prev_days_in_month(years, months);
}

// Get the number of days in the current year
int Time::days_in_year(void) {
  int years, months, days, hours, minutes, seconds, milliseconds;
  this->get_components(years, months, days, hours, minutes, seconds, milliseconds);
  return calendar->days_in_year(years);
}

// Get the number of seconds into the current day
int Time::seconds_into_day(void) {
  int years, months, days, hours, minutes, seconds, milliseconds;
  this->get_components(years, months, days, hours, minutes, seconds, milliseconds);
  return seconds + minutes * (millisecondsInMinute / millisecondsInSecond) + hours * (millisecondsInHour / millisecondsInSecond);
}

//------------------------------------------------------------------------------
// Clock class
//------------------------------------------------------------------------------

// constructors
Clock::Clock(const Time startTime, const TimeInterval clockIntervalIn)
{
  clockTime = startTime;
  clockInterval = clockIntervalIn;
}

// destructor
Clock::~Clock() {};

// over load << operator
std::ostream & operator<<(std::ostream & os, const Clock & clock)
{
  os << "Clock: {" << clock.clockTime << " " << clock.clockInterval << "}";
  return os;
}

void Clock::advance()
{
  clockTime.add_interval(clockInterval);
}

Time Clock::get_time()
{
  return clockTime;
}

//------------------------------------------------------------------------------
// Alarm class
//------------------------------------------------------------------------------

// default constructor
Alarm::Alarm() {};

// constructor
Alarm::Alarm(Clock* clock, Time startTime, TimeInterval alarmIntervalIn)
{
  alarmClock = clock;

  alarmTime = startTime;

  alarmInterval = alarmIntervalIn;
}

// destructor
Alarm::~Alarm() {};

// over load << operator
std::ostream & operator<<(std::ostream & os, const Alarm & alarm)
{
  os << "Alarm: {" << alarm.get_alarm_time() << " {" << alarm.get_alarm_clock() << "} " << alarm.get_alarm_interval() << "}";
  return os;
}

void Alarm::set(Clock* clock, Time startTime, TimeInterval alarmIntervalIn)
{
  alarmClock = clock;

  alarmTime = startTime;

  alarmInterval = alarmIntervalIn;
}

void Alarm::reset()
{
  while (alarmTime <= alarmClock->get_time()) {
    alarmTime.add_interval(alarmInterval);
  }
}

bool Alarm::is_ringing() const
{
  if (alarmClock->get_time() >= alarmTime) {
    return true;
  } else {
    return false;
  }
}

Time Alarm::get_alarm_time() const
{
  return alarmTime;
}

Clock Alarm::get_alarm_clock() const
{
  return (*alarmClock);
}

TimeInterval Alarm::get_alarm_interval() const
{
  return alarmInterval;
}

} // namespace DEMSI
