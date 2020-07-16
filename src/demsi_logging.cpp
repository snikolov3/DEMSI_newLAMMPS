#include "demsi_logging.h"

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <exception>
#include <stdexcept>
#include <chrono>
#include <ctime>

namespace DEMSI {

std::string get_log_filename(const int thisProc, const int nProcs, std::string postfix) {

  // create log filename
  int width;
  if (nProcs < 1E4) {
    width = 4;
  } else if (nProcs < 1E5) {
    width = 5;
  } else if (nProcs < 1E6) {
    width = 6;
  } else if (nProcs < 1E7) {
    width = 7;
  } else if (nProcs < 1E8) {
    width = 8;
  } else {
    width = 9;
  }

  std::stringstream logFilename;
  logFilename.fill('0');
  logFilename << "log.";
  logFilename << std::setfill('0') << std::setw(width) << thisProc;
  logFilename << ".";
  logFilename << postfix;

  return logFilename.str();

} // get_log_filename

std::string system_time(void) {
  auto currentTimeSystem = std::chrono::system_clock::now();
  std::time_t currentTime = std::chrono::system_clock::to_time_t(currentTimeSystem);
  return (std::string) std::ctime(&currentTime);
} // system_time

// Log class constructor.
Log::Log(DEMSI::Partition* partitionIn) {

  logType = DEMSI::LOG::NORMAL;

  partition = partitionIn;

  if (partition->on_master_proc()) {
    writeToLog = true;
    std::string logFilename = get_log_filename(partition->this_proc(), partition->nprocs(), (std::string) "out");
    logFile.open(logFilename);
  } else {
    writeToLog = false;
  }

} // Log::Log

// Log class destructor.
Log::~Log() {
  logFile << "Log file closed at " << system_time() << std::endl << std::flush;
  logFile.close();
  if (errFile.is_open()) {
    errFile << "Err file closed at " << system_time() << std::endl << std::flush;
    errFile.close();
  }
} // Log::~Log

// set up logging with configs
void Log::setup(const bool writeLogAllProcs, const bool writeDebugIn) {

  // set whether log to all procs or not
  if (writeLogAllProcs) {
    writeToLog = true;
    if (not logFile.is_open()) {
      std::string logFilename = get_log_filename(partition->this_proc(), partition->nprocs(), (std::string) "out");
      logFile.open(logFilename);
    }
  }

  // write log header
  if (writeToLog) {
    logFile << std::string(80, '-') << std::endl << std::endl;
    logFile << "  Log file for DEMSI for process " << partition->this_proc() << " of " << partition->nprocs() << " processes" << std::endl;
    logFile << "  Log file started at " << system_time() << std::endl;
    logFile << std::string(80, '-') << std::endl << std::endl << std::flush;
  }

  // whether writing debug messages
  writeDebug = writeDebugIn;

} // Log::setup

// Write the ostream object to the output streams.
void Log::write(std::ostream& (*pf) (std::ostream&)) {

  if (writeToLog and (logType != DEMSI::LOG::DEBUG or writeDebug)) logFile << pf << std::flush;
  if (partition->on_master_proc() and (logType != DEMSI::LOG::DEBUG or writeDebug)) std::cout << pf << std::flush;
  if (logType == DEMSI::LOG::ERROR or logType == DEMSI::LOG::CRITICAL) errFile << pf << std::flush;

} // Log::write

// Overloaded << operator, write ostream object to the output streams.
Log &Log::operator<<(std::ostream& (*pf) (std::ostream&)){
  this->write(pf);
  return *this;
} // Log::operator<<

// Overloaded , operator, write ostream object to the output streams.
Log &Log::operator,(std::ostream& (*pf) (std::ostream&)){
  this->write(pf);
  return *this;
} // Log::operator,

// Parenthesis operator, set the message log type.
Log &Log::operator()(const DEMSI::LOG logTypeIn) {
  logType = logTypeIn;
  if (logType == DEMSI::LOG::WARNING) {
    this->write("Warning: ");
  } else if (logType == DEMSI::LOG::DEBUG) {
    this->write("DEBUG: ");
  } else if (logType == DEMSI::LOG::ERROR or logType == DEMSI::LOG::CRITICAL) {
    if (not errFile.is_open()) {
      errFile.open(get_log_filename(partition->this_proc(), partition->nprocs(), (std::string) "err"));
    }
    if (logType == DEMSI::LOG::ERROR) {
      this->write("Error: ");
    } else {
      this->write("Critical: ");
    }
  }
  return *this;
} // Log::operator()

// Abort the model.
void Log::abort(const std::string message) {
  (*this)(DEMSI::LOG::CRITICAL) << message << std::endl;
  (*this) << "DEMSI aborted at " << system_time() << std::endl << std::flush;
  std::stringstream ss;
  ss << message;
  throw std::runtime_error(ss.str());
} // Log::abort

// Abort the model.
void Log::abort(const std::string message1, const std::string message2) {
  (*this)(DEMSI::LOG::CRITICAL) << message1 << " " << message2 << std::endl;
  (*this) << "DEMSI aborted at " << system_time() << std::endl << std::flush;
  std::stringstream ss;
  ss << message1 << " " << message2;
  throw std::runtime_error(ss.str());
};

// Abort the model.
void Log::abort(const std::string message1, const std::string message2, const std::string message3) {
  (*this)(DEMSI::LOG::CRITICAL) << message1 << " " << message2 << " " << message3 << std::endl;
  (*this) << "DEMSI aborted at " << system_time() << std::endl << std::flush;
  std::stringstream ss;
  ss << message1 << " " << message2 << " " << message3;
  throw std::runtime_error(ss.str());
} // Log::abort

// Abort the model.
void Log::abort(const std::string message1, const std::string message2, const std::string message3, const std::string message4) {
  (*this)(DEMSI::LOG::CRITICAL) << message1 << " " << message2 << " " << message3 << " " << message4 << std::endl;
  (*this) << "DEMSI aborted at " << system_time() << std::endl << std::flush;
  std::stringstream ss;
  ss << message1 << " " << message2 << " " << message3 << " " << message4;
  throw std::runtime_error(ss.str());
} // Log::abort

// Abort the model.
void Log::abort(const std::string message1, const std::string message2, const std::string message3, const std::string message4, const std::string message5) {
  (*this)(DEMSI::LOG::CRITICAL) << message1 << " " << message2 << " " << message3 << " " << message4 << " " << message5 << std::endl;
  (*this) << "DEMSI aborted at " << system_time() << std::endl << std::flush;
  std::stringstream ss;
  ss << message1 << " " << message2 << " " << message3 << " " << message4 << " " << message5;
  throw std::runtime_error(ss.str());
} // Log::abort

// Abort the model.
void Log::abort(const std::string message1, const std::string message2, const std::string message3, const std::string message4, const std::string message5, const std::string message6) {
  (*this)(DEMSI::LOG::CRITICAL) << message1 << " " << message2 << " " << message3 << " " << message4 << " " << message5 << " " << message6 << std::endl;
  (*this) << "DEMSI aborted at " << system_time() << std::endl << std::flush;
  std::stringstream ss;
  ss << message1 << " " << message2 << " " << message3 << " " << message4 << " " << message5 << " " << message6;
  throw std::runtime_error(ss.str());
} // Log::abort

// Assert expression is true, otherwise abort.
void Log::check(const bool expr, const std::string message) {
  if (not expr) {
    this->abort(message);
  }
} // Log::check

// Assert expression is true, otherwise abort.
void Log::check(const bool expr, const std::string message1, const std::string message2) {
  if (not expr) {
    this->abort(message1, message2);
  }
} // Log::check

// Assert expression is true, otherwise abort.
void Log::check(const bool expr, const std::string message1, const std::string message2, const std::string message3) {
  if (not expr) {
    this->abort(message1, message2, message3);
  }
} // Log::check

// Assert expression is true, otherwise abort.
void Log::check(const bool expr, const std::string message1, const std::string message2, const std::string message3, const std::string message4) {
  if (not expr) {
    this->abort(message1, message2, message3, message4);
  }
} // Log::check

// Assert expression is true, otherwise abort.
void Log::check(const bool expr, const std::string message1, const std::string message2, const std::string message3, const std::string message4, const std::string message5) {
  if (not expr) {
    this->abort(message1, message2, message3, message4, message5);
  }
} // Log::check

// Assert expression is true, otherwise abort.
void Log::check(const bool expr, const std::string message1, const std::string message2, const std::string message3, const std::string message4, const std::string message5, const std::string message6) {
  if (not expr) {
    this->abort(message1, message2, message3, message4, message5, message6);
  }
} // Log::check

// Barrier procs then stop
void Log::stop(void) {

  char mpiErrBuffer[MPI_MAX_ERROR_STRING];
  int mpiErrLen;
  int mpiErr;

  mpiErr = MPI_Barrier(partition->comm());
  MPI_Error_string(mpiErr, mpiErrBuffer, &mpiErrLen);
  check(mpiErr == MPI_SUCCESS, "log stop: MPI_Barrier failed: ", (std::string) mpiErrBuffer);

  throw std::runtime_error("Stop");

} // Log::stop

} // DEMSI namespace
