/*
  DEMSI logging class

 */

#ifndef DEMSI_LOGGING_H
#define DEMSI_LOGGING_H

#include "demsi_partition.h"

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <exception>
#include <stdexcept>

/*! \file   demsi_logging.h
    \brief  Header file for the DEMSI::Log class
*/

namespace DEMSI {

  /*! /enum Log message type */
  enum LOG {DEBUG, NORMAL, WARNING, ERROR, CRITICAL};

//------------------------------------------------------------------------------
/*! \class Log
    \brief Class allowing logging messages to be written to stdout and logging
    files.

    Output log messages can be of the following types:
    - DEBUG: Debug messages prefixed with "DEBUG: "
    - NORMAL: Normal messages.
    - WARNING: Warning messages prefixed with "Warning: "
    - ERROR: Error messages prefixed with "Error: "
    - CRITICAL: Fatal error messages prefixed with "Critical: ". These only
    occur when the model aborts.

    Two types of log files are generated: a normal one where all messages are
    written (called "log.XXXX.out", where XXXX is the processor id), and one
    where error and critical error messages are written (called "log.XXXX.err").
    The error log files are only opened if an error or critical message is
    generated on a processor. The normal log file is always generated for the
    master task, but can also be generated for all processes, through a config
    option ("writeLogAllProcs"). Debug messages are not written normally
    although can also be turned on with a config option ("writeLogDebug").

    Messages can be written to the log files with two commands:
    - @code log->write(message); @endcode
    - @code (*log) << message; @endcode

    The second method can be chained.

    The message type is altered with the parenthesis operator:
    - @code (*log)(DEMSI::LOG::ERROR); @endcode

    and can also be chained so that a full log message might be:
    - @code (*log)(DEMSI::LOG::ERROR) << "some" << "message" << std::endl; @endcode

    The default argument for the parenthesis operator is DEMSI::LOG::NORMAL.
    The message type stays fixed until the parenthesis operator is called again.

    The class also supplies two other methods, an abort method thats allows the
    model to be aborted:
    - @code log->abort("Abort message"); @endcode

    and a method that allows assertion of a logical expression, with the model
    aborting if the assertion returns false.
    - @code log->check(1 > 2, "Abort message"); // Abort @endcode

    These two methods can have up to six messages in their argument list but
    they currently must be of type std::string.
*/
//------------------------------------------------------------------------------
class Log {
public:

  /*! \brief Log class constructor.
      \param partitionIn Pointer to the partition object for this logging.
   */
  Log(DEMSI::Partition* partitionIn);

  /*! \brief Log class destructor. */
  ~Log();

  /*! \brief set up logging with configs
      \param writeLogAllProcs Flag whether to write a log on all procs.
      \param writeDebugIn Flag to write out debug messages
   */
  void setup(const bool writeLogAllProcs, const bool writeDebugIn);

  /*! \brief Write the message to the output log streams.
      \param message Message to write to the output streams.
   */
  template <typename T> void write(const T message) {
    if (writeToLog and (logType != DEMSI::LOG::DEBUG or writeDebug)) logFile << message << std::flush;
    if (partition->on_master_proc() and (logType != DEMSI::LOG::DEBUG or writeDebug)) std::cout << message << std::flush;
    if (logType == DEMSI::LOG::ERROR or logType == DEMSI::LOG::CRITICAL) errFile << message << std::flush;
  }

  /*! \brief Write the message to the output log streams (const version).
      \param message Message to write to the output streams.
   */
  template <typename T> void write(const T message) const {
    if (writeToLog and (logType != DEMSI::LOG::DEBUG or writeDebug)) logFile << message << std::flush;
    if (partition->on_master_proc() and (logType != DEMSI::LOG::DEBUG or writeDebug)) std::cout << message << std::flush;
    if (logType == DEMSI::LOG::ERROR or logType == DEMSI::LOG::CRITICAL) errFile << message << std::flush;
  }

  /*! \brief Write the ostream object to the output streams.
      \param pf ostream object to write (e.g. std::endl).
   */
  void write(std::ostream& (*pf) (std::ostream&));

  /*! \brief Overloaded << operator, write message to output streams.
      \param message Message to write to the output streams.
   */
  template <typename T>       Log &operator<<(const T &message)       {this->write(message); return *this;}

  /*! \brief Overloaded << operator, write message to output streams (const version).
      \param message Message to write to the output streams.
   */
  template <typename T> const Log &operator<<(const T &message) const {this->write(message); return *this;}

  /*! \brief Overloaded << operator, write ostream object to the output streams.
      \param pf ostream object to write (e.g. std::endl).
   */
  Log &operator<<(std::ostream& (*pf) (std::ostream&));

  /*! \brief Overloaded , operator, write message to output streams after space.
      \param message Message to write to the output streams.
   */
  template <typename T>       Log &operator,(const T &message)       {this->write(" "); this->write(message); return *this;}

  /*! \brief Overloaded , operator, write message to output streams after space (const version).
      \param message Message to write to the output streams.
   */
  template <typename T> const Log &operator,(const T &message) const {this->write(" "); this->write(message); return *this;}

  /*! \brief Overloaded , operator, write ostream object to the output streams.
      \param pf ostream object to write (e.g. std::endl).
   */
  Log &operator,(std::ostream& (*pf) (std::ostream&));

  /*! \brief Parenthesis operator, set the message log type.
      \param logTypeIn Message log type to set.
   */
  Log &operator()(const DEMSI::LOG logTypeIn = DEMSI::LOG::NORMAL);

  /*! \brief Abort the model.
      \param message Abort message to write to log streams before abort.
   */
  [[noreturn]] void abort(const std::string message="");

  /*! \brief Abort the model.
      \param message1 Abort message (1/2) to write to log streams before abort.
      \param message2 Abort message (2/2) to write to log streams before abort.
   */
  [[noreturn]] void abort(const std::string message1, const std::string message2);

  /*! \brief Abort the model.
      \param message1 Abort message (1/3) to write to log streams before abort.
      \param message2 Abort message (2/3) to write to log streams before abort.
      \param message3 Abort message (3/3) to write to log streams before abort.
   */
  [[noreturn]] void abort(const std::string message1, const std::string message2, const std::string message3);

  /*! \brief Abort the model.
      \param message1 Abort message (1/4) to write to log streams before abort.
      \param message2 Abort message (2/4) to write to log streams before abort.
      \param message3 Abort message (3/4) to write to log streams before abort.
      \param message4 Abort message (4/4) to write to log streams before abort.
   */
  [[noreturn]] void abort(const std::string message1, const std::string message2, const std::string message3, const std::string message4);

  /*! \brief Abort the model.
      \param message1 Abort message (1/5) to write to log streams before abort.
      \param message2 Abort message (2/5) to write to log streams before abort.
      \param message3 Abort message (3/5) to write to log streams before abort.
      \param message4 Abort message (4/5) to write to log streams before abort.
      \param message5 Abort message (5/5) to write to log streams before abort.
   */
  [[noreturn]] void abort(const std::string message1, const std::string message2, const std::string message3, const std::string message4, const std::string message5);

  /*! \brief Abort the model.
      \param message1 Abort message (1/6) to write to log streams before abort.
      \param message2 Abort message (2/6) to write to log streams before abort.
      \param message3 Abort message (3/6) to write to log streams before abort.
      \param message4 Abort message (4/6) to write to log streams before abort.
      \param message5 Abort message (5/6) to write to log streams before abort.
      \param message6 Abort message (6/6) to write to log streams before abort.
   */
  [[noreturn]] void abort(const std::string message1, const std::string message2, const std::string message3, const std::string message4, const std::string message5, const std::string message6);

  /*! \brief Assert expression is true, otherwise abort.
      \param expr Logical expression to check for truth.
      \param message Abort message to write to log streams before abort if assertion fails.
   */
  void check(const bool expr, const std::string message="");

  /*! \brief Assert expression is true, otherwise abort.
      \param expr Logical expression to check for truth.
      \param message1 Abort message (1/2) to write to log streams before abort if assertion fails.
      \param message2 Abort message (2/2) to write to log streams before abort if assertion fails.
   */
  void check(const bool expr, const std::string message1, const std::string message2);

  /*! \brief Assert expression is true, otherwise abort.
      \param expr Logical expression to check for truth.
      \param message1 Abort message (1/3) to write to log streams before abort if assertion fails.
      \param message2 Abort message (2/3) to write to log streams before abort if assertion fails.
      \param message3 Abort message (3/3) to write to log streams before abort if assertion fails.
   */
  void check(const bool expr, const std::string message1, const std::string message2, const std::string message3);

  /*! \brief Assert expression is true, otherwise abort.
      \param expr Logical expression to check for truth.
      \param message1 Abort message (1/4) to write to log streams before abort if assertion fails.
      \param message2 Abort message (2/4) to write to log streams before abort if assertion fails.
      \param message3 Abort message (3/4) to write to log streams before abort if assertion fails.
      \param message4 Abort message (4/4) to write to log streams before abort if assertion fails.
   */
  void check(const bool expr, const std::string message1, const std::string message2, const std::string message3, const std::string message4);

  /*! \brief Assert expression is true, otherwise abort.
      \param expr Logical expression to check for truth.
      \param message1 Abort message (1/5) to write to log streams before abort if assertion fails.
      \param message2 Abort message (2/5) to write to log streams before abort if assertion fails.
      \param message3 Abort message (3/5) to write to log streams before abort if assertion fails.
      \param message4 Abort message (4/5) to write to log streams before abort if assertion fails.
      \param message5 Abort message (5/5) to write to log streams before abort if assertion fails.
   */
  void check(const bool expr, const std::string message1, const std::string message2, const std::string message3, const std::string message4, const std::string message5);

  /*! \brief Assert expression is true, otherwise abort.
      \param expr Logical expression to check for truth.
      \param message1 Abort message (1/6) to write to log streams before abort if assertion fails.
      \param message2 Abort message (2/6) to write to log streams before abort if assertion fails.
      \param message3 Abort message (3/6) to write to log streams before abort if assertion fails.
      \param message4 Abort message (4/6) to write to log streams before abort if assertion fails.
      \param message5 Abort message (5/6) to write to log streams before abort if assertion fails.
      \param message6 Abort message (6/6) to write to log streams before abort if assertion fails.
   */
  void check(const bool expr, const std::string message1, const std::string message2, const std::string message3, const std::string message4, const std::string message5, const std::string message6);

  /*! \brief Barrier procs then stop
   */
  void stop(void);

private:

  /*! Pointer to the partition used for logging. */
  DEMSI::Partition* partition;

  /*! Output stream for all log messages */
  std::ofstream logFile;

  /*! Output stream for error and critical error messages. */
  std::ofstream errFile;

  /*! Flag whether to write to normal log stream on this processor. */
  bool writeToLog;

  /*! Flag whether to write debug messages. */
  bool writeDebug;

  /*! Current log type for messages. */
  DEMSI::LOG logType;

};


} // namespace DEMSI

#endif
