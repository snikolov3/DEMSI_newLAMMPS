#ifndef DEMSI_DIAGNOSTICS_H_
#define DEMSI_DIAGNOSTICS_H_

#include "demsi_logging.h"
#include "demsi_configs.h"
#include "demsi_lmp_instance.h"
#include "demsi_time.h"

/*! \file   demsi_diagnostics.h
    \brief  Header file for the DEMSI::Diagnostics class
*/

namespace DEMSI {

/*! \class Diagnostics
 \brief Class for diagnostics for debugging the code.
*/
class Diagnostics {
public:

  /*! \brief constructor for Diagnostics class. */
  Diagnostics(DEMSI::Log *logIn, DEMSI::Configs *configsIn, DEMSI::LammpsInstance* lammpsInstanceIn, DEMSI::Calendar* calendarIn, DEMSI::Clock* simulationClockIn);

  /*! \brief default destructor for Column class. */
  ~Diagnostics() = default;

  /*! \brief Start a lammps dump output. */
  void start_lammps_dump(void);

private:

  /*! Pointer to the log object */
  DEMSI::Log *log;

  /*! Pointer to the configs object */
  DEMSI::Configs *configs;

  /*! Pointer to lammps instance */
  DEMSI::LammpsInstance* lammpsInstance;

  /*! Pointer to the calendar clock */
  DEMSI::Calendar* calendar;

  /*! Pointer to the simulation clock */
  DEMSI::Clock* simulationClock;

  /*! Do we write a dump file */
  bool writeDump = false;

  /*! Time of start of dump */
  DEMSI::Time* dumpTime = NULL;

  /*! Whether dump has been set yet */
  bool dumpSet = false;

  /*! Frequency of lammps dump output */
  int dumpFrequency = 0;

}; // Diagnostics class

} // namespace DEMSI

#endif /* DIAGNOSTICS_H_ */
