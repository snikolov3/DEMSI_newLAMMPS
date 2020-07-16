#include "demsi_diagnostics.h"

namespace DEMSI {

// constructor for Diagnostics class.
Diagnostics::Diagnostics(DEMSI::Log *logIn, DEMSI::Configs *configsIn, DEMSI::LammpsInstance* lammpsInstanceIn, DEMSI::Calendar* calendarIn, DEMSI::Clock* simulationClockIn) {

  log = logIn;
  configs = configsIn;
  lammpsInstance = lammpsInstanceIn;
  calendar = calendarIn;
  simulationClock = simulationClockIn;

  // LAMMPS dump
  if (configs->is_true({"ConfigGroup:Diagnostics","Config:useDump"})) {
    writeDump = true;
    std::string dumpTimeStr;
    configs->get({"ConfigGroup:Diagnostics","Config:dumpStartTime"},dumpTimeStr);
    dumpTime = new DEMSI::Time(calendar, dumpTimeStr, log);
    configs->get({"ConfigGroup:Diagnostics","Config:dumpFrequency"},dumpFrequency);
  }

}

// Start a lammps dump output.
void Diagnostics::start_lammps_dump(void) {

  if (writeDump and not dumpSet and simulationClock->get_time() >= *dumpTime) {
    (*log)(DEMSI::LOG::DEBUG) << "...Set up dump..." << std::endl;
    std::stringstream lammpsCommand;
    lammpsCommand << "dump DEMSIDump all custom " << dumpFrequency << " demsi.dump id type x y z radius mass fx fy fz" << std::endl;
    lammpsInstance->one(lammpsCommand.str());
    lammpsInstance->one("run 0 pre yes post no");
    dumpSet = true;
  }

}

} // namespace DEMSI
