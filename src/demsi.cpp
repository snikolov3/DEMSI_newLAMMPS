/* ----------------------------------------------------------------------
  Simple driver of lammps with an external force field specified
  on a grid.

 */

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <mpi.h>
#include "demsi_file_utils.h"
#include "demsi_partition.h"
#include "demsi_logging.h"
#include "demsi_kokkosparser.h"
#include "demsi_configs.h"
#include "demsi_lmp_instance.h"
#include "demsi_grid.h"
#include "demsi_grid_io.h"
#include "demsi_external_force.h"
#include "demsi_forcing.h"
#include "demsi_time.h"
#include "demsi_particles.h"
#include "demsi_initial_tessellation.h"
#include "demsi_initial_tessellation_io.h"
#include "demsi_particles_io.h"
#include "demsi_interpolation.h"
#include "demsi_contacts.h"
#include "demsi_column.h"
#include "demsi_remapping.h"
#include "demsi_ocean.h"
#include "demsi_diagnostics.h"

#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>

#include <lammps/lammps.h>
#include <lammps/input.h>
#include <lammps/atom.h>
#include <lammps/library.h>
#include <lammps/comm.h>

#include <Kokkos_Core.hpp>

#include <permonqps.h>

int main(int narg, char **arg)
{
  // Init mpi
  DEMSI::Partition *partition = new DEMSI::Partition(&narg, &arg);

  static char help[] = "";
  int ierrPermon = PermonInitialize(&narg, &arg, (char *)0, help);
  PetscPopSignalHandler();

  // Init logging
  DEMSI::Log *log = new DEMSI::Log(partition);

  // Pass arguments to kokkos
  DEMSI::KokkosParser kokkos_parser(narg, arg, partition->comm(), log);

  // Init kokkos
  auto kokkos_args = kokkos_parser.getKokkosInitArguments();
  Kokkos::initialize(kokkos_args);
  kokkos_parser.logStatus();

  Kokkos::Profiling::pushRegion("Initialize");

  // load run configuration
  log->check(narg >= 2, "Must specify config file as command line argument.");
  (*log)() << "Initialize configs..." << std::endl;
  DEMSI::Configs *configs = new DEMSI::Configs(arg[1], log);

  // set up log from configs
  bool writeLogAllProcs;
  configs->get({"ConfigGroup:logging","Config:writeLogAllProcs"}, writeLogAllProcs);
  bool writeDebug;
  configs->get({"ConfigGroup:logging","Config:writeLogDebug"}, writeDebug);
  log->setup(writeLogAllProcs, writeDebug);

  // init the time manager
  (*log)() << "Init simulation clock..." << std::endl;

  // init simulation calendar
  std::string calendarType;
  configs->get({"ConfigGroup:simulationTiming","Config:calendarType"}, calendarType);
  DEMSI::Calendar calendar(calendarType, log);

  // init simulation start time
  std::string startTimeStr;
  configs->get({"ConfigGroup:simulationTiming","Config:startTime"}, startTimeStr);
  DEMSI::Time startTime(&calendar, startTimeStr, log);

  // init simulation main time step interval
  int timeStep;
  configs->get({"ConfigGroup:simulationTiming","Config:timeStep"}, timeStep);
  std::string timeStepStr = "SECONDS:" + std::to_string(timeStep);
  DEMSI::TimeInterval timeStepInterval(timeStepStr, log);

  // init simulation clock
  DEMSI::Clock clock(startTime, timeStepInterval);

  // init simulation stop time
  std::string simulationStopType;
  configs->get({"ConfigGroup:simulationTiming","Config:simulationStopType"}, simulationStopType);
  DEMSI::Time stopTime(&calendar, log);

  if (simulationStopType == "stopTime") {

    // set stop time explicitly
    std::string stopTimeStr;
    configs->get({"ConfigGroup:simulationTiming","Config:stopTime"}, stopTimeStr);
    stopTime.set(stopTimeStr);

  } else if (simulationStopType == "nSteps") {

    // set stop time from nSteps
    int nSteps;
    configs->get({"ConfigGroup:simulationTiming","Config:nSteps"}, nSteps);
    DEMSI::TimeInterval simulationDuration = timeStepInterval * nSteps;
    stopTime = startTime + simulationDuration;

  } else if (simulationStopType == "simulationDuration") {

    // set stop time from duration
    std::string duration;
    configs->get({"ConfigGroup:simulationTiming","Config:simulationDuration"}, duration);
    DEMSI::TimeInterval simulationDuration(duration, log);
    stopTime = startTime + simulationDuration;

  } else {
    (*log)(DEMSI::LOG::CRITICAL) << "Unknown simulation stop type: " << simulationStopType << std::endl;
    log->abort("Unknown simulation stop type");
  }

  // Initialize LAMMPS
  (*log)() << "Get LAMMPS instances..." << std::endl;
  DEMSI::LammpsInstance *lammpsInstance = new DEMSI::LammpsInstance(partition, 
          kokkos_parser.getNumberOfThreads(), kokkos_parser.getNuma(), kokkos_parser.getNumberOfGPUs());

  // initialize the grid object
  (*log)() << "Create grids objects..." << std::endl;
  std::string gridFileName;
  configs->get({"ConfigGroup:grid","Config:gridFile"}, gridFileName);
  DEMSI::Grid *grid = new DEMSI::Grid(gridFileName, partition, log, lammpsInstance, &clock);

  // create contacts object
  (*log)() << "Create contacts object..." << std::endl;
  DEMSI::Contacts *contacts = new DEMSI::Contacts(log, configs, lammpsInstance);

  // create particle list objects
  (*log)() << "Create demsi particle lists objects..." << std::endl;
  DEMSI::Particles *particles = new DEMSI::Particles(partition, log, configs, lammpsInstance, grid, contacts);

  // create lammps particle files
  (*log)() << "Create lammps particles from input file..." << std::endl;
  particles->create_lammps_particles();

  // setup lammps
  (*log)() << "Setup lammps..." << std::endl;
  particles->setup_lammps();

  // set demsi particles pointers
  (*log)() << "Set DEMSI particle pointers..." << std::endl;
  particles->set_demsi_pointers();

  // set initial particle velocity
  (*log)() << "Set DEMSI particle velocity..." << std::endl;
  particles->initial_velocity();

  // Set the grid partition and read in grid fields
  (*log)() << "Set the grid partition and read in grid fields..." << std::endl;
  grid->read_grid();

  // read in an initial element distribution
  (*log)() << "Read in initial element distribution..." << std::endl;
  DEMSI::Tessellation* tessellation = new DEMSI::Tessellation(partition, log, configs, lammpsInstance, grid, particles);
  // set up tessellation output
  (*log)() << "Setup tessellation write..." << std::endl;
  DEMSI::TessellationOutputStreams* tessellationOutputStreams = new DEMSI::TessellationOutputStreams(partition, log, configs, tessellation, &clock, startTime);

  // initialize forcing object
  (*log)() << "Create forcing object..." << std::endl;
  DEMSI::Forcing *forcing = new DEMSI::Forcing(log, configs, grid, &calendar, &clock);
  forcing->update();

  // column setup
  (*log)() << "Set up column physics..." << std::endl;
  DEMSI::Column *column = new DEMSI::Column(partition, log, configs, lammpsInstance, grid, particles, tessellation, forcing, &clock, &timeStepInterval);

  // set up ocean object
  (*log)() << "Set up ocean..." << std::endl;
  DEMSI::Ocean* ocean = new DEMSI::Ocean(log, configs, grid, tessellation, forcing, column, &clock, &timeStepInterval);

  // initialize column object
  (*log)() << "Initialize column object..." << std::endl;
  column->init();
  column->update_mass();

  // initialize external force object
  (*log)() << "Create external force object..." << std::endl;
  DEMSI::ExternalForce *externalForce = new DEMSI::ExternalForce(log, configs, grid, particles, forcing);

  // set up remapping
  (*log)() << "Set up remapping..." << std::endl;
  DEMSI::Remapping* remapping = new DEMSI::Remapping(partition, log, configs, lammpsInstance, grid, contacts, particles, tessellation, tessellationOutputStreams, column, ocean, &clock);

  //Lazy klooge here - run 0 steps to register the grid forces
  (*log)() << "Run model for zero steps..." << std::endl;
  Kokkos::Profiling::pushRegion("LAMMPS");
  lammpsInstance->run(0);
  Kokkos::Profiling::popRegion();
  particles->set_demsi_pointers();
  column->processor_transfer();

  // Initialize particle contacts
  (*log)() << "Create demsi contacts object..." << std::endl;
  if (contacts->get_contact_type() == DEMSI::CONTACTTYPE::HOPKINS) {
    contacts->initialize_bonded_contacts_hopkins();
    column->compute_mean_min_thickness();
    lammpsInstance->comm_ghosts();
  }
  if (contacts->get_contact_type() == DEMSI::CONTACTTYPE::HOOKE_THICKNESS){
    column->compute_mean_min_thickness();
    lammpsInstance->comm_ghosts();
  }

  // get number of dynamics subcycles
  int nDynamicsSubcycles;
  configs->get({"ConfigGroup:simulationTiming","Config:nDynamicsSubcycles"}, nDynamicsSubcycles);

  // get the dynamics timestep
  double dynamicsTimeStep = (double) timeStep / (double) nDynamicsSubcycles;

  // set lammps timestep
  (*log)() << "Set LAMMPS timestep..." << std::endl;
  lammpsInstance->set_timestep(dynamicsTimeStep);

  // set up particle to grid interpolation
  (*log)() << "Setup interpolation..." << std::endl;
  DEMSI::Interpolation *interp = new DEMSI::Interpolation(partition, log, configs, grid, particles, column);

  // set up particle write
  (*log)() << "Setup particle write..." << std::endl;
  DEMSI::ParticlesOutputStreams* particlesOutputStreams = new DEMSI::ParticlesOutputStreams(partition, log, configs, contacts, particles, column, &clock, startTime);
  (*log)() << "Particle write..." << std::endl;
  particlesOutputStreams->write();

  // set up gridded write
  (*log)() << "Setup grid write..." << std::endl;
  DEMSI::GridOutputStreams* gridOutputStreams = new DEMSI::GridOutputStreams(partition, log, configs, grid, &clock, startTime);

  //interpolate particle values to grid 
  if (interp->interpToGrid)
    interp->particle_fields_to_grid();

  // grid output write
  (*log)() << "Grid write..." << std::endl;
  gridOutputStreams->write();

  // tessellation output write
  (*log)() << "Tessellation write..." << std::endl;
  tessellationOutputStreams->write();

  // use options
  bool useDynamics;
  configs->get({"ConfigGroup:useSections","Config:useDynamics"}, useDynamics);
  bool useLammpsControl = false;
  std::string lammpsControlScript;
  if (configs->exists({"ConfigGroup:lammps", "Config:useLammpsControl"})){
    configs->get({"ConfigGroup:lammps","Config:useLammpsControl"}, useLammpsControl);
    if (useLammpsControl){
      configs->get({"ConfigGroup:lammps","Config:lammpsControlScript"}, lammpsControlScript);
    }
  }

  // diagnostics object
  (*log)() << "Set up diagnostics object..." << std::endl;
  DEMSI::Diagnostics* diagnostics = new DEMSI::Diagnostics(log, configs, lammpsInstance, &calendar, &clock);

  // end initialization
  Kokkos::Profiling::popRegion();

  // advance the simulation clock
  (*log)() << "Advance clock..." << std::endl;
  clock.advance();

  // forcing advance
  (*log)() << "Forcing read..." << std::endl;
  forcing->update();

  // ocean forcing interpolate
  (*log)() << "Ocean forcing interpolate..." << std::endl;
  ocean->interpolate_forcing();

  // initialize ocean sst
  (*log)() << "Initialize ocean SST..." << std::endl;
  ocean->init_sst();

  // starting remap
  (*log)() << "Initial remap..." << std::endl;
  remapping->remap(true, true);

  // column forcing interpolate
  (*log)() << "Column forcing interpolate..." << std::endl;
  column->interpolate_forcing();

  // initialize column sst
  (*log)() << "Initialize column SST..." << std::endl;
  column->init_sst();

  // column shortwave init
  (*log)() << "Init column shortwave..." << std::endl;
  column->init_shortwave();

  // Main loop
  (*log)() << "Main loop..." << std::endl;
  Kokkos::Profiling::pushRegion("Time integration");
  while (clock.get_time() <= stopTime) {

    (*log)() << "Doing timestep " << clock << std::endl << std::flush;

    // pre-dynamics column calculation
    (*log)(DEMSI::LOG::DEBUG) << "...column pre-dynamics..." << std::endl;
    Kokkos::Profiling::pushRegion("Column pre-dynamics");
    column->pre_dynamics();
    Kokkos::Profiling::popRegion();

    // ocean
    (*log)(DEMSI::LOG::DEBUG) << "...ocean update..." << std::endl;
    Kokkos::Profiling::pushRegion("Ocean update");
    ocean->update();
    Kokkos::Profiling::popRegion();

    if (useDynamics) {

      // update particle masses
      (*log)(DEMSI::LOG::DEBUG) << "...column update masses..." << std::endl;
      Kokkos::Profiling::pushRegion("Update mass");
      column->update_mass();
      Kokkos::Profiling::popRegion();

      if (useLammpsControl) {
        lammpsInstance->run_script(lammpsControlScript);
      } else {

	// start lammps dump file
	(*log)(DEMSI::LOG::DEBUG) << "...Start LAMMPS dump file..." << std::endl;
	diagnostics->start_lammps_dump();

        // iterate over dynamics subcycles
        for (int iDynamicsSubcycle=0 ; iDynamicsSubcycle < nDynamicsSubcycles ; iDynamicsSubcycle++) {

          (*log)(DEMSI::LOG::DEBUG) << "...Dynamics subcyle: " << iDynamicsSubcycle << " of " << nDynamicsSubcycles << std::endl;

	  // update particle thicknesses
	  (*log)(DEMSI::LOG::DEBUG) << "...column update thicknesses..." << std::endl;
	  if (contacts->get_contact_type() == DEMSI::CONTACTTYPE::HOPKINS or
	      contacts->get_contact_type() == DEMSI::CONTACTTYPE::HOOKE_THICKNESS) {
	    Kokkos::Profiling::pushRegion("Mean min thickness");
	    column->compute_mean_min_thickness();
	    Kokkos::Profiling::popRegion();
	  }

	  // set ridging variables
	  (*log)(DEMSI::LOG::DEBUG) << "...column ridging variables..." << std::endl;
	  Kokkos::Profiling::pushRegion("Ridging variables");
	  column->set_ridging_variables();
	  Kokkos::Profiling::popRegion();

          //Currently hard-wired for two values on each grid (wind velocity x and y components).
          (*log)(DEMSI::LOG::DEBUG) << "...Add external forces..." << std::endl;
	  Kokkos::Profiling::pushRegion("External forces");
          externalForce->add_forces();
	  Kokkos::Profiling::popRegion();

	  // delete particles
	  (*log)(DEMSI::LOG::DEBUG) << "...Delete particles..." << std::endl;
	  Kokkos::Profiling::pushRegion("Delete particles");
	  column->delete_particles();
	  Kokkos::Profiling::popRegion();

          //Run dynamics with fixed forcing
          (*log)(DEMSI::LOG::DEBUG) << "...LAMMPS run subcycle..." << std::endl;
	  Kokkos::Profiling::pushRegion("Sync data to LAMMPS");
          particles->sync_data_to_lammps();
	  Kokkos::Profiling::popRegion();

          Kokkos::Profiling::pushRegion("LAMMPS");
          lammpsInstance->run_subcycle();
	  Kokkos::Profiling::popRegion();

	  Kokkos::Profiling::pushRegion("LAMMPS reneighbour");
	  lammpsInstance->reneighbor();
          Kokkos::Profiling::popRegion();

          particles->set_demsi_pointers();
	  //particles->check_in_domain();

	  // communicate demsi column particle data
	  (*log)(DEMSI::LOG::DEBUG) << "...column communication across processors..." << std::endl;
	  Kokkos::Profiling::pushRegion("Processor transfer");
	  column->processor_transfer();
	  Kokkos::Profiling::popRegion();

	  // ridging
	  //(*log)(DEMSI::LOG::DEBUG) << "...ridging..." << std::endl;
	  //column->ridge();
	  (*log)(DEMSI::LOG::DEBUG) << "...column dynamics..." << std::endl;
	  Kokkos::Profiling::pushRegion("Column dynamics");
	  column->dynamics();
	  Kokkos::Profiling::popRegion();

        } // dynamics subcycle
      }

    } // use dynamics

    // particle output
    (*log)(DEMSI::LOG::DEBUG) << "...particle write..." << std::endl;
    Kokkos::Profiling::pushRegion("Particle write");
    particlesOutputStreams->write();
    Kokkos::Profiling::popRegion();

    // tessellation output
    (*log)(DEMSI::LOG::DEBUG) << "...tessellation write..." << std::endl;
    Kokkos::Profiling::pushRegion("Tessellation write");
    tessellationOutputStreams->write();
    Kokkos::Profiling::popRegion();

    // remap the particle distribution
    (*log)(DEMSI::LOG::DEBUG) << "...remap particles..." << std::endl;
    Kokkos::Profiling::pushRegion("Remapping");
    remapping->remap();
    Kokkos::Profiling::popRegion();

    // column forcing interpolate - redo because of remapping
    (*log)(DEMSI::LOG::DEBUG) << "Column forcing interpolate..." << std::endl;
    Kokkos::Profiling::pushRegion("Interpolate forcing");
    column->interpolate_forcing();
    Kokkos::Profiling::popRegion();

    // post-dynamics column calculation
    (*log)(DEMSI::LOG::DEBUG) << "...column post-dynamics..." << std::endl;
    Kokkos::Profiling::pushRegion("Column post-dynamics");
    column->post_dynamics();
    Kokkos::Profiling::popRegion();

    // interpolate selected particle fields to grid
    if (interp->interpToGrid) {
       (*log)(DEMSI::LOG::DEBUG) << "...particle-to-grid interpolation..." << std::endl;
       Kokkos::Profiling::pushRegion("Particle-to-grid interpolation");
       interp->particle_fields_to_grid();
       Kokkos::Profiling::popRegion();
    }

    // gridded output
    (*log)(DEMSI::LOG::DEBUG) << "...gridded write..." << std::endl;
    Kokkos::Profiling::pushRegion("Write gridded");
    gridOutputStreams->write();
    Kokkos::Profiling::popRegion();

    // finalize column
    (*log)(DEMSI::LOG::DEBUG) << "...column finalize..." << std::endl;
    Kokkos::Profiling::pushRegion("Column finalize timestep");
    column->finalize_timestep();
    Kokkos::Profiling::popRegion();

    // advance the simulation clock
    (*log)(DEMSI::LOG::DEBUG) << "...advance clock..." << std::endl;
    clock.advance();

    // forcing advance
    (*log)(DEMSI::LOG::DEBUG) << "...forcing read..." << std::endl;
    Kokkos::Profiling::pushRegion("Forcing read");
    forcing->update();
    Kokkos::Profiling::popRegion();

    // forcing intrepolate
    (*log)(DEMSI::LOG::DEBUG) << "...forcing interpolate..." << std::endl;
    Kokkos::Profiling::pushRegion("Interpolate forcing");
    column->interpolate_forcing();
    ocean->interpolate_forcing();
    Kokkos::Profiling::popRegion();

  } // main loop
  Kokkos::Profiling::popRegion();

  // Finalize
  (*log)() << "Finalizing..." << std::endl;
  Kokkos::Profiling::pushRegion("Finalize");

  delete grid;
  delete particles;
  delete particlesOutputStreams;
  delete lammpsInstance;
  delete log;
  delete partition;
  delete forcing;
  delete externalForce;

  ierrPermon = PermonFinalize();

  Kokkos::Profiling::popRegion();

  Kokkos::finalize();

  MPI_Finalize();

}
