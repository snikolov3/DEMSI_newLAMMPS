#include <iostream>
#include <Kokkos_Core.hpp>

#include "../../src/demsi_partition.h"
#include "../../src/demsi_logging.h"
#include "../../src/demsi_kokkosparser.h"
#include "../../src/demsi_configs.h"
#include "../../src/demsi_lmp_instance.h"
#include "../../src/demsi_grid.h"
#include "../../src/demsi_time.h"
#include "../../src/demsi_particles.h"
#include "../../src/demsi_initial_tessellation.h"
#include "../../src/demsi_particles_io.h"
#include "../../src/demsi_column.h"
#include "../../src/demsi_forcing.h"
#include "../../src/demsi_interpolation.h"

int main( int argc, char *argv[] ) {

  DEMSI::Partition *partition = new DEMSI::Partition(&argc, &argv);

  DEMSI::Log *log = new DEMSI::Log(partition);

  // Pass arguments to kokkos
  DEMSI::KokkosParser kokkos_parser(argc, argv, partition->comm(), log);
  auto kokkos_args = kokkos_parser.getKokkosInitArguments();
  Kokkos::initialize(kokkos_args);
  kokkos_parser.logStatus();
  {

  DEMSI::Configs *configs = new DEMSI::Configs("config.xml", log);

  std::string calendarType;
  configs->get({"ConfigGroup:simulationTiming","Config:calendarType"}, calendarType);
  DEMSI::Calendar calendar(calendarType, log);

  std::string startTimeStr;
  configs->get({"ConfigGroup:simulationTiming","Config:startTime"}, startTimeStr);
  DEMSI::Time startTime(&calendar, startTimeStr, log);

  int timeStep;
  configs->get({"ConfigGroup:simulationTiming","Config:timeStep"}, timeStep);
  std::string timeStepStr = "SECONDS:" + std::to_string(timeStep);
  DEMSI::TimeInterval timeStepInterval(timeStepStr, log);

  DEMSI::Clock clock(startTime, timeStepInterval);

  DEMSI::LammpsInstance *lammpsInstance = new DEMSI::LammpsInstance(partition, 
          kokkos_parser.getNumberOfThreads(), kokkos_parser.getNuma(), kokkos_parser.getNumberOfGPUs());

  std::string gridFileName;
  configs->get({"ConfigGroup:grid","Config:gridFile"}, gridFileName);
  DEMSI::Grid *grid = new DEMSI::Grid(gridFileName, partition, log, lammpsInstance, &clock);

  DEMSI::Contacts *contacts = new DEMSI::Contacts(log, configs, lammpsInstance);

  DEMSI::Particles *particles = new DEMSI::Particles(partition, log, configs, lammpsInstance, grid, contacts);

  particles->create_lammps_particles();
  particles->setup_lammps();

  particles->set_demsi_pointers();

   DEMSI::Tessellation* tessellation = new DEMSI::Tessellation(partition, log, configs, lammpsInstance, grid, particles);

  grid->read_grid();

  // initialize forcing object
  DEMSI::Forcing *forcing = new DEMSI::Forcing(log, configs, grid, &calendar, &clock);

  // column setup
  DEMSI::Column *column = new DEMSI::Column(partition, log, configs, lammpsInstance, grid, particles, tessellation, forcing, &clock, &timeStepInterval);
  column->init();
  column->update_mass();

//  DEMSI::Interpolation *interp = new DEMSI::Interpolation(particles, grid, log, partition);
 // set up particle to grid interpolation
  DEMSI::Interpolation *interp = new DEMSI::Interpolation(partition, log, configs, grid, particles, column);

  // Create particle and grid fields for interpolation
  Kokkos::View<double **> gridField = Kokkos::View<double **> ("gridField",grid->nx(),grid->ny());
  Kokkos::View<double *> partField = Kokkos::View<double *> ("partField",*(particles->nParticles));

  // Initialize particle values
  auto particles_x = particles->x;
  Kokkos::parallel_for("ExternalForce::add_forces", 
                        Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, *(particles->nParticles)), 
                        KOKKOS_LAMBDA (const int iParticle) {

    //Get particle position
    double x = particles_x(iParticle,0);
    double y = particles_x(iParticle,1);

    // linear function of position
    partField(iParticle) = x + y;

  });
   
  std::vector<DEMSI::GridPartition*> gridPartitions = grid->getGridPartitions();

  // Interpolate from particles to grid
  interp->particles_to_grid(partField,gridField);

  Kokkos::View<double **>::HostMirror gridField_mirror = Kokkos::create_mirror(gridField);
  Kokkos::deep_copy(gridField_mirror, gridField);
  for (int ix = 0 ; ix < gridPartitions[partition->this_proc()]->size_owned(0) ; ix++) {
    for (int iy = 0 ; iy < gridPartitions[partition->this_proc()]->size_owned(1) ; iy++) {

       int iGlobal = gridPartitions[partition->this_proc()]->iglobal_owned(ix);
       int jGlobal = gridPartitions[partition->this_proc()]->jglobal_owned(iy);

       int iLocalTot = gridPartitions[partition->this_proc()]->ilocal_total(iGlobal);
       int jLocalTot = gridPartitions[partition->this_proc()]->jlocal_total(jGlobal);

       double exactVal = iGlobal*grid->dx() + jGlobal*grid->dy();
       double approxVal = gridField_mirror(iLocalTot,jLocalTot);
       double errorVal = std::abs(exactVal-approxVal);

       double tol = 1.0e-9;

       if (errorVal > tol) {

       //  (*log)() << "Error in function interpolated to grid: \n";
         (*log)() << "  iproc = " << partition->this_proc();
         (*log)() << "  iowned = " << ix;
         (*log)() << "  jowned = " << iy;
         (*log)() << "  iglobal = " << iGlobal;
         (*log)() << "  jglobal = " << jGlobal;
         (*log)() << "  error = " << errorVal;
         (*log)() << "  grid val = " << gridField_mirror(iLocalTot,jLocalTot);
         (*log)() << "  exact val = " << iGlobal*grid->dx() + jGlobal*grid->dy() << "\n";
      }
    }
  }


  delete interp;
  delete grid;
  delete particles;
  delete lammpsInstance;
  delete log;
  delete partition;

  }
  Kokkos::finalize();

  MPI_Finalize();

}
