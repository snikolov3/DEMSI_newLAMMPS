#include "demsi_particles.h"

#include "demsi_communication.h"

#include <netcdf.h>
#include <lammps/atom_kokkos.h>
#include <lammps/atom_vec_demsi_kokkos.h>
#include <lammps/atom_masks.h>
#include <lammps/neighbor.h>
#include <iomanip>
#include <iostream>
#include <fstream>

namespace DEMSI {

Particles::Particles(DEMSI::Partition* partitionIn, DEMSI::Log* logIn, DEMSI::Configs* configsIn, DEMSI::LammpsInstance* lammpsInstanceIn, DEMSI::Grid* gridIn, DEMSI::Contacts* contactsIn) {

  partition = partitionIn;
  log = logIn;
  configs = configsIn;
  lammpsInstance = lammpsInstanceIn;
  grid = gridIn;
  contacts = contactsIn;

  iceDensity = 900;

}

void Particles::sync_data_to_lammps(void) {

  lammpsInstance->lmp->atomKK->k_x.sync<host_execution_space>();
  lammpsInstance->lmp->atomKK->k_v.sync<host_execution_space>();
  lammpsInstance->lmp->atomKK->k_radius.sync<host_execution_space>();
  lammpsInstance->lmp->atomKK->k_rmass.sync<host_execution_space>();
  lammpsInstance->lmp->atomKK->k_forcing.sync<host_execution_space>();
  lammpsInstance->lmp->atomKK->k_min_thickness.sync<host_execution_space>();
  lammpsInstance->lmp->atomKK->k_ridgingIceThickness.sync<host_execution_space>();
  lammpsInstance->lmp->atomKK->k_ridgingIceThicknessWeight.sync<host_execution_space>();
  lammpsInstance->lmp->atomKK->k_netToGrossClosingRatio.sync<host_execution_space>();
  lammpsInstance->lmp->atomKK->k_changeEffectiveElementArea.sync<host_execution_space>();
  lammpsInstance->lmp->atomKK->k_mean_thickness.sync<host_execution_space>();
  lammpsInstance->lmp->atomKK->k_ice_area.sync<host_execution_space>();
  lammpsInstance->lmp->atomKK->k_coriolis.sync<host_execution_space>();
  lammpsInstance->lmp->atomKK->k_ocean_vel.sync<host_execution_space>();
  lammpsInstance->lmp->atomKK->k_bvector.sync<host_execution_space>();

}

void Particles::set_demsi_pointers(void) {

  Kokkos::Profiling::pushRegion("DEMSI pointers");

  nParticles = &(lammpsInstance->lmp->atomKK->nlocal);

  if (*nParticles > 0){

    // gets the host/device view of each variable from lammps
    globalID      = lammpsInstance->lmp->atomKK->k_tag.h_view;
    type          = lammpsInstance->lmp->atomKK->k_type.h_view;
    x             = lammpsInstance->lmp->atomKK->k_x.d_view;
    v             = lammpsInstance->lmp->atomKK->k_v.d_view;
    radius        = lammpsInstance->lmp->atomKK->k_radius.d_view;
    mass          = lammpsInstance->lmp->atomKK->k_rmass.d_view;
    forcing       = lammpsInstance->lmp->atomKK->k_forcing.d_view;
    minThickness  = lammpsInstance->lmp->atomKK->k_min_thickness.d_view;
    ridgingIceThickness        = lammpsInstance->lmp->atomKK->k_ridgingIceThickness.d_view;
    ridgingIceThicknessWeight  = lammpsInstance->lmp->atomKK->k_ridgingIceThicknessWeight.d_view;
    netToGrossClosingRatio     = lammpsInstance->lmp->atomKK->k_netToGrossClosingRatio.d_view;
    changeEffectiveElementArea = lammpsInstance->lmp->atomKK->k_changeEffectiveElementArea.d_view;
    meanThickness = lammpsInstance->lmp->atomKK->k_mean_thickness.d_view;
    ice_area      = lammpsInstance->lmp->atomKK->k_ice_area.d_view;
    coriolis      = lammpsInstance->lmp->atomKK->k_coriolis.d_view;
    ocean_vel     = lammpsInstance->lmp->atomKK->k_ocean_vel.d_view;
    bvector       = lammpsInstance->lmp->atomKK->k_bvector.d_view;

    // ensure that views are up-to-date
    lammpsInstance->lmp->atomKK->k_tag.sync<host_execution_space>();
    lammpsInstance->lmp->atomKK->k_type.sync<host_execution_space>();
    lammpsInstance->lmp->atomKK->k_x.sync<device_execution_space>();
    lammpsInstance->lmp->atomKK->k_v.sync<device_execution_space>();
    lammpsInstance->lmp->atomKK->k_radius.sync<device_execution_space>();
    lammpsInstance->lmp->atomKK->k_rmass.sync<device_execution_space>();
    lammpsInstance->lmp->atomKK->k_forcing.sync<device_execution_space>();
    lammpsInstance->lmp->atomKK->k_min_thickness.sync<device_execution_space>();
    lammpsInstance->lmp->atomKK->k_ridgingIceThickness.sync<device_execution_space>();
    lammpsInstance->lmp->atomKK->k_ridgingIceThicknessWeight.sync<device_execution_space>();
    lammpsInstance->lmp->atomKK->k_netToGrossClosingRatio.sync<device_execution_space>();
    lammpsInstance->lmp->atomKK->k_changeEffectiveElementArea.sync<device_execution_space>();
    lammpsInstance->lmp->atomKK->k_mean_thickness.sync<device_execution_space>();
    lammpsInstance->lmp->atomKK->k_ice_area.sync<device_execution_space>();
    lammpsInstance->lmp->atomKK->k_coriolis.sync<device_execution_space>();
    lammpsInstance->lmp->atomKK->k_ocean_vel.sync<device_execution_space>();
    lammpsInstance->lmp->atomKK->k_bvector.sync<device_execution_space>();

    // indicate to lammps that we are modifying data on the device
    lammpsInstance->lmp->atomKK->k_x.modify<device_execution_space>();
    lammpsInstance->lmp->atomKK->k_v.modify<device_execution_space>();
    lammpsInstance->lmp->atomKK->k_radius.modify<device_execution_space>();
    lammpsInstance->lmp->atomKK->k_rmass.modify<device_execution_space>();
    lammpsInstance->lmp->atomKK->k_forcing.modify<device_execution_space>();
    lammpsInstance->lmp->atomKK->k_min_thickness.modify<device_execution_space>();
    lammpsInstance->lmp->atomKK->k_ridgingIceThickness.modify<device_execution_space>();
    lammpsInstance->lmp->atomKK->k_ridgingIceThicknessWeight.modify<device_execution_space>();
    lammpsInstance->lmp->atomKK->k_netToGrossClosingRatio.modify<device_execution_space>();
    lammpsInstance->lmp->atomKK->k_changeEffectiveElementArea.modify<device_execution_space>();
    lammpsInstance->lmp->atomKK->k_mean_thickness.modify<device_execution_space>();
    lammpsInstance->lmp->atomKK->k_ice_area.modify<device_execution_space>();
    lammpsInstance->lmp->atomKK->k_coriolis.modify<device_execution_space>();
    lammpsInstance->lmp->atomKK->k_ocean_vel.modify<device_execution_space>();
    lammpsInstance->lmp->atomKK->k_bvector.modify<device_execution_space>();

  }

  Kokkos::Profiling::popRegion();

}

// Create a LAMMPS input file from a DEMSI particles netcdf file.
void Particles::create_lammps_particles(void) {

  // input configs
  std::string particleInputFilename;
  configs->get({"ConfigGroup:particleInput","Config:particleInputFile"}, particleInputFilename);

  // check filename isn't none
  std::string filenameUpper = particleInputFilename;
  std::transform(filenameUpper.begin(), filenameUpper.end(), filenameUpper.begin(), toupper);
  if (filenameUpper != "NONE") {

    // master task reads in particles
    if (partition->on_master_proc()) {

      // open file
      int err, ncID;
      err = nc_open(particleInputFilename.c_str(), 0, &ncID);
      log->check(err == NC_NOERR, "Problem opening file: ", nc_strerror(err), ", for file: ", particleInputFilename);

      // get particle number
      int dimidParticles;
      err = nc_inq_dimid(ncID, "nParticles", &dimidParticles);
      log->check(err == NC_NOERR, "Problem getting dimid for nParticles: ", nc_strerror(err), ", for file: ", particleInputFilename);

      size_t nParticlesTmp;
      err = nc_inq_dimlen(ncID, dimidParticles, &nParticlesTmp);
      log->check(err == NC_NOERR, "Problem getting nParticles dimension: ", nc_strerror(err), ", for file: ", particleInputFilename);
      int nParticles = (int) nParticlesTmp;

      // get bond number
      int dimidBonds;
      int nBonds;
      err = nc_inq_dimid(ncID, "nBonds", &dimidBonds);
      if (err == NC_EBADDIM) {
        nBonds = 0;
      } else {
        log->check(err == NC_NOERR, "Problem getting dimid for nBonds: ", nc_strerror(err), ", for file: ", particleInputFilename);

        size_t nBondsTmp;
        err = nc_inq_dimlen(ncID, dimidBonds, &nBondsTmp);
        log->check(err == NC_NOERR, "Problem getting nBonds dimension: ", nc_strerror(err), ", for file: ", particleInputFilename);
        nBonds= (int) nBondsTmp;
      }

      // get global attributes
      double maxRadius;
      err = nc_get_att_double(ncID, NC_GLOBAL, "maxRadius", &maxRadius);
      log->check(err == NC_NOERR, "Problem getting maxRadius attribute: ", nc_strerror(err), ", for file: ", particleInputFilename);

      int nTypes;
      err = nc_get_att_int(ncID, NC_GLOBAL, "nTypes", &nTypes);
      log->check(err == NC_NOERR, "Problem getting nTypes attribute: ", nc_strerror(err), ", for file: ", particleInputFilename);

      // get varids
      int varidGlobalID;
      err = nc_inq_varid(ncID, "globalID", &varidGlobalID);
      log->check(err == NC_NOERR, "Problem getting varid for globalID: ", nc_strerror(err), ", for file: ", particleInputFilename);

      int varidType;
      err = nc_inq_varid(ncID, "type", &varidType);
      log->check(err == NC_NOERR, "Problem getting varid for type: ", nc_strerror(err), ", for file: ", particleInputFilename);

      int varidX;
      err = nc_inq_varid(ncID, "x", &varidX);
      log->check(err == NC_NOERR, "Problem getting varid for x: ", nc_strerror(err), ", for file: ", particleInputFilename);

      int varidRadius;
      err = nc_inq_varid(ncID, "radius", &varidRadius);
      log->check(err == NC_NOERR, "Problem getting varid for radius: ", nc_strerror(err), ", for file: ", particleInputFilename);

      int varidBonds;
      if (nBonds > 0) {
        err = nc_inq_varid(ncID, "bonds", &varidBonds);
        log->check(err == NC_NOERR, "Problem getting varid for bonds: ", nc_strerror(err), ", for file: ", particleInputFilename);
      }

      // get type and calculate nParticlesActive
      int *typeTmp = new int[nParticles];

      size_t startBufferType[1] = {(size_t) 0};
      size_t countBufferType[1] = {(size_t) nParticles};

      err = nc_get_vara_int(ncID, varidType, startBufferType, countBufferType, &typeTmp[0]);
      log->check(err == NC_NOERR, "Problem reading data for type: ", nc_strerror(err), ", for file: ", particleInputFilename);

      int nParticlesActive = 0;
      for (int iParticle=0 ; iParticle < nParticles ; iParticle++) {
	if (typeTmp[iParticle] > 0) nParticlesActive += 1;
      } // iParticle

      delete [] typeTmp;

      // read in buffered input

      // open lammps input file
      std::string filenameOutLammps = "lammps_particle_input.dat";
      std::ofstream lammpsInputFile;
      lammpsInputFile.open(filenameOutLammps);
      lammpsInputFile << std::setprecision(12);

      // lammps file header
      lammpsInputFile << "# particle input file" << std::endl << std::endl;
      lammpsInputFile << nParticlesActive << " atoms" << std::endl;
      if (nBonds > 0) lammpsInputFile << nBonds << " bonds" << std::endl;
      lammpsInputFile << std::endl;
      lammpsInputFile << nTypes << " atom types" << std::endl;
      if (nBonds > 0) lammpsInputFile << "1 bond types" << std::endl;
      lammpsInputFile << std::endl;
      lammpsInputFile << 0.0 << " " << grid->lx() << " xlo xhi" << std::endl;
      lammpsInputFile << 0.0 << " " << grid->ly() << " ylo yhi" << std::endl;
      lammpsInputFile << -maxRadius << " " << maxRadius << " zlo zhi" << std::endl << std::endl;

      // write by particle data
      lammpsInputFile << "Atoms # sphere" << std::endl << std::endl;

      // buffered input
      const int nParticlesPerBuffer = 1000;

      int globalID[nParticlesPerBuffer];
      int type[nParticlesPerBuffer];
      double x[nParticlesPerBuffer][2];
      double radius[nParticlesPerBuffer];

      // min/max diameter
      minimumParticleDiameter =  1.0e30;
      maximumParticleDiameter = -1.0e30;

      int nParticlesLeft = nParticles;
      int iBuffer = 0;
      while (nParticlesLeft > 0) {

        int nParticlesRead = std::min(nParticlesPerBuffer,nParticlesLeft);

        size_t startBuffer[1] = {(size_t) iBuffer*nParticlesPerBuffer};
        size_t countBuffer[1] = {(size_t) nParticlesRead};

        size_t startBuffer2D[2] = {(size_t) iBuffer*nParticlesPerBuffer, 0};
        size_t countBuffer2D[2] = {(size_t) nParticlesRead, 2};

        err = nc_get_vara_int(ncID, varidGlobalID, startBuffer, countBuffer, &globalID[0]);
        log->check(err == NC_NOERR, "Problem reading data for globalID: ", nc_strerror(err), ", for file: ", particleInputFilename);

        err = nc_get_vara_int(ncID, varidType, startBuffer, countBuffer, &type[0]);
        log->check(err == NC_NOERR, "Problem reading data for type: ", nc_strerror(err), ", for file: ", particleInputFilename);

        err = nc_get_vara_double(ncID, varidX, startBuffer2D, countBuffer2D, &x[0][0]);
        log->check(err == NC_NOERR, "Problem reading data for x: ", nc_strerror(err), ", for file: ", particleInputFilename);

        err = nc_get_vara_double(ncID, varidRadius, startBuffer, countBuffer, &radius[0]);
        log->check(err == NC_NOERR, "Problem reading data for radius: ", nc_strerror(err), ", for file: ", particleInputFilename);

        // write to lammps input file
        for (int iParticle=0 ; iParticle < nParticlesRead ; iParticle++) {

	  if (radius[iParticle]*2.0 < minimumParticleDiameter) minimumParticleDiameter = radius[iParticle]*2.0;
	  if (radius[iParticle]*2.0 > maximumParticleDiameter) maximumParticleDiameter = radius[iParticle]*2.0;

	  if (type[iParticle] > 0) {

	    // areal density (kg/m2)
	    double iceThickness = 1.0;
	    double iceFraction = 1.0;
	    double density = iceThickness * iceFraction * 900.0;

	    lammpsInputFile << globalID[iParticle] << " "; // atom id, start at 1
	    lammpsInputFile << type[iParticle] << " "; // atom type
	    lammpsInputFile << radius[iParticle]*2.0 << " "; // diameter
	    lammpsInputFile << density << " "; // density placeholder
	    lammpsInputFile << x[iParticle][0] << " "; // x position
	    lammpsInputFile << x[iParticle][1] << " "; // y position
	    lammpsInputFile << 0.0 << std::endl; // z position

	  } // real particle

        }
        nParticlesLeft = nParticlesLeft - nParticlesPerBuffer;
        iBuffer++;
      }

      // write by bond data
      if (nBonds > 0) {

        lammpsInputFile << std::endl << std::endl;
        lammpsInputFile << "Bonds" << std::endl << std::endl;

        // buffered input
        const int nBondsPerBuffer = 1000;

        int bonds[nBondsPerBuffer][2];

        int nBondsLeft = nBonds;
        int iBuffer = 0;
        while (nBondsLeft > 0) {

          int nBondsRead = std::min(nBondsPerBuffer,nBondsLeft);

          size_t startBuffer[2] = {(size_t) iBuffer*nBondsPerBuffer, 0};
          size_t countBuffer[2] = {(size_t) nBondsRead, 2};

          err = nc_get_vara_int(ncID, varidBonds, startBuffer, countBuffer, &bonds[0][0]);
          log->check(err == NC_NOERR, "Problem reading data for bonds: ", nc_strerror(err), ", for file: ", particleInputFilename);

          // write to lammps input file
          for (int iBond=0 ; iBond < nBondsRead ; iBond++) {

            // bond id, start at 1, 1, atomID1, atomID2
            int bondID = iBuffer * nBondsPerBuffer + iBond + 1;
            lammpsInputFile << bondID << " 1 " << bonds[iBond][0] << " " << bonds[iBond][1] << std::endl;

          }
          nBondsLeft = nBondsLeft - nBondsPerBuffer;
          iBuffer++;
        }
      }
      lammpsInputFile << std::endl << std::endl;

      lammpsInputFile.close();

    } // on master proc

    // mpi barrier
    char mpiErrBuffer[MPI_MAX_ERROR_STRING];
    int mpiErrLen;

    int err = MPI_Barrier(partition->comm());
    MPI_Error_string(err, mpiErrBuffer, &mpiErrLen);
    log->check(err == MPI_SUCCESS, "Problem with particle MPI barrier: ", (std::string) mpiErrBuffer);

  } // filename not "none"

}

// Setup lammps
void Particles::setup_lammps(void) {

  std::string lammpsSetupScript;

  if (configs->exists({"ConfigGroup:lammps","Config:lammpsSetupScript"})) {

    // get input script name
    configs->get({"ConfigGroup:lammps","Config:lammpsSetupScript"}, lammpsSetupScript);

  } else {

    // DEMSI setup lammps
    lammpsSetupScript = "lammps_setup.in";

    // root processor makes the file
    if (partition->on_master_proc()) {

      std::ofstream lammpsSetupFile;
      lammpsSetupFile.open(lammpsSetupScript);

      // general settings
      lammpsSetupFile << "atom_style demsi" << std::endl;
      if (contacts->bonds_active()) lammpsSetupFile << "bond_style zero" << std::endl;
      lammpsSetupFile << "dimension 2" << std::endl;
      lammpsSetupFile << "newton off" << std::endl;
      lammpsSetupFile << "units si" << std::endl;

      // create box
      bool xPeriodic = false;
      if (configs->exists({"ConfigGroup:lammps","Config:xPeriodic"})) {
	configs->get({"ConfigGroup:lammps","Config:xPeriodic"}, xPeriodic);
      }
      bool yPeriodic = false;
      if (configs->exists({"ConfigGroup:lammps","Config:yPeriodic"})) {
	configs->get({"ConfigGroup:lammps","Config:yPeriodic"}, yPeriodic);
      }
      std::string xPeriodicStr = xPeriodic ? "p" : "f";
      std::string yPeriodicStr = yPeriodic ? "p" : "f";
      lammpsSetupFile << "boundary " << xPeriodicStr << " " << yPeriodicStr << " p" << std::endl;

      // map atoms
      if (contacts->bonds_active()) lammpsSetupFile << "atom_modify map yes" << std::endl;

      // add particles
      lammpsSetupFile << "read_data lammps_particle_input.dat" << std::endl;

      // bond coefficients
      if (contacts->bonds_active()) lammpsSetupFile << "bond_coeff * 0.0" << std::endl;

      // time step
      std::ostringstream timestepCommand;
      timestepCommand << std::setprecision(16);
      timestepCommand << "timestep " << contacts->timestep(minimumParticleDiameter);
      lammpsSetupFile << timestepCommand.str() << std::endl;

      // pair style
      if (contacts->get_contact_type() != DEMSI::CONTACTTYPE::NONE) {
	lammpsSetupFile << contacts->pair_style_command() << std::endl;
	lammpsSetupFile << "pair_coeff * *" << std::endl;
      }

      // groups
      if (not configs->exists({"ConfigGroup:lammps","Array:lammpsGroups"})) {

	// mobile elements group
	lammpsSetupFile << "group mobile type 1" << std::endl;

	// coastline
	bool useCoastlineElements;
	configs->get({"ConfigGroup:lammps", "Config:useCoastlineElements"}, useCoastlineElements);
	if (useCoastlineElements) {
	  lammpsSetupFile << "group coastline type 2" << std::endl;
	}

      } else {

	// override groups
	std::vector<std::string> lammpsGroups = configs->get_array({"ConfigGroup:lammps","Array:lammpsGroups"}, "LammpsGroup");
	for (int iGroup = 0 ; iGroup < lammpsGroups.size() ; iGroup++) {
	  lammpsSetupFile << lammpsGroups[iGroup] << std::endl;
	} // iGroup*/

      } // groups

      // fixes
      // enforce 2D motion
      lammpsSetupFile << "fix 1 all enforce2d" << std::endl;

      // integrate the equations of motion for the elements
      bool integrateMotion;
      configs->get({"ConfigGroup:lammps","Config:integrateMotion"}, integrateMotion);
      if (integrateMotion) lammpsSetupFile << "fix 2 mobile nve/sphere/demsi" << std::endl;

      // domain boundaries
      if (configs->exists({"ConfigGroup:lammps","Config:xHasWalls"}) or
	  configs->exists({"ConfigGroup:lammps","Config:yHasWalls"})) {

	bool xHasWalls, yHasWalls;
	configs->get({"ConfigGroup:lammps","Config:xHasWalls"}, xHasWalls);
	configs->get({"ConfigGroup:lammps","Config:yHasWalls"}, yHasWalls);
	std::string wallCommand;
	if (xHasWalls or yHasWalls) wallCommand = contacts->wall_command();
	if (xHasWalls) {
	  std::ostringstream xWallCommand;
	  xWallCommand << "fix wx all " << wallCommand << " xplane " << 0.0 << " " << grid->lx();
	  lammpsSetupFile << xWallCommand.str() << std::endl;
	}
	if (yHasWalls) {
	  std::ostringstream yWallCommand;
	  yWallCommand << "fix wy all " << wallCommand << " yplane " << 0.0 << " " << grid->ly();
	  lammpsSetupFile << yWallCommand.str() << std::endl;
	}

      } // has walls

      // turn off coast to coast interactions for coastline elements
      bool useCoastlineElements;
      configs->get({"ConfigGroup:lammps", "Config:useCoastlineElements"}, useCoastlineElements);
      if (useCoastlineElements) {
	lammpsSetupFile << "neigh_modify exclude group coastline coastline" << std::endl;
      }

      // turn off sorting if no pair wise contact method specified
      if (contacts->get_contact_type() == DEMSI::CONTACTTYPE::NONE) {
	lammpsSetupFile << "atom_modify sort 0 1" << std::endl;
      }

      // extra fixes
      if (configs->exists({"ConfigGroup:lammps","Array:lammpsFixes"})) {

	std::vector<std::string> lammpsFixes = configs->get_array({"ConfigGroup:lammps","Array:lammpsFixes"}, "LammpsFix");
	for (int iFix = 0 ; iFix < lammpsFixes.size() ; iFix++) {
	  lammpsSetupFile << lammpsFixes[iFix] << std::endl;
	} // iFix

      } // fixes

      // skin thickness
      double skinFraction;
      configs->get({"ConfigGroup:lammps","Config:skinFraction"}, skinFraction);
      double skin = skinFraction * maximumParticleDiameter;
      std::ostringstream skinCommand;
      skinCommand << "neighbor " << skin << " bin";
      lammpsSetupFile << skinCommand.str() << std::endl;

      // set neighbor list update settings
      if (contacts->bonds_active()) lammpsSetupFile << "neigh_modify delay 0 every 1 check yes" << std::endl;

      // always communicate particle velocity info
      std::ostringstream commModifyVelCommand;
      commModifyVelCommand << "comm_modify vel yes cutoff " << maximumParticleDiameter;
      lammpsSetupFile << commModifyVelCommand.str() << std::endl;

      // what to do with lost atoms
      std::string lostAtomsType;
      configs->get({"ConfigGroup:lammps", "Config:lostAtomsType"}, lostAtomsType);
      std::ostringstream lostAtomsCommand;
      lostAtomsCommand << "thermo_modify lost " << lostAtomsType;
      lammpsSetupFile << lostAtomsCommand.str() << std::endl;

      // extra commands
      if (configs->exists({"ConfigGroup:lammps","Array:extraCommands"})) {

	std::vector<std::string> extraCommands = configs->get_array({"ConfigGroup:lammps","Array:extraCommands"}, "ExtraCommand");
	for (int iExtra = 0 ; iExtra < extraCommands.size() ; iExtra++) {
	  lammpsSetupFile << extraCommands[iExtra] << std::endl;
	} // iExtra

      } // extra commands

      lammpsSetupFile.close();

    } // on master proc

  }

  // run the setup script
  MPI_Barrier(partition->comm());
  lammpsInstance->run_script(lammpsSetupScript);

} // Particles::setup_lammps

// check all particles lie within the domain
void Particles::check_in_domain(void) {

  // ensure host is up-to-date
  lammpsInstance->lmp->atomKK->k_x.sync<host_execution_space>();

  // get host view of x
  auto x_host_mirror = lammpsInstance->lmp->atomKK->k_x.h_view;

  Kokkos::parallel_for(Kokkos::RangePolicy<host_execution_space>(0,*(nParticles)), KOKKOS_LAMBDA(const int iParticle) {

    log->check(x_host_mirror(iParticle,0) >= 0.0 and x_host_mirror(iParticle,0) <= grid->lx() and x_host_mirror(iParticle,1) >= 0.0 and x_host_mirror(iParticle,1) <= grid->ly(), "Particle ", std::to_string(iParticle), " has left the domain");

  }); // iParticle

}

// initial particle velocity
void Particles::initial_velocity(void) {

  // check if uniform initial velocity option exists
  if (configs->exists({"ConfigGroup:initialization","Config:uVelocityUniform"}) and
      configs->exists({"ConfigGroup:initialization","Config:vVelocityUniform"})) {

    double uVelocityUniform, vVelocityUniform;
    configs->get({"ConfigGroup:initialization","Config:uVelocityUniform"}, uVelocityUniform);
    configs->get({"ConfigGroup:initialization","Config:vVelocityUniform"}, vVelocityUniform);

    auto this_v = this->v;
    Kokkos::parallel_for("initial_velocity", Kokkos::RangePolicy<device_execution_space>(0,*(nParticles)),
        KOKKOS_LAMBDA(const int iParticle) {
      this_v(iParticle,0) = uVelocityUniform;
      this_v(iParticle,1) = vVelocityUniform;
    }); // iParticle

  } // uniform initial velocity

}

} // DEMSI namespace
