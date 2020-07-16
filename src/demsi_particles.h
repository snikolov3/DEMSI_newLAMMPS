/*
  DEMSI particle structure

 */

#ifndef DEMSI_PARTICLES_H
#define DEMSI_PARTICLES_H

#include "demsi_typedefs.h"
#include "demsi_lmp_instance.h"
#include "demsi_configs.h"
#include "demsi_grid.h"
#include "demsi_polygon.h"
#include "demsi_contacts.h"

#include <Kokkos_Core.hpp>
#include <lammps/lammps.h>
#include <lammps/atom_kokkos.h>
#include <lammps/math_const.h>

namespace DEMSI {

//------------------------------------------------------------------------------
/*! \class Particles
    \brief Class defining the DEMSI data structure for particles.

    The class stores data associated with particles in DEMSI.
    It holds a pointer to an existing LAMMPS instance. The particles
    can be created by LAMMPS, or they can be
    appended to an existing LAMMPS instance (including one that does
    not have any particles).
*/
//------------------------------------------------------------------------------
class Particles{
public:

  /*! \brief Constructor for Particles class */
  Particles(DEMSI::Partition* partitionIn, DEMSI::Log* logIn, DEMSI::Configs* configsIn, DEMSI::LammpsInstance* lammpsInstanceIn, DEMSI::Grid* gridIn, DEMSI::Contacts* contactsIn);

  /*! \brief Default destructor for Particles class */
  ~Particles() = default;

  /*! Number of particles - changes over time so must be pointer. */
  int *nParticles;

  /*! Density of ice */
  double iceDensity;

  /*! Unique global ID, as assigned by LAMMPS */
  kokkos_view_type_tagint_host globalID;

  /*! Element type */
  kokkos_view_type_int_host type;

  /*! Positions, where x(i,0) is the x-coordinate of particle i, x(i,1) its y-coordinate */
  kokkos_view_type_x x;

  /*! Velocities */
  kokkos_view_type_v v;

  /*! Particle radii */
  kokkos_view_type_1d_float radius;

  /*! Particle masses */
  kokkos_view_type_1d_float mass;

  /*! Body forces */
  kokkos_view_type_2d_float forcing;

  /*! Minimum ice thickness in element */
  kokkos_view_type_1d_float minThickness;

  /*! Ice thickness for ridging */
  kokkos_view_type_1d_float ridgingIceThickness;

  /*! Weight for ice thickness for ridging */
  kokkos_view_type_1d_float ridgingIceThicknessWeight;

  /*! Net to gross closing ratio */
  kokkos_view_type_1d_float netToGrossClosingRatio;

  /*! change in effective element area */
  kokkos_view_type_1d_float changeEffectiveElementArea;

  /*! Mean ice thickness in element */
  kokkos_view_type_1d_float meanThickness;

  /*! Ice area*/
  kokkos_view_type_1d_float ice_area;

  /*! Coriolis parameter*/
  kokkos_view_type_1d_float coriolis;

  /*! Ocean velocity*/
  kokkos_view_type_2d_float ocean_vel;

  /*! B-vector for semi-implicit integration*/
  kokkos_view_type_2d_float bvector;

  /*! Sync data in LAMMPS to reflect changes in DEMSI */
  void sync_data_to_lammps(void);

  /*! Set the DEMSI pointer to the LAMMPS particle data */
  void set_demsi_pointers(void);

  /*! Placeholder for bond creation */
  void createBond(int iParticle, int jParticle, double elasticModulus, double compressiveStrength, double tensileStrength);

  /*! \brief Create LAMMPS particle input file from DEMSI netcdf particle file and load in LAMMPS.
   */
  void create_lammps_particles(void);

  /*! \brief Setup lammps
   */
  void setup_lammps(void);

  /*! \brief Check if all particles are within the domain
   */
  void check_in_domain(void);

  /*! \brief Set the initial particle velocity if other than zero
   */
  void initial_velocity(void);

  /*! Pointer to lammps instance */
  DEMSI::LammpsInstance* lammpsInstance;

private:

  /*! Pointer to partition object */
  DEMSI::Partition* partition;

  /*! Pointer to Log object */
  DEMSI::Log* log;

  /*! Pointer to configs object */
  DEMSI::Configs* configs;

  /*! Pointer to grid object */
  DEMSI::Grid* grid;

  /*! Pointer to contacts object */
  DEMSI::Contacts* contacts;

  /*! min/max particle diameter */
  double minimumParticleDiameter, maximumParticleDiameter;

};

}

#endif
