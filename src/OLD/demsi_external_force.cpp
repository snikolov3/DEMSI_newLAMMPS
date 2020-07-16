#include "demsi_external_force.h"
#include "demsi_logging.h"
#include "demsi_column.h"
#include "demsi_configs.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include "demsi_forcing.h"

#include <lammps/atom.h>
#include <lammps/math_const.h>
#include <lammps/modify.h>
#include <lammps/fix_nve_sphere_demsi.h>

namespace DEMSI {

  ExternalForce::ExternalForce(DEMSI::Log* logIn, DEMSI::Configs *configsIn, DEMSI::Grid* gridIn, DEMSI::Particles* particlesIn, DEMSI::Forcing* forcingIn) {

    log = logIn;
    configs = configsIn;
    grid = gridIn;
    particles = particlesIn;
    forcing = forcingIn;

    // switches for different forces
    configs->get({"ConfigGroup:useSections","Config:useExternalForces"}, useExternalForces);

    if (useExternalForces) {

      bool useAirStress        = false;
      bool useOceanStress      = false;
      bool useCoriolisForce    = false;
      bool useSurfaceTiltForce = false;

      configs->get({"ConfigGroup:externalForces","Config:useAirStress"}, useAirStress);
      configs->get({"ConfigGroup:externalForces","Config:useOceanStress"}, useOceanStress);
      configs->get({"ConfigGroup:externalForces","Config:useCoriolisForce"}, useCoriolisForce);
      configs->get({"ConfigGroup:externalForces","Config:useSurfaceTiltForce"}, useSurfaceTiltForce);

      // get physical parameters
      if (useAirStress) {
        configs->get({"ConfigGroup:externalForces","Config:atmDensity"}, atmDensity);
        configs->get({"ConfigGroup:externalForces","Config:windDragCoefficient"}, windDragCoefficient);
      } else {
        atmDensity = 0.0;
        windDragCoefficient = 0.0;
      }

      if (useOceanStress) {
        configs->get({"ConfigGroup:externalForces","Config:ocnDensity"}, ocnDensity);
        configs->get({"ConfigGroup:externalForces","Config:oceanDragCoefficient"}, oceanDragCoefficient);
      } else {
        ocnDensity = 0.0;
        oceanDragCoefficient = 0.0;
      }

      //temporary, while testing semi-implicit integration
      //Find integration fix
      class LAMMPS_NS::FixNVESphereDemsi *fix_integrator;
      int index = particles->lammpsInstance->lmp->modify->find_fix_by_style("nve/sphere/demsi");
      fix_integrator = (LAMMPS_NS::FixNVESphereDemsi *) particles->lammpsInstance->lmp->modify->fix[index];
      fix_integrator->ocean_density = ocnDensity;
      fix_integrator->ocean_drag = oceanDragCoefficient;

      if (useCoriolisForce or useSurfaceTiltForce) {
        configs->get({"ConfigGroup:externalForces","Config:earthRotationRate"}, earthRotationRate);
      } else {
        earthRotationRate = 0.0;
      }

      if (useCoriolisForce) {
        coriolisForceMask = 1.0;
      } else {
        coriolisForceMask = 0.0;
      }

      if (useSurfaceTiltForce) {
        surfaceTiltForceMask = 1.0;
      } else {
        surfaceTiltForceMask = 0.0;
      }

    }

  }


  /*==============================================================
   * Extract particle positions and velocities from lammps, extract
   * wind velocity from grid, pass forcing back to lammps.
   *==============================================================*/
  void ExternalForce::add_forces() {

    if (useExternalForces) {

      // prepare variables needed in loop so that they can be
      // capture by KOKKOS_LAMBDA
      auto forcing_xAtmWind = forcing->xAtmWind;
      auto forcing_yAtmWind = forcing->yAtmWind;
      auto forcing_xOcnCurrents = forcing->xOcnCurrents;
      auto forcing_yOcnCurrents = forcing->yOcnCurrents;
      auto grid_latitude = grid->latitude;
      auto particles_forcing = particles->forcing;
      auto particles_x = particles->x;
      auto particles_v = particles->v;
      auto particles_mass = particles->mass;
      auto particles_radius = particles->radius;
      auto particles_coriolis = particles->coriolis;
      auto particles_bvector = particles->bvector;
      auto particles_ocean_vel = particles->ocean_vel;
      auto particles_ice_area = particles->ice_area;

      const double coriolisForceMask = this->coriolisForceMask;
      const double surfaceTiltForceMask = this->surfaceTiltForceMask;
      const double earthRotationRate = this->earthRotationRate;
      const double atmDensity = this->atmDensity;
      const double windDragCoefficient = this->windDragCoefficient;
      const double ocnDensity = this->ocnDensity;
      const double oceanDragCoefficient = this->oceanDragCoefficient;

      Kokkos::View<int*[4]> i("i", particles_x.dimension_0());
      Kokkos::View<int*[4]> j("j", particles_x.dimension_0());
      Kokkos::View<double*[4]> weight("weight", particles_x.dimension_0());

      // get interpolation weights for all particles
      grid->interpolation_weights(particles_x, i, j, weight);

      Kokkos::parallel_for("ExternalForce::add_forces", Kokkos::RangePolicy<device_execution_space>
          (0, *(particles->nParticles)), KOKKOS_LAMBDA (const int iParticle) {

        //Get particle position, velocity
        double x = particles_x(iParticle,0);
        double y = particles_x(iParticle,1);
        double m = particles_mass(iParticle);
        double vx = particles_v(iParticle,0);
        double vy = particles_v(iParticle,1);
        double radius = particles_radius(iParticle);
	double iceArea = particles_ice_area(iParticle);
        double xAtmWind, yAtmWind;
        double xOcnCurrents, yOcnCurrents;
        double elementLatitude;
        double fx, fy;

        //Apply interpolation weights previously calculated

        xAtmWind =
          forcing_xAtmWind(i(iParticle,0),j(iParticle,0)) * weight(iParticle,0) +
          forcing_xAtmWind(i(iParticle,1),j(iParticle,1)) * weight(iParticle,1) +
          forcing_xAtmWind(i(iParticle,2),j(iParticle,2)) * weight(iParticle,2) +
          forcing_xAtmWind(i(iParticle,3),j(iParticle,3)) * weight(iParticle,3);

        yAtmWind =
          forcing_yAtmWind(i(iParticle,0),j(iParticle,0)) * weight(iParticle,0) +
          forcing_yAtmWind(i(iParticle,1),j(iParticle,1)) * weight(iParticle,1) +
          forcing_yAtmWind(i(iParticle,2),j(iParticle,2)) * weight(iParticle,2) +
          forcing_yAtmWind(i(iParticle,3),j(iParticle,3)) * weight(iParticle,3);

        xOcnCurrents =
          forcing_xOcnCurrents(i(iParticle,0),j(iParticle,0)) * weight(iParticle,0) +
          forcing_xOcnCurrents(i(iParticle,1),j(iParticle,1)) * weight(iParticle,1) +
          forcing_xOcnCurrents(i(iParticle,2),j(iParticle,2)) * weight(iParticle,2) +
          forcing_xOcnCurrents(i(iParticle,3),j(iParticle,3)) * weight(iParticle,3);

        yOcnCurrents =
          forcing_yOcnCurrents(i(iParticle,0),j(iParticle,0)) * weight(iParticle,0) +
          forcing_yOcnCurrents(i(iParticle,1),j(iParticle,1)) * weight(iParticle,1) +
          forcing_yOcnCurrents(i(iParticle,2),j(iParticle,2)) * weight(iParticle,2) +
          forcing_yOcnCurrents(i(iParticle,3),j(iParticle,3)) * weight(iParticle,3);

        elementLatitude =
          grid_latitude(i(iParticle,0),j(iParticle,0)) * weight(iParticle,0) +
          grid_latitude(i(iParticle,1),j(iParticle,1)) * weight(iParticle,1) +
          grid_latitude(i(iParticle,2),j(iParticle,2)) * weight(iParticle,2) +
          grid_latitude(i(iParticle,3),j(iParticle,3)) * weight(iParticle,3);

        // Compute force
        fx = 0.0 ; fy = 0.0;
        compute_surface_stress(vx, vy, iceArea, xAtmWind, yAtmWind, atmDensity, windDragCoefficient, fx, fy);
        const double coriolisParameter = 2.0 * earthRotationRate * std::sin(elementLatitude);

        /* Temporary, while testing semi-implicit integration*/
	particles_coriolis(iParticle) = coriolisForceMask*coriolisParameter;

        compute_surface_tilt_force(surfaceTiltForceMask, m, coriolisParameter, xOcnCurrents, yOcnCurrents, fx, fy);

	particles_bvector(iParticle,0) = fx;
	particles_bvector(iParticle,1) = fy;
	particles_ocean_vel(iParticle,0) = xOcnCurrents;
	particles_ocean_vel(iParticle,1) = yOcnCurrents;

      });
    }

    // constant stress gradient
    if (configs->exists({"ConfigGroup:constantStressGradient","Config:xStressGradient"})) {

      double xStressGradient, yStressGradient;
      configs->get({"ConfigGroup:constantStressGradient","Config:xStressGradient"}, xStressGradient);
      configs->get({"ConfigGroup:constantStressGradient","Config:yStressGradient"}, yStressGradient);

      auto grid_lx = grid->lx();
      auto grid_ly = grid->ly();
      auto particles_x = particles->x;
      auto particles_radius = particles->radius;
      auto particles_forcing = particles->forcing;
      Kokkos::parallel_for("ExternalForce::add_constant_stress", Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>
      (0, *(particles->nParticles)), KOKKOS_LAMBDA (const int iParticle) {

        double elementArea = LAMMPS_NS::MathConst::MY_PI * pow(particles_radius(iParticle),2);

        particles_forcing(iParticle,0) = xStressGradient * elementArea * (particles_x(iParticle,0) - 0.5 * grid_lx);
        particles_forcing(iParticle,1) = yStressGradient * elementArea * (particles_x(iParticle,1) - 0.5 * grid_ly);

      }); // iParticle

    } // constant stress gradient

  }

  /*==============================================================
   * Compute surface forces on single particle given its position,
   * velocity and fluid velocity at its position
   *==============================================================*/
  KOKKOS_INLINE_FUNCTION
  void ExternalForce::compute_surface_stress(
    const double xElementVelocity, const double yElementVelocity, // element position
    const double elementArea, // element area
    const double xFluidVelocity, const double yFluidVelocity, // fluid velocity
    const double fluidDensity, const double dragCoefficient,
    double &xForce, double &yForce) const { // wind and ocean force on element

    // speed of wind relative to element
    double relativeSpeed = sqrt(pow(xFluidVelocity-xElementVelocity,2) +
				pow(yFluidVelocity-yElementVelocity,2));

    // wind stress on element
    double xStressForce = fluidDensity * dragCoefficient * relativeSpeed * xFluidVelocity;
    double yStressForce = fluidDensity * dragCoefficient * relativeSpeed * yFluidVelocity;

    // wind and ocean force on element
    xForce = xForce + xStressForce * elementArea;
    yForce = yForce + yStressForce * elementArea;
  }

  /*==============================================================
   * Compute Coriolis force on an element
   *==============================================================*/
  KOKKOS_INLINE_FUNCTION
  void ExternalForce::compute_coriolis_force(
    const double coriolisForceMask,
    const double elementMass, const double coriolisParameter,
    const double xElementVelocity, const double yElementVelocity,
    double &xForce, double &yForce) const {

    xForce = xForce + coriolisForceMask * coriolisParameter * elementMass * yElementVelocity;
    yForce = yForce - coriolisForceMask * coriolisParameter * elementMass * xElementVelocity;
  }

  /*==============================================================
   * Compute surface tilt force on an element
   *==============================================================*/
  KOKKOS_INLINE_FUNCTION
  void ExternalForce::compute_surface_tilt_force(
    const double surfaceTiltForceMask,
    const double elementMass, const double coriolisParameter,
    const double xOcnCurrents, const double yOcnCurrents,
    double &xForce, double &yForce) const {

    xForce = xForce - surfaceTiltForceMask * coriolisParameter * elementMass * yOcnCurrents;
    yForce = yForce + surfaceTiltForceMask * coriolisParameter * elementMass * xOcnCurrents;
  }

} // namespace DEMSI
