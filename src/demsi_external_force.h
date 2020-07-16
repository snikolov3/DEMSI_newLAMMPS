#ifndef DEMSI_EXTERNAL_FORCE_H_
#define DEMSI_EXTERNAL_FORCE_H_

/*! \file   demsi_external_force.h
    \brief  Header file for the DEMSI::ExternalForce class.
*/

#include "demsi_logging.h"
#include "demsi_lmp_instance.h"
#include "demsi_particles.h"
#include "demsi_forcing.h"

#include  <Kokkos_Core.hpp>


namespace DEMSI {
/*! \class ExternalForce
 \brief Class allowing adding of external forces to particles.

 This class allows the addition of external body forces to particles. These
 forces include a force due to drag from atmospheric winds.
*/
class ExternalForce{
public:

  /*! \brief ExternalForce class Constructor
      \param logIn         [in] Pointer to log object.
      \param configs       [in] Pointer to configs object.
      \param gridIn        [in] Pointer to grid object.
      \param particlesIn   [in] Pointer to particles object.
      \param forcingIn     [in] Pointer to forcing object.
   */
  ExternalForce(DEMSI::Log* logIn, DEMSI::Configs *configsIn, DEMSI::Grid* gridIn, DEMSI::Particles* particlesIn, DEMSI::Forcing* forcingIn);

  /*! \brief Default destructor */
  ~ExternalForce() = default;

  /*! \brief Add body forces to LAMMPS particles */
  void add_forces();

  /*! \brief Compute body forces from surface stresses on a single particle
      \param xElementVelocity    [in] - element x velocity
      \param yElementVelocity    [in] - element y velocity
      \param elementArea         [in] - element area
      \param xFluidVelocity      [in] - fluid x velocity
      \param yFluidVelocity      [in] - fluid y velocity
      \param fluidDensity        [in] - fluid density
      \param dragCoefficient     [in] - fluid drag coefficient
      \param xForce             [out] - force in x direction
      \param yForce             [out] - force in y direction
   */
  KOKKOS_INLINE_FUNCTION
  void compute_surface_stress(
    const double xElementVelocity, const double yElementVelocity, // element position
    const double elementArea, // element area
    const double xFluidVelocity, const double yFluidVelocity, // fluid velocity
    const double fluidDensity, const double dragCoefficient,
    double &xForce, double &yForce) const;


  /*! \brief Compute Coriolis forces on a single particle
      \param coriolisForceMask    [in] - coriolis force mask
      \param elementMass          [in] - element mass
      \param coriolisParameter    [in] - element Coriolis parameter
      \param xElementVelocity     [in] - element x velocity
      \param yElementVelocity     [in] - element y velocity
      \param xForce              [out] - force in x direction
      \param yForce              [out] - force in y direction
   */
  KOKKOS_INLINE_FUNCTION
  void compute_coriolis_force(
    const double coriolisForceMask,
    const double elementMass, const double coriolisParameter,
    const double xElementVelocity, const double yElementVelocity,
    double &xForce, double &yForce) const;

  /*! \brief Compute surface tilt forces on a single particle
      \param surfaceTiltForceMask [in] - surface force tilt mask
      \param elementMass          [in] - element mass
      \param coriolisParameter    [in] - element Coriolis parameter
      \param xOcnCurrents         [in] - Ocean surface x velocity
      \param yOcnCurrents         [in] - Ocean surface y velocity
      \param xForce              [out] - force in x direction
      \param yForce              [out] - force in y direction
   */
  KOKKOS_INLINE_FUNCTION
  void compute_surface_tilt_force(
    const double surfaceTiltForceMask,
    const double elementMass, const double coriolisParameter,
    const double xOcnCurrents, const double yOcnCurrents,
    double &xForce, double &yForce) const;

private:

  // Pointer to the logging object.
  DEMSI::Log* log;

  // Pointer to the configs object.
  DEMSI::Configs* configs;

  // Pointer to grid object
  DEMSI::Grid* grid;

  // Pointer to particles object.
  DEMSI::Particles* particles;

  // Pointer to the forcing object.
  DEMSI::Forcing* forcing;

  // Physical parameters for external forces.
  double atmDensity; // atmospheric density density (Kg/m3)
  double ocnDensity; // ocean density (Kg/m3)
  double windDragCoefficient; // Wind drag coefficient
  double oceanDragCoefficient; // Ocean drag coefficient
  double earthRotationRate; // Earth rotation rate (1/s)
  double coriolisForceMask; // 1 if Coriolis force used, 0 otherwise
  double surfaceTiltForceMask; // 1 if Surface tilt force used, 0 otherwise

  /*! If true use external forces, otherwise do not */
  bool useExternalForces;

};

} // namespace DEMSI

#endif /* GRID_H_ */
