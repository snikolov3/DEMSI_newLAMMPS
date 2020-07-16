#ifndef DEMSI_CONTACTS_H_
#define DEMSI_CONTACTS_H_

#include "demsi_typedefs.h"
#include "demsi_lmp_instance.h"
#include "demsi_logging.h"
#include "demsi_configs.h"

// these are LAMMPS include files
#include <lammps/kokkos_type.h>
#include <lammps/atom_kokkos.h>

#include <vector>

/*! \file   demsi_contacts.h
  \brief  Header file for the DEMSI::Contacts class.
*/

namespace DEMSI {

/*! /enum contact type */
enum CONTACTTYPE {HOPKINS, HOOKE, HOOKE_THICKNESS, NONE};

/*! \class Contacts
 \brief Class for particle contacts.
*/
class Contacts {

public:

  /*! Constructor for Contacts class.
    /param logIn Pointer to log object.
    /param configs Pointer to the configs object.
    /param lammpsInstanceIn Pointer to lammps instance object.
   */
  Contacts(DEMSI::Log* logIn, DEMSI::Configs* configs, DEMSI::LammpsInstance* lammpsInstanceIn);

  /*! Default Contacts destructor */
  ~Contacts() = default;

  /*! \brief Get contact time step
      \param minimumParticleDiameter minimum particle diameter
      \return The contact time step
  */
  double timestep(const double minimumParticleDiameter);

  /*! \brief Return the lammps pair style command for the contact.
      \return The lammps pair style command for the contact.
  */
  std::string pair_style_command(void);

  /*! \brief Return the contact part of the lammps wall fix command for the contact.
      \return The contact part of the lammps wall fix command for the contact.
  */
  std::string wall_command(void);

  /*! \brief Initialize bonded contacts for the Hopkins contact model */
  void initialize_bonded_contacts_hopkins(void);

  /*! \brief Get the bond information for the hopkins contact model.
      \param bondGlobalIDs Vector of Global IDs of the two particles in the bond.
      \param bondLength Vector of bond lengths.
      \param bondThickness Vector of bond thicknesses.
      \param bondCrackFraction Vector of the bond crack fractions on either end of the bond.
      \param bondEndPoint1Particle1 Vector of (x,y) positions of the first end of the bond on the first particle side.
      \param bondEndPoint2Particle1 Vector of (x,y) positions of the second end of the bond on the first particle side.
      \param bondEndPoint1Particle2 Vector of (x,y) positions of the first end of the bond on the second particle side.
      \param bondEndPoint2Particle2 Vector of (x,y) positions of the second end of the bond on the second particle side.
   */

  void get_bonds_info_hopkins(Kokkos::View<LAMMPS_NS::tagint*[2]> &bondGlobalIDs, Kokkos::View<double*> &bondLength, Kokkos::View<double*> &bondThickness, Kokkos::View<double*[2]> &bondCrackFraction, Kokkos::View<double*[2]> &bondEndPoint1Particle1, Kokkos::View<double*[2]> &bondEndPoint2Particle1, Kokkos::View<double*[2]> &bondEndPoint1Particle2, Kokkos::View<double*[2]> &bondEndPoint2Particle2);

  /*! \brief Return the current contact type.
      \return The current contact type.
   */
  DEMSI::CONTACTTYPE get_contact_type(void);

  /*! \brief Return whether bonds are currently active.
      \return If true bonds are currently active.
  */
  bool bonds_active(void);

  /*! \brief Return whether newton pair option is on.
      \return True if newton pair option is on.
  */
  bool newton_pair_used(void);

  /*! \brief Return whether half neighbour lists are being used.
      \return True if half neighbour lists used.
  */
  bool half_lists_used(void);

  /*! \brief Return whether elements are on different processors
      \param Local ID of atom 1
      \param Local ID of atom 2
      \return True if elements are on different processors
  */
  KOKKOS_INLINE_FUNCTION
  bool on_different_procs(int, int, int);

private:

  LAMMPS_NS::ExecutionSpace execution_space;

  /*! Pointer to the log object */
  DEMSI::Log* log;

  /*! Pointer to configs object */
  DEMSI::Configs* configs;

  /*! Pointer to lammps instance object */
  DEMSI::LammpsInstance* lammpsInstance;

  /*! contact type: "hopkins", "hooke", "hooke_thickness", "none" */
  std::string contactType;

};

} // namespace DEMSI

#endif /* CONTACTS_H_ */
