#include "demsi_contacts.h"

#include "demsi_lmp_instance.h"
#include "demsi_logging.h"
#include "demsi_configs.h"

#include <mpi.h>
#include <cmath>
#include <array>
#include <stdio.h>
#include <lammps/lammps.h>         // these are LAMMPS include files
#include <lammps/input.h>
#include <lammps/kokkos_type.h>
#include <lammps/atom.h>
#include <lammps/atom_kokkos.h>
#include <lammps/atom_masks.h>
#include <lammps/library.h>
#include <lammps/update.h>
#include <lammps/math_const.h>

#include <lammps/fix_neigh_history.h>
#include <lammps/fix_neigh_history_kokkos.h>
#include <lammps/modify.h>
#include <lammps/neighbor.h>
#include <lammps/neighbor_kokkos.h>
#include <lammps/neigh_list.h>
#include <lammps/neigh_list_kokkos.h>
#include <lammps/pair.h>
#include <lammps/pair_kokkos.h>
#include <lammps/force.h>
#include <lammps/neigh_request.h>

namespace DEMSI {

// Constructor for Contacts class.
Contacts::Contacts(DEMSI::Log* logIn, DEMSI::Configs* configsIn, DEMSI::LammpsInstance* lammpsInstanceIn) {

  log = logIn;
  configs = configsIn;
  lammpsInstance = lammpsInstanceIn;

  // get contact type
  configs->get({"ConfigGroup:contacts","Config:contactType"}, contactType);

}

// Get contact time step
double Contacts::timestep(const double minimumParticleDiameter) {

  if (get_contact_type() != DEMSI::CONTACTTYPE::NONE) {

    double normalStiffness;

    if (get_contact_type() == DEMSI::CONTACTTYPE::HOOKE_THICKNESS or
	get_contact_type() == DEMSI::CONTACTTYPE::HOOKE) {

      // hooke contact model
      configs->get({"ConfigGroup:contacts","Config:normalStiffness"}, normalStiffness);

    } else if (get_contact_type() == DEMSI::CONTACTTYPE::HOPKINS) {

      // hopkins contact model
      double elasticModulus, poissonRatio;

      configs->get({"ConfigGroup:contacts","Config:elasticModulus"}, elasticModulus);
      configs->get({"ConfigGroup:contacts","Config:poissonRatio"}, poissonRatio);

      double shearModulus = elasticModulus / (2.0 * (1.0 + poissonRatio));
      normalStiffness = 4.0 * shearModulus / (3.0 * (1.0 - poissonRatio));

    }

    auto PI = LAMMPS_NS::MathConst::MY_PI;

    double timestepThickness;
    configs->get({"ConfigGroup:lammps","Config:timestepThickness"}, timestepThickness);

    double minimumParticleMass = PI * std::pow(0.5 * minimumParticleDiameter, 2);
    minimumParticleMass = minimumParticleMass * 900.0 * timestepThickness;

    double timeStep = PI / std::sqrt((2.0 * normalStiffness) / minimumParticleMass);
    timeStep = timeStep * 0.02;

    return timeStep;

  } else {

    // no contact model
    double timeStep;
    configs->get({"ConfigGroup:lammps", "Config:timeStep"}, timeStep);
    return timeStep;

  }

}

// Return the lammps pair style command for the contact.
std::string Contacts::pair_style_command(void) {

  std::ostringstream pairStyleCommand;

  if (get_contact_type() == DEMSI::CONTACTTYPE::HOPKINS) {

    // Hopkins contact model
    std::string compressiveBreakingStressType, tensileBreakingStressType;
    double elasticModulus, poissonRatio, compressiveBreakingStress, tensileBreakingStress;
    double frictionAngle, bondedDamping, tangentialFriction;
    double nonbondedNormalDamping, nonbondedTangentialDamping, criticalCrushingThickness;

    configs->get({"ConfigGroup:contacts","Config:elasticModulus"}, elasticModulus);
    configs->get({"ConfigGroup:contacts","Config:poissonRatio"}, poissonRatio);
    configs->get({"ConfigGroup:contacts","Config:compressiveBreakingStressType"}, compressiveBreakingStressType);
    configs->get({"ConfigGroup:contacts","Config:compressiveBreakingStress"}, compressiveBreakingStress);
    configs->get({"ConfigGroup:contacts","Config:tensileBreakingStressType"}, tensileBreakingStressType);
    configs->get({"ConfigGroup:contacts","Config:tensileBreakingStress"}, tensileBreakingStress);
    configs->get({"ConfigGroup:contacts","Config:frictionAngle"}, frictionAngle);
    configs->get({"ConfigGroup:contacts","Config:bondedDamping"}, bondedDamping);
    configs->get({"ConfigGroup:contacts","Config:tangentialFriction"}, tangentialFriction);
    configs->get({"ConfigGroup:contacts","Config:nonbondedNormalDamping"}, nonbondedNormalDamping);
    configs->get({"ConfigGroup:contacts","Config:nonbondedTangentialDamping"}, nonbondedTangentialDamping);
    configs->get({"ConfigGroup:contacts","Config:criticalCrushingThickness"}, criticalCrushingThickness);

    pairStyleCommand << "pair_style gran/hopkins ";
    pairStyleCommand << " " << elasticModulus;
    pairStyleCommand << " " << poissonRatio;
    pairStyleCommand << " " << compressiveBreakingStressType;
    pairStyleCommand << " " << compressiveBreakingStress;
    pairStyleCommand << " " << tensileBreakingStressType;
    pairStyleCommand << " " << tensileBreakingStress;
    pairStyleCommand << " " << frictionAngle;
    pairStyleCommand << " " << bondedDamping;
    pairStyleCommand << " " << tangentialFriction;
    pairStyleCommand << " " << nonbondedNormalDamping;
    pairStyleCommand << " " << nonbondedTangentialDamping;
    pairStyleCommand << " " << criticalCrushingThickness;

    return pairStyleCommand.str();

  } else if (get_contact_type() == DEMSI::CONTACTTYPE::HOOKE_THICKNESS) {

    // Hooke_thickness contact model
    double normalStiffness, tangentialStiffness;
    double dampingCoefficient, frictionCoefficient;

    configs->get({"ConfigGroup:contacts","Config:normalStiffness"}, normalStiffness);
    configs->get({"ConfigGroup:contacts","Config:tangentialStiffness"}, tangentialStiffness);
    configs->get({"ConfigGroup:contacts","Config:dampingCoefficient"}, dampingCoefficient);
    configs->get({"ConfigGroup:contacts","Config:frictionCoefficient"}, frictionCoefficient);

    pairStyleCommand << "pair_style gran/hooke/thickness ";
    pairStyleCommand << " " << normalStiffness;
    pairStyleCommand << " " << tangentialStiffness;
    pairStyleCommand << " " << dampingCoefficient;
    pairStyleCommand << " " << 0;
    pairStyleCommand << " " << frictionCoefficient;
    pairStyleCommand << " " << 1;

    return pairStyleCommand.str();

  } else if (get_contact_type() == DEMSI::CONTACTTYPE::HOOKE) {

    // Hooke contact model
    double normalStiffness, tangentialStiffness;
    double dampingCoefficient, frictionCoefficient;

    configs->get({"ConfigGroup:contacts","Config:normalStiffness"}, normalStiffness);
    configs->get({"ConfigGroup:contacts","Config:tangentialStiffness"}, tangentialStiffness);
    configs->get({"ConfigGroup:contacts","Config:dampingCoefficient"}, dampingCoefficient);
    configs->get({"ConfigGroup:contacts","Config:frictionCoefficient"}, frictionCoefficient);

    pairStyleCommand << "pair_style gran/hooke/thickness ";
    pairStyleCommand << " " << normalStiffness;
    pairStyleCommand << " " << tangentialStiffness;
    pairStyleCommand << " " << dampingCoefficient;
    pairStyleCommand << " " << 0;
    pairStyleCommand << " " << frictionCoefficient;
    pairStyleCommand << " " << 1;

    return pairStyleCommand.str();

  } else {

    log->abort("Contacts: pair_style_command: unsupported contact type.");

  }

}

// Return the contact part of the lammps wall fix command for the contact.
std::string Contacts::wall_command(void) {

  std::ostringstream wallCommand;

  if (get_contact_type() == DEMSI::CONTACTTYPE::HOOKE_THICKNESS or
      get_contact_type() == DEMSI::CONTACTTYPE::HOOKE) {

    double normalStiffness, tangentialStiffness;
    double dampingCoefficient, frictionCoefficient;

    configs->get({"ConfigGroup:contacts","Config:normalStiffness"}, normalStiffness);
    configs->get({"ConfigGroup:contacts","Config:tangentialStiffness"}, tangentialStiffness);
    configs->get({"ConfigGroup:contacts","Config:dampingCoefficient"}, dampingCoefficient);
    configs->get({"ConfigGroup:contacts","Config:frictionCoefficient"}, frictionCoefficient);

    wallCommand << "wall/gran hooke ";
    wallCommand << " " << normalStiffness;
    wallCommand << " " << tangentialStiffness;
    wallCommand << " " << dampingCoefficient;
    wallCommand << " " << 0;
    wallCommand << " " << frictionCoefficient;
    wallCommand << " " << 1;

    return wallCommand.str();

  } else {

    log->abort("Contacts: wall_command: unsupported contact type.");

  }

}

//Initialize bonded contacts for the Hopkins contact model
void Contacts::initialize_bonded_contacts_hopkins(void) {
  //Turns LAMMPS bonds into granular bonded interactions (long story)
  // Initial bond length is based on the criterion in the Hopkins contact model design doc
  // Other parameters are currently hard-coded

  if (lammpsInstance->lmp->atom->nbonds > 0) {

    //Loop over bonds
    // - get 2 bonded partners
    // - get corresponding history arrays
    // - set history array values to correct values
    LAMMPS_NS::NeighborKokkos *neighborKK;
    neighborKK = (LAMMPS_NS::NeighborKokkos *) lammpsInstance->lmp->neighbor;
    int nbondlist = neighborKK->nbondlist;
    int history_ndim = 12;
    class LAMMPS_NS::FixNeighHistory *fix_history;
    int fix_history_index;
    double h0;

    // sync data to device
    lammpsInstance->lmp->atomKK->k_x.sync<device_execution_space>();
    lammpsInstance->lmp->atomKK->k_radius.sync<device_execution_space>();

    // get device views
    auto x = lammpsInstance->lmp->atomKK->k_x.d_view;
    auto R = lammpsInstance->lmp->atomKK->k_radius.d_view;

    LAMMPS_NS::FixNeighHistoryKokkos<device_execution_space> *fix_historyKK;

    //Need to add check to get correct style name
    fix_history_index = lammpsInstance->lmp->modify->find_fix_by_style("NEIGH_HISTORY");
#ifdef KOKKOS_HAVE_CUDA
    fix_history_index = lammpsInstance->lmp->modify->find_fix_by_style("NEIGH_HISTORY/KK/DEVICE");
#else
    fix_history_index = lammpsInstance->lmp->modify->find_fix_by_style("NEIGH_HISTORY/KK/HOST");
#endif

    fix_history = (LAMMPS_NS::FixNeighHistory *) lammpsInstance->lmp->modify->fix[fix_history_index];
    fix_historyKK = (LAMMPS_NS::FixNeighHistoryKokkos<device_execution_space> *)fix_history;

    //see contact model design doc, eventually this will vary between bonds
    configs->get({"ConfigGroup:contacts","Config:bondThickness"}, h0);

    int newton = lammpsInstance->lmp->force->newton_pair;
    if (lammpsInstance->lmp->force->newton_bond != newton){
      log->abort("Newton bond setting must match newton pair");
      //Eventually remove this restriction
    }

    int inum = lammpsInstance->lmp->force->pair->list->inum;
    LAMMPS_NS::NeighListKokkos<device_execution_space>* k_list =
          static_cast<LAMMPS_NS::NeighListKokkos<device_execution_space>*>(lammpsInstance->lmp->force->pair->list);

    // sync data to device
    neighborKK->k_bondlist.sync<device_execution_space>();
    k_list->k_ilist.sync<device_execution_space>();

    // get device views
    auto bondlist = neighborKK->k_bondlist.d_view;
    auto d_numneigh = k_list->d_numneigh;
    auto d_neighbors = k_list->d_neighbors;
    auto d_ilist = k_list->d_ilist;
    auto d_firstflag = fix_historyKK->d_firstflag;
    auto d_firstvalue = fix_historyKK->d_firstvalue;

    int nlocal = lammpsInstance->lmp->atomKK->nlocal;
    Kokkos::parallel_for("Contacts::initialize_bonded_contacts_hopkins", Kokkos::RangePolicy<device_execution_space>
        (0, inum), KOKKOS_LAMBDA (const int ii) {
      double Rmin;
      double delta, dx, dy, rij, L, halfL, nx, ny;
      double tx, ty, cx, cy;
      double s1x, s1y, s2x, s2y;
      bool match;
      int ilocal, jlocal;
      ilocal = d_ilist(ii);
      int jnum = d_numneigh(ilocal);
      for (int jj=0; jj < jnum; jj++){
        jlocal = d_neighbors(ilocal,jj);
        jlocal &= NEIGHMASK;
        for (int k = 0; k < nbondlist; k++){
          match = ((ilocal == bondlist(k,0) && jlocal == bondlist(k,1)) || (ilocal == bondlist(k,1) && jlocal == bondlist(k,0)));
          if (match){
            //Need to check that bondlist ghost id's are handled properly
            //This probably means newton setting *must* be set to
            // the same for pair and bonds (not much of a restriction)
            //May need to handle full lists separately

            Rmin = (R(ilocal) < R(jlocal)) ? R(ilocal) : R(jlocal);
            dx = x(jlocal,0)-x(ilocal,0);
            dy = x(jlocal,1)-x(ilocal,1);
            rij = sqrt(dx*dx + dy*dy);
            delta = R(ilocal) + R(jlocal) - rij;
            nx = dx/rij;
            ny = dy/rij;

            //Goal is to set s1, s2 such that m cross k
            //  is in the direction of n for i-j, opposite for j-i
            //  This way, s1 for i-j is just -s1 for j-i
            tx = -ny;
            ty = nx;
            //if (lammpsInstance->lmp->atom->tag[ilocal] <
            //  lammpsInstance->lmp->atom->tag[jlocal]){
            //tx = -tx;
            //ty = -ty;
            // }

            L = 2*(Rmin + (Rmin - 0.5*delta)*fabs(R(ilocal)-R(jlocal))/rij);
            halfL = 0.5*L;

            cx = x(ilocal,0) + nx*R(ilocal) - 0.5*delta;
            cy = x(ilocal,1) + ny*R(ilocal) - 0.5*delta;
            s1x = cx + tx*halfL;
            s1y = cy + ty*halfL;
            s2x = cx - tx*halfL;
            s2y = cy - ty*halfL;

            d_firstflag(ilocal,jj) = 1;
            d_firstvalue(ilocal,history_ndim*jj) = s1x;
            d_firstvalue(ilocal,history_ndim*jj+1) = s1y;
            d_firstvalue(ilocal,history_ndim*jj+2) = s2x;
            d_firstvalue(ilocal,history_ndim*jj+3) = s2y;
            d_firstvalue(ilocal,history_ndim*jj+4) = s1x;
            d_firstvalue(ilocal,history_ndim*jj+5) = s1y;
            d_firstvalue(ilocal,history_ndim*jj+6) = s2x;
            d_firstvalue(ilocal,history_ndim*jj+7) = s2y;

            d_firstvalue(ilocal,history_ndim*jj+8) = 0;
            d_firstvalue(ilocal,history_ndim*jj+9) = 1;
            d_firstvalue(ilocal,history_ndim*jj+10) = L;
            d_firstvalue(ilocal,history_ndim*jj+11) = h0;

            break;
          }
        }
      }
    });

   //Once bonds are created, reset all forces and torques to 0

    // sync data to device
    lammpsInstance->lmp->atomKK->k_f.sync<device_execution_space>();
    lammpsInstance->lmp->atomKK->k_torque.sync<device_execution_space>();

    // get device views
    auto f = lammpsInstance->lmp->atomKK->k_f.d_view;
    auto tor = lammpsInstance->lmp->atomKK->k_torque.d_view;
    for (int i = 0; i < nlocal; i++){
       f(i,0) = 0.0;
       f(i,1) = 0.0;
       tor(i,0) = 0.0;
       tor(i,1) = 0.0;
    }

   // sync data to host
    lammpsInstance->lmp->atomKK->k_f.sync<host_execution_space>();
    lammpsInstance->lmp->atomKK->k_torque.sync<host_execution_space>();

    fix_history->pre_exchange();
  } // nbonds > 0

  lammpsInstance->lmp->input->one("delete_bonds all multi remove");
  Kokkos::fence();

}

// get the hopkins bond info from lammps
void Contacts::get_bonds_info_hopkins(Kokkos::View<LAMMPS_NS::tagint*[2]> &bondGlobalIDs, Kokkos::View<double*> &bondLength, Kokkos::View<double*> &bondThickness, Kokkos::View<double*[2]> &bondCrackFraction, Kokkos::View<double*[2]> &bondEndPoint1Particle1, Kokkos::View<double*[2]> &bondEndPoint2Particle1, Kokkos::View<double*[2]> &bondEndPoint1Particle2, Kokkos::View<double*[2]> &bondEndPoint2Particle2) {

  LAMMPS_NS::NeighborKokkos *neighborKK;
  neighborKK = (LAMMPS_NS::NeighborKokkos *) lammpsInstance->lmp->neighbor;
  int nbondlist = neighborKK->nbondlist;
  class LAMMPS_NS::FixNeighHistory *fix_history;
  int fix_history_index;
  LAMMPS_NS::FixNeighHistoryKokkos<device_execution_space> *fix_historyKK;

  //Need to add check to get correct style name
  fix_history_index = lammpsInstance->lmp->modify->find_fix_by_style("NEIGH_HISTORY");
#ifdef KOKKOS_HAVE_CUDA
  fix_history_index = lammpsInstance->lmp->modify->find_fix_by_style("NEIGH_HISTORY/KK/DEVICE");
#else
  fix_history_index = lammpsInstance->lmp->modify->find_fix_by_style("NEIGH_HISTORY/KK/HOST");
#endif

  fix_history = (LAMMPS_NS::FixNeighHistory *) lammpsInstance->lmp->modify->fix[fix_history_index];
  fix_historyKK = (LAMMPS_NS::FixNeighHistoryKokkos<device_execution_space> *)fix_history;

  int history_ndim = 12;  //Number of per-contact quantities stored/communicated; currently hard-corded for each contact model


  int inum = lammpsInstance->lmp->force->pair->list->inum;
  LAMMPS_NS::NeighListKokkos<device_execution_space>* k_list =
        static_cast<LAMMPS_NS::NeighListKokkos<device_execution_space>*>(lammpsInstance->lmp->force->pair->list);

  // sync data to device
  neighborKK->k_bondlist.sync<device_execution_space>();
  k_list->k_ilist.sync<device_execution_space>();
  lammpsInstance->lmp->atomKK->k_tag.sync<device_execution_space>();

  // get device views
  auto bondlist = neighborKK->k_bondlist.d_view;
  auto d_numneigh = k_list->d_numneigh;
  auto d_neighbors = k_list->d_neighbors;
  auto d_ilist = k_list->d_ilist;
  auto d_firstvalue = fix_historyKK->d_firstvalue;
  auto d_tag = lammpsInstance->lmp->atomKK->k_tag.d_view;

  int count = 0;
  int nlocal = lammpsInstance->lmp->atomKK->nlocal;
  auto bool_half_lists_used = half_lists_used();
  auto bool_newton_pair_used = newton_pair_used();
  Kokkos::parallel_reduce("Contacts::get_bonds_info_hopkins::count", Kokkos::RangePolicy<device_execution_space>
      (0, inum), KOKKOS_LAMBDA (const int ii, int& t_count) {
    int ilocal, jlocal;
    ilocal = d_ilist(ii);
    int jnum = d_numneigh(ilocal);
    for (int jj = 0; jj < jnum; jj++) {
      jlocal = d_neighbors(ilocal,jj);
      jlocal &= NEIGHMASK;

      //See contact model doc for definitions of various per-contact quantities:

      // For half lists:
      //   for newton on, all contacts are unique.
      //   for newton off, check if they're on different procs; if yes, only store the contact if globalID[i] < globalID[j].
      // For full lists:
      //   all contacts are duplicate, so take only globalID[i]<globalID[j].*/
      LAMMPS_NS::tagint globalID1 = d_tag(ilocal);
      LAMMPS_NS::tagint globalID2 = d_tag(jlocal);

      bool recordBond = false;
      if (bool_half_lists_used){
        if (bool_newton_pair_used) recordBond = true;
        else{
          if (on_different_procs(ilocal, jlocal, nlocal)){
            if (globalID1 < globalID2) recordBond = true;
          }
          else{
            recordBond = true;
          }
        }
      }
      else{
        if (globalID1 < globalID2) recordBond = true;
      }

      if (recordBond) t_count++;

    } // jj
  }, count); // ii

  Kokkos::resize(bondGlobalIDs, count, 2);
  Kokkos::resize(bondEndPoint1Particle1, count, 2);
  Kokkos::resize(bondEndPoint2Particle1, count, 2);
  Kokkos::resize(bondEndPoint1Particle2, count, 2);
  Kokkos::resize(bondEndPoint2Particle2, count, 2);
  Kokkos::resize(bondCrackFraction, count, 2);
  Kokkos::resize(bondLength, count);
  Kokkos::resize(bondThickness, count);

  Kokkos::View<int> insert_index("index being filled", 0);
  Kokkos::parallel_for("Contacts::get_bonds_info_hopkins::fill", Kokkos::RangePolicy<device_execution_space>
      (0, inum), KOKKOS_LAMBDA (const int ii) {
    int ilocal, jlocal;
    ilocal = d_ilist(ii);
    int jnum = d_numneigh(ilocal);
    for (int jj = 0; jj < jnum; jj++) {
      jlocal = d_neighbors(ilocal,jj);
      jlocal &= NEIGHMASK;

      //If needed you can get global IDs now with lmp->atom->tag[ilocal], lmp->atom->tag[jlocal]
      //See contact model doc for definitions of various per-contact quantities:

      // For half lists:
      //   for newton on, all contacts are unique.
      //   for newton off, check if they're on different procs; if yes, only store the contact if globalID[i] < globalID[j].
      // For full lists:
      //   all contacts are duplicate, so take only globalID[i]<globalID[j].*/
      LAMMPS_NS::tagint globalID1 = d_tag(ilocal);
      LAMMPS_NS::tagint globalID2 = d_tag(jlocal);

      bool recordBond = false;
      if (bool_half_lists_used){
        if (bool_newton_pair_used) recordBond = true;
        else{
          if (on_different_procs(ilocal, jlocal, nlocal)){
            if (globalID1 < globalID2) recordBond = true;
          }
          else{
            recordBond = true;
          }
        }
      }
      else{
        if (globalID1 < globalID2) recordBond = true;
      }

      if (recordBond){
        int index_to_fill = Kokkos::atomic_fetch_add(&insert_index(), 1);

        bondGlobalIDs(index_to_fill,0) = globalID1;
        bondGlobalIDs(index_to_fill,1) = globalID2;

        bondEndPoint1Particle1(index_to_fill,0) = d_firstvalue(ii,history_ndim*jj  );
        bondEndPoint1Particle1(index_to_fill,1) = d_firstvalue(ii,history_ndim*jj+1);
        bondEndPoint2Particle1(index_to_fill,0) = d_firstvalue(ii,history_ndim*jj+2);
        bondEndPoint2Particle1(index_to_fill,1) = d_firstvalue(ii,history_ndim*jj+3);
        bondEndPoint1Particle2(index_to_fill,0) = d_firstvalue(ii,history_ndim*jj+4);
        bondEndPoint1Particle2(index_to_fill,1) = d_firstvalue(ii,history_ndim*jj+5);
        bondEndPoint2Particle2(index_to_fill,0) = d_firstvalue(ii,history_ndim*jj+6);
        bondEndPoint2Particle2(index_to_fill,1) = d_firstvalue(ii,history_ndim*jj+7);

        bondCrackFraction(index_to_fill,0) = d_firstvalue(ii,history_ndim*jj+8);
        bondCrackFraction(index_to_fill,1) = d_firstvalue(ii,history_ndim*jj+9);

        bondLength(index_to_fill)    = d_firstvalue(ii,history_ndim*jj+10);
        bondThickness(index_to_fill) = d_firstvalue(ii,history_ndim*jj+11);

      } // avoid duplicates
    } // jj
  }); // ii
  Kokkos::fence();
}

// get the type of contacts currently being used
DEMSI::CONTACTTYPE Contacts::get_contact_type(void) {

  if (contactType == "hopkins") {
    return DEMSI::CONTACTTYPE::HOPKINS;
  } else if (contactType == "hooke") {
    return DEMSI::CONTACTTYPE::HOOKE;
  } else if (contactType == "hooke_thickness") {
    return DEMSI::CONTACTTYPE::HOOKE_THICKNESS;
  } else if (contactType == "none") {
    return DEMSI::CONTACTTYPE::NONE;
  } else {
    log->abort("Unknown contact type: ", contactType);
  }

}

// determine if bonds are being used
bool Contacts::bonds_active(void) {

  if (this->get_contact_type() == DEMSI::CONTACTTYPE::HOPKINS) {
    return true;
  } else {
    return false;
  }

}

// Return whether newton pair option is on
bool Contacts::newton_pair_used(void) {

  if (lammpsInstance->lmp->force->newton_pair) {
    return true;
  } else {
    return false;
  }

}

// Return whether half neighbour lists are being used.
bool Contacts::half_lists_used(void) {

  bool half = false;

  for (int i = 0; i < lammpsInstance->lmp->neighbor->old_nrequest; i++) {
    if (lammpsInstance->lmp->neighbor->old_requests[i]->pair) { //because requests can also come from fixes or computes
      if (lammpsInstance->lmp->neighbor->old_requests[i]->half == 1) {
        half = true;
      }
    }
  }

  return half;

}

// Return whether atoms are on different processor
KOKKOS_INLINE_FUNCTION
bool Contacts::on_different_procs(int ilocal, int jlocal, int nlocal){
  if (ilocal >= nlocal or jlocal >= nlocal) return true;
  return false;
}

} // namespace DEMSI
