/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(gran/hopkins/kk,PairGranHopkinsKokkos<LMPDeviceType>)
PairStyle(gran/hopkins/kk/device,PairGranHopkinsKokkos<LMPDeviceType>)
PairStyle(gran/hopkins/kk/host,PairGranHopkinsKokkos<LMPHostType>)

#else

#ifndef LMP_PAIR_GRAN_HOPKINS_KOKKOS_H
#define LMP_PAIR_GRAN_HOPKINS_KOKKOS_H

#include "pair_gran_hopkins.h"
#include "pair_kokkos.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

template <class DeviceType>
class FixNeighHistoryKokkos;

template<int NEIGHFLAG, int NEWTON_PAIR, int HISTORYUPDATE, int EVFLAG>
struct TagPairGranHopkinsCompute {};

template <class DeviceType>
class PairGranHopkinsKokkos : public PairGranHopkins {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef EV_FLOAT value_type;

  PairGranHopkinsKokkos(class LAMMPS *);
  virtual ~PairGranHopkinsKokkos();
  virtual void compute(int, int);
  void init_style();

  template<int NEIGHFLAG, int NEWTON_PAIR, int HISTORYUPDATE, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairGranHopkinsCompute<NEIGHFLAG,NEWTON_PAIR,HISTORYUPDATE,EVFLAG>, const int, EV_FLOAT &ev) const;

  template<int NEIGHFLAG, int NEWTON_PAIR, int HISTORYUPDATE>
  KOKKOS_INLINE_FUNCTION
  void compute_nonbonded_kokkos(int i, int j, int jj, F_FLOAT &fx, F_FLOAT &fy,
			        F_FLOAT &torque_i, F_FLOAT &torque_j) const;

  template<int NEIGHFLAG, int NEWTON_PAIR, int HISTORYUPDATE>
  KOKKOS_INLINE_FUNCTION
  void compute_bonded_kokkos(int i, int j, int jj, F_FLOAT &fx, F_FLOAT &fy,
	                     F_FLOAT &torque_i, F_FLOAT &torque_j) const;

  KOKKOS_INLINE_FUNCTION
  void update_chi(F_FLOAT kn0, F_FLOAT kt0, F_FLOAT Dn, F_FLOAT Cn,
                  F_FLOAT Dt, F_FLOAT Ct, F_FLOAT hmin, F_FLOAT &chi1,
                  F_FLOAT &chi2) const;

  template<int NEWTON_PAIR>
  KOKKOS_INLINE_FUNCTION
  void ev_tally_xyz(EV_FLOAT &ev, int i, int j,
		    F_FLOAT fx, F_FLOAT fy, F_FLOAT fz,
		    X_FLOAT delx, X_FLOAT dely, X_FLOAT delz) const;

  template<int NEIGHFLAG, int NEWTON_PAIR>
  KOKKOS_INLINE_FUNCTION
  void ev_tally_xyz_atom(EV_FLOAT &ev, int i, int j,
			 F_FLOAT fx, F_FLOAT fy, F_FLOAT fz,
			 X_FLOAT delx, X_FLOAT dely, X_FLOAT delz) const;

  KOKKOS_INLINE_FUNCTION
  void elastic_stiffness(F_FLOAT meanIceThickness1,
			 F_FLOAT meanIceThickness2,
			 F_FLOAT radius1,
			 F_FLOAT radius2,
			 F_FLOAT mass1,
			 F_FLOAT mass2,
			 F_FLOAT bondLength,
			 F_FLOAT &elasticStiffness,
			 F_FLOAT &elasticDamping) const;

  KOKKOS_INLINE_FUNCTION
  void plastic_parameters(F_FLOAT particleRadius,
			  F_FLOAT plasticFrictionCoefficient,
			  F_FLOAT plasticHardeningCoefficient,
			  F_FLOAT ridgingIceThickness1,
			  F_FLOAT ridgingIceThickness2,
			  F_FLOAT ridgingIceThicknessWeight1,
			  F_FLOAT ridgingIceThicknessWeight2,
			  F_FLOAT radius1,
			  F_FLOAT radius2,
			  F_FLOAT &plasticFriction,
			  F_FLOAT &plasticHardeningStiffness) const;

  KOKKOS_INLINE_FUNCTION
  void elastic_plastic_model(F_FLOAT bondLength,
			     F_FLOAT previousForce,
			     F_FLOAT overlap,
			     F_FLOAT convergence,
			     F_FLOAT elasticStiffness,
			     F_FLOAT elasticDamping,
			     F_FLOAT plasticFriction,
			     F_FLOAT plasticHardeningStiffness,
			     F_FLOAT &ridgingForce,
			     F_FLOAT &elasticOverlap,
			     F_FLOAT &plasticOverlap,
			     F_FLOAT &elasticConvergence,
			     F_FLOAT &plasticConvergence) const;

  KOKKOS_INLINE_FUNCTION
  void elastic_model(F_FLOAT bondLength,
		     F_FLOAT elasticOverlap,
		     F_FLOAT elasticConvergence,
		     F_FLOAT elasticStiffness,
		     F_FLOAT elasticDamping,
		     F_FLOAT &elasticForce) const;

  KOKKOS_INLINE_FUNCTION
  void geometry_change(bool modifyOtherElement,
		       F_FLOAT bondLength,
		       F_FLOAT netToGrossClosingRatio1,
		       F_FLOAT netToGrossClosingRatio2,
		       F_FLOAT ridgeSlip,
		       F_FLOAT &ridgeSlipUsed,
		       F_FLOAT &changeEffectiveElementArea1,
		       F_FLOAT &changeEffectiveElementArea2) const;

  KOKKOS_INLINE_FUNCTION
  void hopkins_ridging_model(bool modifyOtherElement,
			     F_FLOAT overlap,
			     F_FLOAT convergence,
			     F_FLOAT meanIceThickness1,
			     F_FLOAT meanIceThickness2,
			     F_FLOAT radius1,
			     F_FLOAT radius2,
			     F_FLOAT mass1,
			     F_FLOAT mass2,
			     F_FLOAT ridgingIceThickness1,
			     F_FLOAT ridgingIceThickness2,
			     F_FLOAT ridgingIceThicknessWeight1,
			     F_FLOAT ridgingIceThicknessWeight2,
			     F_FLOAT netToGrossClosingRatio1,
			     F_FLOAT netToGrossClosingRatio2,
			     F_FLOAT &changeEffectiveElementArea1,
			     F_FLOAT &changeEffectiveElementArea2,
			     F_FLOAT particleRadius,
			     F_FLOAT plasticFrictionCoefficient,
			     F_FLOAT plasticHardeningCoefficient,
			     F_FLOAT bondLength,
			     F_FLOAT &ridgeSlip,
			     F_FLOAT &ridgeSlipUsed,
			     F_FLOAT &previousForce,
			     F_FLOAT &contactForce) const;


 protected:
  typename AT::t_x_array_randomread x;
  typename AT::t_v_array_randomread v;
  typename AT::t_v_array_randomread omega;
  typename AT::t_f_array f;
  typename AT::t_f_array torque;
  typename AT::t_int_1d_randomread type;
  typename AT::t_int_1d_randomread mask;
  typename AT::t_float_1d_randomread rmass;
  typename AT::t_float_1d_randomread radius;
  typename AT::t_float_1d_randomread mean_thickness;
  typename AT::t_float_1d_randomread min_thickness;
  typename AT::t_float_1d_randomread ridgingIceThickness;
  typename AT::t_float_1d_randomread ridgingIceThicknessWeight;
  typename AT::t_float_1d_randomread netToGrossClosingRatio;
  typename AT::t_float_1d changeEffectiveElementArea;

  DAT::tdual_efloat_1d k_eatom;
  DAT::tdual_virial_array k_vatom;
  typename AT::t_efloat_1d d_eatom;
  typename AT::t_virial_array d_vatom;
  typename AT::t_tagint_1d tag;

  typename AT::t_neighbors_2d d_neighbors;
  typename AT::t_int_1d_randomread d_ilist;
  typename AT::t_int_1d_randomread d_numneigh;

  typename Kokkos::View<int**> d_firsttouch;
  typename Kokkos::View<LMP_FLOAT**> d_firsthistory;

  int newton_pair;

  int neighflag;
  int nlocal,nall,eflag,vflag;

  // class storage to avoid strcmp in parallel code execution
  // or update->dt (ptr is not valid when update is copied to device)
  bool strcmp_sig_c0_type_constant;
  bool strcmp_sig_c0_type_KovacsSodhi;
  bool strcmp_sig_t0_type_constant;
  bool strcmp_sig_t0_type_multiply_sig_c0;
  double update_dt;

  FixNeighHistoryKokkos<DeviceType> *fix_historyKK;

  friend void pair_virial_fdotr_compute<PairGranHopkinsKokkos>(PairGranHopkinsKokkos*);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Pair granular requires atom attributes radius, rmass

The atom style defined does not have these attributes.

E: Pair granular requires ghost atoms store velocity

Use the comm_modify vel yes command to enable this.

E: Could not find pair fix neigh history ID

UNDOCUMENTED

U: Pair granular with shear history requires newton pair off

This is a current restriction of the implementation of pair
granular styles with history.

U: Could not find pair fix ID

A fix is created internally by the pair style to store shear
history information.  You cannot delete it.

*/
