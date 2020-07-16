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

#ifdef FIX_CLASS

FixStyle(addforce/kk,FixAddForceKokkos<LMPDeviceType>)
FixStyle(addforce/kk/device,FixAddForceKokkos<LMPDeviceType>)
FixStyle(addforce/kk/host,FixAddForceKokkos<LMPHostType>)

#else

#ifndef LMP_FIX_ADD_FORCE_KOKKOS_H
#define LMP_FIX_ADD_FORCE_KOKKOS_H

#include "fix_addforce.h"
#include "region.h"
#include "kokkos_type.h"
#include "kokkos_few.h"

namespace LAMMPS_NS {
class DomainKokkos;

struct s_double_4 {
  double values[4];
  KOKKOS_INLINE_FUNCTION
  s_double_4() {
    values[0] = 0;
    values[1] = 0;
    values[2] = 0;
    values[3] = 0;
  }
  KOKKOS_INLINE_FUNCTION
  s_double_4& operator+=(const s_double_4 &rhs){
    values[0] += rhs.values[0];
    values[1] += rhs.values[1];
    values[2] += rhs.values[2];
    values[3] += rhs.values[3];
    return *this;
  }

  KOKKOS_INLINE_FUNCTION
  void operator+=(const volatile s_double_4 &rhs) volatile {
    values[0] += rhs.values[0];
    values[1] += rhs.values[1];
    values[2] += rhs.values[2];
    values[3] += rhs.values[3];
  }
};
typedef s_double_4 double_4;

struct TagFixAddForceConstant{};

struct TagFixAddForceNonConstant{};

template<class DeviceType>
class FixAddForceKokkos : public FixAddForce {
 public:
  typedef DeviceType device_type;
  typedef double_4 value_type;
  typedef ArrayTypes<DeviceType> AT;

  FixAddForceKokkos(class LAMMPS *, int, char **);
  ~FixAddForceKokkos();
  void init();
  void post_force(int);

  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixAddForceConstant, const int&, double_4&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixAddForceNonConstant, const int&, double_4&) const;

 private:
  DAT::tdual_ffloat_2d k_sforce;
  typename AT::t_ffloat_2d_randomread d_sforce;
  typename AT::t_int_1d d_match;

  typename AT::t_x_array_randomread x;
  typename AT::t_f_array f;
  typename AT::t_int_1d_randomread mask;
  typename AT::t_imageint_1d_randomread image;

  DomainKokkos *domainKK;
  // member variables of domainKK that need to be captured
  Few<double,3> domainKK_prd;
  Few<double,6> domainKK_h;
  int domainKK_triclinic;

  Region* region;

};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Cannot (yet) use respa with Kokkos

Self-explanatory.

*/
