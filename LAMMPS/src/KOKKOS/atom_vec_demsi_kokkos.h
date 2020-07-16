/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include <cmath>
#include <cstdlib>
#include <cstring>
#include "atom_vec_demsi_kokkos.h"
#include "atom_kokkos.h"
#include "atom_masks.h"
#include "comm_kokkos.h"
#include "domain.h"
#include "modify.h"
#include "force.h"
#include "fix.h"
#include "fix_adapt.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include "memory_kokkos.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define DELTA 10000

/* ---------------------------------------------------------------------- */

AtomVecDemsiKokkos::AtomVecDemsiKokkos(LAMMPS *lmp) : AtomVecKokkos(lmp)
{
  molecular = 1;
  kokkosable = 1;

  comm_x_only = 1;
  comm_f_only = 0;
  size_forward = 3;
  size_reverse = 6;
  size_border = 22;
  size_velocity = 6;
  size_data_atom = 7;
  size_data_vel = 7;
  xcol_data = 5;

  atom->sphere_flag = 1;
  atom->demsi_flag = 1;
  atom->radius_flag = atom->rmass_flag = atom->omega_flag =
    atom->torque_flag = 1;
  bonds_allow = 1;

  k_count = DAT::tdual_int_1d("atom::k_count",1);
  atomKK = (AtomKokkos *) atom;
  commKK = (CommKokkos *) comm;

}

/* ---------------------------------------------------------------------- */

void AtomVecDemsiKokkos::init()
{
  AtomVec::init();

  // set radvary if particle diameters are time-varying due to fix adapt

  radvary = 0;
  comm_x_only = 0;
  size_forward = 3;

  for (int i = 0; i < modify->nfix; i++) {
    if (strcmp(modify->fix[i]->style,"adapt") == 0) {
      FixAdapt *fix = (FixAdapt *) modify->fix[i];
      if (fix->diamflag) {
        radvary = 1;
        comm_x_only = 0;
        size_forward = 5;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   grow atom arrays
   n = 0 grows arrays by DELTA
   n > 0 allocates arrays to size n
------------------------------------------------------------------------- */

void AtomVecDemsiKokkos::grow(int n)
{
  if (n == 0) nmax += DELTA;
  else nmax = n;
  atomKK->nmax = nmax;
  if (nmax < 0 || nmax > MAXSMALLINT)
    error->one(FLERR,"Per-processor system is too big");

  sync(Device,ALL_MASK);
  modified(Device,ALL_MASK);

  memoryKK->grow_kokkos(atomKK->k_tag,atomKK->tag,nmax,"atom:tag");
  memoryKK->grow_kokkos(atomKK->k_type,atomKK->type,nmax,"atom:type");
  memoryKK->grow_kokkos(atomKK->k_mask,atomKK->mask,nmax,"atom:mask");
  memoryKK->grow_kokkos(atomKK->k_image,atomKK->image,nmax,"atom:image");

  memoryKK->grow_kokkos(atomKK->k_x,atomKK->x,nmax,3,"atom:x");
  memoryKK->grow_kokkos(atomKK->k_v,atomKK->v,nmax,3,"atom:v");
  memoryKK->grow_kokkos(atomKK->k_f,atomKK->f,nmax,3,"atom:f");
  memoryKK->grow_kokkos(atomKK->k_radius,atomKK->radius,nmax,"atom:radius");
  memoryKK->grow_kokkos(atomKK->k_rmass,atomKK->rmass,nmax,"atom:rmass");
  memoryKK->grow_kokkos(atomKK->k_omega,atomKK->omega,nmax,3,"atom:omega");
  memoryKK->grow_kokkos(atomKK->k_torque,atomKK->torque,nmax,3,"atom:torque");

  memoryKK->grow_kokkos(atomKK->k_forcing,atomKK->forcing,nmax,2,"atom:forcing");
  memoryKK->grow_kokkos(atomKK->k_mean_thickness,atomKK->mean_thickness,nmax,"atom:mean_thickness");
  memoryKK->grow_kokkos(atomKK->k_min_thickness,atomKK->min_thickness,nmax,"atom:min_thickness");
  memoryKK->grow_kokkos(atomKK->k_ridgingIceThickness,atomKK->ridgingIceThickness,nmax,"atom:ridgingIceThickness");
  memoryKK->grow_kokkos(atomKK->k_ridgingIceThicknessWeight,atomKK->ridgingIceThicknessWeight,nmax,"atom:ridgingIceThicknessWeight");
  memoryKK->grow_kokkos(atomKK->k_netToGrossClosingRatio,atomKK->netToGrossClosingRatio,nmax,"atom:netToGrossClosingRatio");
  memoryKK->grow_kokkos(atomKK->k_changeEffectiveElementArea,atomKK->changeEffectiveElementArea,nmax,"atom:changeEffectiveElementArea");
  memoryKK->grow_kokkos(atomKK->k_ice_area,atomKK->ice_area,nmax,"atom:ice_area");
  memoryKK->grow_kokkos(atomKK->k_coriolis,atomKK->coriolis,nmax,"atom:coriolis");
  memoryKK->grow_kokkos(atomKK->k_ocean_vel,atomKK->ocean_vel,nmax,2,"atom:ocean_vel");
  memoryKK->grow_kokkos(atomKK->k_bvector,atomKK->bvector,nmax,2,"atom:bvector");

  memoryKK->grow_kokkos(atomKK->k_nspecial,atomKK->nspecial,nmax,3,"atom:nspecial");
  memoryKK->grow_kokkos(atomKK->k_special,atomKK->special,nmax,atom->maxspecial,"atom:special");

  memoryKK->grow_kokkos(atomKK->k_num_bond,atomKK->num_bond,nmax,"atom:num_bond");
  memoryKK->grow_kokkos(atomKK->k_bond_type,atomKK->bond_type,nmax,atomKK->bond_per_atom,"atom:bond_type");
  memoryKK->grow_kokkos(atomKK->k_bond_atom,atomKK->bond_atom,nmax,atomKK->bond_per_atom,"atom:bond_atom");

  grow_reset();
  sync(Host,ALL_MASK);

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      modify->fix[atom->extra_grow[iextra]]->grow_arrays(nmax);
}

/* ----------------------------------------------------------------------
   reset local array ptrs
------------------------------------------------------------------------- */

void AtomVecDemsiKokkos::grow_reset()
{
  tag = atomKK->tag;
  d_tag = atomKK->k_tag.d_view;
  h_tag = atomKK->k_tag.h_view;

  type = atomKK->type;
  d_type = atomKK->k_type.d_view;
  h_type = atomKK->k_type.h_view;

  mask = atomKK->mask;
  d_mask = atomKK->k_mask.d_view;
  h_mask = atomKK->k_mask.h_view;

  image = atomKK->image;
  d_image = atomKK->k_image.d_view;
  h_image = atomKK->k_image.h_view;

  x = atomKK->x;
  d_x = atomKK->k_x.d_view;
  h_x = atomKK->k_x.h_view;
  v = atomKK->v;
  d_v = atomKK->k_v.d_view;
  h_v = atomKK->k_v.h_view;
  f = atomKK->f;
  d_f = atomKK->k_f.d_view;
  h_f = atomKK->k_f.h_view;

  radius = atomKK->radius;
  d_radius = atomKK->k_radius.d_view;
  h_radius = atomKK->k_radius.h_view;
  rmass = atomKK->rmass;
  d_rmass = atomKK->k_rmass.d_view;
  h_rmass = atomKK->k_rmass.h_view;

  omega = atomKK->omega;
  d_omega = atomKK->k_omega.d_view;
  h_omega = atomKK->k_omega.h_view;
  torque = atomKK->torque;
  d_torque = atomKK->k_torque.d_view;
  h_torque = atomKK->k_torque.h_view;

  forcing = atomKK->forcing;
  d_forcing = atomKK->k_forcing.d_view;
  h_forcing = atomKK->k_forcing.h_view;
  mean_thickness = atomKK->mean_thickness;
  d_mean_thickness = atomKK->k_mean_thickness.d_view;
  h_mean_thickness = atomKK->k_mean_thickness.h_view;
  min_thickness = atomKK->min_thickness;
  d_min_thickness = atomKK->k_min_thickness.d_view;
  h_min_thickness = atomKK->k_min_thickness.h_view;
  ridgingIceThickness = atomKK->ridgingIceThickness;
  d_ridgingIceThickness = atomKK->k_ridgingIceThickness.d_view;
  h_ridgingIceThickness = atomKK->k_ridgingIceThickness.h_view;
  ridgingIceThicknessWeight = atomKK->ridgingIceThicknessWeight;
  d_ridgingIceThicknessWeight = atomKK->k_ridgingIceThicknessWeight.d_view;
  h_ridgingIceThicknessWeight = atomKK->k_ridgingIceThicknessWeight.h_view;
  netToGrossClosingRatio = atomKK->netToGrossClosingRatio;
  d_netToGrossClosingRatio = atomKK->k_netToGrossClosingRatio.d_view;
  h_netToGrossClosingRatio = atomKK->k_netToGrossClosingRatio.h_view;
  changeEffectiveElementArea = atomKK->changeEffectiveElementArea;
  d_changeEffectiveElementArea = atomKK->k_changeEffectiveElementArea.d_view;
  h_changeEffectiveElementArea = atomKK->k_changeEffectiveElementArea.h_view;
  ice_area = atomKK->ice_area;
  d_ice_area = atomKK->k_ice_area.d_view;
  h_ice_area = atomKK->k_ice_area.h_view;
  coriolis = atomKK->coriolis;
  d_coriolis = atomKK->k_coriolis.d_view;
  h_coriolis = atomKK->k_coriolis.h_view;
  ocean_vel = atomKK->ocean_vel;
  d_ocean_vel = atomKK->k_ocean_vel.d_view;
  h_ocean_vel= atomKK->k_ocean_vel.h_view;
  bvector = atomKK->bvector;
  d_bvector = atomKK->k_bvector.d_view;
  h_bvector= atomKK->k_bvector.h_view;

  nspecial = atomKK->nspecial;
  d_nspecial = atomKK->k_nspecial.d_view;
  h_nspecial = atomKK->k_nspecial.h_view;
  special = atomKK->special;
  d_special = atomKK->k_special.d_view;
  h_special = atomKK->k_special.h_view;

  num_bond = atomKK->num_bond;
  d_num_bond = atomKK->k_num_bond.d_view;
  h_num_bond = atomKK->k_num_bond.h_view;
  bond_type = atomKK->bond_type;
  d_bond_type = atomKK->k_bond_type.d_view;
  h_bond_type = atomKK->k_bond_type.h_view;
  bond_atom = atomKK->bond_atom;
  d_bond_atom = atomKK->k_bond_atom.d_view;
  h_bond_atom = atomKK->k_bond_atom.h_view;

}

/* ----------------------------------------------------------------------
   copy atom I info to atom J
------------------------------------------------------------------------- */

void AtomVecDemsiKokkos::copy(int i, int j, int delflag)
{
/*
  sync(Host,X_MASK | V_MASK | TAG_MASK | TYPE_MASK |
            MASK_MASK | IMAGE_MASK | RADIUS_MASK |
            RMASS_MASK | OMEGA_MASK | THICKNESS_MASK | 
            FORCING_MASK | BOND_MASK | SPECIAL_MASK);
*/

  h_tag[j] = h_tag[i];
  h_type[j] = h_type[i];
  h_mask[j] = h_mask[i];
  h_image[j] = h_image[i];
  h_x(j,0) = h_x(i,0);
  h_x(j,1) = h_x(i,1);
  h_x(j,2) = h_x(i,2);
  h_v(j,0) = h_v(i,0);
  h_v(j,1) = h_v(i,1);
  h_v(j,2) = h_v(i,2);

  h_radius[j] = h_radius[i];
  h_rmass[j] = h_rmass[i];
  h_omega(j,0) = h_omega(i,0);
  h_omega(j,1) = h_omega(i,1);
  h_omega(j,2) = h_omega(i,2);

  h_forcing(j,0) = h_forcing(i,0);
  h_forcing(j,1) = h_forcing(i,1);

  h_mean_thickness(j) = h_mean_thickness(i);
  h_min_thickness(j) = h_min_thickness(i);
  h_ridgingIceThickness(j) = h_ridgingIceThickness(i);
  h_ridgingIceThicknessWeight(j) = h_ridgingIceThicknessWeight(i);
  h_netToGrossClosingRatio(j) = h_netToGrossClosingRatio(i);
  h_changeEffectiveElementArea(j) = h_changeEffectiveElementArea(i);

  h_ice_area(j) = h_ice_area(i);
  h_coriolis(j) = h_coriolis(i);

  h_ocean_vel(j,0) = h_ocean_vel(i,0);
  h_ocean_vel(j,1) = h_ocean_vel(i,1);
  h_bvector(j,0) = h_bvector(i,0);
  h_bvector(j,1) = h_bvector(i,1);

  h_nspecial(j,0) = h_nspecial(i,0);
  h_nspecial(j,1) = h_nspecial(i,1);
  h_nspecial(j,2) = h_nspecial(i,2);

  for (int k = 0; k < h_nspecial(j,2); k++)
     h_special(j,k) = h_special(i,k);

  h_num_bond(j) = h_num_bond(i);
  for (int k = 0; k < h_num_bond(j); k++){
     h_bond_type(j,k) = h_bond_type(i,k);
     h_bond_atom(j,k) = h_bond_atom(i,k);
  }

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      modify->fix[atom->extra_grow[iextra]]->copy_arrays(i,j,delflag);

/*
  modified(Host,X_MASK | V_MASK | TAG_MASK | TYPE_MASK |
                MASK_MASK | IMAGE_MASK | RADIUS_MASK |
                RMASS_MASK | OMEGA_MASK | THICKNESS_MASK | 
                FORCING_MASK | BOND_MASK | SPECIAL_MASK);
*/

}

/* ---------------------------------------------------------------------- */

template<class DeviceType,int PBC_FLAG,int TRICLINIC,int RADVARY>
struct AtomVecDemsiKokkos_PackComm {
  typedef DeviceType device_type;

  typename ArrayTypes<DeviceType>::t_x_array_randomread _x;
  typename ArrayTypes<DeviceType>::t_float_1d _radius,_rmass;
  typename ArrayTypes<DeviceType>::t_xfloat_2d_um _buf;
  typename ArrayTypes<DeviceType>::t_int_2d_const _list;
  const int _iswap;
  X_FLOAT _xprd,_yprd,_zprd,_xy,_xz,_yz;
  X_FLOAT _pbc[6];

  AtomVecDemsiKokkos_PackComm(
    const typename DAT::tdual_x_array &x,
    const typename DAT::tdual_float_1d &radius,
    const typename DAT::tdual_float_1d &rmass,
    const typename DAT::tdual_xfloat_2d &buf,
    const typename DAT::tdual_int_2d &list,
    const int & iswap,
    const X_FLOAT &xprd, const X_FLOAT &yprd, const X_FLOAT &zprd,
    const X_FLOAT &xy, const X_FLOAT &xz, const X_FLOAT &yz, const int* const pbc):
    _x(x.view<DeviceType>()),
    _radius(radius.view<DeviceType>()),
    _rmass(rmass.view<DeviceType>()),
    _list(list.view<DeviceType>()),_iswap(iswap),
    _xprd(xprd),_yprd(yprd),_zprd(zprd),
    _xy(xy),_xz(xz),_yz(yz) {
    const size_t elements = 3 + 2*RADVARY;
    const size_t maxsend = (buf.template view<DeviceType>().extent(0)*buf.template view<DeviceType>().extent(1))/elements;
    _buf = typename ArrayTypes<DeviceType>::t_xfloat_2d_um(buf.view<DeviceType>().data(),maxsend,elements);
    _pbc[0] = pbc[0]; _pbc[1] = pbc[1]; _pbc[2] = pbc[2];
    _pbc[3] = pbc[3]; _pbc[4] = pbc[4]; _pbc[5] = pbc[5];
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const int& i) const {
    const int j = _list(_iswap,i);
    if (PBC_FLAG == 0) {
      _buf(i,0) = _x(j,0);
      _buf(i,1) = _x(j,1);
      _buf(i,2) = _x(j,2);
    } else {
      if (TRICLINIC == 0) {
	_buf(i,0) = _x(j,0) + _pbc[0]*_xprd;
	_buf(i,1) = _x(j,1) + _pbc[1]*_yprd;
	_buf(i,2) = _x(j,2) + _pbc[2]*_zprd;
      } else {
	_buf(i,0) = _x(j,0) + _pbc[0]*_xprd + _pbc[5]*_xy + _pbc[4]*_xz;
	_buf(i,1) = _x(j,1) + _pbc[1]*_yprd + _pbc[3]*_yz;
	_buf(i,2) = _x(j,2) + _pbc[2]*_zprd;
      }
    }
    if (RADVARY) {
       _buf(i,3) = _radius(j);
       _buf(i,4) = _rmass(j);
    }
  }
};

/* ---------------------------------------------------------------------- */

int AtomVecDemsiKokkos::pack_comm_kokkos(
  const int &n,
  const DAT::tdual_int_2d &list,
  const int & iswap,
  const DAT::tdual_xfloat_2d &buf,
  const int &pbc_flag,
  const int* const pbc)
{
  
  // Check whether to always run forward communication on the host
  // Choose correct forward PackComm kernel
  if(commKK->forward_comm_on_host) {
    sync(Host,X_MASK|RADIUS_MASK|RMASS_MASK);
    if (radvary) {
      if(pbc_flag) {
        if(domain->triclinic) {
    	  struct AtomVecDemsiKokkos_PackComm<LMPHostType,1,1,1> f(
           atomKK->k_x,
    	   atomKK->k_radius,atomKK->k_rmass,
	   buf,list,iswap,
	   domain->xprd,domain->yprd,domain->zprd,
	   domain->xy,domain->xz,domain->yz,pbc);
	  Kokkos::parallel_for(n,f);
        } else {
  	  struct AtomVecDemsiKokkos_PackComm<LMPHostType,1,0,1> f(
          atomKK->k_x,
	  atomKK->k_radius,atomKK->k_rmass,
	  buf,list,iswap,
	  domain->xprd,domain->yprd,domain->zprd,
	  domain->xy,domain->xz,domain->yz,pbc);
	  Kokkos::parallel_for(n,f);
        }
      } else {
        if(domain->triclinic) {
 	  struct AtomVecDemsiKokkos_PackComm<LMPHostType,0,1,1> f(
            atomKK->k_x,
	    atomKK->k_radius,atomKK->k_rmass,
	    buf,list,iswap,
	    domain->xprd,domain->yprd,domain->zprd,
	    domain->xy,domain->xz,domain->yz,pbc);
	  Kokkos::parallel_for(n,f);
        } else {
	  struct AtomVecDemsiKokkos_PackComm<LMPHostType,0,0,1> f(
            atomKK->k_x,
	    atomKK->k_radius,atomKK->k_rmass,
	    buf,list,iswap,
	    domain->xprd,domain->yprd,domain->zprd,
	    domain->xy,domain->xz,domain->yz,pbc);
	  Kokkos::parallel_for(n,f);
        }
      }
    } else {
      if(pbc_flag) {
        if(domain->triclinic) {
    	  struct AtomVecDemsiKokkos_PackComm<LMPHostType,1,1,0> f(
           atomKK->k_x,
    	   atomKK->k_radius,atomKK->k_rmass,
	   buf,list,iswap,
	   domain->xprd,domain->yprd,domain->zprd,
	   domain->xy,domain->xz,domain->yz,pbc);
	  Kokkos::parallel_for(n,f);
        } else {
  	  struct AtomVecDemsiKokkos_PackComm<LMPHostType,1,0,0> f(
          atomKK->k_x,
	  atomKK->k_radius,atomKK->k_rmass,
	  buf,list,iswap,
	  domain->xprd,domain->yprd,domain->zprd,
	  domain->xy,domain->xz,domain->yz,pbc);
	  Kokkos::parallel_for(n,f);
        }
      } else {
        if(domain->triclinic) {
 	  struct AtomVecDemsiKokkos_PackComm<LMPHostType,0,1,0> f(
            atomKK->k_x,
	    atomKK->k_radius,atomKK->k_rmass,
	    buf,list,iswap,
	    domain->xprd,domain->yprd,domain->zprd,
	    domain->xy,domain->xz,domain->yz,pbc);
	  Kokkos::parallel_for(n,f);
        } else {
	  struct AtomVecDemsiKokkos_PackComm<LMPHostType,0,0,0> f(
            atomKK->k_x,
	    atomKK->k_radius,atomKK->k_rmass,
	    buf,list,iswap,
	    domain->xprd,domain->yprd,domain->zprd,
	    domain->xy,domain->xz,domain->yz,pbc);
	  Kokkos::parallel_for(n,f);
        }
      }
    }
  } else {
    sync(Device,X_MASK|RADIUS_MASK|RMASS_MASK);
    if(radvary) {
      if(pbc_flag) {
        if(domain->triclinic) {
	  struct AtomVecDemsiKokkos_PackComm<LMPDeviceType,1,1,1> f(
            atomKK->k_x,
	    atomKK->k_radius,atomKK->k_rmass,
	    buf,list,iswap,
	    domain->xprd,domain->yprd,domain->zprd,
	    domain->xy,domain->xz,domain->yz,pbc);
	  Kokkos::parallel_for(n,f);
        } else {
	  struct AtomVecDemsiKokkos_PackComm<LMPDeviceType,1,0,1> f(
            atomKK->k_x,
	    atomKK->k_radius,atomKK->k_rmass,
	    buf,list,iswap,
	    domain->xprd,domain->yprd,domain->zprd,
	    domain->xy,domain->xz,domain->yz,pbc);
	  Kokkos::parallel_for(n,f);
        }
      } else {
        if(domain->triclinic) {
	  struct AtomVecDemsiKokkos_PackComm<LMPDeviceType,0,1,1> f(
            atomKK->k_x,
	    atomKK->k_radius,atomKK->k_rmass,
	    buf,list,iswap,
	    domain->xprd,domain->yprd,domain->zprd,
	    domain->xy,domain->xz,domain->yz,pbc);
	  Kokkos::parallel_for(n,f);
        } else {
	  struct AtomVecDemsiKokkos_PackComm<LMPDeviceType,0,0,1> f(
            atomKK->k_x,
	    atomKK->k_radius,atomKK->k_rmass,
	    buf,list,iswap,
	    domain->xprd,domain->yprd,domain->zprd,
	    domain->xy,domain->xz,domain->yz,pbc);
	  Kokkos::parallel_for(n,f);
        }
      }
    } else {
      if(pbc_flag) {
        if(domain->triclinic) {
	  struct AtomVecDemsiKokkos_PackComm<LMPDeviceType,1,1,0> f(
            atomKK->k_x,
	    atomKK->k_radius,atomKK->k_rmass,
	    buf,list,iswap,
	    domain->xprd,domain->yprd,domain->zprd,
	    domain->xy,domain->xz,domain->yz,pbc);
	  Kokkos::parallel_for(n,f);
        } else {
	  struct AtomVecDemsiKokkos_PackComm<LMPDeviceType,1,0,0> f(
            atomKK->k_x,
	    atomKK->k_radius,atomKK->k_rmass,
	    buf,list,iswap,
	    domain->xprd,domain->yprd,domain->zprd,
	    domain->xy,domain->xz,domain->yz,pbc);
	  Kokkos::parallel_for(n,f);
        }
      } else {
        if(domain->triclinic) {
	  struct AtomVecDemsiKokkos_PackComm<LMPDeviceType,0,1,0> f(
            atomKK->k_x,
	    atomKK->k_radius,atomKK->k_rmass,
	    buf,list,iswap,
	    domain->xprd,domain->yprd,domain->zprd,
	    domain->xy,domain->xz,domain->yz,pbc);
	  Kokkos::parallel_for(n,f);
        } else {
	  struct AtomVecDemsiKokkos_PackComm<LMPDeviceType,0,0,0> f(
            atomKK->k_x,
	    atomKK->k_radius,atomKK->k_rmass,
	    buf,list,iswap,
	    domain->xprd,domain->yprd,domain->zprd,
	    domain->xy,domain->xz,domain->yz,pbc);
	  Kokkos::parallel_for(n,f);
        }
      }
    }
  }
  return n*size_forward;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType,int RADVARY,int PBC_FLAG,int TRICLINIC,int DEFORM_VREMAP>
struct AtomVecDemsiKokkos_PackCommVel {
  typedef DeviceType device_type;

  typename ArrayTypes<DeviceType>::t_x_array_randomread _x;
  typename ArrayTypes<DeviceType>::t_int_1d _mask;
  typename ArrayTypes<DeviceType>::t_float_1d _radius,_rmass;
  typename ArrayTypes<DeviceType>::t_v_array _v, _omega;
  typename ArrayTypes<DeviceType>::t_xfloat_2d_um _buf;
  typename ArrayTypes<DeviceType>::t_int_2d_const _list;
  const int _iswap;
  X_FLOAT _xprd,_yprd,_zprd,_xy,_xz,_yz;
  X_FLOAT _pbc[6];
  X_FLOAT _h_rate[6];
  const int _deform_vremap;

 AtomVecDemsiKokkos_PackCommVel(
    const typename DAT::tdual_x_array &x,
    const typename DAT::tdual_int_1d &mask,
    const typename DAT::tdual_float_1d &radius,
    const typename DAT::tdual_float_1d &rmass,
    const typename DAT::tdual_v_array &v,
    const typename DAT::tdual_v_array &omega,
    const typename DAT::tdual_xfloat_2d &buf,
    const typename DAT::tdual_int_2d &list,
    const int &iswap,
    const X_FLOAT &xprd, const X_FLOAT &yprd, const X_FLOAT &zprd,
    const X_FLOAT &xy, const X_FLOAT &xz, const X_FLOAT &yz, const int* const pbc,
    const double * const h_rate,
    const int &deform_vremap):
    _x(x.view<DeviceType>()),
    _mask(mask.view<DeviceType>()),
    _radius(radius.view<DeviceType>()),
    _rmass(rmass.view<DeviceType>()),
    _v(v.view<DeviceType>()),
    _omega(omega.view<DeviceType>()),
    _list(list.view<DeviceType>()),_iswap(iswap),
    _xprd(xprd),_yprd(yprd),_zprd(zprd),
    _xy(xy),_xz(xz),_yz(yz),
    _deform_vremap(deform_vremap)
  {
    const size_t elements = 9 + 2 * RADVARY;
    const int maxsend = (buf.template view<DeviceType>().extent(0)*buf.template view<DeviceType>().extent(1))/elements;
    _buf = typename ArrayTypes<DeviceType>::t_xfloat_2d_um(buf.view<DeviceType>().data(),maxsend,elements);
    _pbc[0] = pbc[0]; _pbc[1] = pbc[1]; _pbc[2] = pbc[2];
    _pbc[3] = pbc[3]; _pbc[4] = pbc[4]; _pbc[5] = pbc[5];
    _h_rate[0] = h_rate[0]; _h_rate[1] = h_rate[1]; _h_rate[2] = h_rate[2];
    _h_rate[3] = h_rate[3]; _h_rate[4] = h_rate[4]; _h_rate[5] = h_rate[5];
  }

 KOKKOS_INLINE_FUNCTION
  void operator() (const int& i) const {
    const int j = _list(_iswap,i);
    if (PBC_FLAG == 0) {
      _buf(i,0) = _x(j,0);
      _buf(i,1) = _x(j,1);
      _buf(i,2) = _x(j,2);
    } else {
      if (TRICLINIC == 0) {
        _buf(i,0) = _x(j,0) + _pbc[0]*_xprd;
        _buf(i,1) = _x(j,1) + _pbc[1]*_yprd;
        _buf(i,2) = _x(j,2) + _pbc[2]*_zprd;
      } else {
        _buf(i,0) = _x(j,0) + _pbc[0]*_xprd + _pbc[5]*_xy + _pbc[4]*_xz;
        _buf(i,1) = _x(j,1) + _pbc[1]*_yprd + _pbc[3]*_yz;
        _buf(i,2) = _x(j,2) + _pbc[2]*_zprd;
      }
    }
    if (DEFORM_VREMAP == 0) {
      _buf(i,3) = _v(j,0);
      _buf(i,4) = _v(j,1);
      _buf(i,5) = _v(j,2);
    } else {
      if (_mask(i) & _deform_vremap) {
        _buf(i,3) = _v(j,0) + _pbc[0]*_h_rate[0] + _pbc[5]*_h_rate[5] + _pbc[4]*_h_rate[4];
        _buf(i,4) = _v(j,1) + _pbc[1]*_h_rate[1] + _pbc[3]*_h_rate[3];
        _buf(i,5) = _v(j,2) + _pbc[2]*_h_rate[2];
      } else {
        _buf(i,3) = _v(j,0);
        _buf(i,4) = _v(j,1);
        _buf(i,5) = _v(j,2);
      }
    }
    _buf(i,6) = _omega(j,0);
    _buf(i,7) = _omega(j,1);
    _buf(i,8) = _omega(j,2);
    if (RADVARY) {
      _buf(i,9) = _radius(j);
      _buf(i,10) = _rmass(j);
    }
  }
};

/* ---------------------------------------------------------------------- */

int AtomVecDemsiKokkos::pack_comm_vel_kokkos(
  const int &n,
  const DAT::tdual_int_2d &list,
  const int & iswap,
  const DAT::tdual_xfloat_2d &buf,
  const int &pbc_flag,
  const int* const pbc)
{
 if(commKK->forward_comm_on_host) {
    sync(Host,X_MASK|RADIUS_MASK|RMASS_MASK|V_MASK|OMEGA_MASK);
    if(pbc_flag) {
      if(deform_vremap) {
        if(domain->triclinic) {
          if (radvary == 0) {
            struct AtomVecDemsiKokkos_PackCommVel<LMPHostType,0,1,1,1> f(
              atomKK->k_x,atomKK->k_mask,
              atomKK->k_radius,atomKK->k_rmass,
              atomKK->k_v,atomKK->k_omega,
              buf,list,iswap,
              domain->xprd,domain->yprd,domain->zprd,
              domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap);
            Kokkos::parallel_for(n,f);
          } else {
            struct AtomVecDemsiKokkos_PackCommVel<LMPHostType,1,1,1,1> f(
              atomKK->k_x,atomKK->k_mask,
              atomKK->k_radius,atomKK->k_rmass,
              atomKK->k_v,atomKK->k_omega,
              buf,list,iswap,
              domain->xprd,domain->yprd,domain->zprd,
              domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap);
            Kokkos::parallel_for(n,f);
          }
        } else {
          if (radvary == 0) {
            struct AtomVecDemsiKokkos_PackCommVel<LMPHostType,0,1,0,1> f(
              atomKK->k_x,atomKK->k_mask,
              atomKK->k_radius,atomKK->k_rmass,
              atomKK->k_v,atomKK->k_omega,
              buf,list,iswap,
              domain->xprd,domain->yprd,domain->zprd,
              domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap);
            Kokkos::parallel_for(n,f);
          } else {
            struct AtomVecDemsiKokkos_PackCommVel<LMPHostType,1,1,0,1> f(
              atomKK->k_x,atomKK->k_mask,
              atomKK->k_radius,atomKK->k_rmass,
              atomKK->k_v,atomKK->k_omega,
              buf,list,iswap,
              domain->xprd,domain->yprd,domain->zprd,
              domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap);
            Kokkos::parallel_for(n,f);
          }
        }
      } else {
       if(domain->triclinic) {
          if (radvary == 0) {
            struct AtomVecDemsiKokkos_PackCommVel<LMPHostType,0,1,1,0> f(
              atomKK->k_x,atomKK->k_mask,
              atomKK->k_radius,atomKK->k_rmass,
              atomKK->k_v,atomKK->k_omega,
              buf,list,iswap,
              domain->xprd,domain->yprd,domain->zprd,
              domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap);
            Kokkos::parallel_for(n,f);
          } else {
            struct AtomVecDemsiKokkos_PackCommVel<LMPHostType,1,1,1,0> f(
              atomKK->k_x,atomKK->k_mask,
              atomKK->k_radius,atomKK->k_rmass,
              atomKK->k_v,atomKK->k_omega,
              buf,list,iswap,
              domain->xprd,domain->yprd,domain->zprd,
              domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap);
            Kokkos::parallel_for(n,f);
          }
        } else {
          if (radvary == 0) {
            struct AtomVecDemsiKokkos_PackCommVel<LMPHostType,0,1,0,0> f(
              atomKK->k_x,atomKK->k_mask,
              atomKK->k_radius,atomKK->k_rmass,
              atomKK->k_v,atomKK->k_omega,
              buf,list,iswap,
              domain->xprd,domain->yprd,domain->zprd,
              domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap);
            Kokkos::parallel_for(n,f);
          } else {
            struct AtomVecDemsiKokkos_PackCommVel<LMPHostType,1,1,0,0> f(
              atomKK->k_x,atomKK->k_mask,
              atomKK->k_radius,atomKK->k_rmass,
              atomKK->k_v,atomKK->k_omega,
              buf,list,iswap,
              domain->xprd,domain->yprd,domain->zprd,
              domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap);
            Kokkos::parallel_for(n,f);
          }
        }
      }
    } else {
      if(domain->triclinic) {
        if (radvary == 0) {
          struct AtomVecDemsiKokkos_PackCommVel<LMPHostType,0,0,1,0> f(
            atomKK->k_x,atomKK->k_mask,
            atomKK->k_radius,atomKK->k_rmass,
            atomKK->k_v,atomKK->k_omega,
            buf,list,iswap,
            domain->xprd,domain->yprd,domain->zprd,
            domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap);
          Kokkos::parallel_for(n,f);
        } else {
          struct AtomVecDemsiKokkos_PackCommVel<LMPHostType,1,0,1,0> f(
            atomKK->k_x,atomKK->k_mask,
            atomKK->k_radius,atomKK->k_rmass,
            atomKK->k_v,atomKK->k_omega,
            buf,list,iswap,
            domain->xprd,domain->yprd,domain->zprd,
            domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap);
          Kokkos::parallel_for(n,f);
        }
      } else {
       if (radvary == 0) {
          struct AtomVecDemsiKokkos_PackCommVel<LMPHostType,0,0,0,0> f(
            atomKK->k_x,atomKK->k_mask,
            atomKK->k_radius,atomKK->k_rmass,
            atomKK->k_v,atomKK->k_omega,
            buf,list,iswap,
            domain->xprd,domain->yprd,domain->zprd,
            domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap);
          Kokkos::parallel_for(n,f);
        } else {
          struct AtomVecDemsiKokkos_PackCommVel<LMPHostType,1,0,0,0> f(
            atomKK->k_x,atomKK->k_mask,
            atomKK->k_radius,atomKK->k_rmass,
            atomKK->k_v,atomKK->k_omega,
            buf,list,iswap,
            domain->xprd,domain->yprd,domain->zprd,
            domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap);
          Kokkos::parallel_for(n,f);
        }
      }
    }
  } else {
    sync(Device,X_MASK|RADIUS_MASK|RMASS_MASK|V_MASK|OMEGA_MASK);
    if(pbc_flag) {
      if(deform_vremap) {
        if(domain->triclinic) {
          if (radvary == 0) {
            struct AtomVecDemsiKokkos_PackCommVel<LMPDeviceType,0,1,1,1> f(
              atomKK->k_x,atomKK->k_mask,
              atomKK->k_radius,atomKK->k_rmass,
              atomKK->k_v,atomKK->k_omega,
              buf,list,iswap,
              domain->xprd,domain->yprd,domain->zprd,
              domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap);
            Kokkos::parallel_for(n,f);
          } else {
            struct AtomVecDemsiKokkos_PackCommVel<LMPDeviceType,1,1,1,1> f(
              atomKK->k_x,atomKK->k_mask,
              atomKK->k_radius,atomKK->k_rmass,
              atomKK->k_v,atomKK->k_omega,
              buf,list,iswap,
              domain->xprd,domain->yprd,domain->zprd,
              domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap);
            Kokkos::parallel_for(n,f);
          }
       } else {
          if (radvary == 0) {
            struct AtomVecDemsiKokkos_PackCommVel<LMPDeviceType,0,1,0,1> f(
              atomKK->k_x,atomKK->k_mask,
              atomKK->k_radius,atomKK->k_rmass,
              atomKK->k_v,atomKK->k_omega,
              buf,list,iswap,
              domain->xprd,domain->yprd,domain->zprd,
              domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap);
            Kokkos::parallel_for(n,f);
          } else {
            struct AtomVecDemsiKokkos_PackCommVel<LMPDeviceType,1,1,0,1> f(
              atomKK->k_x,atomKK->k_mask,
              atomKK->k_radius,atomKK->k_rmass,
              atomKK->k_v,atomKK->k_omega,
              buf,list,iswap,
              domain->xprd,domain->yprd,domain->zprd,
              domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap);
            Kokkos::parallel_for(n,f);
          }
        }
      } else {
        if(domain->triclinic) {
          if (radvary == 0) {
            struct AtomVecDemsiKokkos_PackCommVel<LMPDeviceType,0,1,1,0> f(
              atomKK->k_x,atomKK->k_mask,
              atomKK->k_radius,atomKK->k_rmass,
              atomKK->k_v,atomKK->k_omega,
              buf,list,iswap,
              domain->xprd,domain->yprd,domain->zprd,
              domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap);
            Kokkos::parallel_for(n,f);
          } else {
            struct AtomVecDemsiKokkos_PackCommVel<LMPDeviceType,1,1,1,0> f(
              atomKK->k_x,atomKK->k_mask,
              atomKK->k_radius,atomKK->k_rmass,
              atomKK->k_v,atomKK->k_omega,
              buf,list,iswap,
              domain->xprd,domain->yprd,domain->zprd,
              domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap);
            Kokkos::parallel_for(n,f);
          }
        } else {
          if (radvary == 0) {
            struct AtomVecDemsiKokkos_PackCommVel<LMPDeviceType,0,1,0,0> f(
              atomKK->k_x,atomKK->k_mask,
              atomKK->k_radius,atomKK->k_rmass,
              atomKK->k_v,atomKK->k_omega,
              buf,list,iswap,
              domain->xprd,domain->yprd,domain->zprd,
              domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap);
            Kokkos::parallel_for(n,f);
          } else {
            struct AtomVecDemsiKokkos_PackCommVel<LMPDeviceType,1,1,0,0> f(
              atomKK->k_x,atomKK->k_mask,
              atomKK->k_radius,atomKK->k_rmass,
              atomKK->k_v,atomKK->k_omega,
              buf,list,iswap,
              domain->xprd,domain->yprd,domain->zprd,
              domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap);
            Kokkos::parallel_for(n,f);
          }
        }
      }
    } else {
      if(domain->triclinic) {
        if (radvary == 0) {
          struct AtomVecDemsiKokkos_PackCommVel<LMPDeviceType,0,0,1,0> f(
            atomKK->k_x,atomKK->k_mask,
            atomKK->k_radius,atomKK->k_rmass,
            atomKK->k_v,atomKK->k_omega,
            buf,list,iswap,
            domain->xprd,domain->yprd,domain->zprd,
            domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap);
          Kokkos::parallel_for(n,f);
        } else {
          struct AtomVecDemsiKokkos_PackCommVel<LMPDeviceType,1,0,1,0> f(
            atomKK->k_x,atomKK->k_mask,
            atomKK->k_radius,atomKK->k_rmass,
            atomKK->k_v,atomKK->k_omega,
            buf,list,iswap,
            domain->xprd,domain->yprd,domain->zprd,
            domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap);
          Kokkos::parallel_for(n,f);
        }
      } else {
       if (radvary == 0) {
          struct AtomVecDemsiKokkos_PackCommVel<LMPDeviceType,0,0,0,0> f(
            atomKK->k_x,atomKK->k_mask,
            atomKK->k_radius,atomKK->k_rmass,
            atomKK->k_v,atomKK->k_omega,
            buf,list,iswap,
            domain->xprd,domain->yprd,domain->zprd,
            domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap);
          Kokkos::parallel_for(n,f);
        } else {
          struct AtomVecDemsiKokkos_PackCommVel<LMPDeviceType,1,0,0,0> f(
            atomKK->k_x,atomKK->k_mask,
            atomKK->k_radius,atomKK->k_rmass,
            atomKK->k_v,atomKK->k_omega,
            buf,list,iswap,
            domain->xprd,domain->yprd,domain->zprd,
            domain->xy,domain->xz,domain->yz,pbc,h_rate,deform_vremap);
          Kokkos::parallel_for(n,f);
        }
      }
    }
  }
  return n*(size_forward+size_velocity);
}


/* ---------------------------------------------------------------------- */

template<class DeviceType,int PBC_FLAG,int TRICLINIC,int RADVARY>
struct AtomVecDemsiKokkos_PackCommSelf {
  typedef DeviceType device_type;

  typename ArrayTypes<DeviceType>::t_x_array_randomread _x;
  typename ArrayTypes<DeviceType>::t_x_array _xw;
  typename ArrayTypes<DeviceType>::t_float_1d _radius,_rmass;
  int _nfirst;
  typename ArrayTypes<DeviceType>::t_int_2d_const _list;
  const int _iswap;
  X_FLOAT _xprd,_yprd,_zprd,_xy,_xz,_yz;
  X_FLOAT _pbc[6];

  AtomVecDemsiKokkos_PackCommSelf(
    const typename DAT::tdual_x_array &x,
    const typename DAT::tdual_float_1d &radius,
    const typename DAT::tdual_float_1d &rmass,
    const int &nfirst,
    const typename DAT::tdual_int_2d &list,
    const int & iswap,
    const X_FLOAT &xprd, const X_FLOAT &yprd, const X_FLOAT &zprd,
    const X_FLOAT &xy, const X_FLOAT &xz, const X_FLOAT &yz, const int* const pbc):
    _x(x.view<DeviceType>()),_xw(x.view<DeviceType>()),
    _radius(radius.view<DeviceType>()),
    _rmass(rmass.view<DeviceType>()),
    _nfirst(nfirst),_list(list.view<DeviceType>()),_iswap(iswap),
    _xprd(xprd),_yprd(yprd),_zprd(zprd),
    _xy(xy),_xz(xz),_yz(yz) {
    _pbc[0] = pbc[0]; _pbc[1] = pbc[1]; _pbc[2] = pbc[2];
    _pbc[3] = pbc[3]; _pbc[4] = pbc[4]; _pbc[5] = pbc[5];
  };

  KOKKOS_INLINE_FUNCTION
  void operator() (const int& i) const {
    const int j = _list(_iswap,i);
    if (PBC_FLAG == 0) {
      _xw(i+_nfirst,0) = _x(j,0);
      _xw(i+_nfirst,1) = _x(j,1);
      _xw(i+_nfirst,2) = _x(j,2);
    } else {
      if (TRICLINIC == 0) {
	_xw(i+_nfirst,0) = _x(j,0) + _pbc[0]*_xprd;
	_xw(i+_nfirst,1) = _x(j,1) + _pbc[1]*_yprd;
	_xw(i+_nfirst,2) = _x(j,2) + _pbc[2]*_zprd;
      } else {
	_xw(i+_nfirst,0) = _x(j,0) + _pbc[0]*_xprd + _pbc[5]*_xy + _pbc[4]*_xz;
	_xw(i+_nfirst,1) = _x(j,1) + _pbc[1]*_yprd + _pbc[3]*_yz;
	_xw(i+_nfirst,2) = _x(j,2) + _pbc[2]*_zprd;
      }
    }
    if (RADVARY) {
      _radius(i+_nfirst) = _radius(j);
      _rmass(i+_nfirst) = _rmass(j);
    }
  }
};

/* ---------------------------------------------------------------------- */

int AtomVecDemsiKokkos::pack_comm_self(
  const int &n, const DAT::tdual_int_2d &list, const int & iswap,
  const int nfirst, const int &pbc_flag, const int* const pbc) {

  if(commKK->forward_comm_on_host) {
    sync(Host,X_MASK|RADIUS_MASK|RMASS_MASK);
    modified(Host,X_MASK|RADIUS_MASK|RMASS_MASK);
    if(radvary) {
      if(pbc_flag) {
        if(domain->triclinic) {
	  struct AtomVecDemsiKokkos_PackCommSelf<LMPHostType,1,1,1> f(
            atomKK->k_x,
	    atomKK->k_radius,atomKK->k_rmass,
	    nfirst,list,iswap,
	    domain->xprd,domain->yprd,domain->zprd,
	    domain->xy,domain->xz,domain->yz,pbc);
	  Kokkos::parallel_for(n,f);
        } else {
	  struct AtomVecDemsiKokkos_PackCommSelf<LMPHostType,1,0,1> f(
            atomKK->k_x,
	    atomKK->k_radius,atomKK->k_rmass,
	    nfirst,list,iswap,
	    domain->xprd,domain->yprd,domain->zprd,
	    domain->xy,domain->xz,domain->yz,pbc);
	  Kokkos::parallel_for(n,f);
        }
      } else {
        if(domain->triclinic) {
	  struct AtomVecDemsiKokkos_PackCommSelf<LMPHostType,0,1,1> f(
            atomKK->k_x,
	    atomKK->k_radius,atomKK->k_rmass,
	    nfirst,list,iswap,
	    domain->xprd,domain->yprd,domain->zprd,
	    domain->xy,domain->xz,domain->yz,pbc);
	  Kokkos::parallel_for(n,f);
        } else {
	  struct AtomVecDemsiKokkos_PackCommSelf<LMPHostType,0,0,1> f(
            atomKK->k_x,
	    atomKK->k_radius,atomKK->k_rmass,
	    nfirst,list,iswap,
	    domain->xprd,domain->yprd,domain->zprd,
	    domain->xy,domain->xz,domain->yz,pbc);
	  Kokkos::parallel_for(n,f);
        }
      }
    } else {
      if(pbc_flag) {
        if(domain->triclinic) {
	  struct AtomVecDemsiKokkos_PackCommSelf<LMPHostType,1,1,0> f(
            atomKK->k_x,
	    atomKK->k_radius,atomKK->k_rmass,
	    nfirst,list,iswap,
	    domain->xprd,domain->yprd,domain->zprd,
	    domain->xy,domain->xz,domain->yz,pbc);
	  Kokkos::parallel_for(n,f);
        } else {
	  struct AtomVecDemsiKokkos_PackCommSelf<LMPHostType,1,0,0> f(
            atomKK->k_x,
	    atomKK->k_radius,atomKK->k_rmass,
	    nfirst,list,iswap,
	    domain->xprd,domain->yprd,domain->zprd,
	    domain->xy,domain->xz,domain->yz,pbc);
	  Kokkos::parallel_for(n,f);
        }
      } else {
        if(domain->triclinic) {
	  struct AtomVecDemsiKokkos_PackCommSelf<LMPHostType,0,1,0> f(
            atomKK->k_x,
	    atomKK->k_radius,atomKK->k_rmass,
	    nfirst,list,iswap,
	    domain->xprd,domain->yprd,domain->zprd,
	    domain->xy,domain->xz,domain->yz,pbc);
	  Kokkos::parallel_for(n,f);
        } else {
	  struct AtomVecDemsiKokkos_PackCommSelf<LMPHostType,0,0,0> f(
            atomKK->k_x,
	    atomKK->k_radius,atomKK->k_rmass,
	    nfirst,list,iswap,
	    domain->xprd,domain->yprd,domain->zprd,
	    domain->xy,domain->xz,domain->yz,pbc);
	  Kokkos::parallel_for(n,f);
        }
      }
    }
  } else {
    sync(Device,X_MASK|RADIUS_MASK|RMASS_MASK|THICKNESS_MASK);
    modified(Device,X_MASK|RADIUS_MASK|RMASS_MASK|THICKNESS_MASK);
    if(radvary) {
      if(pbc_flag) {
        if(domain->triclinic) {
	  struct AtomVecDemsiKokkos_PackCommSelf<LMPDeviceType,1,1,1> f(
            atomKK->k_x,
	    atomKK->k_radius,atomKK->k_rmass,
	    nfirst,list,iswap,
	    domain->xprd,domain->yprd,domain->zprd,
	    domain->xy,domain->xz,domain->yz,pbc);
	  Kokkos::parallel_for(n,f);
        } else {
	  struct AtomVecDemsiKokkos_PackCommSelf<LMPDeviceType,1,0,1> f(
            atomKK->k_x,
	    atomKK->k_radius,atomKK->k_rmass,
	    nfirst,list,iswap,
	    domain->xprd,domain->yprd,domain->zprd,
	    domain->xy,domain->xz,domain->yz,pbc);
	  Kokkos::parallel_for(n,f);
        }
      } else {
        if(domain->triclinic) {
	  struct AtomVecDemsiKokkos_PackCommSelf<LMPDeviceType,0,1,1> f(
            atomKK->k_x,
	    atomKK->k_radius,atomKK->k_rmass,
	    nfirst,list,iswap,
	    domain->xprd,domain->yprd,domain->zprd,
	    domain->xy,domain->xz,domain->yz,pbc);
	  Kokkos::parallel_for(n,f);
        } else {
	  struct AtomVecDemsiKokkos_PackCommSelf<LMPDeviceType,0,0,1> f(
            atomKK->k_x,
	    atomKK->k_radius,atomKK->k_rmass,
	    nfirst,list,iswap,
	    domain->xprd,domain->yprd,domain->zprd,
	    domain->xy,domain->xz,domain->yz,pbc);
	  Kokkos::parallel_for(n,f);
        }
      }
    } else {
      if(pbc_flag) {
        if(domain->triclinic) {
	  struct AtomVecDemsiKokkos_PackCommSelf<LMPDeviceType,1,1,0> f(
            atomKK->k_x,
	    atomKK->k_radius,atomKK->k_rmass,
	    nfirst,list,iswap,
	    domain->xprd,domain->yprd,domain->zprd,
	    domain->xy,domain->xz,domain->yz,pbc);
	  Kokkos::parallel_for(n,f);
        } else {
	  struct AtomVecDemsiKokkos_PackCommSelf<LMPDeviceType,1,0,0> f(
            atomKK->k_x,
	    atomKK->k_radius,atomKK->k_rmass,
	    nfirst,list,iswap,
	    domain->xprd,domain->yprd,domain->zprd,
	    domain->xy,domain->xz,domain->yz,pbc);
	  Kokkos::parallel_for(n,f);
        }
      } else {
        if(domain->triclinic) {
	  struct AtomVecDemsiKokkos_PackCommSelf<LMPDeviceType,0,1,0> f(
            atomKK->k_x,
	    atomKK->k_radius,atomKK->k_rmass,
	    nfirst,list,iswap,
	    domain->xprd,domain->yprd,domain->zprd,
	    domain->xy,domain->xz,domain->yz,pbc);
	  Kokkos::parallel_for(n,f);
        } else {
	  struct AtomVecDemsiKokkos_PackCommSelf<LMPDeviceType,0,0,0> f(
            atomKK->k_x,
	    atomKK->k_radius,atomKK->k_rmass,
	    nfirst,list,iswap,
	    domain->xprd,domain->yprd,domain->zprd,
	    domain->xy,domain->xz,domain->yz,pbc);
	  Kokkos::parallel_for(n,f);
        }
      }
    }
  }
  return n*size_forward;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType,int RADVARY>
struct AtomVecDemsiKokkos_UnpackComm {
  typedef DeviceType device_type;

  typename ArrayTypes<DeviceType>::t_x_array _x;
  typename ArrayTypes<DeviceType>::t_float_1d _radius,_rmass;
  typename ArrayTypes<DeviceType>::t_xfloat_2d_const _buf;
  int _first;

  AtomVecDemsiKokkos_UnpackComm(
    const typename DAT::tdual_x_array &x,
    const typename DAT::tdual_float_1d &radius,
    const typename DAT::tdual_float_1d &rmass,
    const typename DAT::tdual_xfloat_2d &buf,
    const int& first):
    _x(x.view<DeviceType>()),
    _radius(radius.view<DeviceType>()),
    _rmass(rmass.view<DeviceType>()),
    _buf(buf.view<DeviceType>()),
    _first(first) {
    const size_t elements = 3 + 2 * RADVARY;
    const int maxsend = (buf.template view<DeviceType>().extent(0)*buf.template view<DeviceType>().extent(1))/elements;
    buffer_view<DeviceType>(_buf,buf,maxsend,elements);
   };

  KOKKOS_INLINE_FUNCTION
  void operator() (const int& i) const {
    _x(i+_first,0) = _buf(i,0);
    _x(i+_first,1) = _buf(i,1);
    _x(i+_first,2) = _buf(i,2);
    if (RADVARY){
      _radius(i+_first) = _buf(i,3);
      _rmass(i+_first) = _buf(i,4);
    }
  }
};

/* ---------------------------------------------------------------------- */

void AtomVecDemsiKokkos::unpack_comm_kokkos(
  const int &n, const int &first,
  const DAT::tdual_xfloat_2d &buf ) {
  if(commKK->forward_comm_on_host) {
    if (radvary) {
      modified(Host,X_MASK|RADIUS_MASK|RMASS_MASK);
      struct AtomVecDemsiKokkos_UnpackComm<LMPHostType,1> f(
        atomKK->k_x,
        atomKK->k_radius,atomKK->k_rmass,
        buf,first);
      Kokkos::parallel_for(n,f);
    } else {
      modified(Host,X_MASK|RADIUS_MASK|RMASS_MASK);
      struct AtomVecDemsiKokkos_UnpackComm<LMPHostType,0> f(
        atomKK->k_x,
        atomKK->k_radius,atomKK->k_rmass,
        buf,first);
      Kokkos::parallel_for(n,f);
    }
  } else {
    if (radvary) {
      modified(Device,X_MASK|RADIUS_MASK|RMASS_MASK);
      struct AtomVecDemsiKokkos_UnpackComm<LMPDeviceType,1> f(
        atomKK->k_x,
        atomKK->k_radius,atomKK->k_rmass,
        buf,first);
      Kokkos::parallel_for(n,f);
    } else {
      modified(Device,X_MASK|RADIUS_MASK|RMASS_MASK);
      struct AtomVecDemsiKokkos_UnpackComm<LMPDeviceType,0> f(
        atomKK->k_x,
        atomKK->k_radius,atomKK->k_rmass,
        buf,first);
      Kokkos::parallel_for(n,f);
    }
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType,int RADVARY>
struct AtomVecDemsiKokkos_UnpackCommVel {
  typedef DeviceType device_type;

  typename ArrayTypes<DeviceType>::t_x_array _x;
  typename ArrayTypes<DeviceType>::t_float_1d _radius,_rmass;
  typename ArrayTypes<DeviceType>::t_v_array _v, _omega;
  typename ArrayTypes<DeviceType>::t_xfloat_2d_const _buf;
  int _first;

  AtomVecDemsiKokkos_UnpackCommVel(
    const typename DAT::tdual_x_array &x,
    const typename DAT::tdual_float_1d &radius,
    const typename DAT::tdual_float_1d &rmass,
    const typename DAT::tdual_v_array &v,
    const typename DAT::tdual_v_array &omega,
    const typename DAT::tdual_xfloat_2d &buf,
    const int& first):
    _x(x.view<DeviceType>()),
    _radius(radius.view<DeviceType>()),
    _rmass(rmass.view<DeviceType>()),
    _v(v.view<DeviceType>()),
    _omega(omega.view<DeviceType>()),
    _first(first)
  {
    const size_t elements = 9 + 2 * RADVARY;
    const int maxsend = (buf.template view<DeviceType>().extent(0)*buf.template view<DeviceType>().extent(1))/elements;
    buffer_view<DeviceType>(_buf,buf,maxsend,elements);
  };

 KOKKOS_INLINE_FUNCTION
  void operator() (const int& i) const {
    _x(i+_first,0) = _buf(i,0);
    _x(i+_first,1) = _buf(i,1);
    _x(i+_first,2) = _buf(i,2);
    _v(i+_first,0) = _buf(i,3);
    _v(i+_first,1) = _buf(i,4);
    _v(i+_first,2) = _buf(i,5);
    _omega(i+_first,0) = _buf(i,6);
    _omega(i+_first,1) = _buf(i,7);
    _omega(i+_first,2) = _buf(i,8);
    if (RADVARY) {
      _radius(i+_first) = _buf(i,9);
      _rmass(i+_first) = _buf(i,10);
    }
  }
};

/* ---------------------------------------------------------------------- */

void AtomVecDemsiKokkos::unpack_comm_vel_kokkos(
  const int &n, const int &first,
  const DAT::tdual_xfloat_2d &buf ) {
  if(commKK->forward_comm_on_host) {
    modified(Host,X_MASK|RADIUS_MASK|RMASS_MASK|V_MASK|OMEGA_MASK);
    if (radvary == 0) {
      struct AtomVecDemsiKokkos_UnpackCommVel<LMPHostType,0> f(
        atomKK->k_x,
        atomKK->k_radius,atomKK->k_rmass,
        atomKK->k_v,atomKK->k_omega,
        buf,first);
      Kokkos::parallel_for(n,f);
    } else {
      struct AtomVecDemsiKokkos_UnpackCommVel<LMPHostType,1> f(
        atomKK->k_x,
        atomKK->k_radius,atomKK->k_rmass,
        atomKK->k_v,atomKK->k_omega,
        buf,first);
      Kokkos::parallel_for(n,f);
    }
  } else {
    modified(Device,X_MASK|RADIUS_MASK|RMASS_MASK|V_MASK|OMEGA_MASK);
    if (radvary == 0) {
      struct AtomVecDemsiKokkos_UnpackCommVel<LMPDeviceType,0> f(
        atomKK->k_x,
        atomKK->k_radius,atomKK->k_rmass,
        atomKK->k_v,atomKK->k_omega,
        buf,first);
      Kokkos::parallel_for(n,f);
    } else {
      struct AtomVecDemsiKokkos_UnpackCommVel<LMPDeviceType,1> f(
        atomKK->k_x,
        atomKK->k_radius,atomKK->k_rmass,
        atomKK->k_v,atomKK->k_omega,
        buf,first);
      Kokkos::parallel_for(n,f);
    }
  }
}


/* ---------------------------------------------------------------------- */

int AtomVecDemsiKokkos::pack_comm(int n, int *list, double *buf,
			           int pbc_flag, int *pbc)
{
  int i,j,m;
  double dx,dy,dz;

  if (radvary == 0) {
    // Not sure if we need to call sync for X here
    sync(Host,X_MASK);
    m = 0;
    if (pbc_flag == 0) {
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = h_x(j,0);
        buf[m++] = h_x(j,1);
        buf[m++] = h_x(j,2);
      }
    } else {
      if (domain->triclinic == 0) {
        dx = pbc[0]*domain->xprd;
        dy = pbc[1]*domain->yprd;
        dz = pbc[2]*domain->zprd;
      } else {
        dx = pbc[0]*domain->xprd + pbc[5]*domain->xy + pbc[4]*domain->xz;
        dy = pbc[1]*domain->yprd + pbc[3]*domain->yz;
        dz = pbc[2]*domain->zprd;
      }
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = h_x(j,0) + dx;
        buf[m++] = h_x(j,1) + dy;
        buf[m++] = h_x(j,2) + dz;
      }
    }
  } else {
    sync(Host,X_MASK|RADIUS_MASK|RMASS_MASK);
    m = 0;
    if (pbc_flag == 0) {
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = h_x(j,0);
        buf[m++] = h_x(j,1);
        buf[m++] = h_x(j,2);
        buf[m++] = h_radius[j];
        buf[m++] = h_rmass[j];
      }
    } else {
      if (domain->triclinic == 0) {
        dx = pbc[0]*domain->xprd;
        dy = pbc[1]*domain->yprd;
        dz = pbc[2]*domain->zprd;
      } else {
        dx = pbc[0]*domain->xprd + pbc[5]*domain->xy + pbc[4]*domain->xz;
        dy = pbc[1]*domain->yprd + pbc[3]*domain->yz;
        dz = pbc[2]*domain->zprd;
      }
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = h_x(j,0) + dx;
        buf[m++] = h_x(j,1) + dy;
        buf[m++] = h_x(j,2) + dz;
        buf[m++] = h_radius[j];
        buf[m++] = h_rmass[j];
      }
    }
  }

  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecDemsiKokkos::pack_comm_vel(int n, int *list, double *buf,
				       int pbc_flag, int *pbc)
{
  int i,j,m;
  double dx,dy,dz,dvx,dvy,dvz;

  if (radvary == 0) {
    sync(Host,X_MASK|V_MASK|OMEGA_MASK);
    m = 0;
    if (pbc_flag == 0) {
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = h_x(j,0);
        buf[m++] = h_x(j,1);
        buf[m++] = h_x(j,2);
        buf[m++] = h_v(j,0);
        buf[m++] = h_v(j,1);
        buf[m++] = h_v(j,2);
        buf[m++] = h_omega(j,0);
        buf[m++] = h_omega(j,1);
        buf[m++] = h_omega(j,2);
      }
    } else {
      if (domain->triclinic == 0) {
        dx = pbc[0]*domain->xprd;
        dy = pbc[1]*domain->yprd;
        dz = pbc[2]*domain->zprd;
      } else {
        dx = pbc[0]*domain->xprd + pbc[5]*domain->xy + pbc[4]*domain->xz;
        dy = pbc[1]*domain->yprd + pbc[3]*domain->yz;
        dz = pbc[2]*domain->zprd;
      }
      if (!deform_vremap) {
        for (i = 0; i < n; i++) {
          j = list[i];
          buf[m++] = h_x(j,0) + dx;
          buf[m++] = h_x(j,1) + dy;
          buf[m++] = h_x(j,2) + dz;
          buf[m++] = h_v(j,0);
          buf[m++] = h_v(j,1);
          buf[m++] = h_v(j,2);
          buf[m++] = h_omega(j,0);
          buf[m++] = h_omega(j,1);
          buf[m++] = h_omega(j,2);
        }
      } else {
        dvx = pbc[0]*h_rate[0] + pbc[5]*h_rate[5] + pbc[4]*h_rate[4];
        dvy = pbc[1]*h_rate[1] + pbc[3]*h_rate[3];
        dvz = pbc[2]*h_rate[2];
        for (i = 0; i < n; i++) {
          j = list[i];
          buf[m++] = h_x(j,0) + dx;
          buf[m++] = h_x(j,1) + dy;
          buf[m++] = h_x(j,2) + dz;
         if (mask[i] & deform_groupbit) {
            buf[m++] = h_v(j,0) + dvx;
            buf[m++] = h_v(j,1) + dvy;
            buf[m++] = h_v(j,2) + dvz;
          } else {
            buf[m++] = h_v(j,0);
            buf[m++] = h_v(j,1);
            buf[m++] = h_v(j,2);
          }
          buf[m++] = h_omega(j,0);
          buf[m++] = h_omega(j,1);
          buf[m++] = h_omega(j,2);
        }
      }
    }

  } else {
    sync(Host,X_MASK|RADIUS_MASK|RMASS_MASK|V_MASK|OMEGA_MASK);
    m = 0;
    if (pbc_flag == 0) {
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = h_x(j,0);
        buf[m++] = h_x(j,1);
        buf[m++] = h_x(j,2);
        buf[m++] = h_radius[j];
        buf[m++] = h_rmass[j];
        buf[m++] = h_v(j,0);
        buf[m++] = h_v(j,1);
        buf[m++] = h_v(j,2);
        buf[m++] = h_omega(j,0);
        buf[m++] = h_omega(j,1);
        buf[m++] = h_omega(j,2);
      }
    } else {
      if (domain->triclinic == 0) {
        dx = pbc[0]*domain->xprd;
        dy = pbc[1]*domain->yprd;
        dz = pbc[2]*domain->zprd;
      } else {
        dx = pbc[0]*domain->xprd + pbc[5]*domain->xy + pbc[4]*domain->xz;
        dy = pbc[1]*domain->yprd + pbc[3]*domain->yz;
        dz = pbc[2]*domain->zprd;
      }
      if (!deform_vremap) {
        for (i = 0; i < n; i++) {
          j = list[i];
          buf[m++] = h_x(j,0) + dx;
          buf[m++] = h_x(j,1) + dy;
          buf[m++] = h_x(j,2) + dz;
          buf[m++] = h_radius[j];
          buf[m++] = h_rmass[j];
          buf[m++] = h_v(j,0);
          buf[m++] = h_v(j,1);
          buf[m++] = h_v(j,2);
          buf[m++] = h_omega(j,0);
          buf[m++] = h_omega(j,1);
          buf[m++] = h_omega(j,2);
        }
      } else {
        dvx = pbc[0]*h_rate[0] + pbc[5]*h_rate[5] + pbc[4]*h_rate[4];
        dvy = pbc[1]*h_rate[1] + pbc[3]*h_rate[3];
        dvz = pbc[2]*h_rate[2];
        for (i = 0; i < n; i++) {
          j = list[i];
          buf[m++] = h_x(j,0) + dx;
          buf[m++] = h_x(j,1) + dy;
          buf[m++] = h_x(j,2) + dz;
          buf[m++] = h_radius[j];
          buf[m++] = h_rmass[j];
          if (mask[i] & deform_groupbit) {
            buf[m++] = h_v(j,0) + dvx;
            buf[m++] = h_v(j,1) + dvy;
            buf[m++] = h_v(j,2) + dvz;
          } else {
            buf[m++] = h_v(j,0);
            buf[m++] = h_v(j,1);
            buf[m++] = h_v(j,2);
          }
          buf[m++] = h_omega(j,0);
          buf[m++] = h_omega(j,1);
          buf[m++] = h_omega(j,2);
        }
      }
    }
  }

  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecDemsiKokkos::pack_comm_hybrid(int n, int *list, double *buf)
{
  int m = 0;
  if (radvary == 0) {
 
     return 0;
     
  } else {
    sync(Host,RADIUS_MASK|RMASS_MASK);
    for (int i = 0; i < n; i++) {
      const int j = list[i];
      buf[m++] = h_radius[j];
      buf[m++] = h_rmass[j];
    }
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecDemsiKokkos::unpack_comm(int n, int first, double *buf)
{
  if (radvary == 0) {
    int m = 0;
    const int last = first + n;
    for (int i = first; i < last; i++) {
      h_x(i,0) = buf[m++];
      h_x(i,1) = buf[m++];
      h_x(i,2) = buf[m++];
    }
    modified(Host,X_MASK);
  } else {
    int m = 0;
    const int last = first + n;
    for (int i = first; i < last; i++) {
      h_x(i,0) = buf[m++];
      h_x(i,1) = buf[m++];
      h_x(i,2) = buf[m++];
      h_radius(i) = buf[m++];
      h_rmass(i) = buf[m++];
    }
    modified(Host,X_MASK|RADIUS_MASK|RMASS_MASK);
  }
}

/* ---------------------------------------------------------------------- */

void AtomVecDemsiKokkos::unpack_comm_vel(int n, int first, double *buf)
{
  if (radvary == 0) {
    int m = 0;
    const int last = first + n;
    for (int i = first; i < last; i++) {
      h_x(i,0) = buf[m++];
      h_x(i,1) = buf[m++];
      h_x(i,2) = buf[m++];
      h_v(i,0) = buf[m++];
      h_v(i,1) = buf[m++];
      h_v(i,2) = buf[m++];
      h_omega(i,0) = buf[m++];
      h_omega(i,1) = buf[m++];
      h_omega(i,2) = buf[m++];
    }
    modified(Host,X_MASK|V_MASK|OMEGA_MASK);
  } else {
    int m = 0;
    const int last = first + n;
    for (int i = first; i < last; i++) {
      h_x(i,0) = buf[m++];
      h_x(i,1) = buf[m++];
      h_x(i,2) = buf[m++];
      h_radius(i) = buf[m++];
      h_rmass(i) = buf[m++];
      h_v(i,0) = buf[m++];
      h_v(i,1) = buf[m++];
      h_v(i,2) = buf[m++];
      h_omega(i,0) = buf[m++];
      h_omega(i,1) = buf[m++];
      h_omega(i,2) = buf[m++];
    }
    modified(Host,X_MASK|RADIUS_MASK|RMASS_MASK|V_MASK|OMEGA_MASK);
  }
}

/* ---------------------------------------------------------------------- */

int AtomVecDemsiKokkos::unpack_comm_hybrid(int n, int first, double *buf)
{
  if (radvary == 0) return 0;

  int m = 0;
  const int last = first + n;

  if (radvary == 0) {
     return 0;
  }
  else {
     for (int i = first; i < last; i++) {
        h_radius(i) = buf[m++];
        h_rmass(i) = buf[m++];
     }
     modified(Host,RADIUS_MASK|RMASS_MASK);
  }
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecDemsiKokkos::pack_reverse(int n, int first, double *buf)
{
  if(n > 0)
    sync(Host,F_MASK|TORQUE_MASK);

  int m = 0;
  const int last = first + n;
  for (int i = first; i < last; i++) {
    buf[m++] = h_f(i,0);
    buf[m++] = h_f(i,1);
    buf[m++] = h_f(i,2);
    buf[m++] = h_torque(i,0);
    buf[m++] = h_torque(i,1);
    buf[m++] = h_torque(i,2);
  }
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecDemsiKokkos::pack_reverse_hybrid(int n, int first, double *buf)
{
  if(n > 0)
    sync(Host,TORQUE_MASK);

  int m = 0;
  const int last = first + n;
  for (int i = first; i < last; i++) {
    buf[m++] = h_torque(i,0);
    buf[m++] = h_torque(i,1);
    buf[m++] = h_torque(i,2);
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecDemsiKokkos::unpack_reverse(int n, int *list, double *buf)
{
  if(n > 0) {
    modified(Host,F_MASK|TORQUE_MASK);
  }

  int m = 0;
  for (int i = 0; i < n; i++) {
    const int j = list[i];
    h_f(j,0) += buf[m++];
    h_f(j,1) += buf[m++];
    h_f(j,2) += buf[m++];
    h_torque(j,0) += buf[m++];
    h_torque(j,1) += buf[m++];
    h_torque(j,2) += buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

int AtomVecDemsiKokkos::unpack_reverse_hybrid(int n, int *list, double *buf)
{
  if(n > 0) {
    modified(Host,TORQUE_MASK);
  }

  int m = 0;
  for (int i = 0; i < n; i++) {
    const int j = list[i];
    h_torque(j,0) += buf[m++];
    h_torque(j,1) += buf[m++];
    h_torque(j,2) += buf[m++];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType,int PBC_FLAG>
struct AtomVecDemsiKokkos_PackBorder {
  typedef DeviceType device_type;

  typename ArrayTypes<DeviceType>::t_xfloat_2d _buf;
  const typename ArrayTypes<DeviceType>::t_int_2d_const _list;
  const int _iswap;
  const typename ArrayTypes<DeviceType>::t_x_array_randomread _x;
  const typename ArrayTypes<DeviceType>::t_tagint_1d _tag;
  const typename ArrayTypes<DeviceType>::t_int_1d _type;
  const typename ArrayTypes<DeviceType>::t_int_1d _mask;
  const typename ArrayTypes<DeviceType>::t_float_2d _forcing;
  typename ArrayTypes<DeviceType>::t_float_1d _radius,_rmass;
  typename ArrayTypes<DeviceType>::t_float_1d _mean_thickness,_min_thickness;
  typename ArrayTypes<DeviceType>::t_float_1d _ridgingIceThickness;
  typename ArrayTypes<DeviceType>::t_float_1d _ridgingIceThicknessWeight;
  typename ArrayTypes<DeviceType>::t_float_1d _netToGrossClosingRatio;
  typename ArrayTypes<DeviceType>::t_float_1d _changeEffectiveElementArea;
  typename ArrayTypes<DeviceType>::t_float_1d _ice_area,_coriolis;
  typename ArrayTypes<DeviceType>::t_float_2d _ocean_vel,_bvector;
  X_FLOAT _dx,_dy,_dz;

  AtomVecDemsiKokkos_PackBorder(const typename ArrayTypes<DeviceType>::t_xfloat_2d &buf,
			      const typename ArrayTypes<DeviceType>::t_int_2d_const &list,
			      const int & iswap,
			      const typename ArrayTypes<DeviceType>::t_x_array &x,
			      const typename ArrayTypes<DeviceType>::t_tagint_1d &tag,
			      const typename ArrayTypes<DeviceType>::t_int_1d &type,
			      const typename ArrayTypes<DeviceType>::t_int_1d &mask,
			      const typename ArrayTypes<DeviceType>::t_float_1d &radius,
			      const typename ArrayTypes<DeviceType>::t_float_1d &rmass,
			      const typename ArrayTypes<DeviceType>::t_float_2d &forcing,
			      const typename ArrayTypes<DeviceType>::t_float_1d &mean_thickness,
			      const typename ArrayTypes<DeviceType>::t_float_1d &min_thickness,
			      const typename ArrayTypes<DeviceType>::t_float_1d &ridgingIceThickness,
			      const typename ArrayTypes<DeviceType>::t_float_1d &ridgingIceThicknessWeight,
			      const typename ArrayTypes<DeviceType>::t_float_1d &netToGrossClosingRatio,
			      const typename ArrayTypes<DeviceType>::t_float_1d &changeEffectiveElementArea,
			      const typename ArrayTypes<DeviceType>::t_float_1d &ice_area,
			      const typename ArrayTypes<DeviceType>::t_float_1d &coriolis,
			      const typename ArrayTypes<DeviceType>::t_float_2d &ocean_vel,
			      const typename ArrayTypes<DeviceType>::t_float_2d &bvector,
			      const X_FLOAT &dx, const X_FLOAT &dy, const X_FLOAT &dz):
    _buf(buf),_list(list),_iswap(iswap),
    _x(x),_tag(tag),_type(type),_mask(mask),
    _radius(radius),
    _rmass(rmass),
    _forcing(forcing),
    _mean_thickness(mean_thickness), _min_thickness(min_thickness),
    _ridgingIceThickness(ridgingIceThickness),
    _ridgingIceThicknessWeight(ridgingIceThicknessWeight),
    _netToGrossClosingRatio(netToGrossClosingRatio),
    _changeEffectiveElementArea(changeEffectiveElementArea),
    _ice_area(ice_area), _coriolis(coriolis),
    _ocean_vel(ocean_vel), _bvector(bvector),
    _dx(dx),_dy(dy),_dz(dz) {}
  
  KOKKOS_INLINE_FUNCTION
  void operator() (const int& i) const {
    const int j = _list(_iswap,i);
    if (PBC_FLAG == 0) {
      _buf(i,0) = _x(j,0);
      _buf(i,1) = _x(j,1);
      _buf(i,2) = _x(j,2);
    } else {
      _buf(i,0) = _x(j,0) + _dx;
      _buf(i,1) = _x(j,1) + _dy;
      _buf(i,2) = _x(j,2) + _dz;
    }
    _buf(i,3)  = d_ubuf(_tag(j)).d;
    _buf(i,4)  = d_ubuf(_type(j)).d;
    _buf(i,5)  = d_ubuf(_mask(j)).d;
    _buf(i,6)  = _radius(j);
    _buf(i,7)  = _rmass(j);
    _buf(i,8)  = _forcing(j,0);
    _buf(i,9)  = _forcing(j,1);
    _buf(i,10) = _mean_thickness(j);
    _buf(i,11) = _min_thickness(j);
    _buf(i,12) = _ridgingIceThickness(j);
    _buf(i,13) = _ridgingIceThicknessWeight(j);
    _buf(i,14) = _netToGrossClosingRatio(j);
    _buf(i,15) = _changeEffectiveElementArea(j);
    _buf(i,16) = _ice_area(j);
    _buf(i,17) = _coriolis(j);
    _buf(i,18) = _ocean_vel(j,0);
    _buf(i,19) = _ocean_vel(j,1);
    _buf(i,20) = _bvector(j,0);
    _buf(i,21) = _bvector(j,1);
  }
};

/* ---------------------------------------------------------------------- */

int AtomVecDemsiKokkos::pack_border_kokkos(int n, DAT::tdual_int_2d k_sendlist, DAT::tdual_xfloat_2d buf,int iswap,
					    int pbc_flag, int *pbc, ExecutionSpace space)
{
  X_FLOAT dx,dy,dz;

  //sync(space,ALL_MASK);

  if (pbc_flag != 0) {
    if (domain->triclinic == 0) {
      dx = pbc[0]*domain->xprd;
      dy = pbc[1]*domain->yprd;
      dz = pbc[2]*domain->zprd;
    } else {
      dx = pbc[0];
      dy = pbc[1];
      dz = pbc[2];
    }
    if(space==Host) {
      AtomVecDemsiKokkos_PackBorder<LMPHostType,1> f(
        buf.view<LMPHostType>(), k_sendlist.view<LMPHostType>(),
        iswap,h_x,h_tag,h_type,h_mask,
        h_radius,h_rmass,h_forcing,h_mean_thickness,h_min_thickness,
	h_ridgingIceThickness,h_ridgingIceThicknessWeight,h_netToGrossClosingRatio,h_changeEffectiveElementArea,
        h_ice_area,h_coriolis,h_ocean_vel,h_bvector,dx,dy,dz);
      Kokkos::parallel_for(n,f);
    } else {
      AtomVecDemsiKokkos_PackBorder<LMPDeviceType,1> f(
        buf.view<LMPDeviceType>(), k_sendlist.view<LMPDeviceType>(),
        iswap,d_x,d_tag,d_type,d_mask,
        d_radius,d_rmass,d_forcing,d_mean_thickness,d_min_thickness,
	d_ridgingIceThickness,d_ridgingIceThicknessWeight,d_netToGrossClosingRatio,d_changeEffectiveElementArea,
        d_ice_area,d_coriolis,d_ocean_vel,d_bvector,dx,dy,dz);
      Kokkos::parallel_for(n,f);
    }
  } else {
    dx = dy = dz = 0;
    if(space==Host) {
      AtomVecDemsiKokkos_PackBorder<LMPHostType,0> f(
        buf.view<LMPHostType>(), k_sendlist.view<LMPHostType>(),
        iswap,h_x,h_tag,h_type,h_mask,
        h_radius,h_rmass,h_forcing,h_mean_thickness,h_min_thickness,
	h_ridgingIceThickness,h_ridgingIceThicknessWeight,h_netToGrossClosingRatio,h_changeEffectiveElementArea,
        h_ice_area,h_coriolis,h_ocean_vel,h_bvector,dx,dy,dz);
      Kokkos::parallel_for(n,f);
    } else {
      AtomVecDemsiKokkos_PackBorder<LMPDeviceType,0> f(
        buf.view<LMPDeviceType>(), k_sendlist.view<LMPDeviceType>(),
        iswap,d_x,d_tag,d_type,d_mask,
        d_radius,d_rmass,d_forcing,d_mean_thickness,d_min_thickness,
	d_ridgingIceThickness,d_ridgingIceThicknessWeight,d_netToGrossClosingRatio,d_changeEffectiveElementArea,
        d_ice_area,d_coriolis,d_ocean_vel,d_bvector,dx,dy,dz);
      Kokkos::parallel_for(n,f);
    }
  }
  return n*size_border;
}

/* ---------------------------------------------------------------------- */

int AtomVecDemsiKokkos::pack_border(int n, int *list, double *buf,
                               int pbc_flag, int *pbc)
{
  int i,j,m;
  double dx,dy,dz;


   // need this?
  //sync(Host,ALL_MASK);

  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = h_x(j,0);
      buf[m++] = h_x(j,1);
      buf[m++] = h_x(j,2);
      buf[m++] = ubuf(h_tag(j)).d;
      buf[m++] = ubuf(h_type(j)).d;
      buf[m++] = ubuf(h_mask(j)).d;
      buf[m++] = h_radius(j);
      buf[m++] = h_rmass(j);
      buf[m++] = h_forcing(j,0);
      buf[m++] = h_forcing(j,1);
      buf[m++] = h_mean_thickness(j);
      buf[m++] = h_min_thickness(j);
      buf[m++] = h_ridgingIceThickness(j);
      buf[m++] = h_ridgingIceThicknessWeight(j);
      buf[m++] = h_netToGrossClosingRatio(j);
      buf[m++] = h_changeEffectiveElementArea(j);
      buf[m++] = h_ice_area(j);
      buf[m++] = h_coriolis(j);
      buf[m++] = h_ocean_vel(j,0);
      buf[m++] = h_ocean_vel(j,1);
      buf[m++] = h_bvector(j,0);
      buf[m++] = h_bvector(j,1);
    }
  } else {
    if (domain->triclinic == 0) {
      dx = pbc[0]*domain->xprd;
      dy = pbc[1]*domain->yprd;
      dz = pbc[2]*domain->zprd;
    } else {
      dx = pbc[0];
      dy = pbc[1];
      dz = pbc[2];
    }
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = h_x(j,0) + dx;
      buf[m++] = h_x(j,1) + dy;
      buf[m++] = h_x(j,2) + dz;
      buf[m++] = ubuf(h_tag(j)).d;
      buf[m++] = ubuf(h_type(j)).d;
      buf[m++] = ubuf(h_mask(j)).d;
      buf[m++] = h_radius(j);
      buf[m++] = h_rmass(j);
      buf[m++] = h_forcing(j,0);
      buf[m++] = h_forcing(j,1);
      buf[m++] = h_mean_thickness(j);
      buf[m++] = h_min_thickness(j);
      buf[m++] = h_ridgingIceThickness(j);
      buf[m++] = h_ridgingIceThicknessWeight(j);
      buf[m++] = h_netToGrossClosingRatio(j);
      buf[m++] = h_changeEffectiveElementArea(j);
      buf[m++] = h_ice_area(j);
      buf[m++] = h_coriolis(j);
      buf[m++] = h_ocean_vel(j,0);
      buf[m++] = h_ocean_vel(j,1);
      buf[m++] = h_bvector(j,0);
      buf[m++] = h_bvector(j,1);
    }
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->pack_border(n,list,&buf[m]);

  return m;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType,int PBC_FLAG,int DEFORM_VREMAP>
struct AtomVecDemsiKokkos_PackBorderVel {
  typedef DeviceType device_type;

  typename ArrayTypes<DeviceType>::t_xfloat_2d_um _buf;
  const typename ArrayTypes<DeviceType>::t_int_2d_const _list;
  const int _iswap;
  const typename ArrayTypes<DeviceType>::t_x_array_randomread _x;
  const typename ArrayTypes<DeviceType>::t_tagint_1d _tag;
  const typename ArrayTypes<DeviceType>::t_int_1d _type;
  const typename ArrayTypes<DeviceType>::t_int_1d _mask;
  typename ArrayTypes<DeviceType>::t_float_1d _radius,_rmass;
  typename ArrayTypes<DeviceType>::t_float_1d _mean_thickness,_min_thickness;
  typename ArrayTypes<DeviceType>::t_float_1d _ridgingIceThickness;
  typename ArrayTypes<DeviceType>::t_float_1d _ridgingIceThicknessWeight;
  typename ArrayTypes<DeviceType>::t_float_1d _netToGrossClosingRatio;
  typename ArrayTypes<DeviceType>::t_float_1d _changeEffectiveElementArea;
  typename ArrayTypes<DeviceType>::t_float_2d _forcing;
  typename ArrayTypes<DeviceType>::t_float_1d _ice_area,_coriolis;
  typename ArrayTypes<DeviceType>::t_float_2d _ocean_vel,_bvector;
  typename ArrayTypes<DeviceType>::t_v_array _v, _omega;
  X_FLOAT _dx,_dy,_dz, _dvx, _dvy, _dvz;
  const int _deform_groupbit;

  AtomVecDemsiKokkos_PackBorderVel(
    const typename ArrayTypes<DeviceType>::t_xfloat_2d &buf,
    const typename ArrayTypes<DeviceType>::t_int_2d_const &list,
    const int &iswap,
    const typename ArrayTypes<DeviceType>::t_x_array &x,
    const typename ArrayTypes<DeviceType>::t_tagint_1d &tag,
    const typename ArrayTypes<DeviceType>::t_int_1d &type,
    const typename ArrayTypes<DeviceType>::t_int_1d &mask,
    const typename ArrayTypes<DeviceType>::t_float_1d &radius,
    const typename ArrayTypes<DeviceType>::t_float_1d &rmass,
    const typename ArrayTypes<DeviceType>::t_float_2d &forcing,
    const typename ArrayTypes<DeviceType>::t_float_1d &mean_thickness,
    const typename ArrayTypes<DeviceType>::t_float_1d &min_thickness,
    const typename ArrayTypes<DeviceType>::t_float_1d &ridgingIceThickness,
    const typename ArrayTypes<DeviceType>::t_float_1d &ridgingIceThicknessWeight,
    const typename ArrayTypes<DeviceType>::t_float_1d &netToGrossClosingRatio,
    const typename ArrayTypes<DeviceType>::t_float_1d &changeEffectiveElementArea,
    const typename ArrayTypes<DeviceType>::t_float_1d &ice_area,
    const typename ArrayTypes<DeviceType>::t_float_1d &coriolis,
    const typename ArrayTypes<DeviceType>::t_float_2d &ocean_vel,
    const typename ArrayTypes<DeviceType>::t_float_2d &bvector,
    const typename ArrayTypes<DeviceType>::t_v_array &v,
    const typename ArrayTypes<DeviceType>::t_v_array &omega,
    const X_FLOAT &dx, const X_FLOAT &dy, const X_FLOAT &dz,
    const X_FLOAT &dvx, const X_FLOAT &dvy, const X_FLOAT &dvz,
    const int &deform_groupbit):
    _buf(buf),_list(list),_iswap(iswap),
    _x(x),_tag(tag),_type(type),_mask(mask),
    _radius(radius),
    _rmass(rmass),
    _forcing(forcing),
    _mean_thickness(mean_thickness),
    _min_thickness(min_thickness),
    _ridgingIceThickness(ridgingIceThickness),
    _ridgingIceThicknessWeight(ridgingIceThicknessWeight),
    _netToGrossClosingRatio(netToGrossClosingRatio),
    _changeEffectiveElementArea(changeEffectiveElementArea),
    _ice_area(ice_area), 
    _coriolis(coriolis),
    _ocean_vel(ocean_vel), 
    _bvector(bvector),
    _v(v), _omega(omega),
    _dx(dx),_dy(dy),_dz(dz),
    _dvx(dvx),_dvy(dvy),_dvz(dvz),
    _deform_groupbit(deform_groupbit)
  {
    const size_t elements = 28;
    const int maxsend = (buf.extent(0)*buf.extent(1))/elements;
    _buf = typename ArrayTypes<DeviceType>::t_xfloat_2d_um(buf.data(),maxsend,elements);
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const int& i) const {
    const int j = _list(_iswap,i);
    if (PBC_FLAG == 0) {
      _buf(i,0) = _x(j,0);
      _buf(i,1) = _x(j,1);
      _buf(i,2) = _x(j,2);
    } else {
      _buf(i,0) = _x(j,0) + _dx;
      _buf(i,1) = _x(j,1) + _dy;
      _buf(i,2) = _x(j,2) + _dz;
    }
    _buf(i,3) = d_ubuf(_tag(j)).d;
    _buf(i,4) = d_ubuf(_type(j)).d;
    _buf(i,5) = d_ubuf(_mask(j)).d;
    _buf(i,6) = _radius(j);
    _buf(i,7) = _rmass(j);
    _buf(i,8) = _forcing(j,0);
    _buf(i,9) = _forcing(j,1);
    _buf(i,10) = _mean_thickness(j);
    _buf(i,11) = _min_thickness(j);
    _buf(i,12) = _ridgingIceThickness(j);
    _buf(i,13) = _ridgingIceThicknessWeight(j);
    _buf(i,14) = _netToGrossClosingRatio(j);
    _buf(i,15) = _changeEffectiveElementArea(j);
    _buf(i,16) = _ice_area(j);
    _buf(i,17) = _coriolis(j);
    _buf(i,18) = _ocean_vel(j,0);
    _buf(i,19) = _ocean_vel(j,1);
    _buf(i,20) = _bvector(j,0);
    _buf(i,21) = _bvector(j,1);
    if (DEFORM_VREMAP) {
      if (_mask(i) & _deform_groupbit) {
        _buf(i,22) = _v(j,0) + _dvx;
        _buf(i,23) = _v(j,1) + _dvy;
        _buf(i,24) = _v(j,2) + _dvz;
      }
    }
    else {
      _buf(i,22) = _v(j,0);
      _buf(i,23) = _v(j,1);
      _buf(i,24) = _v(j,2);
    }
    _buf(i,25) = _omega(j,0);
    _buf(i,26) = _omega(j,1);
    _buf(i,27) = _omega(j,2);
  }
};

/* ---------------------------------------------------------------------- */

int AtomVecDemsiKokkos::pack_border_vel_kokkos(
  int n, DAT::tdual_int_2d k_sendlist, DAT::tdual_xfloat_2d buf,int iswap,
  int pbc_flag, int *pbc, ExecutionSpace space)
{
  X_FLOAT dx=0,dy=0,dz=0;
  X_FLOAT dvx=0,dvy=0,dvz=0;

  //sync(space,ALL_MASK);

  if (pbc_flag != 0) {
    if (domain->triclinic == 0) {
      dx = pbc[0]*domain->xprd;
      dy = pbc[1]*domain->yprd;
      dz = pbc[2]*domain->zprd;
    } else {
      dx = pbc[0];
      dy = pbc[1];
      dz = pbc[2];
    }
   if (!deform_vremap) {
      if(space==Host) {
        AtomVecDemsiKokkos_PackBorderVel<LMPHostType,1,0> f(
          buf.view<LMPHostType>(), k_sendlist.view<LMPHostType>(),
          iswap,h_x,h_tag,h_type,h_mask,
          h_radius,h_rmass,
          h_forcing,
          h_mean_thickness, h_min_thickness,
	  h_ridgingIceThickness,
	  h_ridgingIceThicknessWeight,
	  h_netToGrossClosingRatio,
	  h_changeEffectiveElementArea,
          h_ice_area, h_coriolis,
          h_ocean_vel, h_bvector,
          h_v, h_omega,
          dx,dy,dz,dvx,dvy,dvz,
          deform_groupbit);
        Kokkos::parallel_for(n,f);
      } else {
        AtomVecDemsiKokkos_PackBorderVel<LMPDeviceType,1,0> f(
          buf.view<LMPDeviceType>(), k_sendlist.view<LMPDeviceType>(),
          iswap,d_x,d_tag,d_type,d_mask,
          d_radius,d_rmass,
          d_forcing,
          d_mean_thickness, d_min_thickness,
	  d_ridgingIceThickness,
	  d_ridgingIceThicknessWeight,
	  d_netToGrossClosingRatio,
	  d_changeEffectiveElementArea,
          d_ice_area, d_coriolis,
          d_ocean_vel, d_bvector,
          d_v, d_omega,
          dx,dy,dz,dvx,dvy,dvz,
          deform_groupbit);
        Kokkos::parallel_for(n,f);
      }
    }
   else {
      dvx = pbc[0]*h_rate[0] + pbc[5]*h_rate[5] + pbc[4]*h_rate[4];
      dvy = pbc[1]*h_rate[1] + pbc[3]*h_rate[3];
      dvz = pbc[2]*h_rate[2];
      if(space==Host) {
        AtomVecDemsiKokkos_PackBorderVel<LMPHostType,1,1> f(
          buf.view<LMPHostType>(), k_sendlist.view<LMPHostType>(),
          iswap,h_x,h_tag,h_type,h_mask,
          h_radius,h_rmass,
          h_forcing,
          h_mean_thickness, h_min_thickness,
	  h_ridgingIceThickness,
	  h_ridgingIceThicknessWeight,
	  h_netToGrossClosingRatio,
	  h_changeEffectiveElementArea,
          h_ice_area, h_coriolis,
          h_ocean_vel, h_bvector,
          h_v, h_omega,
          dx,dy,dz,dvx,dvy,dvz,
          deform_groupbit);
        Kokkos::parallel_for(n,f);
      } else {
        AtomVecDemsiKokkos_PackBorderVel<LMPDeviceType,1,1> f(
          buf.view<LMPDeviceType>(), k_sendlist.view<LMPDeviceType>(),
          iswap,d_x,d_tag,d_type,d_mask,
          d_radius,d_rmass,
          d_forcing,
          d_mean_thickness, d_min_thickness,
	  d_ridgingIceThickness,
	  d_ridgingIceThicknessWeight,
	  d_netToGrossClosingRatio,
	  d_changeEffectiveElementArea,
          d_ice_area, d_coriolis,
          d_ocean_vel, d_bvector,
          d_v, d_omega,
          dx,dy,dz,dvx,dvy,dvz,
          deform_groupbit);
        Kokkos::parallel_for(n,f);
      }
    }
  } else {
    if(space==Host) {
      AtomVecDemsiKokkos_PackBorderVel<LMPHostType,0,0> f(
        buf.view<LMPHostType>(), k_sendlist.view<LMPHostType>(),
        iswap,h_x,h_tag,h_type,h_mask,
        h_radius,h_rmass,
        h_forcing,
        h_mean_thickness, h_min_thickness,
	h_ridgingIceThickness,
	h_ridgingIceThicknessWeight,
	h_netToGrossClosingRatio,
	h_changeEffectiveElementArea,
        h_ice_area, h_coriolis,
        h_ocean_vel, h_bvector,
        h_v, h_omega,
        dx,dy,dz,dvx,dvy,dvz,
        deform_groupbit);
      Kokkos::parallel_for(n,f);
    } else {
      AtomVecDemsiKokkos_PackBorderVel<LMPDeviceType,0,0> f(
        buf.view<LMPDeviceType>(), k_sendlist.view<LMPDeviceType>(),
        iswap,d_x,d_tag,d_type,d_mask,
        d_radius,d_rmass,
        d_forcing,
        d_mean_thickness, d_min_thickness,
	d_ridgingIceThickness,
	d_ridgingIceThicknessWeight,
	d_netToGrossClosingRatio,
	d_changeEffectiveElementArea,
        d_ice_area, d_coriolis,
        d_ocean_vel, d_bvector,
        d_v, d_omega,
        dx,dy,dz,dvx,dvy,dvz,
        deform_groupbit);
      Kokkos::parallel_for(n,f);
    }
  }

  return n*(size_border + size_velocity);
}

/* ---------------------------------------------------------------------- */

int AtomVecDemsiKokkos::pack_border_vel(int n, int *list, double *buf,
                                   int pbc_flag, int *pbc)
{
  int i,j,m;
  double dx,dy,dz,dvx,dvy,dvz;

  sync(Host,ALL_MASK);

  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = h_x(j,0);
      buf[m++] = h_x(j,1);
      buf[m++] = h_x(j,2);
      buf[m++] = ubuf(h_tag(j)).d;
      buf[m++] = ubuf(h_type(j)).d;
      buf[m++] = ubuf(h_mask(j)).d;
      buf[m++] = h_radius(j);
      buf[m++] = h_rmass(j);
      buf[m++] = h_forcing(j,0);
      buf[m++] = h_forcing(j,1);
      buf[m++] = h_mean_thickness(j);
      buf[m++] = h_min_thickness(j);
      buf[m++] = h_ridgingIceThickness(j);
      buf[m++] = h_ridgingIceThicknessWeight(j);
      buf[m++] = h_netToGrossClosingRatio(j);
      buf[m++] = h_changeEffectiveElementArea(j);
      buf[m++] = h_ice_area(j);
      buf[m++] = h_coriolis(j);
      buf[m++] = h_ocean_vel(j,0);
      buf[m++] = h_ocean_vel(j,1);
      buf[m++] = h_bvector(j,0);
      buf[m++] = h_bvector(j,1);
      buf[m++] = h_v(j,0);
      buf[m++] = h_v(j,1);
      buf[m++] = h_v(j,2);
      buf[m++] = h_omega(j,0);
      buf[m++] = h_omega(j,1);
      buf[m++] = h_omega(j,2);
    }
  } else {
    if (domain->triclinic == 0) {
      dx = pbc[0]*domain->xprd;
      dy = pbc[1]*domain->yprd;
      dz = pbc[2]*domain->zprd;
    } else {
      dx = pbc[0];
      dy = pbc[1];
      dz = pbc[2];
    }
    if (!deform_vremap) {
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = h_x(j,0) + dx;
        buf[m++] = h_x(j,1) + dy;
        buf[m++] = h_x(j,2) + dz;
        buf[m++] = ubuf(h_tag(j)).d;
        buf[m++] = ubuf(h_type(j)).d;
        buf[m++] = ubuf(h_mask(j)).d;
        buf[m++] = h_radius(j);
        buf[m++] = h_rmass(j);
        buf[m++] = h_forcing(j,0);
        buf[m++] = h_forcing(j,1);
        buf[m++] = h_mean_thickness(j);
        buf[m++] = h_min_thickness(j);
	buf[m++] = h_ridgingIceThickness(j);
	buf[m++] = h_ridgingIceThicknessWeight(j);
	buf[m++] = h_netToGrossClosingRatio(j);
	buf[m++] = h_changeEffectiveElementArea(j);
        buf[m++] = h_ice_area(j);
        buf[m++] = h_coriolis(j);
        buf[m++] = h_ocean_vel(j,0);
        buf[m++] = h_ocean_vel(j,1);
        buf[m++] = h_bvector(j,0);
        buf[m++] = h_bvector(j,1);
        buf[m++] = h_v(j,0);
        buf[m++] = h_v(j,1);
        buf[m++] = h_v(j,2);
        buf[m++] = h_omega(j,0);
        buf[m++] = h_omega(j,1);
        buf[m++] = h_omega(j,2);
      }
    } else {
      dvx = pbc[0]*h_rate[0] + pbc[5]*h_rate[5] + pbc[4]*h_rate[4];
      dvy = pbc[1]*h_rate[1] + pbc[3]*h_rate[3];
      dvz = pbc[2]*h_rate[2];
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = h_x(j,0) + dx;
        buf[m++] = h_x(j,1) + dy;
        buf[m++] = h_x(j,2) + dz;
        buf[m++] = ubuf(h_tag[j]).d;
        buf[m++] = ubuf(h_type[j]).d;
        buf[m++] = ubuf(h_mask[j]).d;
        buf[m++] = h_radius[j];
        buf[m++] = h_rmass[j];
        buf[m++] = h_forcing(j,0);
        buf[m++] = h_forcing(j,1);
        buf[m++] = h_mean_thickness(j);
        buf[m++] = h_min_thickness(j);
	buf[m++] = h_ridgingIceThickness(j);
	buf[m++] = h_ridgingIceThicknessWeight(j);
	buf[m++] = h_netToGrossClosingRatio(j);
	buf[m++] = h_changeEffectiveElementArea(j);
        buf[m++] = h_ice_area(j);
        buf[m++] = h_coriolis(j);
        buf[m++] = h_ocean_vel(j,0);
        buf[m++] = h_ocean_vel(j,1);
        buf[m++] = h_bvector(j,0);
        buf[m++] = h_bvector(j,1);
        if (mask[i] & deform_groupbit) {
          buf[m++] = h_v(j,0) + dvx;
          buf[m++] = h_v(j,1) + dvy;
          buf[m++] = h_v(j,2) + dvz;
        } else {
          buf[m++] = h_v(j,0);
          buf[m++] = h_v(j,1);
          buf[m++] = h_v(j,2);
        }
        buf[m++] = h_omega(j,0);
        buf[m++] = h_omega(j,1);
        buf[m++] = h_omega(j,2);
      }
    }
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->pack_border(n,list,&buf[m]);

  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecDemsiKokkos::pack_border_hybrid(int n, int *list, double *buf)
{
  sync(Host,RADIUS_MASK|RMASS_MASK);

  int m = 0;
  for (int i = 0; i < n; i++) {
    const int j = list[i];
    buf[m++] = h_radius(j);
    buf[m++] = h_rmass(j);
  }
  return m;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
struct AtomVecDemsiKokkos_UnpackBorder {
  typedef DeviceType device_type;

  const typename ArrayTypes<DeviceType>::t_xfloat_2d_const _buf;
  typename ArrayTypes<DeviceType>::t_x_array _x;
  typename ArrayTypes<DeviceType>::t_tagint_1d _tag;
  typename ArrayTypes<DeviceType>::t_int_1d _type;
  typename ArrayTypes<DeviceType>::t_int_1d _mask;
  typename ArrayTypes<DeviceType>::t_float_1d _radius,_rmass;
  typename ArrayTypes<DeviceType>::t_float_2d _forcing;
  typename ArrayTypes<DeviceType>::t_float_1d _mean_thickness,_min_thickness;
  typename ArrayTypes<DeviceType>::t_float_1d _ridgingIceThickness;
  typename ArrayTypes<DeviceType>::t_float_1d _ridgingIceThicknessWeight;
  typename ArrayTypes<DeviceType>::t_float_1d _netToGrossClosingRatio;
  typename ArrayTypes<DeviceType>::t_float_1d _changeEffectiveElementArea;
  typename ArrayTypes<DeviceType>::t_float_1d _ice_area,_coriolis;
  typename ArrayTypes<DeviceType>::t_float_2d _ocean_vel,_bvector;
  int _first;

  AtomVecDemsiKokkos_UnpackBorder(
    const typename ArrayTypes<DeviceType>::t_xfloat_2d_const &buf,
    typename ArrayTypes<DeviceType>::t_x_array &x,
    typename ArrayTypes<DeviceType>::t_tagint_1d &tag,
    typename ArrayTypes<DeviceType>::t_int_1d &type,
    typename ArrayTypes<DeviceType>::t_int_1d &mask,
    const typename ArrayTypes<DeviceType>::t_float_1d &radius,
    const typename ArrayTypes<DeviceType>::t_float_1d &rmass,
    typename ArrayTypes<DeviceType>::t_float_2d &forcing,
    typename ArrayTypes<DeviceType>::t_float_1d &mean_thickness,
    typename ArrayTypes<DeviceType>::t_float_1d &min_thickness,
    typename ArrayTypes<DeviceType>::t_float_1d &ridgingIceThickness,
    typename ArrayTypes<DeviceType>::t_float_1d &ridgingIceThicknessWeight,
    typename ArrayTypes<DeviceType>::t_float_1d &netToGrossClosingRatio,
    typename ArrayTypes<DeviceType>::t_float_1d &changeEffectiveElementArea,
    typename ArrayTypes<DeviceType>::t_float_1d &ice_area,
    typename ArrayTypes<DeviceType>::t_float_1d &coriolis,
    typename ArrayTypes<DeviceType>::t_float_2d &ocean_vel,
    typename ArrayTypes<DeviceType>::t_float_2d &bvector,
    const int& first):
    _buf(buf),_x(x),_tag(tag),_type(type),_mask(mask),
    _radius(radius),
    _rmass(rmass),_forcing(forcing),
    _mean_thickness(mean_thickness),_min_thickness(min_thickness),
    _ridgingIceThickness(ridgingIceThickness),
    _ridgingIceThicknessWeight(ridgingIceThicknessWeight),
    _netToGrossClosingRatio(netToGrossClosingRatio),
    _changeEffectiveElementArea(changeEffectiveElementArea),
    _ice_area(ice_area),_coriolis(coriolis),
    _ocean_vel(ocean_vel),_bvector(bvector),
    _first(first) {};

  KOKKOS_INLINE_FUNCTION
  void operator() (const int& i) const {
    _x(i+_first,0) = _buf(i,0);
    _x(i+_first,1) = _buf(i,1);
    _x(i+_first,2) = _buf(i,2);
    _tag(i+_first) = static_cast<tagint> (d_ubuf(_buf(i,3)).i);
    _type(i+_first) = static_cast<int>  (d_ubuf(_buf(i,4)).i);
    _mask(i+_first) = static_cast<int>  (d_ubuf(_buf(i,5)).i);
    _radius(i+_first) = _buf(i,6);
    _rmass(i+_first) = _buf(i,7);
    _forcing(i+_first,0) = _buf(i,8);
    _forcing(i+_first,1) = _buf(i,9);
    _mean_thickness(i+_first) = _buf(i,10);
    _min_thickness(i+_first) = _buf(i,11);
    _ridgingIceThickness(i+_first) = _buf(i,12);
    _ridgingIceThicknessWeight(i+_first) = _buf(i,13);
    _netToGrossClosingRatio(i+_first) = _buf(i,14);
    _changeEffectiveElementArea(i+_first) = _buf(i,15);
    _ice_area(i+_first) = _buf(i,16);
    _coriolis(i+_first) = _buf(i,17);
    _ocean_vel(i+_first,0) = _buf(i,18);
    _ocean_vel(i+_first,1) = _buf(i,19);
    _bvector(i+_first,0) = _buf(i,20);
    _bvector(i+_first,1) = _buf(i,21);
  }
};

/* ---------------------------------------------------------------------- */

void AtomVecDemsiKokkos::unpack_border_kokkos(const int &n, const int &first,
					       const DAT::tdual_xfloat_2d &buf,
                                               ExecutionSpace space) {

  while (first+n >= nmax) grow(0);
  if(space==Host) {
    struct AtomVecDemsiKokkos_UnpackBorder<LMPHostType> f(buf.view<LMPHostType>(),
      h_x,h_tag,h_type,h_mask,
      h_radius,h_rmass,h_forcing,h_mean_thickness,h_min_thickness,
      h_ridgingIceThickness,h_ridgingIceThicknessWeight,h_netToGrossClosingRatio,h_changeEffectiveElementArea,
      h_ice_area,h_coriolis,h_ocean_vel,h_bvector,first);
    Kokkos::parallel_for(n,f);
  } else {
    struct AtomVecDemsiKokkos_UnpackBorder<LMPDeviceType> f(buf.view<LMPDeviceType>(),
      d_x,d_tag,d_type,d_mask,
      d_radius,d_rmass,d_forcing,d_mean_thickness,d_min_thickness,
      d_ridgingIceThickness,d_ridgingIceThicknessWeight,d_netToGrossClosingRatio,d_changeEffectiveElementArea,
      d_ice_area,d_coriolis,d_ocean_vel,d_bvector,first);
    Kokkos::parallel_for(n,f);
  }

  modified(space,X_MASK|TAG_MASK|TYPE_MASK|MASK_MASK|
	         RADIUS_MASK|RMASS_MASK|FORCING_MASK|
                 THICKNESS_MASK);
}

/* ---------------------------------------------------------------------- */

void AtomVecDemsiKokkos::unpack_border(int n, int first, double *buf)
{
  int m = 0;
  const int last = first + n;
  for (int i = first; i < last; i++) {
    if (i == nmax) grow(0);
    h_x(i,0) = buf[m++];
    h_x(i,1) = buf[m++];
    h_x(i,2) = buf[m++];
    h_tag[i] = (tagint) ubuf(buf[m++]).i;
    h_type[i] = (int) ubuf(buf[m++]).i;
    h_mask[i] = (int) ubuf(buf[m++]).i;
    h_radius[i] = buf[m++];
    h_rmass[i] = buf[m++];
    h_forcing(i,0) = buf[m++];
    h_forcing(i,1) = buf[m++];
    h_mean_thickness(i) = buf[m++];
    h_min_thickness(i) = buf[m++];
    h_ridgingIceThickness(i) = buf[m++];
    h_ridgingIceThicknessWeight(i) = buf[m++];
    h_netToGrossClosingRatio(i) = buf[m++];
    h_changeEffectiveElementArea(i) = buf[m++];
    h_ice_area(i) = buf[m++];
    h_coriolis(i) = buf[m++];
    h_ocean_vel(i,0) = buf[m++];
    h_ocean_vel(i,1) = buf[m++];
    h_bvector(i,0) = buf[m++];
    h_bvector(i,1) = buf[m++];
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->
        unpack_border(n,first,&buf[m]);

  modified(Host,X_MASK|TAG_MASK|TYPE_MASK|MASK_MASK|RADIUS_MASK|RMASS_MASK|
                FORCING_MASK|THICKNESS_MASK);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
struct AtomVecDemsiKokkos_UnpackBorderVel {
  typedef DeviceType device_type;

  typename ArrayTypes<DeviceType>::t_xfloat_2d_const_um _buf;
  typename ArrayTypes<DeviceType>::t_x_array _x;
  typename ArrayTypes<DeviceType>::t_tagint_1d _tag;
  typename ArrayTypes<DeviceType>::t_int_1d _type;
  typename ArrayTypes<DeviceType>::t_int_1d _mask;
  typename ArrayTypes<DeviceType>::t_float_1d _radius,_rmass;
  typename ArrayTypes<DeviceType>::t_float_2d _forcing;
  typename ArrayTypes<DeviceType>::t_float_1d _mean_thickness,_min_thickness;
  typename ArrayTypes<DeviceType>::t_float_1d _ridgingIceThickness;
  typename ArrayTypes<DeviceType>::t_float_1d _ridgingIceThicknessWeight;
  typename ArrayTypes<DeviceType>::t_float_1d _netToGrossClosingRatio;
  typename ArrayTypes<DeviceType>::t_float_1d _changeEffectiveElementArea;
  typename ArrayTypes<DeviceType>::t_float_1d _ice_area,_coriolis;
  typename ArrayTypes<DeviceType>::t_float_2d _ocean_vel,_bvector;
  typename ArrayTypes<DeviceType>::t_v_array _v;
  typename ArrayTypes<DeviceType>::t_v_array _omega;
  int _first;

  AtomVecDemsiKokkos_UnpackBorderVel(
    const typename ArrayTypes<DeviceType>::t_xfloat_2d_const &buf,
    const typename ArrayTypes<DeviceType>::t_x_array &x,
    const typename ArrayTypes<DeviceType>::t_tagint_1d &tag,
    const typename ArrayTypes<DeviceType>::t_int_1d &type,
    const typename ArrayTypes<DeviceType>::t_int_1d &mask,
    const typename ArrayTypes<DeviceType>::t_float_1d &radius,
    const typename ArrayTypes<DeviceType>::t_float_1d &rmass,
    const typename ArrayTypes<DeviceType>::t_float_2d &forcing,
    const typename ArrayTypes<DeviceType>::t_float_1d &mean_thickness,
    const typename ArrayTypes<DeviceType>::t_float_1d &min_thickness,
    const typename ArrayTypes<DeviceType>::t_float_1d &ridgingIceThickness,
    const typename ArrayTypes<DeviceType>::t_float_1d &ridgingIceThicknessWeight,
    const typename ArrayTypes<DeviceType>::t_float_1d &netToGrossClosingRatio,
    const typename ArrayTypes<DeviceType>::t_float_1d &changeEffectiveElementArea,
    const typename ArrayTypes<DeviceType>::t_float_1d &ice_area,
    const typename ArrayTypes<DeviceType>::t_float_1d &coriolis,
    const typename ArrayTypes<DeviceType>::t_float_2d &ocean_vel,
    const typename ArrayTypes<DeviceType>::t_float_2d &bvector,
    const typename ArrayTypes<DeviceType>::t_v_array &v,
    const typename ArrayTypes<DeviceType>::t_v_array &omega,
    const int& first):
    _buf(buf),_x(x),_tag(tag),_type(type),_mask(mask),
    _radius(radius),
    _rmass(rmass),
    _forcing(forcing),
    _mean_thickness(mean_thickness),
    _min_thickness(min_thickness),
    _ridgingIceThickness(ridgingIceThickness),
    _ridgingIceThicknessWeight(ridgingIceThicknessWeight),
    _netToGrossClosingRatio(netToGrossClosingRatio),
    _changeEffectiveElementArea(changeEffectiveElementArea),
    _ice_area(ice_area),
    _coriolis(coriolis),
    _ocean_vel(ocean_vel),
    _bvector(bvector),
    _v(v), _omega(omega),
    _first(first)
  {
    const size_t elements = 28;
    const int maxsend = (buf.extent(0)*buf.extent(1))/elements;
    _buf = typename ArrayTypes<DeviceType>::t_xfloat_2d_const_um(buf.data(),maxsend,elements);
  };

  KOKKOS_INLINE_FUNCTION
  void operator() (const int& i) const {
    _x(i+_first,0) = _buf(i,0);
    _x(i+_first,1) = _buf(i,1);
    _x(i+_first,2) = _buf(i,2);
    _tag(i+_first) = static_cast<tagint> (d_ubuf(_buf(i,3)).i);
    _type(i+_first) = static_cast<int>  (d_ubuf(_buf(i,4)).i);
    _mask(i+_first) = static_cast<int>  (d_ubuf(_buf(i,5)).i);
    _radius(i+_first) = _buf(i,6);
    _rmass(i+_first) = _buf(i,7);
    _forcing(i+_first,0) = _buf(i,8);
    _forcing(i+_first,1) = _buf(i,9);
    _mean_thickness(i+_first) = _buf(i,10);
    _min_thickness(i+_first) = _buf(i,11);
    _ridgingIceThickness(i+_first) = _buf(i,12);
    _ridgingIceThicknessWeight(i+_first) = _buf(i,13);
    _netToGrossClosingRatio(i+_first) = _buf(i,14);
    _changeEffectiveElementArea(i+_first) = _buf(i,15);
    _ice_area(i+_first) = _buf(i,16);
    _coriolis(i+_first) = _buf(i,17);
    _ocean_vel(i+_first,0) = _buf(i,18);
    _ocean_vel(i+_first,1) = _buf(i,19);
    _bvector(i+_first,0) = _buf(i,20);
    _bvector(i+_first,1) = _buf(i,21);
    _v(i+_first,0) = _buf(i,22);
    _v(i+_first,1) = _buf(i,23);
    _v(i+_first,2) = _buf(i,24);
    _omega(i+_first,0) = _buf(i,25);
    _omega(i+_first,1) = _buf(i,26);
    _omega(i+_first,2) = _buf(i,27);
  }
};

/* ---------------------------------------------------------------------- */

void AtomVecDemsiKokkos::unpack_border_vel_kokkos(
  const int &n, const int &first,
  const DAT::tdual_xfloat_2d &buf,ExecutionSpace space) {
  while (first+n >= nmax) grow(0);
  if(space==Host) {
    struct AtomVecDemsiKokkos_UnpackBorderVel<LMPHostType> f(buf.view<LMPHostType>(),
      h_x,h_tag,h_type,h_mask,
      h_radius,h_rmass,
      h_forcing,
      h_mean_thickness,h_min_thickness,
      h_ridgingIceThickness,h_ridgingIceThicknessWeight,h_netToGrossClosingRatio,h_changeEffectiveElementArea,
      h_ice_area,h_coriolis,
      h_ocean_vel,h_bvector,
      h_v, h_omega,
      first);
    Kokkos::parallel_for(n,f);
  } else {
    struct AtomVecDemsiKokkos_UnpackBorderVel<LMPDeviceType> f(buf.view<LMPDeviceType>(),
      d_x,d_tag,d_type,d_mask,
      d_radius,d_rmass,
      d_forcing,
      d_mean_thickness,d_min_thickness,
      d_ridgingIceThickness,d_ridgingIceThicknessWeight,d_netToGrossClosingRatio,d_changeEffectiveElementArea,
      d_ice_area,d_coriolis,
      d_ocean_vel,d_bvector,
      d_v, d_omega,
      first);
    Kokkos::parallel_for(n,f);
  }

  modified(space,X_MASK|TAG_MASK|TYPE_MASK|MASK_MASK|
                 RADIUS_MASK|RMASS_MASK|V_MASK|OMEGA_MASK|
                 THICKNESS_MASK|FORCING_MASK);
}

/* ---------------------------------------------------------------------- */

void AtomVecDemsiKokkos::unpack_border_vel(int n, int first, double *buf)
{
  int m = 0;
  const int last = first + n;
  for (int i = first; i < last; i++) {
    if (i == nmax) grow(0);
    h_x(i,0) = buf[m++];
    h_x(i,1) = buf[m++];
    h_x(i,2) = buf[m++];
    h_tag(i) = (tagint) ubuf(buf[m++]).i;
    h_type(i) = (int) ubuf(buf[m++]).i;
    h_mask(i) = (int) ubuf(buf[m++]).i;
    h_radius(i) = buf[m++];
    h_rmass(i) = buf[m++];
    h_forcing(i,0) = buf[m++];
    h_forcing(i,1) = buf[m++];
    h_mean_thickness(i) = buf[m++];
    h_min_thickness(i) = buf[m++];
    h_ridgingIceThickness(i) = buf[m++];
    h_ridgingIceThicknessWeight(i) = buf[m++];
    h_netToGrossClosingRatio(i) = buf[m++];
    h_changeEffectiveElementArea(i) = buf[m++];
    h_ice_area(i) = buf[m++];
    h_coriolis(i) = buf[m++];
    h_ocean_vel(i,0) = buf[m++];
    h_ocean_vel(i,1) = buf[m++];
    h_bvector(i,0) = buf[m++];
    h_bvector(i,1) = buf[m++];
    h_v(i,0) = buf[m++];
    h_v(i,1) = buf[m++];
    h_v(i,2) = buf[m++];
    h_omega(i,0) = buf[m++];
    h_omega(i,1) = buf[m++];
    h_omega(i,2) = buf[m++];
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->
        unpack_border(n,first,&buf[m]);

  modified(Host,X_MASK|TAG_MASK|TYPE_MASK|MASK_MASK|RADIUS_MASK|RMASS_MASK|
                V_MASK|OMEGA_MASK|FORCING_MASK|THICKNESS_MASK);
}

/* ---------------------------------------------------------------------- */

int AtomVecDemsiKokkos::unpack_border_hybrid(int n, int first, double *buf)
{

  int i,m,last;

  m = 0;
  last = first + n;
  for (int i = first; i < last; i++) {
    h_radius(i) = buf[m++];
    h_rmass(i) = buf[m++];
  }
  modified(Host,RADIUS_MASK|RMASS_MASK);
  return m;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
struct AtomVecDemsiKokkos_PackExchangeFunctor {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typename AT::t_x_array_randomread _x;
  typename AT::t_v_array_randomread _v;
  typename AT::t_tagint_1d_randomread _tag;
  typename AT::t_int_1d_randomread _type;
  typename AT::t_int_1d_randomread _mask;
  typename AT::t_imageint_1d_randomread _image;
  typename AT::t_float_1d_randomread _radius,_rmass;
  typename AT::t_v_array_randomread _omega;
  typename AT::t_float_2d_randomread _forcing;
  typename AT::t_float_1d_randomread _mean_thickness,_min_thickness;
  typename AT::t_float_1d_randomread _ridgingIceThickness;
  typename AT::t_float_1d_randomread _ridgingIceThicknessWeight;
  typename AT::t_float_1d_randomread _netToGrossClosingRatio;
  typename AT::t_float_1d_randomread _changeEffectiveElementArea;
  typename AT::t_float_1d_randomread _ice_area,_coriolis;
  typename AT::t_float_2d_randomread _ocean_vel,_bvector;
  typename AT::t_int_1d_randomread _num_bond;
  typename AT::t_int_2d_randomread _bond_type;
  typename AT::t_tagint_2d_randomread _bond_atom;
  typename AT::t_int_2d_randomread _nspecial;
  typename AT::t_tagint_2d_randomread _special;
  typename AT::t_x_array _xw;
  typename AT::t_v_array _vw;
  typename AT::t_tagint_1d _tagw;
  typename AT::t_int_1d _typew;
  typename AT::t_int_1d _maskw;
  typename AT::t_imageint_1d _imagew;
  typename AT::t_float_1d _radiusw,_rmassw;
  typename AT::t_v_array _omegaw;
  typename AT::t_float_2d _forcingw;
  typename AT::t_float_1d _mean_thicknessw,_min_thicknessw;
  typename AT::t_float_1d _ridgingIceThicknessw;
  typename AT::t_float_1d _ridgingIceThicknessWeightw;
  typename AT::t_float_1d _netToGrossClosingRatiow;
  typename AT::t_float_1d _changeEffectiveElementAreaw;
  typename AT::t_float_1d _ice_areaw,_coriolisw;
  typename AT::t_float_2d _ocean_velw,_bvectorw;
  typename AT::t_int_1d _num_bondw;
  typename AT::t_int_2d _bond_typew;
  typename AT::t_tagint_2d _bond_atomw;
  typename AT::t_int_2d _nspecialw;
  typename AT::t_tagint_2d _specialw;
  typename AT::t_xfloat_2d_um _buf;
  typename AT::t_int_1d_const _sendlist;
  typename AT::t_int_1d_const _copylist;
  int _nlocal,_dim;
  X_FLOAT _lo,_hi;
  size_t elements;

  AtomVecDemsiKokkos_PackExchangeFunctor(
    const AtomKokkos* atom,
    const typename AT::tdual_xfloat_2d buf,
    typename AT::tdual_int_1d sendlist,
    typename AT::tdual_int_1d copylist,int nlocal, int dim,X_FLOAT lo, X_FLOAT hi):
    _x(atom->k_x.view<DeviceType>()),
    _v(atom->k_v.view<DeviceType>()),
    _tag(atom->k_tag.view<DeviceType>()),
    _type(atom->k_type.view<DeviceType>()),
    _mask(atom->k_mask.view<DeviceType>()),
    _image(atom->k_image.view<DeviceType>()),
    _radius(atom->k_radius.view<DeviceType>()),
    _rmass(atom->k_rmass.view<DeviceType>()),
    _omega(atom->k_omega.view<DeviceType>()),
    _forcing(atom->k_forcing.view<DeviceType>()),
    _mean_thickness(atom->k_mean_thickness.view<DeviceType>()),
    _min_thickness(atom->k_min_thickness.view<DeviceType>()),
    _ridgingIceThickness(atom->k_ridgingIceThickness.view<DeviceType>()),
    _ridgingIceThicknessWeight(atom->k_ridgingIceThicknessWeight.view<DeviceType>()),
    _netToGrossClosingRatio(atom->k_netToGrossClosingRatio.view<DeviceType>()),
    _changeEffectiveElementArea(atom->k_changeEffectiveElementArea.view<DeviceType>()),
    _ice_area(atom->k_ice_area.view<DeviceType>()),
    _coriolis(atom->k_coriolis.view<DeviceType>()),
    _ocean_vel(atom->k_ocean_vel.view<DeviceType>()),
    _bvector(atom->k_bvector.view<DeviceType>()),
    _num_bond(atom->k_num_bond.view<DeviceType>()),
    _bond_type(atom->k_bond_type.view<DeviceType>()),
    _bond_atom(atom->k_bond_atom.view<DeviceType>()),
    _nspecial(atom->k_nspecial.view<DeviceType>()),
    _special(atom->k_special.view<DeviceType>()),
    _xw(atom->k_x.view<DeviceType>()),
    _vw(atom->k_v.view<DeviceType>()),
    _tagw(atom->k_tag.view<DeviceType>()),
    _typew(atom->k_type.view<DeviceType>()),
    _maskw(atom->k_mask.view<DeviceType>()),
    _imagew(atom->k_image.view<DeviceType>()),
    _radiusw(atom->k_radius.view<DeviceType>()),
    _rmassw(atom->k_rmass.view<DeviceType>()),
    _omegaw(atom->k_omega.view<DeviceType>()),
    _forcingw(atom->k_forcing.view<DeviceType>()),
    _mean_thicknessw(atom->k_mean_thickness.view<DeviceType>()),
    _min_thicknessw(atom->k_min_thickness.view<DeviceType>()),
    _ridgingIceThicknessw(atom->k_ridgingIceThickness.view<DeviceType>()),
    _ridgingIceThicknessWeightw(atom->k_ridgingIceThicknessWeight.view<DeviceType>()),
    _netToGrossClosingRatiow(atom->k_netToGrossClosingRatio.view<DeviceType>()),
    _changeEffectiveElementAreaw(atom->k_changeEffectiveElementArea.view<DeviceType>()),
    _ice_areaw(atom->k_ice_area.view<DeviceType>()),
    _coriolisw(atom->k_coriolis.view<DeviceType>()),
    _ocean_velw(atom->k_ocean_vel.view<DeviceType>()),
    _bvectorw(atom->k_bvector.view<DeviceType>()),
    _num_bondw(atom->k_num_bond.view<DeviceType>()),
    _bond_typew(atom->k_bond_type.view<DeviceType>()),
    _bond_atomw(atom->k_bond_atom.view<DeviceType>()),
    _nspecialw(atom->k_nspecial.view<DeviceType>()),
    _specialw(atom->k_special.view<DeviceType>()),
    _sendlist(sendlist.template view<DeviceType>()),
    _copylist(copylist.template view<DeviceType>()),
    _nlocal(nlocal),_dim(dim),
    _lo(lo),_hi(hi){
    elements = 33+atom->maxspecial+2*atom->bond_per_atom;
    const int maxsend = (buf.template view<DeviceType>().extent(0)*buf.template view<DeviceType>().extent(1))/elements;
    _buf = typename AT::t_xfloat_2d_um(buf.template view<DeviceType>().data(),maxsend,elements);
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const int &mysend) const {
    const int i = _sendlist(mysend);
    int m = 1;
    _buf(mysend,0) = elements;
    _buf(mysend,m++) = _x(i,0);
    _buf(mysend,m++) = _x(i,1);
    _buf(mysend,m++) = _x(i,2);
    _buf(mysend,m++) = _v(i,0);
    _buf(mysend,m++) = _v(i,1);
    _buf(mysend,m++) = _v(i,2);
    _buf(mysend,m++) = d_ubuf(_tag(i)).d;
    _buf(mysend,m++) = d_ubuf(_type(i)).d;
    _buf(mysend,m++) = d_ubuf(_mask(i)).d;
    _buf(mysend,m++) = d_ubuf(_image(i)).d;
    _buf(mysend,m++) = _radius(i);
    _buf(mysend,m++) = _rmass(i);
    _buf(mysend,m++) = _omega(i,0);
    _buf(mysend,m++) = _omega(i,1);
    _buf(mysend,m++) = _omega(i,2);
    _buf(mysend,m++) = _forcing(i,0);
    _buf(mysend,m++) = _forcing(i,1);
    _buf(mysend,m++) = _mean_thickness(i);
    _buf(mysend,m++) = _min_thickness(i);
    _buf(mysend,m++) = _ridgingIceThickness(i);
    _buf(mysend,m++) = _ridgingIceThicknessWeight(i);
    _buf(mysend,m++) = _netToGrossClosingRatio(i);
    _buf(mysend,m++) = _changeEffectiveElementArea(i);
    _buf(mysend,m++) = _ice_area(i);
    _buf(mysend,m++) = _coriolis(i);
    _buf(mysend,m++) = _ocean_vel(i,0);
    _buf(mysend,m++) = _ocean_vel(i,1);
    _buf(mysend,m++) = _bvector(i,0);
    _buf(mysend,m++) = _bvector(i,1);
    _buf(mysend,m++) = d_ubuf(_num_bond(i)).d;
    for (int k = 0; k < _num_bond(i); k++) {
      _buf(mysend,m++) = d_ubuf(_bond_type(i,k)).d;
      _buf(mysend,m++) = d_ubuf(_bond_atom(i,k)).d;
    }
    _buf(mysend,m++) = d_ubuf(_nspecial(i,0)).d;
    _buf(mysend,m++) = d_ubuf(_nspecial(i,1)).d;
    _buf(mysend,m++) = d_ubuf(_nspecial(i,2)).d;
    for (int k = 0; k < _nspecial(i,2); k++)
      _buf(mysend,m++) = d_ubuf(_special(i,k)).d;

    const int j = _copylist(mysend);

    if (j>-1) {
      _xw(i,0) = _x(j,0);
      _xw(i,1) = _x(j,1);
      _xw(i,2) = _x(j,2);
      _vw(i,0) = _v(j,0);
      _vw(i,1) = _v(j,1);
      _vw(i,2) = _v(j,2);
      _tagw(i) = _tag(j);
      _typew(i) = _type(j);
      _maskw(i) = _mask(j);
      _imagew(i) = _image(j);
      _radiusw(i) = _radius(j);
      _rmassw(i) = _rmass(j);
      _omegaw(i,0) = _omega(j,0);
      _omegaw(i,1) = _omega(j,1);
      _omegaw(i,2) = _omega(j,2);
      _forcingw(i,0) = _forcing(j,0);
      _forcingw(i,1) = _forcing(j,1);
      _mean_thicknessw(i) = _mean_thickness(j);
      _min_thicknessw(i) = _min_thickness(j);
      _ridgingIceThicknessw(i) = _ridgingIceThickness(j);
      _ridgingIceThicknessWeightw(i) = _ridgingIceThicknessWeight(j);
      _netToGrossClosingRatiow(i) = _netToGrossClosingRatio(j);
      _changeEffectiveElementAreaw(i) = _changeEffectiveElementArea(j);
      _ice_areaw(i) = _ice_area(j);
      _coriolisw(i) = _coriolis(j);
      _ocean_velw(i,0) = _ocean_vel(j,0);
      _ocean_velw(i,1) = _ocean_vel(j,1);
      _bvectorw(i,0) = _bvector(j,0);
      _bvectorw(i,1) = _bvector(j,1);
      _num_bondw(i) = _num_bond(j);
      for (int k = 0; k < _num_bond(j); k++) {
        _bond_typew(i,k) = _bond_type(j,k);
        _bond_atomw(i,k) = _bond_atom(j,k);
      }
      _nspecialw(i,0) = _nspecial(j,0);
      _nspecialw(i,1) = _nspecial(j,1);
      _nspecialw(i,2) = _nspecial(j,2);
      for (int k = 0; k < _nspecial(j,2); k++)
        _specialw(i,k) = _special(j,k);
    }
  }
};

/* ---------------------------------------------------------------------- */

int AtomVecDemsiKokkos::pack_exchange_kokkos(const int &nsend,DAT::tdual_xfloat_2d &k_buf, 
                                             DAT::tdual_int_1d k_sendlist,DAT::tdual_int_1d k_copylist,
                                             ExecutionSpace space,int dim,X_FLOAT lo,X_FLOAT hi) 
{


  const int elements = 33+atom->maxspecial+2*atom->bond_per_atom;

  if(nsend > (int) (k_buf.view<LMPHostType>().extent(0)*k_buf.view<LMPHostType>().extent(1))/elements) {
    int newsize = nsend*elements/k_buf.view<LMPHostType>().extent(1)+1;
    k_buf.resize(newsize,k_buf.view<LMPHostType>().extent(1));
  }

  sync(space,X_MASK | V_MASK | TAG_MASK | TYPE_MASK |
             MASK_MASK | IMAGE_MASK| RADIUS_MASK | RMASS_MASK |
             OMEGA_MASK | THICKNESS_MASK | FORCING_MASK | 
             BOND_MASK | SPECIAL_MASK);

  if(space == Host) {
    AtomVecDemsiKokkos_PackExchangeFunctor<LMPHostType> f(atomKK,k_buf,k_sendlist,k_copylist,atom->nlocal,dim,lo,hi);
    Kokkos::parallel_for(nsend,f);
  } else {
    AtomVecDemsiKokkos_PackExchangeFunctor<LMPDeviceType> f(atomKK,k_buf,k_sendlist,k_copylist,atom->nlocal,dim,lo,hi);
    Kokkos::parallel_for(nsend,f);
  }
  return nsend*elements;
}

/* ----------------------------------------------------------------------
   pack data for atom I for sending to another proc
   xyz must be 1st 3 values, so comm::exchange() can test on them
------------------------------------------------------------------------- */

int AtomVecDemsiKokkos::pack_exchange(int i, double *buf)
{

 //  Need this? 
  sync(Host,X_MASK | V_MASK | TAG_MASK | TYPE_MASK |
            MASK_MASK | IMAGE_MASK| RADIUS_MASK | RMASS_MASK |
            OMEGA_MASK | FORCING_MASK | THICKNESS_MASK | 
            BOND_MASK | SPECIAL_MASK);

  int m = 1;
  buf[m++] = h_x(i,0);
  buf[m++] = h_x(i,1);
  buf[m++] = h_x(i,2);
  buf[m++] = h_v(i,0);
  buf[m++] = h_v(i,1);
  buf[m++] = h_v(i,2);
  buf[m++] = ubuf(h_tag(i)).d;
  buf[m++] = ubuf(h_type(i)).d;
  buf[m++] = ubuf(h_mask(i)).d;
  buf[m++] = ubuf(h_image(i)).d;

  buf[m++] = h_radius(i);
  buf[m++] = h_rmass(i);
  buf[m++] = h_omega(i,0);
  buf[m++] = h_omega(i,1);
  buf[m++] = h_omega(i,2);

  buf[m++] = h_forcing(i,0);
  buf[m++] = h_forcing(i,1);
  buf[m++] = h_mean_thickness(i);
  buf[m++] = h_min_thickness(i);
  buf[m++] = h_ridgingIceThickness(i);
  buf[m++] = h_ridgingIceThicknessWeight(i);
  buf[m++] = h_netToGrossClosingRatio(i);
  buf[m++] = h_changeEffectiveElementArea(i);
  buf[m++] = h_ice_area(i);
  buf[m++] = h_coriolis(i);
  buf[m++] = h_ocean_vel(i,0);
  buf[m++] = h_ocean_vel(i,1);
  buf[m++] = h_bvector(i,0);
  buf[m++] = h_bvector(i,1);

  buf[m++] = ubuf(h_num_bond(i)).d;
  for (int k = 0; k < h_num_bond(i); k++) {
    buf[m++] = ubuf(h_bond_type(i,k)).d;
    buf[m++] = ubuf(h_bond_atom(i,k)).d;
  }

  buf[m++] = ubuf(h_nspecial(i,0)).d;
  buf[m++] = ubuf(h_nspecial(i,1)).d;
  buf[m++] = ubuf(h_nspecial(i,2)).d;
  for (int k = 0; k < h_nspecial(i,2); k++)
    buf[m++] = ubuf(h_special(i,k)).d;

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      m += modify->fix[atom->extra_grow[iextra]]->pack_exchange(i,&buf[m]);

  buf[0] = m;
  return m;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
struct AtomVecDemsiKokkos_UnpackExchangeFunctor {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typename AT::t_x_array _x;
  typename AT::t_v_array _v;
  typename AT::t_tagint_1d _tag;
  typename AT::t_int_1d _type;
  typename AT::t_int_1d _mask;
  typename AT::t_imageint_1d _image;
  typename AT::t_float_1d _radius;
  typename AT::t_float_1d _rmass;
  typename AT::t_v_array _omega;
  typename AT::t_float_2d _forcing;
  typename AT::t_float_1d _mean_thickness;
  typename AT::t_float_1d _min_thickness;
  typename AT::t_float_1d _ridgingIceThickness;
  typename AT::t_float_1d _ridgingIceThicknessWeight;
  typename AT::t_float_1d _netToGrossClosingRatio;
  typename AT::t_float_1d _changeEffectiveElementArea;
  typename AT::t_float_1d _ice_area;
  typename AT::t_float_1d _coriolis;
  typename AT::t_float_2d _ocean_vel;
  typename AT::t_float_2d _bvector;
  typename AT::t_int_1d _num_bond;
  typename AT::t_int_2d _bond_type;
  typename AT::t_tagint_2d _bond_atom;
  typename AT::t_int_2d _nspecial;
  typename AT::t_tagint_2d _special;

  typename AT::t_xfloat_2d_um _buf;
  typename AT::t_int_1d _nlocal;
  int _dim;
  X_FLOAT _lo,_hi;
  size_t elements;

  AtomVecDemsiKokkos_UnpackExchangeFunctor(
    const AtomKokkos* atom,
    const typename AT::tdual_xfloat_2d buf,
    typename AT::tdual_int_1d nlocal,
    int dim, X_FLOAT lo, X_FLOAT hi):
    _x(atom->k_x.view<DeviceType>()),
    _v(atom->k_v.view<DeviceType>()),
    _tag(atom->k_tag.view<DeviceType>()),
    _type(atom->k_type.view<DeviceType>()),
    _mask(atom->k_mask.view<DeviceType>()),
    _image(atom->k_image.view<DeviceType>()),
    _radius(atom->k_radius.view<DeviceType>()),
    _rmass(atom->k_rmass.view<DeviceType>()),
    _omega(atom->k_omega.view<DeviceType>()),
    _forcing(atom->k_forcing.view<DeviceType>()),
    _mean_thickness(atom->k_mean_thickness.view<DeviceType>()),
    _min_thickness(atom->k_min_thickness.view<DeviceType>()),
    _ridgingIceThickness(atom->k_ridgingIceThickness.view<DeviceType>()),
    _ridgingIceThicknessWeight(atom->k_ridgingIceThicknessWeight.view<DeviceType>()),
    _netToGrossClosingRatio(atom->k_netToGrossClosingRatio.view<DeviceType>()),
    _changeEffectiveElementArea(atom->k_changeEffectiveElementArea.view<DeviceType>()),
    _ice_area(atom->k_ice_area.view<DeviceType>()),
    _coriolis(atom->k_coriolis.view<DeviceType>()),
    _ocean_vel(atom->k_ocean_vel.view<DeviceType>()),
    _bvector(atom->k_bvector.view<DeviceType>()),
    _num_bond(atom->k_num_bond.view<DeviceType>()),
    _bond_type(atom->k_bond_type.view<DeviceType>()),
    _bond_atom(atom->k_bond_atom.view<DeviceType>()),
    _nspecial(atom->k_nspecial.view<DeviceType>()),
    _special(atom->k_special.view<DeviceType>()),
    _nlocal(nlocal.template view<DeviceType>()),_dim(dim),
    _lo(lo),_hi(hi){
    elements = 33+atom->maxspecial+2*atom->bond_per_atom;
    const int maxsendlist = (buf.template view<DeviceType>().extent(0)*buf.template view<DeviceType>().extent(1))/elements;

    buffer_view<DeviceType>(_buf,buf,maxsendlist,elements);
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const int &myrecv) const {
    X_FLOAT x = _buf(myrecv,_dim+1);
    if (x >= _lo && x < _hi) {
      int i = Kokkos::atomic_fetch_add(&_nlocal(0),1);
      int m = 1;
      _x(i,0) = _buf(myrecv,m++);
      _x(i,1) = _buf(myrecv,m++);
      _x(i,2) = _buf(myrecv,m++);
      _v(i,0) = _buf(myrecv,m++);
      _v(i,1) = _buf(myrecv,m++);
      _v(i,2) = _buf(myrecv,m++);
      _tag(i) = (tagint) d_ubuf(_buf(myrecv,m++)).i;
      _type(i) = (int) d_ubuf(_buf(myrecv,m++)).i;
      _mask(i) = (int) d_ubuf(_buf(myrecv,m++)).i;
      _image(i) = (imageint) d_ubuf(_buf(myrecv,m++)).i;
      _radius(i) = _buf(myrecv,m++);
      _rmass(i) = _buf(myrecv,m++);
      _omega(i,0) = _buf(myrecv,m++);
      _omega(i,1) = _buf(myrecv,m++);
      _omega(i,2) = _buf(myrecv,m++);
      _forcing(i,0) = _buf(myrecv,m++);
      _forcing(i,1) = _buf(myrecv,m++);
      _mean_thickness(i) = _buf(myrecv,m++);
      _min_thickness(i) = _buf(myrecv,m++);
      _ridgingIceThickness(i) = _buf(myrecv,m++);
      _ridgingIceThicknessWeight(i) = _buf(myrecv,m++);
      _netToGrossClosingRatio(i) = _buf(myrecv,m++);
      _changeEffectiveElementArea(i) = _buf(myrecv,m++);
      _ice_area(i) = _buf(myrecv,m++);
      _coriolis(i) = _buf(myrecv,m++);
      _ocean_vel(i,0) = _buf(myrecv,m++);
      _ocean_vel(i,1) = _buf(myrecv,m++);
      _bvector(i,0) = _buf(myrecv,m++);
      _bvector(i,1) = _buf(myrecv,m++);
      _num_bond(i) = (int) d_ubuf(_buf(myrecv,m++)).i;
      int k;
      for (k = 0; k < _num_bond(i); k++) {
        _bond_type(i,k) = (int) d_ubuf(_buf(myrecv,m++)).i;
        _bond_atom(i,k) = (tagint) d_ubuf(_buf(myrecv,m++)).i;
      }
      _nspecial(i,0) = (int) d_ubuf(_buf(myrecv,m++)).i;
      _nspecial(i,1) = (int) d_ubuf(_buf(myrecv,m++)).i;
      _nspecial(i,2) = (int) d_ubuf(_buf(myrecv,m++)).i;
      for (k = 0; k < _nspecial(i,2); k++)
        _special(i,k) = (tagint) d_ubuf(_buf(myrecv,m++)).i;
    }
  }
};

/* ---------------------------------------------------------------------- */

int AtomVecDemsiKokkos::unpack_exchange_kokkos(DAT::tdual_xfloat_2d &k_buf,int nrecv,
                                               int nlocal,int dim,X_FLOAT lo,X_FLOAT hi,
                                               ExecutionSpace space) {


  const size_t elements = 33+atom->maxspecial+2*atom->bond_per_atom;
  if(space == Host) {
    k_count.h_view(0) = nlocal;
    AtomVecDemsiKokkos_UnpackExchangeFunctor<LMPHostType> f(atomKK,k_buf,k_count,dim,lo,hi);
    Kokkos::parallel_for(nrecv/elements,f); 
  } else {
    k_count.h_view(0) = nlocal;
    k_count.modify<LMPHostType>();
    k_count.sync<LMPDeviceType>();
    AtomVecDemsiKokkos_UnpackExchangeFunctor<LMPDeviceType> f(atomKK,k_buf,k_count,dim,lo,hi);
    Kokkos::parallel_for(nrecv/elements,f);
    k_count.modify<LMPDeviceType>();
    k_count.sync<LMPHostType>();
  }

//  Need this?
  modified(space,X_MASK | V_MASK | TAG_MASK | TYPE_MASK |
                 MASK_MASK | IMAGE_MASK| RADIUS_MASK | RMASS_MASK |
                 OMEGA_MASK | THICKNESS_MASK | FORCING_MASK |
                 BOND_MASK | SPECIAL_MASK);


  return k_count.h_view(0);
}

/* ---------------------------------------------------------------------- */

int AtomVecDemsiKokkos::unpack_exchange(double *buf)
{
 

  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);

  int m = 1;
  h_x(nlocal,0) = buf[m++];
  h_x(nlocal,1) = buf[m++];
  h_x(nlocal,2) = buf[m++];
  h_v(nlocal,0) = buf[m++];
  h_v(nlocal,1) = buf[m++];
  h_v(nlocal,2) = buf[m++];
  h_tag(nlocal) = (tagint) ubuf(buf[m++]).i;
  h_type(nlocal) = (int) ubuf(buf[m++]).i;
  h_mask(nlocal) = (int) ubuf(buf[m++]).i;
  h_image(nlocal) = (imageint) ubuf(buf[m++]).i;

  h_radius(nlocal) = buf[m++];
  h_rmass(nlocal) = buf[m++];
  h_omega(nlocal,0) = buf[m++];
  h_omega(nlocal,1) = buf[m++];
  h_omega(nlocal,2) = buf[m++];

  h_forcing(nlocal,0) = buf[m++];
  h_forcing(nlocal,1) = buf[m++];
  h_mean_thickness(nlocal) = buf[m++];
  h_min_thickness(nlocal) = buf[m++];
  h_ridgingIceThickness(nlocal) = buf[m++];
  h_ridgingIceThicknessWeight(nlocal) = buf[m++];
  h_netToGrossClosingRatio(nlocal) = buf[m++];
  h_changeEffectiveElementArea(nlocal) = buf[m++];
  h_ice_area(nlocal) = buf[m++];
  h_coriolis(nlocal) = buf[m++];
  h_ocean_vel(nlocal,0) = buf[m++];
  h_ocean_vel(nlocal,1) = buf[m++];
  h_bvector(nlocal,0) = buf[m++];
  h_bvector(nlocal,1) = buf[m++];

  h_num_bond(nlocal) = (int) ubuf(buf[m++]).i;
  for (int k = 0; k < h_num_bond(nlocal); k++) {
    h_bond_type(nlocal,k) = (int) ubuf(buf[m++]).i;
    h_bond_atom(nlocal,k) = (tagint) ubuf(buf[m++]).i;
  }

  h_nspecial(nlocal,0) = (int) ubuf(buf[m++]).i;
  h_nspecial(nlocal,1) = (int) ubuf(buf[m++]).i;
  h_nspecial(nlocal,2) = (int) ubuf(buf[m++]).i;
  for (int k = 0; k < h_nspecial(nlocal,2); k++)
     h_special(nlocal,k) = (tagint) ubuf(buf[m++]).i;

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      m += modify->fix[atom->extra_grow[iextra]]->
        unpack_exchange(nlocal,&buf[m]);

  modified(Host,X_MASK | V_MASK | TAG_MASK | TYPE_MASK |
           MASK_MASK | IMAGE_MASK | RADIUS_MASK | RMASS_MASK |
           OMEGA_MASK | FORCING_MASK | THICKNESS_MASK | 
           BOND_MASK | SPECIAL_MASK);

  atom->nlocal++;
  return m;
}

/* ----------------------------------------------------------------------
   size of restart data for all atoms owned by this proc
   include extra data stored by fixes
------------------------------------------------------------------------- */

int AtomVecDemsiKokkos::size_restart()
{
  int i;

  int nlocal = atom->nlocal;
  int n = 0;
  for (int i = 0; i < nlocal; i++)
     n += 20 + 2*num_bond[i];

  if (atom->nextra_restart)
    for (int iextra = 0; iextra < atom->nextra_restart; iextra++)
      for (i = 0; i < nlocal; i++)
        n += modify->fix[atom->extra_restart[iextra]]->size_restart(i);

  return n;
}

/* ----------------------------------------------------------------------
   pack atom I's data for restart file including extra quantities
   xyz must be 1st 3 values, so that read_restart can test on them
   molecular types may be negative, but write as positive
------------------------------------------------------------------------- */

int AtomVecDemsiKokkos::pack_restart(int i, double *buf)
{
  sync(Host,X_MASK | TAG_MASK | TYPE_MASK |
            MASK_MASK | IMAGE_MASK | V_MASK |
            RADIUS_MASK | RMASS_MASK | OMEGA_MASK |
            FORCING_MASK | THICKNESS_MASK | BOND_MASK);

  int m = 1;
  buf[m++] = h_x(i,0);
  buf[m++] = h_x(i,1);
  buf[m++] = h_x(i,2);
  buf[m++] = ubuf(h_tag(i)).d;
  buf[m++] = ubuf(h_type(i)).d;
  buf[m++] = ubuf(h_mask(i)).d;
  buf[m++] = ubuf(h_image(i)).d;
  buf[m++] = h_v(i,0);
  buf[m++] = h_v(i,1);
  buf[m++] = h_v(i,2);

  buf[m++] = h_radius(i);
  buf[m++] = h_rmass(i);
  buf[m++] = h_omega(i,0);
  buf[m++] = h_omega(i,1);
  buf[m++] = h_omega(i,2);

  buf[m++] = h_forcing(i,0);
  buf[m++] = h_forcing(i,1);
  buf[m++] = h_mean_thickness(i);
  buf[m++] = h_min_thickness(i);
  buf[m++] = h_ridgingIceThickness(i);
  buf[m++] = h_ridgingIceThicknessWeight(i);
  buf[m++] = h_netToGrossClosingRatio(i);
  buf[m++] = h_changeEffectiveElementArea(i);
  buf[m++] = h_ice_area(i);
  buf[m++] = h_coriolis(i);
  buf[m++] = h_ocean_vel(i,0);
  buf[m++] = h_ocean_vel(i,1);
  buf[m++] = h_bvector(i,0);
  buf[m++] = h_bvector(i,1);

  buf[m++] = ubuf(h_num_bond(i)).d;
  for (int k = 0; k < h_num_bond(i); k++) {
    buf[m++] = ubuf(MAX(h_bond_type(i,k),-h_bond_type(i,k))).d;
    buf[m++] = ubuf(h_bond_atom(i,k)).d;
  }

  if (atom->nextra_restart)
    for (int iextra = 0; iextra < atom->nextra_restart; iextra++)
      m += modify->fix[atom->extra_restart[iextra]]->pack_restart(i,&buf[m]);

  buf[0] = m;
  return m;
}

/* ----------------------------------------------------------------------
   unpack data for one atom from restart file including extra quantities
------------------------------------------------------------------------- */

int AtomVecDemsiKokkos::unpack_restart(double *buf)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) {
    grow(0);
    if (atom->nextra_store)
      memory->grow(atom->extra,nmax,atom->nextra_store,"atom:extra");
  }

  modified(Host,X_MASK | TAG_MASK | TYPE_MASK |
                MASK_MASK | IMAGE_MASK | V_MASK |
	        RADIUS_MASK | RMASS_MASK | OMEGA_MASK |
                FORCING_MASK | THICKNESS_MASK |
                BOND_MASK | SPECIAL_MASK);

  int m = 1;
  h_x(nlocal,0) = buf[m++];
  h_x(nlocal,1) = buf[m++];
  h_x(nlocal,2) = buf[m++];
  h_tag(nlocal) = (tagint) ubuf(buf[m++]).i;
  h_type(nlocal) = (int) ubuf(buf[m++]).i;
  h_mask(nlocal) = (int) ubuf(buf[m++]).i;
  h_image(nlocal) = (imageint) ubuf(buf[m++]).i;
  h_v(nlocal,0) = buf[m++];
  h_v(nlocal,1) = buf[m++];
  h_v(nlocal,2) = buf[m++];

  h_radius(nlocal) = buf[m++];
  h_rmass(nlocal) = buf[m++];
  h_omega(nlocal,0) = buf[m++];
  h_omega(nlocal,1) = buf[m++];
  h_omega(nlocal,2) = buf[m++];

  h_forcing(nlocal,0) = buf[m++];
  h_forcing(nlocal,1) = buf[m++];
  h_mean_thickness(nlocal) = buf[m++];
  h_min_thickness(nlocal) = buf[m++];
  h_ridgingIceThickness(nlocal) = buf[m++];
  h_ridgingIceThicknessWeight(nlocal) = buf[m++];
  h_netToGrossClosingRatio(nlocal) = buf[m++];
  h_changeEffectiveElementArea(nlocal) = buf[m++];
  h_ice_area(nlocal) = buf[m++];
  h_coriolis(nlocal) = buf[m++];
  h_ocean_vel(nlocal,0) = buf[m++];
  h_ocean_vel(nlocal,1) = buf[m++];
  h_bvector(nlocal,0) = buf[m++];
  h_bvector(nlocal,1) = buf[m++];

  h_num_bond(nlocal) = (int) ubuf(buf[m++]).i;
  for (int k = 0; k < h_num_bond(nlocal); k++) {
    h_bond_type(nlocal,k) = (int) ubuf(buf[m++]).i;
    h_bond_atom(nlocal,k) = (tagint) ubuf(buf[m++]).i;
  }

  h_nspecial(nlocal,0) = h_nspecial(nlocal,1) = h_nspecial(nlocal,2) = 0;

  double **extra = atom->extra;
  if (atom->nextra_store) {
    int size = static_cast<int> (buf[0]) - m;
    for (int i = 0; i < size; i++) extra[nlocal][i] = buf[m++];
  }

  atom->nlocal++;
  return m;
}

/* ----------------------------------------------------------------------
   create one atom of itype at coord
   set other values to defaults
------------------------------------------------------------------------- */

void AtomVecDemsiKokkos::create_atom(int itype, double *coord)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) {
    atomKK->modified(Host,ALL_MASK);
    grow(0);
  }
  atomKK->modified(Host,ALL_MASK);

  h_tag[nlocal] = 0;
  h_type[nlocal] = itype;
  h_x(nlocal,0) = coord[0];
  h_x(nlocal,1) = coord[1];
  h_x(nlocal,2) = coord[2];
  h_mask(nlocal) = 1;
  h_image(nlocal) = ((imageint) IMGMAX << IMG2BITS) |
    ((imageint) IMGMAX << IMGBITS) | IMGMAX;
  h_v(nlocal,0) = 0.0;
  h_v(nlocal,1) = 0.0;
  h_v(nlocal,2) = 0.0;

  h_radius(nlocal) = 0.5;
  h_rmass(nlocal) = MY_PI*h_radius(nlocal)*h_radius(nlocal);
  h_omega(nlocal,0) = 0.0;
  h_omega(nlocal,1) = 0.0;
  h_omega(nlocal,2) = 0.0;

  h_forcing(nlocal,0) = 0.0;
  h_forcing(nlocal,1) = 0.0;
  h_mean_thickness(nlocal) = 0.0;
  h_min_thickness(nlocal) = 0.0;
  h_ridgingIceThickness(nlocal) = 0.0;
  h_ridgingIceThicknessWeight(nlocal) = 0.0;
  h_netToGrossClosingRatio(nlocal) = 0.0;
  h_changeEffectiveElementArea(nlocal) = 0.0;
  h_ice_area(nlocal) = 0.0;
  h_coriolis(nlocal) = 0.0;
  h_ocean_vel(nlocal,0) = 0.0;
  h_ocean_vel(nlocal,1) = 0.0;
  h_bvector(nlocal,0) = 0.0;
  h_bvector(nlocal,1) = 0.0;

  h_num_bond[nlocal] = 0;
  h_nspecial(nlocal,0) = h_nspecial(nlocal,1) = h_nspecial(nlocal,2) = 0;

  atom->nlocal++;
}

/* ----------------------------------------------------------------------
   unpack one line from Atoms section of data file
   initialize other atom quantities
------------------------------------------------------------------------- */

void AtomVecDemsiKokkos::data_atom(double *coord, imageint imagetmp, char **values)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);
  atomKK->modified(Host,ALL_MASK);

  h_tag(nlocal) = ATOTAGINT(values[0]);
  h_type(nlocal) = atoi(values[1]);
  if (h_type(nlocal) <= 0 || h_type(nlocal) > atom->ntypes)
    error->one(FLERR,"Invalid atom type in Atoms section of data file");

  h_radius(nlocal) = 0.5 * atof(values[2]);
  if (h_radius(nlocal) < 0.0)
    error->one(FLERR,"Invalid radius in Atoms section of data file");

  double density = atof(values[3]);
  if (density <= 0.0)
    error->one(FLERR,"Invalid density in Atoms section of data file");

  if (h_radius(nlocal) == 0.0) h_rmass(nlocal) = density;
  else
    h_rmass(nlocal) = MY_PI * h_radius(nlocal)*h_radius(nlocal);

  h_x(nlocal,0) = coord[0];
  h_x(nlocal,1) = coord[1];
  h_x(nlocal,2) = coord[2];

  h_image(nlocal) = imagetmp;

  h_mask(nlocal) = 1;
  h_v(nlocal,0) = 0.0;
  h_v(nlocal,1) = 0.0;
  h_v(nlocal,2) = 0.0;
  h_omega(nlocal,0) = 0.0;
  h_omega(nlocal,1) = 0.0;
  h_omega(nlocal,2) = 0.0;

  h_forcing(nlocal,0) = 0.0;
  h_forcing(nlocal,1) = 0.0;
  h_mean_thickness(nlocal) = 0.0;
  h_min_thickness(nlocal) = 0.0;
  h_ridgingIceThickness(nlocal) = 0.0;
  h_ridgingIceThicknessWeight(nlocal) = 0.0;
  h_netToGrossClosingRatio(nlocal) = 0.0;
  h_changeEffectiveElementArea(nlocal) = 0.0;

  h_num_bond(nlocal) = 0.0;

  atom->nlocal++;
}

/* ----------------------------------------------------------------------
   unpack hybrid quantities from one line in Atoms section of data file
   initialize other atom quantities for this sub-style
------------------------------------------------------------------------- */

int AtomVecDemsiKokkos::data_atom_hybrid(int nlocal, char **values)
{
  radius[nlocal] = 0.5 * atof(values[0]);
  if (radius[nlocal] < 0.0)
    error->one(FLERR,"Invalid radius in Atoms section of data file");

  double density = atof(values[1]);
  if (density <= 0.0)
    error->one(FLERR,"Invalid density in Atoms section of data file");

  if (radius[nlocal] == 0.0) rmass[nlocal] = density;
  else
    rmass[nlocal] = MY_PI * radius[nlocal]*radius[nlocal];

  atomKK->modified(Host,RADIUS_MASK|RMASS_MASK);

  return 2;
}

/* ----------------------------------------------------------------------
   unpack one line from Velocities section of data file
------------------------------------------------------------------------- */

void AtomVecDemsiKokkos::data_vel(int m, char **values)
{
  sync(Host,V_MASK|OMEGA_MASK);
  h_v(m,0) = atof(values[0]);
  h_v(m,1) = atof(values[1]);
  h_v(m,2) = atof(values[2]);
  h_omega(m,0) = atof(values[3]);
  h_omega(m,1) = atof(values[4]);
  h_omega(m,2) = atof(values[5]);
  modified(Host,V_MASK|OMEGA_MASK);
}

/* ----------------------------------------------------------------------
   unpack hybrid quantities from one line in Velocities section of data file
------------------------------------------------------------------------- */

int AtomVecDemsiKokkos::data_vel_hybrid(int m, char **values)
{
  sync(Host,OMEGA_MASK);
  h_omega(m,0) = atof(values[0]);
  h_omega(m,1) = atof(values[1]);
  h_omega(m,2) = atof(values[2]);
  modified(Host,OMEGA_MASK);
  return 3;
}

/* ----------------------------------------------------------------------
   pack atom info for data file including 3 image flags
------------------------------------------------------------------------- */

void AtomVecDemsiKokkos::pack_data(double **buf)
{
  atomKK->sync(Host,TAG_MASK|TYPE_MASK|RADIUS_MASK|RMASS_MASK|X_MASK|IMAGE_MASK);

  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    buf[i][0] = ubuf(h_tag[i]).d;
    buf[i][1] = ubuf(h_type[i]).d;
    buf[i][2] = 2.0*h_radius[i];
    if (h_radius[i] == 0.0) buf[i][3] = h_rmass[i];
    else
      buf[i][3] = h_rmass[i] / (MY_PI * h_radius[i]*h_radius[i]);
    buf[i][4] = h_x(i,0);
    buf[i][5] = h_x(i,1);
    buf[i][6] = h_x(i,2);
    buf[i][7] = ubuf((h_image[i] & IMGMASK) - IMGMAX).d;
    buf[i][8] = ubuf((h_image[i] >> IMGBITS & IMGMASK) - IMGMAX).d;
    buf[i][9] = ubuf((h_image[i] >> IMG2BITS) - IMGMAX).d;
  }
}

/* ----------------------------------------------------------------------
   pack hybrid atom info for data file
------------------------------------------------------------------------- */

int AtomVecDemsiKokkos::pack_data_hybrid(int i, double *buf)
{
  atomKK->sync(Host,RADIUS_MASK|RMASS_MASK);

  buf[0] = 2.0*h_radius[i];
  if (h_radius[i] == 0.0) buf[1] = h_rmass[i];
  else buf[1] = h_rmass[i] / (MY_PI * h_radius[i]*h_radius[i]);
  return 2;
}

/* ----------------------------------------------------------------------
   write atom info to data file including 3 image flags
------------------------------------------------------------------------- */

void AtomVecDemsiKokkos::write_data(FILE *fp, int n, double **buf)
{
  for (int i = 0; i < n; i++)
    fprintf(fp,TAGINT_FORMAT
            " %d %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e %d %d %d\n",
            (tagint) ubuf(buf[i][0]).i,(int) ubuf(buf[i][1]).i,
            buf[i][2],buf[i][3],
            buf[i][4],buf[i][5],buf[i][6],
            (int) ubuf(buf[i][7]).i,(int) ubuf(buf[i][8]).i,
            (int) ubuf(buf[i][9]).i);
}

/* ----------------------------------------------------------------------
   write hybrid atom info to data file
------------------------------------------------------------------------- */

int AtomVecDemsiKokkos::write_data_hybrid(FILE *fp, double *buf)
{
  fprintf(fp," %-1.16e %-1.16e",buf[0],buf[1]);
  return 2;
}

/* ----------------------------------------------------------------------
   pack velocity info for data file
------------------------------------------------------------------------- */

void AtomVecDemsiKokkos::pack_vel(double **buf)
{
  sync(Host,TAG_MASK|V_MASK|OMEGA_MASK);

  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    buf[i][0] = ubuf(h_tag[i]).d;
    buf[i][1] = h_v(i,0);
    buf[i][2] = h_v(i,1);
    buf[i][3] = h_v(i,2);
    buf[i][4] = h_omega(i,0);
    buf[i][5] = h_omega(i,1);
    buf[i][6] = h_omega(i,2);
  }
}

/* ----------------------------------------------------------------------
   pack hybrid velocity info for data file
------------------------------------------------------------------------- */

int AtomVecDemsiKokkos::pack_vel_hybrid(int i, double *buf)
{
  sync(Host,OMEGA_MASK);

  buf[0] = h_omega(i,0);
  buf[1] = h_omega(i,1);
  buf[2] = h_omega(i,2);
  return 3;
}

/* ----------------------------------------------------------------------
   write velocity info to data file
------------------------------------------------------------------------- */

void AtomVecDemsiKokkos::write_vel(FILE *fp, int n, double **buf)
{
  for (int i = 0; i < n; i++)
    fprintf(fp,TAGINT_FORMAT
            " %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e\n",
            (tagint) ubuf(buf[i][0]).i,buf[i][1],buf[i][2],buf[i][3],
            buf[i][4],buf[i][5],buf[i][6]);
}

/* ----------------------------------------------------------------------
   write hybrid velocity info to data file
------------------------------------------------------------------------- */

int AtomVecDemsiKokkos::write_vel_hybrid(FILE *fp, double *buf)
{
  fprintf(fp," %-1.16e %-1.16e %-1.16e",buf[0],buf[1],buf[2]);
  return 3;
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
------------------------------------------------------------------------- */

bigint AtomVecDemsiKokkos::memory_usage()
{
  bigint bytes = 0;

  if (atom->memcheck("tag")) bytes += memory->usage(tag,nmax);
  if (atom->memcheck("type")) bytes += memory->usage(type,nmax);
  if (atom->memcheck("mask")) bytes += memory->usage(mask,nmax);
  if (atom->memcheck("image")) bytes += memory->usage(image,nmax);
  if (atom->memcheck("x")) bytes += memory->usage(x,nmax,3);
  if (atom->memcheck("v")) bytes += memory->usage(v,nmax,3);
  if (atom->memcheck("f")) bytes += memory->usage(f,nmax*comm->nthreads,3);

  if (atom->memcheck("radius")) bytes += memory->usage(radius,nmax);
  if (atom->memcheck("rmass")) bytes += memory->usage(rmass,nmax);
  if (atom->memcheck("omega")) bytes += memory->usage(omega,nmax,3);
  if (atom->memcheck("torque"))
    bytes += memory->usage(torque,nmax*comm->nthreads,3);

  if (atom->memcheck("forcing")) bytes += memory->usage(forcing,nmax,2);
  if (atom->memcheck("mean_thickness")) bytes += memory->usage(mean_thickness,nmax);
  if (atom->memcheck("min_thickness")) bytes += memory->usage(min_thickness,nmax);
  if (atom->memcheck("ridgingIceThickness")) bytes += memory->usage(ridgingIceThickness,nmax);
  if (atom->memcheck("ridgingIceThicknessWeight")) bytes += memory->usage(ridgingIceThicknessWeight,nmax);
  if (atom->memcheck("netToGrossClosingRatio")) bytes += memory->usage(netToGrossClosingRatio,nmax);
  if (atom->memcheck("changeEffectiveElementArea")) bytes += memory->usage(changeEffectiveElementArea,nmax);
  if (atom->memcheck("ice_area")) bytes += memory->usage(ice_area,nmax);
  if (atom->memcheck("coriolis")) bytes += memory->usage(coriolis,nmax);
  if (atom->memcheck("ocean_vel")) bytes += memory->usage(ocean_vel,nmax,2);
  if (atom->memcheck("bvector")) bytes += memory->usage(bvector,nmax,2);

  if (atom->memcheck("num_bond")) bytes += memory->usage(num_bond,nmax);
  if (atom->memcheck("bond_type"))
    bytes += memory->usage(bond_type,nmax,atom->bond_per_atom);
  if (atom->memcheck("bond_atom"))
    bytes += memory->usage(bond_atom,nmax,atom->bond_per_atom);

  if (atom->memcheck("nspecial")) bytes += memory->usage(nspecial,nmax,3);
  if (atom->memcheck("special"))
     bytes += memory->usage(special,nmax,atom->maxspecial);

  return bytes;
}

/* ---------------------------------------------------------------------- */

void AtomVecDemsiKokkos::sync(ExecutionSpace space, unsigned int mask)
{
  if (space == Device) {
    if (mask & X_MASK) atomKK->k_x.sync<LMPDeviceType>();
    if (mask & V_MASK) atomKK->k_v.sync<LMPDeviceType>();
    if (mask & F_MASK) atomKK->k_f.sync<LMPDeviceType>();
    if (mask & TAG_MASK) atomKK->k_tag.sync<LMPDeviceType>();
    if (mask & TYPE_MASK) atomKK->k_type.sync<LMPDeviceType>();
    if (mask & MASK_MASK) atomKK->k_mask.sync<LMPDeviceType>();
    if (mask & IMAGE_MASK) atomKK->k_image.sync<LMPDeviceType>();
    if (mask & RADIUS_MASK) atomKK->k_radius.sync<LMPDeviceType>();
    if (mask & RMASS_MASK) atomKK->k_rmass.sync<LMPDeviceType>();
    if (mask & OMEGA_MASK) atomKK->k_omega.sync<LMPDeviceType>();
    if (mask & TORQUE_MASK) atomKK->k_torque.sync<LMPDeviceType>();
    if (mask & FORCING_MASK) atomKK->k_forcing.sync<LMPDeviceType>();
    if (mask & THICKNESS_MASK) {
        atomKK->k_mean_thickness.sync<LMPDeviceType>();
        atomKK->k_min_thickness.sync<LMPDeviceType>();
	atomKK->k_ridgingIceThickness.sync<LMPDeviceType>();
	atomKK->k_ridgingIceThicknessWeight.sync<LMPDeviceType>();
	atomKK->k_netToGrossClosingRatio.sync<LMPDeviceType>();
	atomKK->k_changeEffectiveElementArea.sync<LMPDeviceType>();
        atomKK->k_ice_area.sync<LMPDeviceType>();
        atomKK->k_coriolis.sync<LMPDeviceType>();
        atomKK->k_ocean_vel.sync<LMPDeviceType>();
        atomKK->k_bvector.sync<LMPDeviceType>();
    }
    if (mask & BOND_MASK) {
        atomKK->k_num_bond.sync<LMPDeviceType>();
        atomKK->k_bond_type.sync<LMPDeviceType>();
        atomKK->k_bond_atom.sync<LMPDeviceType>();
    }
    if (mask & SPECIAL_MASK) {
        atomKK->k_nspecial.sync<LMPDeviceType>();
        atomKK->k_special.sync<LMPDeviceType>();
    }
  } else {
    if (mask & X_MASK) atomKK->k_x.sync<LMPHostType>();
    if (mask & V_MASK) atomKK->k_v.sync<LMPHostType>();
    if (mask & F_MASK) atomKK->k_f.sync<LMPHostType>();
    if (mask & TAG_MASK) atomKK->k_tag.sync<LMPHostType>();
    if (mask & TYPE_MASK) atomKK->k_type.sync<LMPHostType>();
    if (mask & MASK_MASK) atomKK->k_mask.sync<LMPHostType>();
    if (mask & IMAGE_MASK) atomKK->k_image.sync<LMPHostType>();
    if (mask & RADIUS_MASK) atomKK->k_radius.sync<LMPHostType>();
    if (mask & RMASS_MASK) atomKK->k_rmass.sync<LMPHostType>();
    if (mask & OMEGA_MASK) atomKK->k_omega.sync<LMPHostType>();
    if (mask & TORQUE_MASK) atomKK->k_torque.sync<LMPHostType>();
    if (mask & FORCING_MASK) atomKK->k_forcing.sync<LMPHostType>();
    if (mask & THICKNESS_MASK) {
        atomKK->k_mean_thickness.sync<LMPHostType>();
        atomKK->k_min_thickness.sync<LMPHostType>();
	atomKK->k_ridgingIceThickness.sync<LMPHostType>();
	atomKK->k_ridgingIceThicknessWeight.sync<LMPHostType>();
	atomKK->k_netToGrossClosingRatio.sync<LMPHostType>();
	atomKK->k_changeEffectiveElementArea.sync<LMPHostType>();
        atomKK->k_ice_area.sync<LMPHostType>();
        atomKK->k_coriolis.sync<LMPHostType>();
        atomKK->k_ocean_vel.sync<LMPHostType>();
        atomKK->k_bvector.sync<LMPHostType>();
    }
    if (mask & BOND_MASK) {
        atomKK->k_num_bond.sync<LMPHostType>();
        atomKK->k_bond_type.sync<LMPHostType>();
        atomKK->k_bond_atom.sync<LMPHostType>();
    }
    if (mask & SPECIAL_MASK) {
        atomKK->k_nspecial.sync<LMPHostType>();
        atomKK->k_special.sync<LMPHostType>();
    }
  }
}

/* ---------------------------------------------------------------------- */

void AtomVecDemsiKokkos::sync_overlapping_device(ExecutionSpace space, unsigned int mask)
{
  if (space == Device) {
    if ((mask & X_MASK) && atomKK->k_x.need_sync<LMPDeviceType>())
      perform_async_copy<DAT::tdual_x_array>(atomKK->k_x,space);
    if ((mask & V_MASK) && atomKK->k_v.need_sync<LMPDeviceType>())
      perform_async_copy<DAT::tdual_v_array>(atomKK->k_v,space);
    if ((mask & F_MASK) && atomKK->k_f.need_sync<LMPDeviceType>())
      perform_async_copy<DAT::tdual_f_array>(atomKK->k_f,space);
    if ((mask & TAG_MASK) && atomKK->k_tag.need_sync<LMPDeviceType>())
      perform_async_copy<DAT::tdual_tagint_1d>(atomKK->k_tag,space);
    if ((mask & TYPE_MASK) && atomKK->k_type.need_sync<LMPDeviceType>())
      perform_async_copy<DAT::tdual_int_1d>(atomKK->k_type,space);
    if ((mask & MASK_MASK) && atomKK->k_mask.need_sync<LMPDeviceType>())
      perform_async_copy<DAT::tdual_int_1d>(atomKK->k_mask,space);
    if ((mask & IMAGE_MASK) && atomKK->k_image.need_sync<LMPDeviceType>())
      perform_async_copy<DAT::tdual_imageint_1d>(atomKK->k_image,space);
    if ((mask & RADIUS_MASK) && atomKK->k_radius.need_sync<LMPDeviceType>())
      perform_async_copy<DAT::tdual_float_1d>(atomKK->k_radius,space);
    if ((mask & RMASS_MASK) && atomKK->k_rmass.need_sync<LMPDeviceType>())
      perform_async_copy<DAT::tdual_float_1d>(atomKK->k_rmass,space);
    if ((mask & OMEGA_MASK) && atomKK->k_omega.need_sync<LMPDeviceType>())
      perform_async_copy<DAT::tdual_v_array>(atomKK->k_omega,space);
    if ((mask & TORQUE_MASK) && atomKK->k_torque.need_sync<LMPDeviceType>())
      perform_async_copy<DAT::tdual_f_array>(atomKK->k_torque,space);
    if ((mask & FORCING_MASK) && atomKK->k_forcing.need_sync<LMPDeviceType>())
      perform_async_copy<DAT::tdual_float_2d>(atomKK->k_forcing,space);
    if (mask & THICKNESS_MASK) {
        if (atomKK->k_mean_thickness.need_sync<LMPDeviceType>())
           perform_async_copy<DAT::tdual_float_1d>(atomKK->k_mean_thickness,space);
        if (atomKK->k_min_thickness.need_sync<LMPDeviceType>())
           perform_async_copy<DAT::tdual_float_1d>(atomKK->k_min_thickness,space);
	if (atomKK->k_ridgingIceThickness.need_sync<LMPDeviceType>())
	   perform_async_copy<DAT::tdual_float_1d>(atomKK->k_ridgingIceThickness,space);
	if (atomKK->k_ridgingIceThicknessWeight.need_sync<LMPDeviceType>())
	   perform_async_copy<DAT::tdual_float_1d>(atomKK->k_ridgingIceThicknessWeight,space);
	if (atomKK->k_netToGrossClosingRatio.need_sync<LMPDeviceType>())
	   perform_async_copy<DAT::tdual_float_1d>(atomKK->k_netToGrossClosingRatio,space);
	if (atomKK->k_changeEffectiveElementArea.need_sync<LMPDeviceType>())
	   perform_async_copy<DAT::tdual_float_1d>(atomKK->k_changeEffectiveElementArea,space);
        if (atomKK->k_ice_area.need_sync<LMPDeviceType>())
           perform_async_copy<DAT::tdual_float_1d>(atomKK->k_ice_area,space);
        if (atomKK->k_coriolis.need_sync<LMPDeviceType>())
           perform_async_copy<DAT::tdual_float_1d>(atomKK->k_coriolis,space);
        if (atomKK->k_ocean_vel.need_sync<LMPDeviceType>())
           perform_async_copy<DAT::tdual_float_2d>(atomKK->k_ocean_vel,space);
        if (atomKK->k_bvector.need_sync<LMPDeviceType>())
           perform_async_copy<DAT::tdual_float_2d>(atomKK->k_bvector,space);
    }
    if (mask & SPECIAL_MASK) {
      if (atomKK->k_nspecial.need_sync<LMPDeviceType>())
        perform_async_copy<DAT::tdual_int_2d>(atomKK->k_nspecial,space);
      if (atomKK->k_special.need_sync<LMPDeviceType>())
        perform_async_copy<DAT::tdual_tagint_2d>(atomKK->k_special,space);
    }
    if (mask & BOND_MASK) {
      if (atomKK->k_num_bond.need_sync<LMPDeviceType>())
        perform_async_copy<DAT::tdual_int_1d>(atomKK->k_num_bond,space);
      if (atomKK->k_bond_type.need_sync<LMPDeviceType>())
        perform_async_copy<DAT::tdual_int_2d>(atomKK->k_bond_type,space);
      if (atomKK->k_bond_atom.need_sync<LMPDeviceType>())
        perform_async_copy<DAT::tdual_tagint_2d>(atomKK->k_bond_atom,space);
    }
  } else {
    if ((mask & X_MASK) && atomKK->k_x.need_sync<LMPHostType>())
      perform_async_copy<DAT::tdual_x_array>(atomKK->k_x,space);
    if ((mask & V_MASK) && atomKK->k_v.need_sync<LMPHostType>())
      perform_async_copy<DAT::tdual_v_array>(atomKK->k_v,space);
    if ((mask & F_MASK) && atomKK->k_f.need_sync<LMPHostType>())
      perform_async_copy<DAT::tdual_f_array>(atomKK->k_f,space);
    if ((mask & TAG_MASK) && atomKK->k_tag.need_sync<LMPHostType>())
      perform_async_copy<DAT::tdual_tagint_1d>(atomKK->k_tag,space);
    if ((mask & TYPE_MASK) && atomKK->k_type.need_sync<LMPHostType>())
      perform_async_copy<DAT::tdual_int_1d>(atomKK->k_type,space);
    if ((mask & MASK_MASK) && atomKK->k_mask.need_sync<LMPHostType>())
      perform_async_copy<DAT::tdual_int_1d>(atomKK->k_mask,space);
    if ((mask & IMAGE_MASK) && atomKK->k_image.need_sync<LMPHostType>())
      perform_async_copy<DAT::tdual_imageint_1d>(atomKK->k_image,space);
    if ((mask & RADIUS_MASK) && atomKK->k_radius.need_sync<LMPHostType>())
      perform_async_copy<DAT::tdual_float_1d>(atomKK->k_radius,space);
    if ((mask & RMASS_MASK) && atomKK->k_rmass.need_sync<LMPHostType>())
      perform_async_copy<DAT::tdual_float_1d>(atomKK->k_rmass,space);
    if ((mask & OMEGA_MASK) && atomKK->k_omega.need_sync<LMPHostType>())
      perform_async_copy<DAT::tdual_v_array>(atomKK->k_omega,space);
    if ((mask & TORQUE_MASK) && atomKK->k_torque.need_sync<LMPHostType>())
      perform_async_copy<DAT::tdual_f_array>(atomKK->k_torque,space);
    if ((mask & FORCING_MASK) && atomKK->k_forcing.need_sync<LMPHostType>())
      perform_async_copy<DAT::tdual_float_2d>(atomKK->k_forcing,space);
    if (mask & THICKNESS_MASK) {
        if (atomKK->k_mean_thickness.need_sync<LMPHostType>())
           perform_async_copy<DAT::tdual_float_1d>(atomKK->k_mean_thickness,space);
        if (atomKK->k_min_thickness.need_sync<LMPHostType>())
           perform_async_copy<DAT::tdual_float_1d>(atomKK->k_min_thickness,space);
	if (atomKK->k_ridgingIceThickness.need_sync<LMPHostType>())
	   perform_async_copy<DAT::tdual_float_1d>(atomKK->k_ridgingIceThickness,space);
	if (atomKK->k_ridgingIceThicknessWeight.need_sync<LMPHostType>())
	   perform_async_copy<DAT::tdual_float_1d>(atomKK->k_ridgingIceThicknessWeight,space);
	if (atomKK->k_netToGrossClosingRatio.need_sync<LMPHostType>())
	   perform_async_copy<DAT::tdual_float_1d>(atomKK->k_netToGrossClosingRatio,space);
	if (atomKK->k_changeEffectiveElementArea.need_sync<LMPHostType>())
	   perform_async_copy<DAT::tdual_float_1d>(atomKK->k_changeEffectiveElementArea,space);
        if (atomKK->k_ice_area.need_sync<LMPHostType>())
           perform_async_copy<DAT::tdual_float_1d>(atomKK->k_ice_area,space);
        if (atomKK->k_coriolis.need_sync<LMPHostType>())
           perform_async_copy<DAT::tdual_float_1d>(atomKK->k_coriolis,space);
        if (atomKK->k_ocean_vel.need_sync<LMPHostType>())
           perform_async_copy<DAT::tdual_float_2d>(atomKK->k_ocean_vel,space);
        if (atomKK->k_bvector.need_sync<LMPHostType>())
           perform_async_copy<DAT::tdual_float_2d>(atomKK->k_bvector,space);
    }
    if (mask & SPECIAL_MASK) {
      if (atomKK->k_nspecial.need_sync<LMPHostType>())
        perform_async_copy<DAT::tdual_int_2d>(atomKK->k_nspecial,space);
      if (atomKK->k_special.need_sync<LMPHostType>())
        perform_async_copy<DAT::tdual_tagint_2d>(atomKK->k_special,space);
    }
    if (mask & BOND_MASK) {
      if (atomKK->k_num_bond.need_sync<LMPHostType>())
        perform_async_copy<DAT::tdual_int_1d>(atomKK->k_num_bond,space);
      if (atomKK->k_bond_type.need_sync<LMPHostType>())
        perform_async_copy<DAT::tdual_int_2d>(atomKK->k_bond_type,space);
      if (atomKK->k_bond_atom.need_sync<LMPHostType>())
        perform_async_copy<DAT::tdual_tagint_2d>(atomKK->k_bond_atom,space);
    }
  }
}

/* ---------------------------------------------------------------------- */

void AtomVecDemsiKokkos::modified(ExecutionSpace space, unsigned int mask)
{
  if (space == Device) {
    if (mask & X_MASK) atomKK->k_x.modify<LMPDeviceType>();
    if (mask & V_MASK) atomKK->k_v.modify<LMPDeviceType>();
    if (mask & F_MASK) atomKK->k_f.modify<LMPDeviceType>();
    if (mask & TAG_MASK) atomKK->k_tag.modify<LMPDeviceType>();
    if (mask & TYPE_MASK) atomKK->k_type.modify<LMPDeviceType>();
    if (mask & MASK_MASK) atomKK->k_mask.modify<LMPDeviceType>();
    if (mask & IMAGE_MASK) atomKK->k_image.modify<LMPDeviceType>();
    if (mask & RADIUS_MASK) atomKK->k_radius.modify<LMPDeviceType>();
    if (mask & RMASS_MASK) atomKK->k_rmass.modify<LMPDeviceType>();
    if (mask & OMEGA_MASK) atomKK->k_omega.modify<LMPDeviceType>();
    if (mask & TORQUE_MASK) atomKK->k_torque.modify<LMPDeviceType>();
    if (mask & FORCING_MASK) atomKK->k_forcing.modify<LMPDeviceType>();
    if (mask & THICKNESS_MASK) {
        atomKK->k_mean_thickness.modify<LMPDeviceType>();
        atomKK->k_min_thickness.modify<LMPDeviceType>();
	atomKK->k_ridgingIceThickness.modify<LMPDeviceType>();
	atomKK->k_ridgingIceThicknessWeight.modify<LMPDeviceType>();
	atomKK->k_netToGrossClosingRatio.modify<LMPDeviceType>();
	atomKK->k_changeEffectiveElementArea.modify<LMPDeviceType>();
	atomKK->k_ice_area.modify<LMPDeviceType>();
        atomKK->k_coriolis.modify<LMPDeviceType>();
        atomKK->k_ocean_vel.modify<LMPDeviceType>();
        atomKK->k_bvector.modify<LMPDeviceType>();
    }
    if (mask & BOND_MASK) {
        atomKK->k_num_bond.modify<LMPDeviceType>();
        atomKK->k_bond_type.modify<LMPDeviceType>();
        atomKK->k_bond_atom.modify<LMPDeviceType>();
    }
    if (mask & SPECIAL_MASK) {
        atomKK->k_nspecial.modify<LMPDeviceType>();
        atomKK->k_special.modify<LMPDeviceType>();
    }
  } else {
    if (mask & X_MASK) atomKK->k_x.modify<LMPHostType>();
    if (mask & V_MASK) atomKK->k_v.modify<LMPHostType>();
    if (mask & F_MASK) atomKK->k_f.modify<LMPHostType>();
    if (mask & TAG_MASK) atomKK->k_tag.modify<LMPHostType>();
    if (mask & TYPE_MASK) atomKK->k_type.modify<LMPHostType>();
    if (mask & MASK_MASK) atomKK->k_mask.modify<LMPHostType>();
    if (mask & IMAGE_MASK) atomKK->k_image.modify<LMPHostType>();
    if (mask & RADIUS_MASK) atomKK->k_radius.modify<LMPHostType>();
    if (mask & RMASS_MASK) atomKK->k_rmass.modify<LMPHostType>();
    if (mask & OMEGA_MASK) atomKK->k_omega.modify<LMPHostType>();
    if (mask & TORQUE_MASK) atomKK->k_torque.modify<LMPHostType>();
    if (mask & FORCING_MASK) atomKK->k_forcing.modify<LMPHostType>();
    if (mask & THICKNESS_MASK) {
        atomKK->k_mean_thickness.modify<LMPHostType>();
        atomKK->k_min_thickness.modify<LMPHostType>();
	atomKK->k_ridgingIceThickness.modify<LMPHostType>();
	atomKK->k_ridgingIceThicknessWeight.modify<LMPHostType>();
	atomKK->k_netToGrossClosingRatio.modify<LMPHostType>();
	atomKK->k_changeEffectiveElementArea.modify<LMPHostType>();
        atomKK->k_ice_area.modify<LMPHostType>();
        atomKK->k_coriolis.modify<LMPHostType>();
        atomKK->k_ocean_vel.modify<LMPHostType>();
        atomKK->k_bvector.modify<LMPHostType>();
    }
    if (mask & BOND_MASK) {
        atomKK->k_num_bond.modify<LMPHostType>();
        atomKK->k_bond_type.modify<LMPHostType>();
        atomKK->k_bond_atom.modify<LMPHostType>();
    }
    if (mask & SPECIAL_MASK) {
        atomKK->k_nspecial.modify<LMPHostType>();
        atomKK->k_special.modify<LMPHostType>();
    }
  }
}
