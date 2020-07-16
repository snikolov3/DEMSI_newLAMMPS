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

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "atom_vec_demsi.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "modify.h"
#include "force.h"
#include "fix.h"
#include "fix_adapt.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

AtomVecDemsi::AtomVecDemsi(LAMMPS *lmp) : AtomVec(lmp)
{
molecular = 1;

  comm_x_only = 1;
  comm_f_only = 0;
size_forward = 3;
size_reverse = 6;
size_border = 18;
  size_velocity = 6;
  size_data_atom = 7;
  size_data_vel = 7;
  xcol_data = 5;

  atom->sphere_flag = 1;
atom->demsi_flag = 1;
  atom->radius_flag = atom->rmass_flag = atom->omega_flag =
    atom->torque_flag = 1;
bonds_allow = 1;
}

/* ---------------------------------------------------------------------- */

void AtomVecDemsi::init()
{
  AtomVec::init();

  // set radvary if particle diameters are time-varying due to fix adapt

  radvary = 0;
comm_x_only = 0;
size_forward = 3;

  for (int i = 0; i < modify->nfix; i++)
    if (strcmp(modify->fix[i]->style,"adapt") == 0) {
      FixAdapt *fix = (FixAdapt *) modify->fix[i];
      if (fix->diamflag) {
        radvary = 1;
        comm_x_only = 0;
size_forward = 3;
      }
    }
}

/* ----------------------------------------------------------------------
   grow atom arrays
   n = 0 grows arrays by a chunk
   n > 0 allocates arrays to size n
------------------------------------------------------------------------- */

void AtomVecDemsi::grow(int n)
{
  if (n == 0) grow_nmax();
  else nmax = n;
  atom->nmax = nmax;
  if (nmax < 0 || nmax > MAXSMALLINT)
    error->one(FLERR,"Per-processor system is too big");

  tag = memory->grow(atom->tag,nmax,"atom:tag");
  type = memory->grow(atom->type,nmax,"atom:type");
  mask = memory->grow(atom->mask,nmax,"atom:mask");
  image = memory->grow(atom->image,nmax,"atom:image");
  x = memory->grow(atom->x,nmax,3,"atom:x");
  v = memory->grow(atom->v,nmax,3,"atom:v");
  f = memory->grow(atom->f,nmax*comm->nthreads,3,"atom:f");

  radius = memory->grow(atom->radius,nmax,"atom:radius");
  rmass = memory->grow(atom->rmass,nmax,"atom:rmass");
  omega = memory->grow(atom->omega,nmax,3,"atom:omega");
  torque = memory->grow(atom->torque,nmax*comm->nthreads,3,"atom:torque");

  forcing = memory->grow(atom->forcing,nmax,2,"atom:forcing");
  mean_thickness = memory->grow(atom->mean_thickness,nmax,"atom:mean_thickness");
  min_thickness = memory->grow(atom->min_thickness,nmax,"atom:min_thickness");
  ice_area = memory->grow(atom->ice_area,nmax,"atom:ice_area");
  coriolis = memory->grow(atom->coriolis,nmax,"atom:coriolis");
  ocean_vel = memory->grow(atom->ocean_vel,nmax,2,"atom:ocean_vel");
  bvector = memory->grow(atom->bvector,nmax,2,"atom:bvector");

  nspecial = memory->grow(atom->nspecial,nmax,3,"atom:nspecial");
  special = memory->grow(atom->special,nmax,atom->maxspecial,"atom:special");
  num_bond = memory->grow(atom->num_bond,nmax,"atom:num_bond");
  bond_type = memory->grow(atom->bond_type,nmax,atom->bond_per_atom,"atom:bond_type");
  bond_atom = memory->grow(atom->bond_atom,nmax,atom->bond_per_atom,"atom:bond_atom");

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      modify->fix[atom->extra_grow[iextra]]->grow_arrays(nmax);
}

/* ----------------------------------------------------------------------
   reset local array ptrs
------------------------------------------------------------------------- */

void AtomVecDemsi::grow_reset()
{
  tag = atom->tag; type = atom->type;
  mask = atom->mask; image = atom->image;
  x = atom->x; v = atom->v; f = atom->f;
  radius = atom->radius; rmass = atom->rmass;
  omega = atom->omega; torque = atom->torque;

  forcing = atom->forcing;
  mean_thickness = atom->mean_thickness;
  min_thickness = atom->min_thickness;
  ice_area = atom->ice_area;
  coriolis = atom->coriolis;
  ocean_vel = atom->ocean_vel;
  bvector = atom->bvector;

  nspecial = atom->nspecial; special = atom->special;
  num_bond = atom->num_bond; bond_type = atom->bond_type;
  bond_atom = atom->bond_atom;
}

/* ----------------------------------------------------------------------
   copy atom I info to atom J
------------------------------------------------------------------------- */

void AtomVecDemsi::copy(int i, int j, int delflag)
{
  tag[j] = tag[i];
  type[j] = type[i];
  mask[j] = mask[i];
  image[j] = image[i];
  x[j][0] = x[i][0];
  x[j][1] = x[i][1];
  x[j][2] = x[i][2];
  v[j][0] = v[i][0];
  v[j][1] = v[i][1];
  v[j][2] = v[i][2];

  radius[j] = radius[i];
  rmass[j] = rmass[i];
  omega[j][0] = omega[i][0];
  omega[j][1] = omega[i][1];
  omega[j][2] = omega[i][2];

  forcing[j][0] = forcing[i][0];
  forcing[j][1] = forcing[i][1];
  mean_thickness[j] = mean_thickness[i];
  min_thickness[j] = min_thickness[i];
  ice_area[j] = ice_area[i];
  coriolis[j] = coriolis[i];
  ocean_vel[j][0] = ocean_vel[i][0];
  ocean_vel[j][1] = ocean_vel[i][1];
  bvector[j][0] = bvector[i][0];
  bvector[j][1] = bvector[i][1];
  nspecial[j][0] = nspecial[i][0];
  nspecial[j][1] = nspecial[i][1];
  nspecial[j][2] = nspecial[i][2];
  for (int k = 0; k < nspecial[j][2]; k++) special[j][k] = special[i][k];
    num_bond[j] = num_bond[i];
  for (int k = 0; k < num_bond[j]; k++) {
    bond_type[j][k] = bond_type[i][k];
    bond_atom[j][k] = bond_atom[i][k];
}

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      modify->fix[atom->extra_grow[iextra]]->copy_arrays(i,j,delflag);
}

/* ---------------------------------------------------------------------- */

int AtomVecDemsi::pack_comm(int n, int *list, double *buf,
                               int pbc_flag, int *pbc)
{
  int i,j,m;
  double dx,dy,dz;

  if (radvary == 0) {
    m = 0;
    if (pbc_flag == 0) {
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = x[j][0];
        buf[m++] = x[j][1];
        buf[m++] = x[j][2];
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
        buf[m++] = x[j][0] + dx;
        buf[m++] = x[j][1] + dy;
        buf[m++] = x[j][2] + dz;
      }
    }

  } else {
    m = 0;
    if (pbc_flag == 0) {
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = x[j][0];
        buf[m++] = x[j][1];
        buf[m++] = x[j][2];
        buf[m++] = radius[j];
        buf[m++] = rmass[j];
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
        buf[m++] = x[j][0] + dx;
        buf[m++] = x[j][1] + dy;
        buf[m++] = x[j][2] + dz;
        buf[m++] = radius[j];
        buf[m++] = rmass[j];
      }
    }
  }

  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecDemsi::pack_comm_vel(int n, int *list, double *buf,
                                   int pbc_flag, int *pbc)
{
  int i,j,m;
  double dx,dy,dz,dvx,dvy,dvz;

  if (radvary == 0) {
    m = 0;
    if (pbc_flag == 0) {
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = x[j][0];
        buf[m++] = x[j][1];
        buf[m++] = x[j][2];
        buf[m++] = v[j][0];
        buf[m++] = v[j][1];
        buf[m++] = v[j][2];
        buf[m++] = omega[j][0];
        buf[m++] = omega[j][1];
        buf[m++] = omega[j][2];
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
          buf[m++] = x[j][0] + dx;
          buf[m++] = x[j][1] + dy;
          buf[m++] = x[j][2] + dz;
          buf[m++] = v[j][0];
          buf[m++] = v[j][1];
          buf[m++] = v[j][2];
          buf[m++] = omega[j][0];
          buf[m++] = omega[j][1];
          buf[m++] = omega[j][2];
        }
      } else {
        dvx = pbc[0]*h_rate[0] + pbc[5]*h_rate[5] + pbc[4]*h_rate[4];
        dvy = pbc[1]*h_rate[1] + pbc[3]*h_rate[3];
        dvz = pbc[2]*h_rate[2];
        for (i = 0; i < n; i++) {
          j = list[i];
          buf[m++] = x[j][0] + dx;
          buf[m++] = x[j][1] + dy;
          buf[m++] = x[j][2] + dz;
          if (mask[i] & deform_groupbit) {
            buf[m++] = v[j][0] + dvx;
            buf[m++] = v[j][1] + dvy;
            buf[m++] = v[j][2] + dvz;
          } else {
            buf[m++] = v[j][0];
            buf[m++] = v[j][1];
            buf[m++] = v[j][2];
          }
          buf[m++] = omega[j][0];
          buf[m++] = omega[j][1];
          buf[m++] = omega[j][2];
        }
      }
    }

  } else {
    m = 0;
    if (pbc_flag == 0) {
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = x[j][0];
        buf[m++] = x[j][1];
        buf[m++] = x[j][2];
        buf[m++] = radius[j];
        buf[m++] = rmass[j];
        buf[m++] = v[j][0];
        buf[m++] = v[j][1];
        buf[m++] = v[j][2];
        buf[m++] = omega[j][0];
        buf[m++] = omega[j][1];
        buf[m++] = omega[j][2];
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
          buf[m++] = x[j][0] + dx;
          buf[m++] = x[j][1] + dy;
          buf[m++] = x[j][2] + dz;
          buf[m++] = radius[j];
          buf[m++] = rmass[j];
          buf[m++] = v[j][0];
          buf[m++] = v[j][1];
          buf[m++] = v[j][2];
          buf[m++] = omega[j][0];
          buf[m++] = omega[j][1];
          buf[m++] = omega[j][2];
        }
      } else {
        dvx = pbc[0]*h_rate[0] + pbc[5]*h_rate[5] + pbc[4]*h_rate[4];
        dvy = pbc[1]*h_rate[1] + pbc[3]*h_rate[3];
        dvz = pbc[2]*h_rate[2];
        for (i = 0; i < n; i++) {
          j = list[i];
          buf[m++] = x[j][0] + dx;
          buf[m++] = x[j][1] + dy;
          buf[m++] = x[j][2] + dz;
          buf[m++] = radius[j];
          buf[m++] = rmass[j];
          if (mask[i] & deform_groupbit) {
            buf[m++] = v[j][0] + dvx;
            buf[m++] = v[j][1] + dvy;
            buf[m++] = v[j][2] + dvz;
          } else {
            buf[m++] = v[j][0];
            buf[m++] = v[j][1];
            buf[m++] = v[j][2];
          }
          buf[m++] = omega[j][0];
          buf[m++] = omega[j][1];
          buf[m++] = omega[j][2];
        }
      }
    }
  }

  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecDemsi::pack_comm_hybrid(int n, int *list, double *buf)
{
  int i,j,m;

  if (radvary == 0) return 0;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = radius[j];
    buf[m++] = rmass[j];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecDemsi::unpack_comm(int n, int first, double *buf)
{
  int i,m,last;

  if (radvary == 0) {
    m = 0;
    last = first + n;
    for (i = first; i < last; i++) {
      x[i][0] = buf[m++];
      x[i][1] = buf[m++];
      x[i][2] = buf[m++];
    }
  } else {
    m = 0;
    last = first + n;
    for (i = first; i < last; i++) {
      x[i][0] = buf[m++];
      x[i][1] = buf[m++];
      x[i][2] = buf[m++];
      radius[i] = buf[m++];
      rmass[i] = buf[m++];
    }
  }
}

/* ---------------------------------------------------------------------- */

void AtomVecDemsi::unpack_comm_vel(int n, int first, double *buf)
{
  int i,m,last;

  if (radvary == 0) {
    m = 0;
    last = first + n;
    for (i = first; i < last; i++) {
      x[i][0] = buf[m++];
      x[i][1] = buf[m++];
      x[i][2] = buf[m++];
      v[i][0] = buf[m++];
      v[i][1] = buf[m++];
      v[i][2] = buf[m++];
      omega[i][0] = buf[m++];
      omega[i][1] = buf[m++];
      omega[i][2] = buf[m++];
    }
  } else {
    m = 0;
    last = first + n;
    for (i = first; i < last; i++) {
      x[i][0] = buf[m++];
      x[i][1] = buf[m++];
      x[i][2] = buf[m++];
      radius[i] = buf[m++];
      rmass[i] = buf[m++];
      v[i][0] = buf[m++];
      v[i][1] = buf[m++];
      v[i][2] = buf[m++];
      omega[i][0] = buf[m++];
      omega[i][1] = buf[m++];
      omega[i][2] = buf[m++];
    }
  }
}

/* ---------------------------------------------------------------------- */

int AtomVecDemsi::unpack_comm_hybrid(int n, int first, double *buf)
{
  int i,m,last;

  if (radvary == 0) return 0;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    radius[i] = buf[m++];
    rmass[i] = buf[m++];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecDemsi::pack_reverse(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = f[i][0];
    buf[m++] = f[i][1];
    buf[m++] = f[i][2];
    buf[m++] = torque[i][0];
    buf[m++] = torque[i][1];
    buf[m++] = torque[i][2];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecDemsi::pack_reverse_hybrid(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = torque[i][0];
    buf[m++] = torque[i][1];
    buf[m++] = torque[i][2];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecDemsi::unpack_reverse(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    f[j][0] += buf[m++];
    f[j][1] += buf[m++];
    f[j][2] += buf[m++];
    torque[j][0] += buf[m++];
    torque[j][1] += buf[m++];
    torque[j][2] += buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

int AtomVecDemsi::unpack_reverse_hybrid(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    torque[j][0] += buf[m++];
    torque[j][1] += buf[m++];
    torque[j][2] += buf[m++];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecDemsi::pack_border(int n, int *list, double *buf,
                               int pbc_flag, int *pbc)
{
  int i,j,m;
  double dx,dy,dz;

  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];
  buf[m++] = forcing[j][0];
  buf[m++] = forcing[j][1];
  buf[m++] = mean_thickness[j];
  buf[m++] = min_thickness[j];
  buf[m++] = ice_area[j];
  buf[m++] = coriolis[j];
  buf[m++] = ocean_vel[j][0];
  buf[m++] = ocean_vel[j][1];
  buf[m++] = bvector[j][0];
  buf[m++] = bvector[j][1];
      buf[m++] = ubuf(tag[j]).d;
      buf[m++] = ubuf(type[j]).d;
      buf[m++] = ubuf(mask[j]).d;
      buf[m++] = radius[j];
      buf[m++] = rmass[j];
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
      buf[m++] = x[j][0] + dx;
      buf[m++] = x[j][1] + dy;
      buf[m++] = x[j][2] + dz;
  buf[m++] = forcing[j][0];
  buf[m++] = forcing[j][1];
  buf[m++] = mean_thickness[j];
  buf[m++] = min_thickness[j];
  buf[m++] = ice_area[j];
  buf[m++] = coriolis[j];
  buf[m++] = ocean_vel[j][0];
  buf[m++] = ocean_vel[j][1];
  buf[m++] = bvector[j][0];
  buf[m++] = bvector[j][1];
      buf[m++] = ubuf(tag[j]).d;
      buf[m++] = ubuf(type[j]).d;
      buf[m++] = ubuf(mask[j]).d;
      buf[m++] = radius[j];
      buf[m++] = rmass[j];
    }
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->pack_border(n,list,&buf[m]);

  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecDemsi::pack_border_vel(int n, int *list, double *buf,
                                   int pbc_flag, int *pbc)
{
  int i,j,m;
  double dx,dy,dz,dvx,dvy,dvz;

  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];
  buf[m++] = forcing[j][0];
  buf[m++] = forcing[j][1];
  buf[m++] = mean_thickness[j];
  buf[m++] = min_thickness[j];
  buf[m++] = ice_area[j];
  buf[m++] = coriolis[j];
  buf[m++] = ocean_vel[j][0];
  buf[m++] = ocean_vel[j][1];
  buf[m++] = bvector[j][0];
  buf[m++] = bvector[j][1];
      buf[m++] = ubuf(tag[j]).d;
      buf[m++] = ubuf(type[j]).d;
      buf[m++] = ubuf(mask[j]).d;
      buf[m++] = radius[j];
      buf[m++] = rmass[j];
      buf[m++] = v[j][0];
      buf[m++] = v[j][1];
      buf[m++] = v[j][2];
      buf[m++] = omega[j][0];
      buf[m++] = omega[j][1];
      buf[m++] = omega[j][2];
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
        buf[m++] = x[j][0] + dx;
        buf[m++] = x[j][1] + dy;
        buf[m++] = x[j][2] + dz;
  buf[m++] = forcing[j][0];
  buf[m++] = forcing[j][1];
  buf[m++] = mean_thickness[j];
  buf[m++] = min_thickness[j];
  buf[m++] = ice_area[j];
  buf[m++] = coriolis[j];
  buf[m++] = ocean_vel[j][0];
  buf[m++] = ocean_vel[j][1];
  buf[m++] = bvector[j][0];
  buf[m++] = bvector[j][1];
        buf[m++] = ubuf(tag[j]).d;
        buf[m++] = ubuf(type[j]).d;
        buf[m++] = ubuf(mask[j]).d;
        buf[m++] = radius[j];
        buf[m++] = rmass[j];
        buf[m++] = v[j][0];
        buf[m++] = v[j][1];
        buf[m++] = v[j][2];
        buf[m++] = omega[j][0];
        buf[m++] = omega[j][1];
        buf[m++] = omega[j][2];
      }
    } else {
      dvx = pbc[0]*h_rate[0] + pbc[5]*h_rate[5] + pbc[4]*h_rate[4];
      dvy = pbc[1]*h_rate[1] + pbc[3]*h_rate[3];
      dvz = pbc[2]*h_rate[2];
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = x[j][0] + dx;
        buf[m++] = x[j][1] + dy;
        buf[m++] = x[j][2] + dz;
  buf[m++] = forcing[j][0];
  buf[m++] = forcing[j][1];
  buf[m++] = mean_thickness[j];
  buf[m++] = min_thickness[j];
  buf[m++] = ice_area[j];
  buf[m++] = coriolis[j];
  buf[m++] = ocean_vel[j][0];
  buf[m++] = ocean_vel[j][1];
  buf[m++] = bvector[j][0];
  buf[m++] = bvector[j][1];
        buf[m++] = ubuf(tag[j]).d;
        buf[m++] = ubuf(type[j]).d;
        buf[m++] = ubuf(mask[j]).d;
        buf[m++] = radius[j];
        buf[m++] = rmass[j];
        if (mask[i] & deform_groupbit) {
          buf[m++] = v[j][0] + dvx;
          buf[m++] = v[j][1] + dvy;
          buf[m++] = v[j][2] + dvz;
        } else {
          buf[m++] = v[j][0];
          buf[m++] = v[j][1];
          buf[m++] = v[j][2];
        }
        buf[m++] = omega[j][0];
        buf[m++] = omega[j][1];
        buf[m++] = omega[j][2];
      }
    }
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->pack_border(n,list,&buf[m]);

  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecDemsi::pack_border_hybrid(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = radius[j];
    buf[m++] = rmass[j];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecDemsi::unpack_border(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    if (i == nmax) grow(0);
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
forcing[i][0] = buf[m++];
forcing[i][1] = buf[m++];
mean_thickness[i] = buf[m++];
min_thickness[i] = buf[m++];
ice_area[i] = buf[m++];
coriolis[i] = buf[m++];
ocean_vel[i][0] = buf[m++];
ocean_vel[i][1] = buf[m++];
bvector[i][0] = buf[m++];
bvector[i][1] = buf[m++];
    tag[i] = (tagint) ubuf(buf[m++]).i;
    type[i] = (int) ubuf(buf[m++]).i;
    mask[i] = (int) ubuf(buf[m++]).i;
    radius[i] = buf[m++];
    rmass[i] = buf[m++];
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->
        unpack_border(n,first,&buf[m]);
}


/* ---------------------------------------------------------------------- */

void AtomVecDemsi::unpack_border_vel(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    if (i == nmax) grow(0);
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
forcing[i][0] = buf[m++];
forcing[i][1] = buf[m++];
mean_thickness[i] = buf[m++];
min_thickness[i] = buf[m++];
ice_area[i] = buf[m++];
coriolis[i] = buf[m++];
ocean_vel[i][0] = buf[m++];
ocean_vel[i][1] = buf[m++];
bvector[i][0] = buf[m++];
bvector[i][1] = buf[m++];
    tag[i] = (tagint) ubuf(buf[m++]).i;
    type[i] = (int) ubuf(buf[m++]).i;
    mask[i] = (int) ubuf(buf[m++]).i;
    radius[i] = buf[m++];
    rmass[i] = buf[m++];
    v[i][0] = buf[m++];
    v[i][1] = buf[m++];
    v[i][2] = buf[m++];
    omega[i][0] = buf[m++];
    omega[i][1] = buf[m++];
    omega[i][2] = buf[m++];
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->
        unpack_border(n,first,&buf[m]);
}

/* ---------------------------------------------------------------------- */

int AtomVecDemsi::unpack_border_hybrid(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    radius[i] = buf[m++];
    rmass[i] = buf[m++];
  }
  return m;
}

/* ----------------------------------------------------------------------
   pack data for atom I for sending to another proc
   xyz must be 1st 3 values, so comm::exchange() can test on them
------------------------------------------------------------------------- */

int AtomVecDemsi::pack_exchange(int i, double *buf)
{
  int m = 1;
  buf[m++] = x[i][0];
  buf[m++] = x[i][1];
  buf[m++] = x[i][2];
  buf[m++] = forcing[i][0];
  buf[m++] = forcing[i][1];
  buf[m++] = mean_thickness[i];
  buf[m++] = min_thickness[i];
  buf[m++] = ice_area[i];
  buf[m++] = coriolis[i];
  buf[m++] = ocean_vel[i][0];
  buf[m++] = ocean_vel[i][1];
  buf[m++] = bvector[i][0];
  buf[m++] = bvector[i][1];
  buf[m++] = v[i][0];
  buf[m++] = v[i][1];
  buf[m++] = v[i][2];
  buf[m++] = ubuf(tag[i]).d;
  buf[m++] = ubuf(type[i]).d;
  buf[m++] = ubuf(mask[i]).d;
  buf[m++] = ubuf(image[i]).d;

  buf[m++] = radius[i];
  buf[m++] = rmass[i];
  buf[m++] = omega[i][0];
  buf[m++] = omega[i][1];
  buf[m++] = omega[i][2];


  buf[m++] = ubuf(num_bond[i]).d;
  for (int k = 0; k < num_bond[i]; k++) {
    buf[m++] = ubuf(bond_type[i][k]).d;
    buf[m++] = ubuf(bond_atom[i][k]).d;
  }
  buf[m++] = ubuf(nspecial[i][0]).d;
  buf[m++] = ubuf(nspecial[i][1]).d;
  buf[m++] = ubuf(nspecial[i][2]).d;
  for (int k = 0; k < nspecial[i][2]; k++) buf[m++] = ubuf(special[i][k]).d;

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      m += modify->fix[atom->extra_grow[iextra]]->pack_exchange(i,&buf[m]);

  buf[0] = m;
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecDemsi::unpack_exchange(double *buf)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);

  int m = 1;
  x[nlocal][0] = buf[m++];
  x[nlocal][1] = buf[m++];
  x[nlocal][2] = buf[m++];
forcing[nlocal][0] = buf[m++];
forcing[nlocal][1] = buf[m++];
mean_thickness[nlocal] = buf[m++];
min_thickness[nlocal] = buf[m++];
ice_area[nlocal] = buf[m++];
coriolis[nlocal] = buf[m++];
ocean_vel[nlocal][0] = buf[m++];
ocean_vel[nlocal][1] = buf[m++];
bvector[nlocal][0] = buf[m++];
bvector[nlocal][1] = buf[m++];
  v[nlocal][0] = buf[m++];
  v[nlocal][1] = buf[m++];
  v[nlocal][2] = buf[m++];
  tag[nlocal] = (tagint) ubuf(buf[m++]).i;
  type[nlocal] = (int) ubuf(buf[m++]).i;
  mask[nlocal] = (int) ubuf(buf[m++]).i;
  image[nlocal] = (imageint) ubuf(buf[m++]).i;

  radius[nlocal] = buf[m++];
  rmass[nlocal] = buf[m++];
  omega[nlocal][0] = buf[m++];
  omega[nlocal][1] = buf[m++];
  omega[nlocal][2] = buf[m++];

num_bond[nlocal] = (int) ubuf(buf[m++]).i;
for (int k = 0; k < num_bond[nlocal]; k++) {
  bond_type[nlocal][k] = (int) ubuf(buf[m++]).i;
  bond_atom[nlocal][k] = (tagint) ubuf(buf[m++]).i;
}
nspecial[nlocal][0] = (int) ubuf(buf[m++]).i;
nspecial[nlocal][1] = (int) ubuf(buf[m++]).i;
nspecial[nlocal][2] = (int) ubuf(buf[m++]).i;
for (int k = 0; k < nspecial[nlocal][2]; k++)
  special[nlocal][k] = (tagint) ubuf(buf[m++]).i;

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      m += modify->fix[atom->extra_grow[iextra]]->
        unpack_exchange(nlocal,&buf[m]);

  atom->nlocal++;
  return m;
}

/* ----------------------------------------------------------------------
   size of restart data for all atoms owned by this proc
   include extra data stored by fixes
------------------------------------------------------------------------- */

int AtomVecDemsi::size_restart()
{
  int i;

  int nlocal = atom->nlocal;
  int n = 16 * nlocal;

  if (atom->nextra_restart)
    for (int iextra = 0; iextra < atom->nextra_restart; iextra++)
      for (i = 0; i < nlocal; i++)
        n += modify->fix[atom->extra_restart[iextra]]->size_restart(i);

  return n;
}

/* ----------------------------------------------------------------------
   pack atom I's data for restart file including extra quantities
   xyz must be 1st 3 values, so that read_restart can test on them
molecular = 1;
------------------------------------------------------------------------- */

int AtomVecDemsi::pack_restart(int i, double *buf)
{
  int m = 1;
  buf[m++] = x[i][0];
  buf[m++] = x[i][1];
  buf[m++] = x[i][2];
  buf[m++] = forcing[i][0];
  buf[m++] = forcing[i][1];
  buf[m++] = mean_thickness[i];
  buf[m++] = min_thickness[i];
  buf[m++] = ice_area[i];
  buf[m++] = coriolis[i];
  buf[m++] = ocean_vel[i][0];
  buf[m++] = ocean_vel[i][1];
  buf[m++] = bvector[i][0];
  buf[m++] = bvector[i][1];
  buf[m++] = ubuf(tag[i]).d;
  buf[m++] = ubuf(type[i]).d;
  buf[m++] = ubuf(mask[i]).d;
  buf[m++] = ubuf(image[i]).d;
  buf[m++] = v[i][0];
  buf[m++] = v[i][1];
  buf[m++] = v[i][2];

  buf[m++] = radius[i];
  buf[m++] = rmass[i];
  buf[m++] = omega[i][0];
  buf[m++] = omega[i][1];
  buf[m++] = omega[i][2];


  buf[m++] = ubuf(num_bond[i]).d;
  for (int k = 0; k < num_bond[i]; k++) {
    buf[m++] = ubuf(bond_type[i][k]).d;
    buf[m++] = ubuf(bond_atom[i][k]).d;
  }
  buf[m++] = ubuf(nspecial[i][0]).d;
  buf[m++] = ubuf(nspecial[i][1]).d;
  buf[m++] = ubuf(nspecial[i][2]).d;
  for (int k = 0; k < nspecial[i][2]; k++) buf[m++] = ubuf(special[i][k]).d;

  if (atom->nextra_restart)
    for (int iextra = 0; iextra < atom->nextra_restart; iextra++)
      m += modify->fix[atom->extra_restart[iextra]]->pack_restart(i,&buf[m]);

  buf[0] = m;
  return m;
}

/* ----------------------------------------------------------------------
   unpack data for one atom from restart file including extra quantities
------------------------------------------------------------------------- */

int AtomVecDemsi::unpack_restart(double *buf)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) {
    grow(0);
    if (atom->nextra_store)
      memory->grow(atom->extra,nmax,atom->nextra_store,"atom:extra");
  }

  int m = 1;
  x[nlocal][0] = buf[m++];
  x[nlocal][1] = buf[m++];
  x[nlocal][2] = buf[m++];
forcing[nlocal][0] = buf[m++];
forcing[nlocal][1] = buf[m++];
mean_thickness[nlocal] = buf[m++];
min_thickness[nlocal] = buf[m++];
ice_area[nlocal] = buf[m++];
coriolis[nlocal] = buf[m++];
ocean_vel[nlocal][0] = buf[m++];
ocean_vel[nlocal][1] = buf[m++];
bvector[nlocal][0] = buf[m++];
bvector[nlocal][1] = buf[m++];
  tag[nlocal] = (tagint) ubuf(buf[m++]).i;
  type[nlocal] = (int) ubuf(buf[m++]).i;
  mask[nlocal] = (int) ubuf(buf[m++]).i;
  image[nlocal] = (imageint) ubuf(buf[m++]).i;
  v[nlocal][0] = buf[m++];
  v[nlocal][1] = buf[m++];
  v[nlocal][2] = buf[m++];

  radius[nlocal] = buf[m++];
  rmass[nlocal] = buf[m++];
  omega[nlocal][0] = buf[m++];
  omega[nlocal][1] = buf[m++];
  omega[nlocal][2] = buf[m++];

num_bond[nlocal] = (int) ubuf(buf[m++]).i;
for (int k = 0; k < num_bond[nlocal]; k++) {
  bond_type[nlocal][k] = (int) ubuf(buf[m++]).i;
  bond_atom[nlocal][k] = (tagint) ubuf(buf[m++]).i;
}
nspecial[nlocal][0] = (int) ubuf(buf[m++]).i;
nspecial[nlocal][1] = (int) ubuf(buf[m++]).i;
nspecial[nlocal][2] = (int) ubuf(buf[m++]).i;
for (int k = 0; k < nspecial[nlocal][2]; k++)
  special[nlocal][k] = (tagint) ubuf(buf[m++]).i;

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

void AtomVecDemsi::create_atom(int itype, double *coord)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);

  tag[nlocal] = 0;
  type[nlocal] = itype;
  x[nlocal][0] = coord[0];
  x[nlocal][1] = coord[1];
  x[nlocal][2] = coord[2];
  mask[nlocal] = 1;
  image[nlocal] = ((imageint) IMGMAX << IMG2BITS) |
    ((imageint) IMGMAX << IMGBITS) | IMGMAX;
  v[nlocal][0] = 0.0;
  v[nlocal][1] = 0.0;
  v[nlocal][2] = 0.0;

  radius[nlocal] = 0.5;
rmass[nlocal] = MY_PI*radius[nlocal]*radius[nlocal];
  omega[nlocal][1] = 0.0;
  omega[nlocal][2] = 0.0;

  forcing[nlocal][0] = 0.0;
  forcing[nlocal][1] = 0.0;
  mean_thickness[nlocal] = 0.0;
  min_thickness[nlocal] = 0.0;
  ice_area[nlocal] = 0.0;
  coriolis[nlocal] = 0.0;
  ocean_vel[nlocal][0] = 0.0;
  ocean_vel[nlocal][1] = 0.0;
  bvector[nlocal][0] = 0.0;
  bvector[nlocal][1] = 0.0;

num_bond[nlocal] = 0;
nspecial[nlocal][0] = nspecial[nlocal][1] = nspecial[nlocal][2] = 0;

  atom->nlocal++;
}

/* ----------------------------------------------------------------------
   unpack one line from Atoms section of data file
   initialize other atom quantities
------------------------------------------------------------------------- */

void AtomVecDemsi::data_atom(double *coord, imageint imagetmp, char **values)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);

  tag[nlocal] = ATOTAGINT(values[0]);
  type[nlocal] = atoi(values[1]);
  if (type[nlocal] <= 0 || type[nlocal] > atom->ntypes)
    error->one(FLERR,"Invalid atom type in Atoms section of data file");

  radius[nlocal] = 0.5 * atof(values[2]);
  if (radius[nlocal] < 0.0)
    error->one(FLERR,"Invalid radius in Atoms section of data file");

  double density = atof(values[3]);
  if (density <= 0.0)
    error->one(FLERR,"Invalid density in Atoms section of data file");

  if (radius[nlocal] == 0.0) rmass[nlocal] = density;
  else
rmass[nlocal] = MY_PI*radius[nlocal]*radius[nlocal];

  x[nlocal][0] = coord[0];
  x[nlocal][1] = coord[1];
  x[nlocal][2] = coord[2];

  image[nlocal] = imagetmp;

  mask[nlocal] = 1;
  v[nlocal][0] = 0.0;
  v[nlocal][1] = 0.0;
  v[nlocal][2] = 0.0;
  omega[nlocal][0] = 0.0;
  omega[nlocal][1] = 0.0;
  omega[nlocal][2] = 0.0;

  forcing[nlocal][0] = 0.0;
  forcing[nlocal][1] = 0.0;
  mean_thickness[nlocal] = 0.0;
  min_thickness[nlocal] = 0.0;
  ice_area[nlocal] = 0.0;
  coriolis[nlocal] = 0.0;
  ocean_vel[nlocal][0] = 0.0;
  ocean_vel[nlocal][1] = 0.0;
  bvector[nlocal][0] = 0.0;
  bvector[nlocal][1] = 0.0;

num_bond[nlocal] = 0;
nspecial[nlocal][0] = nspecial[nlocal][1] = nspecial[nlocal][2] = 0;

  atom->nlocal++;
}

/* ----------------------------------------------------------------------
   unpack hybrid quantities from one line in Atoms section of data file
   initialize other atom quantities for this sub-style
------------------------------------------------------------------------- */

int AtomVecDemsi::data_atom_hybrid(int nlocal, char **values)
{
  radius[nlocal] = 0.5 * atof(values[0]);
  if (radius[nlocal] < 0.0)
    error->one(FLERR,"Invalid radius in Atoms section of data file");

  double density = atof(values[1]);
  if (density <= 0.0)
    error->one(FLERR,"Invalid density in Atoms section of data file");

  if (radius[nlocal] == 0.0) rmass[nlocal] = density;
  else
rmass[nlocal] = MY_PI*radius[nlocal]*radius[nlocal];

  return 2;
}

/* ----------------------------------------------------------------------
   unpack one line from Velocities section of data file
------------------------------------------------------------------------- */

void AtomVecDemsi::data_vel(int m, char **values)
{
  v[m][0] = atof(values[0]);
  v[m][1] = atof(values[1]);
  v[m][2] = atof(values[2]);
  omega[m][0] = atof(values[3]);
  omega[m][1] = atof(values[4]);
  omega[m][2] = atof(values[5]);
}

/* ----------------------------------------------------------------------
   unpack hybrid quantities from one line in Velocities section of data file
------------------------------------------------------------------------- */

int AtomVecDemsi::data_vel_hybrid(int m, char **values)
{
  omega[m][0] = atof(values[0]);
  omega[m][1] = atof(values[1]);
  omega[m][2] = atof(values[2]);
  return 3;
}

/* ----------------------------------------------------------------------
   pack atom info for data file including 3 image flags
------------------------------------------------------------------------- */

void AtomVecDemsi::pack_data(double **buf)
{
  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    buf[i][0] = ubuf(tag[i]).d;
    buf[i][1] = ubuf(type[i]).d;
    buf[i][2] = 2.0*radius[i];
    if (radius[i] == 0.0) buf[i][3] = rmass[i];
    else
      buf[i][3] = rmass[i] / (MY_PI*radius[i]*radius[i]);
    buf[i][4] = x[i][0];
    buf[i][5] = x[i][1];
    buf[i][6] = x[i][2];
    buf[i][7] = ubuf((image[i] & IMGMASK) - IMGMAX).d;
    buf[i][8] = ubuf((image[i] >> IMGBITS & IMGMASK) - IMGMAX).d;
    buf[i][9] = ubuf((image[i] >> IMG2BITS) - IMGMAX).d;
  }
}

/* ----------------------------------------------------------------------
   pack hybrid atom info for data file
------------------------------------------------------------------------- */

int AtomVecDemsi::pack_data_hybrid(int i, double *buf)
{
  buf[0] = 2.0*radius[i];
  if (radius[i] == 0.0) buf[1] = rmass[i];
  else buf[1] = rmass[i] / (MY_PI*radius[i]*radius[i]);
  return 2;
}

/* ----------------------------------------------------------------------
   write atom info to data file including 3 image flags
------------------------------------------------------------------------- */

void AtomVecDemsi::write_data(FILE *fp, int n, double **buf)
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

int AtomVecDemsi::write_data_hybrid(FILE *fp, double *buf)
{
  fprintf(fp," %-1.16e %-1.16e",buf[0],buf[1]);
  return 2;
}

/* ----------------------------------------------------------------------
   pack velocity info for data file
------------------------------------------------------------------------- */

void AtomVecDemsi::pack_vel(double **buf)
{
  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    buf[i][0] = ubuf(tag[i]).d;
    buf[i][1] = v[i][0];
    buf[i][2] = v[i][1];
    buf[i][3] = v[i][2];
    buf[i][4] = omega[i][0];
    buf[i][5] = omega[i][1];
    buf[i][6] = omega[i][2];
  }
}

/* ----------------------------------------------------------------------
   pack hybrid velocity info for data file
------------------------------------------------------------------------- */

int AtomVecDemsi::pack_vel_hybrid(int i, double *buf)
{
  buf[0] = omega[i][0];
  buf[1] = omega[i][1];
  buf[2] = omega[i][2];
  return 3;
}

/* ----------------------------------------------------------------------
   write velocity info to data file
------------------------------------------------------------------------- */

void AtomVecDemsi::write_vel(FILE *fp, int n, double **buf)
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

int AtomVecDemsi::write_vel_hybrid(FILE *fp, double *buf)
{
  fprintf(fp," %-1.16e %-1.16e %-1.16e",buf[0],buf[1],buf[2]);
  return 3;
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
------------------------------------------------------------------------- */

bigint AtomVecDemsi::memory_usage()
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
