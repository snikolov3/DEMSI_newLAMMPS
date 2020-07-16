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
#include <cstdio>
#include <cstring>
#include "fix_nve_sphere_demsi.h"
#include "atom.h"
#include "domain.h"
#include "atom_vec.h"
#include "update.h"
#include "respa.h"
#include "force.h"
#include "error.h"
#include "math_vector.h"
#include "math_extra.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathExtra;

enum{NONE,DIPOLE};
enum{NODLM,DLM};

/* ---------------------------------------------------------------------- */

FixNVESphereDemsi::FixNVESphereDemsi(LAMMPS *lmp, int narg, char **arg) :
  FixNVE(lmp, narg, arg)
{
  if (narg < 3) error->all(FLERR,"Illegal fix nve/sphere/demsi command");

  time_integrate = 1;

  // process extra keywords
  // inertia = moment of inertia prefactor for sphere or disc
  inertia = 0.5;

  ocean_density = ocean_drag = 0;

  if (domain->dimension != 2)
    error->all(FLERR,"Fix nve/sphere demsi requires 2d simulation");
  if (!atom->demsi_flag)
    error->all(FLERR,"Fix nve/sphere requires atom style demsi");
}

/* ---------------------------------------------------------------------- */

void FixNVESphereDemsi::init()
{
  FixNVE::init();

  // check that all particles are finite-size spheres
  // no point particles allowed

  double *radius = atom->radius;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit)
      if (radius[i] == 0.0)
        error->one(FLERR,"Fix nve/sphere requires extended particles");
}

/* ---------------------------------------------------------------------- */

void FixNVESphereDemsi::initial_integrate(int /*vflag*/)
{
  double dtfm;

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **omega = atom->omega;
  double **torque = atom->torque;
  double *radius = atom->radius;
  double *rmass = atom->rmass;
  double *ice_area = atom->ice_area;
  double *coriolis = atom->coriolis;
  double **ocean_vel = atom->ocean_vel;
  double **bvector = atom->bvector;
  double **forcing = atom->forcing;

  double D, vel_diff, m_prime;
  double a00, a01, a10, a11;
  double detinv;
  double b0, b1;
  double det;

  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  // set timestep here since dt may have changed or come via rRESPA

  double dtfrotate = dtf / inertia;
  double dtirotate;
  // update v,x,omega for all particles
  // d_omega/dt = torque / inertia

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      vel_diff = sqrt((ocean_vel[i][0]-v[i][0])*(ocean_vel[i][0]-v[i][0]) +
          (ocean_vel[i][1]-v[i][1])*(ocean_vel[i][1]-v[i][1]));
      D = ice_area[i]*ocean_drag*ocean_density*vel_diff;
      m_prime = rmass[i]/dtf;
      a00 = a11 = m_prime+D;
      a10 = rmass[i]*coriolis[i];
      a01 = -a10;

      b0 = m_prime*v[i][0] + f[i][0] + bvector[i][0] + forcing[i][0] + D*ocean_vel[i][0];
      b1 = m_prime*v[i][1] + f[i][1] + bvector[i][1] + forcing[i][1] + D*ocean_vel[i][1];

      detinv = 1.0/(a00*a11 - a01*a10);
      v[i][0] = detinv*( a11*b0 - a01*b1);
      v[i][1] = detinv*(-a10*b0 + a00*b1);

      x[i][0] += dtv * v[i][0];
      x[i][1] += dtv * v[i][1];

      dtirotate = dtfrotate / (radius[i]*radius[i]*rmass[i]);
      omega[i][2] += dtirotate * torque[i][2];
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixNVESphereDemsi::final_integrate()
{
  double dtfm,dtirotate;

  double **v = atom->v;
  double **f = atom->f;
  double **omega = atom->omega;
  double **torque = atom->torque;
  double *rmass = atom->rmass;
  double *ice_area = atom->ice_area;
  double *coriolis = atom->coriolis;
  double **ocean_vel = atom->ocean_vel;
  double **bvector = atom->bvector;
  double **forcing = atom->forcing;

  double *radius = atom->radius;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  double D, vel_diff, m_prime;
  double a00, a01, a10, a11;
  double detinv;
  double b0, b1;
  double det;

  // set timestep here since dt may have changed or come via rRESPA

  double dtfrotate = dtf / inertia;

  // update v,omega for all particles
  // d_omega/dt = torque / inertia

  double rke = 0.0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      vel_diff = sqrt((ocean_vel[i][0]-v[i][0])*(ocean_vel[i][0]-v[i][0]) +
          (ocean_vel[i][1]-v[i][1])*(ocean_vel[i][1]-v[i][1]));
      D = ice_area[i]*ocean_drag*ocean_density*vel_diff;
      m_prime = rmass[i]/dtf;
      a00 = a11 = m_prime+D;
      a10 = rmass[i]*coriolis[i];
      a01 = -a10;

      b0 = m_prime*v[i][0] + f[i][0] + bvector[i][0] + forcing[i][0] + D*ocean_vel[i][0];
      b1 = m_prime*v[i][1] + f[i][1] + bvector[i][1] + forcing[i][1] + D*ocean_vel[i][1];

      detinv = 1.0/(a00*a11 - a01*a10);
      v[i][0] = detinv*( a11*b0 - a01*b1);
      v[i][1] = detinv*(-a10*b0 + a00*b1);

      dtirotate = dtfrotate / (radius[i]*radius[i]*rmass[i]);
      omega[i][2] += dtirotate * torque[i][2];
      rke += (omega[i][0]*omega[i][0] + omega[i][1]*omega[i][1] +
              omega[i][2]*omega[i][2])*radius[i]*radius[i]*rmass[i];
    }

}
