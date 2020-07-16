/* ----------------------------------------------------------------------
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
   ------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Dan S. Bolintineanu (SNL), Adrian K. Turner (LANL)
   ------------------------------------------------------------------------- */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "pair_gran_hopkins.h"
#include "atom.h"
#include "update.h"
#include "force.h"
#include "fix.h"
#include "fix_neigh_history.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "comm.h"
#include "memory.h"
#include "error.h"
#include "math_const.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define EPSILON 1e-10

#define DEBUGID_1 35
#define DEBUGID_2 55
#define DEBUG_TIMESTEP 32696
/* ---------------------------------------------------------------------- */

PairGranHopkins::PairGranHopkins(LAMMPS *lmp) :
              PairGranHookeHistory(lmp, 12)
{
  nondefault_history_transfer = 1;
  beyond_contact = 1;
}

/* ---------------------------------------------------------------------- */
PairGranHopkins::~PairGranHopkins()
{
}
/* ---------------------------------------------------------------------- */

void PairGranHopkins::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum;
  int itype,jtype;

  int *ilist,*jlist,*numneigh,**firstneigh;
  int *touch,**firsttouch;
  double *history,*allhistory,**firsthistory;

  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  // update rigid body info for owned & ghost atoms if using FixRigid masses
  // body[i] = which body atom I is in, -1 if none
  // mass_body = mass of each rigid body
  // Not yet applicable for DEMSI, but may be something to look into

  if (fix_rigid && neighbor->ago == 0){
    int tmp;
    int *body = (int *) fix_rigid->extract("body",tmp);
    double *mass_body = (double *) fix_rigid->extract("masstotal",tmp);
    if (atom->nmax > nmax) {
      memory->destroy(mass_rigid);
      nmax = atom->nmax;
      memory->create(mass_rigid,nmax,"pair:mass_rigid");
    }
    int nlocal = atom->nlocal;
    for (i = 0; i < nlocal; i++)
      if (body[i] >= 0) mass_rigid[i] = mass_body[body[i]];
      else mass_rigid[i] = 0.0;
    comm->forward_comm_pair(this);
  }

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  firsttouch = fix_history->firstflag;
  firsthistory = fix_history->firstvalue;

  // loop over my elements
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itype = atom->type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    allhistory = firsthistory[i];

    // loop over neighbors of each element
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      jtype = atom->type[j];
      j &= NEIGHMASK;

      history = &allhistory[size_history*jj];
      //'history' now points to the ii-jj array that stores
      //all the history associated with pair ii-jj
      //For bonded pairs:
      // history[0-7]: x,y components of s1i, s2i, s1j, s2j
      // history[8,9]: chi1, chi2
      // history[10]: bond length
      // history[11]: bond thickness, h
      //For unbonded pairs:
      // history[0,1]: normal force at previous time step, x and y components
      // history[2,3]: accumulated tangential displacement at contact, x and y
      // history[4] : contact thickness h
      // history[5] : delta_0
      // history[6] : flag for whether plastic deformation has occurred

      if (history[8] >= history[9]){ // Un-bonded, chi1 >= chi2
        compute_nonbonded(history, &firsttouch[i][jj], i, j);
      }
      else { //Bonded
        compute_bonded(history, &firsttouch[i][jj], i, j);
      }
    }
  }
  if (vflag_fdotr) virial_fdotr_compute();
}

void PairGranHopkins::compute_nonbonded(double *history, int* touch, int i, int j){

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **omega = atom->omega;
  double **torque = atom->torque;
  double *radius = atom->radius;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double delx, dely;
  double r, rsq, rinv, radsum, radmin;
  double nx, ny, vrx, vry, vnnr, vnx, vny;
  double wrz, vtrx, vtry, vtx, vty, vrel;
  double delta, delta_dot;

  double fnx, fny;
  double sig_c, hmin;
  double hprime, kp, kr, ke, L, hstar;
  double num, denom, fnmag_plastic, fnmag_elastic, fnmag;
  double fx, fy;
  double ncrossF;

  double ndisp, dispmag, scalefac;
  double ftx, fty, ftmag, ftcrit;

  int historyupdate = 1;
  if (update->setupflag) historyupdate = 0;

  //Hard-coded, ask Adrian what this should actually be??
  hstar = 0.3;

  //This is different than in LammpsInstance::initialize_bonded_contacts,
  // to conform to the sign convention in comparable granular pair styles
  delx = x[i][0] - x[j][0];
  dely = x[i][1] - x[j][1];
  rsq = delx*delx + dely*dely;
  radsum = radius[i] + radius[j];

  if (rsq >= radsum*radsum){
    *touch = 0;
    for (int k = 0; k < size_history; k++){
      history[k] = 0;
    }
    delx = dely = fx = fy = 0; //for virial calculation
    //history[0] = history[4] = 1234; //for debugging
  }
  else{
    if (!*touch){ //If this is first contact
      *touch = 1;
      history[0] = history[1] = history[2] = history[3] = history[5] = 0;
      history[4] = hprime_0;
    }
    r = sqrt(rsq);
    rinv = 1.0/r;
    nx = delx/r;
    ny = dely/r;

    radmin = MIN(radius[i], radius[j]);
    L = 2*radmin*(1+(abs(radius[i] - radius[j])/r));

    // relative translational velocity
    vrx = v[i][0] - v[j][0];
    vry = v[i][1] - v[j][1];


    delta = radsum - r - history[5];
    if (delta < 0) return ; //delta = 0;

    // Compute tangential force
    // normal component of relative translational velocity
    vnnr = vrx*nx + vry*ny;
    vnx = nx*vnnr;
    vny = ny*vnnr;

    // subtract to compute tangential component of relative translational velocity
    vtrx = vrx - vnx;
    vtry = vry - vny;

    // total relative tangential velocities at contact
    wrz = radius[i]*omega[i][2] + radius[j]*omega[j][2];
    vtx = vtrx + ny*wrz;
    vty = vtry - nx*wrz;

    vrel = vtx*vtx + vty*vty;
    vrel = sqrt(vrel);

    delta_dot = -vnnr;

    // Compute plastic normal force
    hprime = history[4];
    ke = Emod/L*(1/(1/atom->mean_thickness[i] + 1/atom->mean_thickness[j]));
    hmin = MIN(atom->min_thickness[i], atom->min_thickness[j]);
    if (hprime < hstar){
      kr = 26126*hprime;
      kp = 928*hprime*hprime;
    }
    else{
      kr = kp = hprime*sig_c;
    }

    num = history[0]/(kp*update->dt) + delta_dot*L + delta*L*ke/damp_normal + ke*kr/(damp_normal*kp);
    denom = 1/(kp*update->dt) + 1/damp_normal*(1+ke/kp);
    fnmag_plastic = num/denom;

    // Elastic normal force
    fnmag_elastic = ke*delta*L + damp_normal*delta_dot*L;

    if (fabs(fnmag_elastic) < fabs(fnmag_plastic)){
      fnmag = fnmag_elastic;
      history[6] = 0;
    }
    else{
      fnmag = fnmag_plastic;
      history[6] = 1;
    }

    //fnmag = fnmag_elastic;

    fnx = fnmag*nx;
    fny = fnmag*ny;

    // update tangential displacement, rotate if needed
    if (historyupdate){
      history[0] = fnx;
      history[1] = fny;
      ndisp = nx*history[2] + ny*history[3];
      dispmag =sqrt( history[2]*history[2] + history[3]*history[3]);
      denom = dispmag - ndisp;
      if (ndisp > EPSILON && denom != 0){
        scalefac = dispmag/denom;
        history[2] -= ndisp*nx;
        history[3] -= ndisp*ny;
        history[2] *= scalefac;
        history[3] *= scalefac;
      }
      history[2] += vtx*update->dt;
      history[3] += vty*update->dt;
    }

    dispmag = sqrt( history[2]*history[2] + history[3]*history[3]);

    // total tangential force
    kt = Gmod/L*(1/(1/atom->mean_thickness[i] + 1/atom->mean_thickness[j]));
    ftx = - (kt*history[2] + damp_tangential*vtx);
    fty = - (kt*history[3] + damp_tangential*vty);

    ftmag = sqrt(ftx*ftx + fty*fty);
    ftcrit = friction_tangential*fabs(fnmag);
    if (ftmag > ftcrit){
      if (dispmag != 0){
        ftx *= ftcrit/ftmag;
        fty *= ftcrit/ftmag;
        history[2] = -(ftcrit + damp_tangential*vtx)/kt;
        history[3] = -(ftcrit + damp_tangential*vty)/kt;
      }
      else ftx = fty = 0;
    }

    //Apply forces
    fx = fnx + ftx;
    fy = fny + fty;

    f[i][0] += fx;
    f[i][1] += fy;

    // torque induced by tangential force
    ncrossF = nx*fty - ny*ftx;
    torque[i][2] += -radius[i]*ncrossF;

    if (force->newton_pair || j < atom->nlocal){
      f[j][0] -= fx;
      f[j][1] -= fy;
      torque[j][2] += -radius[j]*ncrossF;
    }
  }
  if (evflag) ev_tally_xyz(i,j,atom->nlocal, force->newton_pair,
                                   0.0,0.0,fx,fy,0,delx,dely,0);
}

void PairGranHopkins::compute_bonded(double *history, int* touch, int i, int j){
  //See design document for definitions of these variables
  double s1x, s1y, s2x, s2y, mx, my, mmag, mex, mey, bex, bey, rx, ry, rxj, ryj;
  double An, Bn, Cn, Dn, Bt, Ct, Dt, Bnj, Cnj, Dnj, Btj;
  double Fnmag, Ftmag, Nn, Nt, Nnj, Ntj;

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **omega = atom->omega;
  double **torque = atom->torque;

  double chi1 = history[8];
  double chi2 = history[9];
  double chidiff, chidiff2, chidiff3;

  double sig_n1, sig_s1;
  double chi_c, chi_t, chi_s1, chi_s2;
  double nprefac, sprefac;
  double damp_prefac, area_bond, fdampx, fdampy, torquedamp;
  double fx, fy;
  double hmin;

  double dvx, dvy;

  int historyupdate = 1;
  if (update->setupflag) historyupdate = 0;

  if (historyupdate){
    // Update bond end points based on particle translations
    history[0] += dt*v[i][0];
    history[1] += dt*v[i][1];
    history[2] += dt*v[i][0];
    history[3] += dt*v[i][1];

    history[4] += dt*v[j][0];
    history[5] += dt*v[j][1];
    history[6] += dt*v[j][0];
    history[7] += dt*v[j][1];

    // Update bond end points based on particle rotations
    history[0] += -dt*omega[i][2]*(history[1]-x[i][1]);
    history[1] +=  dt*omega[i][2]*(history[0]-x[i][0]);
    history[2] += -dt*omega[i][2]*(history[3]-x[i][1]);
    history[3] +=  dt*omega[i][2]*(history[2]-x[i][0]);

    history[4] += -dt*omega[j][2]*(history[5]-x[j][1]);
    history[5] +=  dt*omega[j][2]*(history[4]-x[j][0]);
    history[6] += -dt*omega[j][2]*(history[7]-x[j][1]);
    history[7] +=  dt*omega[j][2]*(history[6]-x[j][0]);
  }

  //Compute s_1, s_2, m, m_e, b_e
  s1x = history[4] - history[0];
  s1y = history[5] - history[1];
  s2x = history[6] - history[2];
  s2y = history[7] - history[3];

  mx = history[2] + 0.5*s2x - history[0] - 0.5*s1x;
  my = history[3] + 0.5*s2y - history[1] - 0.5*s1y;
  mmag = sqrt(mx*mx + my*my);
  mex = mx/mmag;
  mey = my/mmag;

  bex = mey;
  bey = -mex;

  rx = history[0] + 0.5*s1x - x[i][0];
  ry = history[1] + 0.5*s1y - x[i][1];

  //Compute forces and torques
  Dn = s1x*bex + s1y*bey;
  Cn = s2x*bex + s2y*bey - Dn;
  Bn = rx*bey - ry*bex;
  An = mx*bey - my*bex;

  Dt = s1x*mex + s1y*mey;
  Ct = s2x*mex + s2y*mey - Dt;
  Bt = rx*mey - ry*mex;

  chidiff = chi2-chi1;
  chidiff2 = 0.5*(chi2*chi2 - chi1*chi1);
  chidiff3 = MathConst::THIRD*(chi2*chi2*chi2 - chi1*chi1*chi1);
  kn = history[11]*Emod/history[10];
  kt = history[11]*Gmod/history[10];
  nprefac = history[10]*kn;
  sprefac = history[10]*kt;

  Fnmag = nprefac*(Dn*chidiff + Cn*chidiff2);
  Ftmag = sprefac*(Dt*chidiff + Ct*chidiff2);
  Nn = nprefac*(An*Cn*chidiff3 + (Bn*Cn+An*Dn)*chidiff2 + Bn*Dn*chidiff);
  Nt = Bt*Ftmag;

  //Damping force
  area_bond = history[10]*history[11]*chidiff;
  damp_prefac = damp_bonded*area_bond;

  dvx = v[j][0] - v[i][0];
  dvy = v[j][1] - v[i][1];
  fdampx = damp_prefac*dvx;
  fdampy = damp_prefac*dvy;

  //Damping torque
  torquedamp = -damp_bonded*area_bond*area_bond*(omega[i][2]-omega[j][2]);

  //Update forces and torque
  fx = Fnmag*bex + Ftmag*mex;
  fy = Fnmag*bey + Ftmag*mey;

  //Ensure that damping does not cause force to change sign
 /* if (fx < 0 && fdampx > 0){
    fx = MIN(fx+fdampx, 0);
  }
  else if (fx > 0 && fdampx < 0){
    fx = MAX(fx+fdampx, 0);
  }
  else */
  fx = fx+fdampx;

  /*if (fy < 0 && fdampy > 0){
    fy = MIN(fy+fdampy, 0);
  }
  else if (fy > 0 && fdampy < 0){
    fy = MAX(fy+fdampy, 0);
  }
  else*/
  fy = fy+fdampy;

  f[i][0] += fx;
  f[i][1] += fy;

  torque[i][2] += Nn + Nt + torquedamp;

  if (force->newton_pair || j < atom->nlocal){
    f[j][0] -= fx;
    f[j][1] -= fy;

    rxj = history[0] + 0.5*s1x - x[j][0];
    ryj = history[1] + 0.5*s1y - x[j][1];
    Bnj = rxj*bey - ryj*bex;
    Cnj = -Cn;
    Dnj = -Dn;
    Btj = rxj*mey - ryj*mex;

    Nnj = nprefac*(An*Cnj*chidiff3 + (Bnj*Cnj+An*Dnj)*chidiff2 + Bnj*Dnj*chidiff);
    Ntj = -Btj*Ftmag;
    torque[j][2] += Nnj + Ntj - torquedamp;
  }

  //Update chi1, chi2
  hmin = MIN(atom->min_thickness[i], atom->min_thickness[j]);
  if (historyupdate){
    double c1, c2;
    c1 = history[8];
    c2 = history[9];
    update_chi(kn, kt, Dn, Cn, Dt, Ct, hmin, history[8], history[9]);
    *touch = 1;
    if (history[8] >= history[9]){ //Bond just broke
        double dx = x[i][0] - x[j][0];
        double dy = x[i][1] - x[j][1];
        double rij = sqrt(dx*dx + dy*dy);
        double delta_0 = atom->radius[i] + atom->radius[j] - rij;
        if (delta_0 < 0) delta_0 = 0;
        for (int k = 0; k < size_history; ++k)
          history[k] = 0;
        history[4] = hprime_0;
        history[5] = delta_0;
    }
  }
  if (evflag) ev_tally_xyz(i,j,atom->nlocal, force->newton_pair,
                                     0.0,0.0,fx,fy,0,x[i][0]-x[j][0],x[i][1]-x[j][1],0);
}


void PairGranHopkins::update_chi(double kn, double kt, double Dn, double Cn, double Dt, double Ct, double hmin, double &chi1, double &chi2){
  double sig_n1 = kn*(Dn + Cn*chi1);
  double sig_s1 = kt*(Dt + Ct*chi1);
  double sig_n2 = kn*(Dn + Cn*chi2);
  double sig_s2 = kt*(Dt + Ct*chi2);

  // function pointers for greater efficiency?
  double sig_c;
  if (strcmp(sig_c0_type,"constant") == 0) {
    sig_c = sig_c0;
  } else if (strcmp(sig_c0_type,"KovacsSodhi") == 0) {
    sig_c = sig_c0*pow(hmin,(2.0/3.0)) * 1000.0;
  } else {
    error->all(FLERR,"Unknown sig_c0_type");
  }
  double sig_t;
  if (strcmp(sig_t0_type,"constant") == 0) {
    sig_t = sig_t0;
  } else if (strcmp(sig_t0_type,"multiply_sig_c0") == 0) {
    sig_t = sig_t0 * sig_c;
  } else {
    error->all(FLERR,"Unknown sig_t0_type");
  }

  double denom;
  sig_t = -sig_t; //Somewhat against convention, tensile load is taken to be negative

  double c1, c2;
  c1 = chi1;
  c2 = chi2;
  //Check for purely tensile/compressive failure
  if (sig_n1 < sig_t || sig_n2 < sig_t || sig_n1 > sig_c || sig_n2 > sig_c){
    if (sig_n1 < sig_t){
      if (Cn != 0) chi1 = (sig_t/kn - Dn)/Cn;
      else {
        chi1 = chi2;
        return;
      }
    }
    if (sig_n2 < sig_t) chi2 = (sig_t/kn - Dn)/Cn;
    //if Cn==0, then sig_n1 = sig_n2, and function would've returned above

    if (sig_n1 > sig_c){
      if (Cn != 0) chi1 = (sig_c/kn - Dn)/Cn;
      else{
        chi1 = chi2;
        return;
      }
    }
    if (sig_n2 > sig_c) chi2 = (sig_c/kn - Dn)/Cn;
    //if Cn==0, function would've returned above
  }

  //Re-compute, since stress state could still be outside of failure envelope
  sig_n1 = kn*(Dn + Cn*chi1);
  sig_s1 = kt*(Dt + Ct*chi1);
  sig_n2 = kn*(Dn + Cn*chi2);
  sig_s2 = kt*(Dt + Ct*chi2);

  //Check for 'cohesion' shear failure at chi1
  //Top branch of envelope
  if (sig_s1 > tanphi*sig_n1 - tanphi*sig_t){
    denom = tanphi*kn*Cn - kt*Ct;
    if (denom != 0) chi1 = (kt*Dt+tanphi*(sig_t-kn*Dn))/denom;
    else {
      if (Cn != 0 && Ct != 0){
        //Rare case of kt*Ct = -kn*tanphi*Cn, yield criterion is independent of chi,
        // therefore bond breaks.
        chi1 = chi2;
        return;
      }
      else if (Cn == 0 && Ct == 0){
        //Rare case of s1 = s2, again yield criterion is independent of chi,
        // and there is a possibility of shear failure of entire bond
        if (kt*Dt > tanphi*(kn*Dn - sig_t)){
          chi1 = chi2;
          return;
        }
      }
    }
  }

  //Bottom branch of envelope
  if (sig_s1 < -tanphi*sig_n1 + tanphi*sig_t){
    denom = -kn*tanphi*Cn - kt*Ct;
    if (denom != 0) chi1 = (kt*Dt+tanphi*(kn*Dn-sig_t))/denom;
    else {
      if (Cn != 0 && Ct != 0){
        chi1 = chi2;
        return;
      }
      else if (Cn == 0 && Ct == 0){
        if (kt*Dt < -tanphi*(kn*Dn - sig_t)){
          chi1 = chi2;
          return;
        }
      }
    }
  }

  //Check for 'cohesion' shear failure at chi2
  //Top branch of envelope
  if (sig_s2 > tanphi*sig_n2 - tanphi*sig_t){
    denom = tanphi*kn*Cn - kt*Ct;
    if (denom != 0) chi2 = (kt*Dt+tanphi*(sig_t-kn*Dn))/denom;
    else {
      if (Cn != 0 && Ct != 0){
        //Rare case of kt*Ct = -kn*tanphi*Cn, yield criterion is independent of chi,
        // therefore bond breaks.
        chi2 = chi1;
        return;
      }
      else if (Cn == 0 && Ct == 0){
        if (kt*Dt > tanphi*(kn*Dn - sig_t)){
          chi1 = chi2;
          return;
        }
      }
    }
  }

  //Bottom branch of envelope
  if (sig_s2 < -tanphi*sig_n2 + tanphi*sig_t){
    denom = -kn*tanphi*Cn - kt*Ct;
    if (denom != 0) chi2 = (kt*Dt+tanphi*(kn*Dn-sig_t))/denom;
    else {
      if (Cn != 0 && Ct != 0){
        chi1 = chi2;
        return;
      }
      else if (Cn == 0 && Ct == 0){
        if (kt*Dt < -tanphi*(kn*Dn - sig_t)){
          chi1 = chi2;
          return;
        }
      }
    }
  }  
  if (chi1 < 0 || chi1 > 1 || chi2 < 0 || chi2 > 1)
    chi1 = chi2 = 0;
}

/* ----------------------------------------------------------------------
   global settings
   ------------------------------------------------------------------------- */

void PairGranHopkins::settings(int narg, char **arg)
{
  if (narg != 12) error->all(FLERR,"Illegal pair_style command");

  Emod = force->numeric(FLERR, arg[0]);
  poiss = force->numeric(FLERR, arg[1]);
  strcpy(sig_c0_type,arg[2]);
  sig_c0 = force->numeric(FLERR, arg[3]);
  strcpy(sig_t0_type,arg[4]);
  sig_t0 = force->numeric(FLERR, arg[5]);
  phi = force->numeric(FLERR, arg[6]);
  damp_bonded = force->numeric(FLERR, arg[7]);
  friction_tangential = force->numeric(FLERR, arg[8]);
  damp_normal = force->numeric(FLERR, arg[9]);
  damp_tangential = force->numeric(FLERR, arg[10]);
  hprime_0 = force->numeric(FLERR, arg[11]);

  tanphi = tan(phi*MathConst::MY_PI/180.0);
  Gmod = Emod/(2*(1+poiss));
}

/* ---------------------------------------------------------------------- */
double PairGranHopkins::init_one(int i, int j)
{
  double cutoff;
  cutoff = PairGranHookeHistory::init_one(i, j);
  cutoff += maxrad_dynamic[i]*0.1; //This could be an input parameter?
  return cutoff;
}

/* ---------------------------------------------------------------------- */

double PairGranHopkins::single(int i, int j, int itype, int jtype,
    double rsq,
    double factor_coul, double factor_lj,
    double &fforce)
{
  return 0.0;
}

/* ---------------------------------------------------------------------- */
void PairGranHopkins::transfer_history(double* sourcevalues, double* targetvalues){
  if (sourcevalues[8] < sourcevalues[9]){
    targetvalues[0] = sourcevalues[4];
    targetvalues[1] = sourcevalues[5];
    targetvalues[2] = sourcevalues[6];
    targetvalues[3] = sourcevalues[7];

    targetvalues[4] = sourcevalues[0];
    targetvalues[5] = sourcevalues[1];
    targetvalues[6] = sourcevalues[2];
    targetvalues[7] = sourcevalues[3];

    targetvalues[8] = sourcevalues[8];
    targetvalues[9] = sourcevalues[9];
    targetvalues[10] = sourcevalues[10];
    targetvalues[11] = sourcevalues[11];
  }
  else{
    for (int i = 0; i < 4; i++){
      targetvalues[i] = -sourcevalues[i];
    }
    targetvalues[4] = sourcevalues[4];
    targetvalues[5] = sourcevalues[5];
    targetvalues[8] = sourcevalues[8];
    targetvalues[9] = sourcevalues[9];
  }
};
