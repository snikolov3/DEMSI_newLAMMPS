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

#include "mliap_model_quadratic.h"
#include "pair_mliap.h"
#include <cmath>
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define MAXLINE 1024
#define MAXWORD 3

/* ---------------------------------------------------------------------- */

MLIAPModelQuadratic::MLIAPModelQuadratic(LAMMPS* lmp, char* coefffilename) : 
  MLIAPModel(lmp, coefffilename)
{
  nonlinearflag = 1;
  ndescriptors = sqrt(2*nparams)-1;
}

/* ---------------------------------------------------------------------- */

MLIAPModelQuadratic::~MLIAPModelQuadratic(){}

/* ----------------------------------------------------------------------
   Calculate model gradients w.r.t descriptors for each atom dE(B_i)/dB_i
   ---------------------------------------------------------------------- */

void MLIAPModelQuadratic::gradient(PairMLIAP* pairmliap, NeighList* list, double **descriptors, double **beta, int eflag)
{
  int i;
  int *type = atom->type;

  for (int ii = 0; ii < list->inum; ii++) {
    i = list->ilist[ii];
    const int itype = type[i];
    const int ielem = pairmliap->map[itype];
    double* coeffi = coeffelem[ielem];

    for (int icoeff = 0; icoeff < ndescriptors; icoeff++)
      beta[ii][icoeff] = coeffi[icoeff+1];

    int k = ndescriptors+1;
    for (int icoeff = 0; icoeff < ndescriptors; icoeff++) {
      double bveci = descriptors[ii][icoeff];
      beta[ii][icoeff] += coeffi[k]*bveci;
      k++;
      for (int jcoeff = icoeff+1; jcoeff < ndescriptors; jcoeff++) {
        double bvecj = descriptors[ii][jcoeff];
        beta[ii][icoeff] += coeffi[k]*bvecj;
        beta[ii][jcoeff] += coeffi[k]*bveci;
        k++;
      }
    }

    // add in contributions to global and per-atom energy
    // this is optional and has no effect on force calculation
 
    if (eflag) {

      // energy of atom I

      double* coeffi = coeffelem[ielem];
      double etmp = coeffi[0];

      // E_i = beta.B_i + 0.5*B_i^t.alpha.B_i

      for (int icoeff = 0; icoeff < ndescriptors; icoeff++)
        etmp += coeffi[icoeff+1]*descriptors[ii][icoeff];

      // quadratic contributions

      int k = ndescriptors+1;
      for (int icoeff = 0; icoeff < ndescriptors; icoeff++) {
        double bveci = descriptors[ii][icoeff];
        etmp += 0.5*coeffi[k++]*bveci*bveci;
        for (int jcoeff = icoeff+1; jcoeff < ndescriptors; jcoeff++) {
          double bvecj = descriptors[ii][jcoeff];
          etmp += coeffi[k++]*bveci*bvecj;
        }
      }
      pairmliap->e_tally(i,etmp);
    }
  }
}

