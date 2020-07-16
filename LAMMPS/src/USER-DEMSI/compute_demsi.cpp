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

#include <string.h>
#include "compute_demsi.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeDEMSIAtom::ComputeDEMSIAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  demsi_info(NULL)
{
  if (narg != 3) error->all(FLERR,"Illegal compute demsi/atom command");

  peratom_flag = 1;
  size_peratom_cols = 4;

  nmax = 0;
}

/* ---------------------------------------------------------------------- */

ComputeDEMSIAtom::~ComputeDEMSIAtom()
{
  memory->destroy(demsi_info);
}

/* ---------------------------------------------------------------------- */

void ComputeDEMSIAtom::init()
{
  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"demsi/atom") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute demsi/atom");
}

/* ---------------------------------------------------------------------- */

void ComputeDEMSIAtom::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow demsi info array if necessary

  if (atom->nmax > nmax) {
    memory->destroy(demsi_info);
    nmax = atom->nmax;
    memory->create(demsi_info,nmax,4,"demsi_info/atom:demsi_info");
    array_atom = demsi_info;
  }

  // compute demsi info for each atom in group

  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
	  demsi_info[i][0] = atom->forcing[i][0];
	  demsi_info[i][1] = atom->forcing[i][1];
	  demsi_info[i][2] = atom->mean_thickness[i];
	  demsi_info[i][3] = atom->min_thickness[i];
      }
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeDEMSIAtom::memory_usage()
{
  double bytes = nmax * sizeof(double) * 4;
  return bytes;
}
