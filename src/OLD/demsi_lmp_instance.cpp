/*
  Lightweight wrapper to lammps library. Possibly useful in managing
  multiple instances of lammps.
 */

#include "demsi_lmp_instance.h"

#include <mpi.h>
#include <math.h>
#include <lammps/lammps.h>         // these are LAMMPS include files
#include <lammps/input.h>
#include <lammps/atom.h>
#include <lammps/library.h>
#include <lammps/update.h>
#include <lammps/modify.h>
#include <lammps/neighbor.h>
#include <lammps/domain.h>
#include <lammps/comm.h>
#include <lammps/fix.h>

#include <Kokkos_Core.hpp>

namespace DEMSI {

LammpsInstance::LammpsInstance(DEMSI::Partition* partitionIn, const int nthreads, const int nnuma, const int ngpus) {
  partition = partitionIn;

  auto s_nthreads = std::to_string(nthreads);
  auto s_nnuma = std::to_string(nnuma);
  auto s_ngpus = std::to_string(ngpus);

#ifdef KOKKOS_HAVE_CUDA
  int narg = 17;
  char** arg = new char*[17];

  arg[0] = (char *) "lmp";
  arg[1] = (char *) "-k";
  arg[2] = (char *) "on";
  arg[3] = (char *) "threads";
  arg[4] = (char *) s_nthreads.c_str();
  arg[5] = (char *) "numa";
  arg[6] = (char *) s_nnuma.c_str();
  arg[7] = (char *) "gpus";
  arg[8] = (char *) s_ngpus.c_str();
  arg[9] = (char *) "-sf";
  arg[10] = (char *) "kk";
  arg[11] = (char *) "-pk";
  arg[12] = (char *) "kokkos";
  arg[13] = (char *) "neigh";
  arg[14] = (char *) "half";
  arg[15] = (char *) "gpu/direct";
  arg[16] = (char *) "off";
  lmp = new LAMMPS_NS::LAMMPS(narg, arg, partition->comm());
#else
  int narg = 13;
  char** arg = new char*[13];

  arg[0] = (char *) "lmp";
  arg[1] = (char *) "-k";
  arg[2] = (char *) "on";
  arg[3] = (char *) "threads";
  arg[4] = (char *) s_nthreads.c_str();
  arg[5] = (char *) "numa";
  arg[6] = (char *) s_nnuma.c_str();
  arg[7] = (char *) "-sf";
  arg[8] = (char *) "kk";
  arg[9] = (char *) "-pk";
  arg[10] = (char *) "kokkos";
  arg[11] = (char *) "neigh";
  arg[12] = (char *) "half";
  lmp = new LAMMPS_NS::LAMMPS(narg, arg, partition->comm());
#endif

}

void LammpsInstance::run_script(const std::string filename) {
  lammps_file(lmp, (char *) filename.c_str());
}

void LammpsInstance::run(int nsteps){
  char str[1024];
  sprintf(str, "run %d pre no post no", nsteps);
  //Depending on what we end up doing between runs,
  // pre and/or post may need to be set to 'yes'.
  // See LAMMPS doc page of run command for more info.
  lmp->input->one(str);
}

void LammpsInstance::reneighbor(){
  if (lmp->modify->n_pre_exchange) lmp->modify->pre_exchange();
  lmp->domain->pbc();
  if (lmp->domain->box_change){
    lmp->domain->reset_box();
    lmp->comm->setup();
    if (lmp->neighbor->style) lmp->neighbor->setup_bins();
  }
  lmp->comm->exchange();
  lmp->comm->borders();
  if (lmp->modify->n_pre_neighbor) lmp->modify->pre_neighbor();
  lmp->neighbor->build(1);
  if (lmp->modify->n_post_neighbor) lmp->modify->post_neighbor();
}

void LammpsInstance::comm_ghosts(){
  lmp->comm->exchange();
  lmp->comm->borders();
}

void LammpsInstance::run_subcycle(){
  this->run(nSteps);
}

void LammpsInstance::set_timestep(const double dynamicsTimeStep){
  // make dynamics timestep integer multiple of lammps timestep
  nSteps = ceil(dynamicsTimeStep / lmp->update->dt);
  lmp->update->dt = dynamicsTimeStep / (double) nSteps;
  // re-initialize fixes, since some of them set internal parameters based on timestep
  for (int i = 0; i < lmp->modify->nfix; ++i)
    lmp->modify->fix[i]->init();
}

void LammpsInstance::one(const std::string str) {
  lmp->input->one(str.c_str());
}

LammpsInstance::~LammpsInstance(){
  delete lmp;
}

} // namespace DEMSI
