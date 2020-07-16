#ifndef DEMSI_LAMMPS_INSTANCE_H
#define DEMSI_LAMMPS_INSTANCE_H

#include <string>
#include <mpi.h>
#include <lammps/lammps.h>
#include "demsi_partition.h"

/*! \file   demsi_lmp_instance.h
    \brief  Header file for the DEMSI::LammpsInstance classes
*/

namespace DEMSI {

//------------------------------------------------------------------------------
/*! \class LammpsInstance
    \brief Wrapper class for LAMMPS object
*/
//------------------------------------------------------------------------------
class LammpsInstance{
public:

  /*! \brief LammpsInstance class constructor.
      \param partitionIn Pointer to partition object.
      \param nthreads    number of threads to execute on
      \param nnuma       number of numa regions
      \param ngpus       number of gpus available
   */
  LammpsInstance(DEMSI::Partition* partitionIn, const int nthreads, const int nnuma=1, const int ngpus=0);

  /*! \brief LammpsInstance class destructor.*/
  ~LammpsInstance();

  /*! \brief Run a LAMMPS input script.
      \param filename
  */
  void run_script(const std::string filename);

  /*! \brief Run LAMMPS for a number of time steps.
      \param nsteps Number of timesteps to run LAMMPS for.
  */
  void run(int nsteps);

  /*! \brief Initialize bonded contacts.
  */
  void initialize_bonded_contacts();

  /*! \brief Run LAMMPS for nSteps number of time steps.
  */
  void run_subcycle();

  /*! \brief Set the LAMMPS timestep and number of timesteps from the dynamics subcycle timestep.
      \param dynamicsTimeStep
  */
  void set_timestep(const double dynamicsTimeStep);

  /*! \brief Send one command to lammps
      \param str command to send
  */
  void one(const std::string str);

  /*! \brief Carry out LAMMPS reneighboring, including remapping based on periodic boundary conditions
  */
  void reneighbor();

  /*! \brief Create ghost particles for parallel runs, communicate relevant per-particle data to them
  */
  void comm_ghosts();

  /*! Pointer to LAMMPS object. */
  LAMMPS_NS::LAMMPS *lmp;

private:

  /*! Pointer to partition object. */
  DEMSI::Partition* partition;

  /*! Number of LAMMPS timesteps per dynamics subcycle.*/
  int nSteps;

};

} // namespace DEMSI

#endif
