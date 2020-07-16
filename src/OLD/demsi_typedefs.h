/*
  DEMSI typedefs

 */

#ifndef DEMSI_TYPEDEFS_H
#define DEMSI_TYPEDEFS_H

#include <Kokkos_Core.hpp>
#include <lammps/atom_kokkos.h>

namespace DEMSI {

typedef Kokkos::DefaultExecutionSpace device_execution_space;
typedef Kokkos::DefaultHostExecutionSpace host_execution_space;
typedef device_execution_space::memory_space device_memory_space;
typedef host_execution_space::memory_space host_memory_space;

typedef DAT::tdual_tagint_1d::t_dev kokkos_view_type_tagint;
typedef DAT::tdual_tagint_1d::t_host kokkos_view_type_tagint_host;
typedef Kokkos::
  DualView<int*, LMPDeviceType::array_layout, LMPDeviceType> tdual_int_1d;
typedef tdual_int_1d::t_host kokkos_view_type_int_host;
typedef DAT::tdual_x_array::t_dev kokkos_view_type_x;
typedef DAT::tdual_v_array::t_dev kokkos_view_type_v;
typedef DAT::tdual_float_1d::t_dev kokkos_view_type_1d_float;
typedef DAT::tdual_float_2d::t_dev kokkos_view_type_2d_float;

} // DEMSI namespace

#endif
