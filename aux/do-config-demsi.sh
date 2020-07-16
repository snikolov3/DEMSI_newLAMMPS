#!/bin/bash -ex
cmake $HOME/src/DEMSI \
-DCMAKE_INSTALL_PREFIX:PATH=$HOME/install/gcc/DEMSI \
-DCMAKE_CXX_COMPILER:FILEPATH=$HOME/install/gcc/mpich/bin/mpicxx \
-DCMAKE_C_COMPILER:FILEPATH=$HOME/install/gcc/mpich/bin/mpicc \
-DnetCDF_PREFIX:PATH=$HOME/install/gcc/netcdf-demsi \
2>&1 | tee config_log
