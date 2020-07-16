Model prequisites
=================

Compiling DEMSI
---------------

The DEMSI code is written in C++ and Fortran. The following is needed
to compile DEMSI:

- Git
- CMake
- A modern C++ compiler
- A compatable Fortran compiler
- An MPI installation
- `The Netcdf library <https://www.unidata.ucar.edu/software/netcdf/>`_

Supported versions:

+------------------+--------------------+
| Software         | Supported versions |
+==================+====================+
| C++ Compiler     | gcc7               |
+------------------+--------------------+
| Fortran Compiler | gcc7               |
+------------------+--------------------+
| MPI              | Open MPI 3.0.0     |
+------------------+--------------------+
| Netcdf           | 4.1.2              |
+------------------+--------------------+

The following environment variables need to be set:

- The ``NETCDF`` environment variable needs to be set to the location of the installed Netcdf library.

Pre- and Post-processing
------------------------

Pre- and Post-processing is performed with python. The DEMSI project
supports anaconda with the following environment:

.. code::

   > conda create -n pydemsi python=3.5
   > conda install -n pydemsi six numpy netCDF4 matplotlib scipy
   > conda install -n pydemsi -c anaconda sphinx
   > conda install -n pydemsi -c conda-forge doxygen breathe shapely


To activate and deactivate the DEMSI python environment:

.. code::

   > conda activate pydemsi
   > conda deactivate pydemsi

Some test cases also use the `Voro++ <http://math.lbl.gov/voro++/>`_ library to calculate Voronoi tessellations. To use these test cases, download the Voro++ library and build the executable. Then set the ``DEMSI_VORO_EXE`` environment variable to the Voro++ executable location.
