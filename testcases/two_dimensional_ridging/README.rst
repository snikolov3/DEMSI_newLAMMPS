Two Dimensional Convergence
===========================

Description
-----------

Purpose
-------

Sub-tests
---------

*None*


Initial condition
-----------------

To build the input particle distribution:

.. code::

   > python ../../utils/initialize_particles/lloyds_disc_packing/lloyds_disc_packing.py -v ${DEMSI_VORO_EXE} --x0 0.0 --x1 100000.0 --y0 0.0 --y1 100000.0 --r0 6000.0 --r1 4000.0
   > python ../../utils/initialize_particles/lloyds_disc_packing/lloyds_disc_packing.py -v ${DEMSI_VORO_EXE} --x0 0.0 --x1 100000.0 --y0 0.0 --y1 100000.0 --r0 11000.0 --r1 9000.0  

   > python ../../utils/initialize_particles/lloyds_disc_packing/expand_particle_distribution.py -i particles_in.nc -o particles_in_expand.nc --dx 1030000 --dy 1000000
   
With the input particle distribution created we can create the test case input data:

.. code::

   > python make_testcase.py
