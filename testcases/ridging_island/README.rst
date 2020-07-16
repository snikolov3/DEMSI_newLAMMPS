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

   > python ../../utils/initialize_particles/lloyds_disc_packing/lloyds_disc_packing.py -o particles_in_voro.nc -v ${DEMSI_VORO_EXE} --x0 0.0 --x1 200000.0 --y0 0.0 --y1 200000.0 --r0 24000.0 --r1 16000.0

   > python ../../utils/initialize_particles/lloyds_disc_packing/expand_particle_distribution.py -i particles_in_voro.nc -o particles_in_expand.nc --dx 1000000 --dy 1000000

With the input particle distribution created we can create the test case input data:

.. code::

   > python make_testcase.py -i particles_in_expand.nc

.. code::

   > python ../../utils/visualization/make_particle_plot.py -g grid.nc -i particles_in_init.nc -o ridging_island_initial_condition.png --removeticks -t
