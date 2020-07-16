Arctic Basin
============

Description
-----------

Arctic basin test case simulates the Arctic basin and surrounding seas.

Purpose
-------

This test case tests the infrastructure for simulating the Arctic basin.

Sub-tests
---------

*None*

Initial condition
-----------------

To build the input particle distribution use a variant of Lloyds algorithm:

.. code::

   > python -o particles_in_voro.nc ../../utils/initialize_particles/lloyds_disc_packing/lloyds_disc_packing.py -v ${DEMSI_VORO_EXE} --x0 0.0 --x1 1000000.0 --y0 0.0 --y1 1000000.0 --r0 20000.0 --r1 30000.0

   > python ../../utils/initialize_particles/lloyds_disc_packing/expand_particle_distribution.py -i particles_in_voro.nc -o particles_in_expand.nc --dx 11000000 --dy 11000000

With the input particle distribution created we can create the test case input data:

.. code::

   > python make_testcase_init.py -i particles_in_expand.nc -w 11000000.0

   > python make_testcase_forcing.py

Create an image of the test case initial condition:

.. code::

   > python ../../utils/visualization/make_particle_plot.py -g grid.nc -i particles_in.nc -o arctic_basin_initial_condition.png --removeticks -t

.. figure:: https://i.imgur.com/KRPGbcU.png
   :width: 80%

   **arctic_basin_initial_condition.png**: Initial condition for Arctic Basin test case. Yellow elements show coastlines, purple elements are sea ice.

Results
-------

Run the simulation:

.. code::

   > ../../demsi config.xml

Create an image of the output:

.. code::

   > python ../../utils/visualization/make_particle_plot.py -g grid.nc -i output/particles_arctic_basin.2000-01-13_00.nc -o arctic_basin.png -t --removeticks

.. figure:: https://i.imgur.com/0p8shLn.png
   :width: 80%

   **arctic_basin.png**: Arctic Basin elements after 12 days of simulation.

Create a movie of the simulation:

.. code::

   > mkdir tmp

   > python ../../utils/visualization/make_particle_movie.py -g grid.nc -f "./output/particles_arctic_basin.*" -o arctic.mov -v iceAreaCell
