Single Particle
===============

Description
-----------

The single particle test case consists of a single particle with an initial x velocity of 1 m/s and that experiences no forces. The particle moves with thsi fixed velocity in the x direction for the duration of the test case.

Purpose
-------

This test case tests that the time step is correctly translated from DEMSI to LAMMPS.

Sub-tests
---------

*None*

Initial condition
-----------------

.. code::

   > ../../utils/visualization/make_particle_plot.py -g grid.nc -i particles_in.nc -o single_particle_initial_condition.png --removeticks

.. figure:: https://i.imgur.com/SY9Ztf5.png
   :width: 50%

   **single_particle_initial_condition.png**: Single particle is initially at the left most boundary of the domain mid way up.

Results
-------

.. code::

   > ../../demsi config.xml
   > python check_testcase_output.py output/
   x position is correct at 10010.000000
   y position is correct at 500000.000000
