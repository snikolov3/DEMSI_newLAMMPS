Cantilever Lattice
==================

Description
-----------

The cantilever lattice test case consists of a bonded lattice of elements. The left most elements are fixed, while the right most elements feel a force in the negative y direction.

Purpose
-------

The cantilever lattice test case tests the bonded Hopkins contact model.

Sub-tests
---------

+---------------------+-------------------------------------------------------------------------------------------+
| Sub test name       | Description                                                                               |
+=====================+===========================================================================================+
| break               | Breaking stresses set so cantilever breaks                                                |
+---------------------+-------------------------------------------------------------------------------------------+
| nobreak             | Breaking stresses set so cantilever does not break and undergoes oscillations             |
+---------------------+-------------------------------------------------------------------------------------------+
| nobreak_equilibrium | As "break" but sufficiently long that oscillations damp and equilibrium position achieved |
+---------------------+-------------------------------------------------------------------------------------------+

Initial condition
-----------------

.. code::

   > ../../utils/visualization/make_particle_plot.py -g grid.nc -i particles_in.nc -o cantilever_lattice_initial_condition.png -t --removeticks

.. figure:: https://i.imgur.com/3lebj93.png
   :width: 50%

   **cantilever_lattice_initial_condition.png**: Cantilever lattice initial condition. Green elements are fixed, yellow elements feel a force in the negative y direction. All elements are bonded.

Results
-------

.. code::

   > ../../demsi config_nobreak.xml
   > ../../utils/visualization/make_particle_plot.py -g grid.nc -i ./output/particles_out.0001-01-01_01\:00\:00.nc -o cantilever_lattice_nobreak_1.png --removeticks

.. figure:: https://i.imgur.com/MEmRuMr.png
   :width: 50%

   **cantilever_lattice_nobreak_1.png**: No break cantilever after 1 hour

.. code::

   > ../../utils/visualization/make_particle_plot.py -g grid.nc -i ./output/particles_out.0002-01-01_01\:00\:00.nc -o cantilever_lattice_nobreak_2.png --removeticks

.. figure:: https://i.imgur.com/UxXjLiK.png
   :width: 50%

   **cantilever_lattice_nobreak_2.png**: No break cantilever after 2 hours
