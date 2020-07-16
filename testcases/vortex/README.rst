Vortex
======

Description
-----------

The vortex test case consists of a 500 km by 500 km square domain with non-periodic boundaries. There are 5000 circular elements arranged in a regular array in the upper half of the domain. A vortex wind field causes the elements to undergo a circular motion around the domain. The test case is based on *Flato (1993)* [#Flato93]_.

Purpose
-------

The vortex test case is a general purpose test case for testing general aspects of the DEMSI model.

Sub-tests
---------

+-------------------------+---------------------------------------------------------------------------------------+
| Sub test name           | Description                                                                           |
+=========================+=======================================================================================+
| fixed_forcing           | Fixed wind forcing, Ice fraction: 1, Ice thickness: 1m                                |
+-------------------------+---------------------------------------------------------------------------------------+
| fixed_forcing_random    | Fixed wind forcing, Ice fraction: 1, Ice thickness: Uniform distribution 0.8-1m       |
+-------------------------+---------------------------------------------------------------------------------------+
| fixed_forcing_stability | Fixed wind forcing, Ice fraction: 1, Ice thickness: Uniform distribution 0-1m         |
+-------------------------+---------------------------------------------------------------------------------------+
| varying_forcing         | Time varying wind forcing that reverses direction, Ice fraction: 1, Ice thickness: 1m |
+-------------------------+---------------------------------------------------------------------------------------+

Initial condition
-----------------

.. code::

   > ../../utils/visualization/make_particle_plot.py -g grid.nc -i particles_in.nc -o vortex_initial_condition.png --removeticks

.. figure:: https://i.imgur.com/v1FOEOn.png
   :width: 50%

   **vortex_initial_condition.png**: Vortex test case initial condition.

Results
-------

.. code::

   > ../../demsi config_fixed_forcing.xml
   > ../../utils/visualization/make_particle_plot.py -g grid.nc -i ./output/particles_out.0001-01-03_00.nc -o vortex_1.png --removeticks

.. figure:: https://i.imgur.com/ziOdIQo.png
   :width: 50%

   **vortex_1.png**: Fixed forcing vortex test case after 2 days.

.. code::

   > ../../utils/visualization/make_particle_plot.py -g grid.nc -i ./output/particles_out.0001-01-05_00.nc -o vortex_2.png --removeticks

.. figure:: https://i.imgur.com/hSNRiPR.png
   :width: 50%

   **vortex_2.png**: Fixed forcing vortex test case after 4 days.

References
----------

.. [#Flato93] **Gregory M. Flato (1993)** A particle‐in‐cell sea‐ice model, *Atmosphere-Ocean*, 31:3, 339-358, DOI: 10.1080/07055900.1993.9649475 (`link <https://doi.org/10.1080/07055900.1993.9649475>`_)
