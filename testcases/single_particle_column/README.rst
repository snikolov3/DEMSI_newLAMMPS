Single Particle Column
======================

Description
-----------

This test case consists of a single fixed particle, subject to atmospheric and oceanic forcing from the MPAS-Seaice single cell test case integrated for a single year.

Purpose
-------

The single particle column test case is used to test the column physics implementation in the model.

Sub-tests
---------

*None*

Initial condition
-----------------

.. code::

   > ../../utils/visualization/make_particle_plot.py -g grid.nc -i particles_in.nc -o single_particle_column_initial_condition.png --removeticks

.. figure:: https://i.imgur.com/Tmluzzi.png
   :width: 50%

   **single_particle_column_initial_condition.png**: The test case consists of a single fixed element in a domain that perfectly fits the element.

Results
-------

.. code::

   > ../../demsi config.xml
   > python plot_column_output.py

.. figure:: https://i.imgur.com/hRwQqBZ.png
   :width: 50%

   **single_particle_column_output.png**: Element ice (*red*) and snow (*blue*) thicknesses and surface temperature (*green*) over time for a single year.