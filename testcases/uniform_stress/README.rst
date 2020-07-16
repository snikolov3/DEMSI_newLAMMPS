Uniform stress
==============

Description
-----------

This test case simulates a square region of bonded elements undergoing
uniform stress. The test case is based on *Hopkins, M. A., S. Frankenstein, and A. S. Thorndike (2004)* [#Hopkins04]_.

Purpose
-------

This test case demonstrates the bonded contact model fracture.

Sub-tests
---------

*None*

Initial condition
-----------------

To build the initial condition using the standalone LAMMPS executable that is built with DEMSI:

Build the initial particle placement code:

.. code::

   > cd DEMSI/utils/initialize_particles/place_nonoverlapping

   > make

Create an initial particle distribution:

.. code::

   > python ../../utils/initialize_particles/lloyds_disc_packing/lloyds_disc_packing.py -o particles_in_voro.nc -v ${DEMSI_VORO_EXE} --x0 0.0 --x1 100000.0 --y0 0.0 --y1 100000.0 --r0 3000.0 --r1 5000.0

   > python ../../utils/initialize_particles/lloyds_disc_packing/expand_particle_distribution.py -i particles_in_voro.nc -o particles_in_expand.nc --dx 1000000 --dy 1000000

Create the test case input data:

.. code::

   > python make_testcase.py

Produce a figure of the initial condition:

.. code::

   > python ../../utils/visualization/make_particle_plot.py -g grid.nc -i particles_in.nc -o uniform_stress_initial_condition.png --removeticks -t

.. figure:: https://i.imgur.com/NyERn1x.png
   :width: 80%

   **uniform_stress_initial_condition.png**: Initial condition for
   uniform stress test case.

Results
-------

Run the test case:

.. code::

   > ../../demsi config.xml

Produce a figure of the final state:

.. code::

    > python ../../utils/visualization/make_particle_plot.py -g grid.nc -i ./output/particles_out.0001-01-01_04:10:00.nc -o uniform_stress.png --removeticks

.. figure:: https://i.imgur.com/YDdSBpw.png
   :width: 80%

   **uniform_stress.png**: Uniform stress test case after 15,000 seconds.

Other analysis scripts:

.. code::

   > python plot_bond_timeseries.py

   > python plot_bonds.py

   > python plot_connected_region_timeseries.py

   > python plot_connected_regions.py

References
----------

.. [#Hopkins04] **Hopkins, M. A., S. Frankenstein, and A. S. Thorndike (2004)** Formation of an aggregate scale in Arctic sea ice, *J. Geophys.Res.*, 109, C01032, DOI:10.1029/2003JC001855 (`link <https://doi.org/10.1029/2003JC001855>`_)
