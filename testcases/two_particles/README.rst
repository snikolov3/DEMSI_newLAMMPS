Two Particles
=============

Description
-----------

The two particle test case consists of two elements initially in contact. One element is fixed while the motion of the other is prescribed. This relative motion produces forces in the contact model and potentially causes fracture of any bond between the particles. For most bonded test cases fracture of the bond should occur at 5000m relative motion. Exceptions are the rotate and compression force test cases, where the later breaks at 1250m relative motion.

Purpose
-------

The two particle test case tests the bonded and unbonded element contact model.

Sub-tests
---------

+----------------------------+-------------------------------------------------------------------------------------------+
| Sub test name              | Description                                                                               |
+============================+===========================================================================================+
| bonded_compression         | Bonded elements undergo pure compression                                                  |
+----------------------------+-------------------------------------------------------------------------------------------+
| bonded_compression_force   | Bonded elements undergo pure compression. Bond force is output.                           |
+----------------------------+-------------------------------------------------------------------------------------------+
| bonded_compression_shear   | Bonded elements undergo both compression and shear                                        |
+----------------------------+-------------------------------------------------------------------------------------------+
| bonded_rotate              | Bonded elements; one element is fixed, while the other rotates and translates.            |
+----------------------------+-------------------------------------------------------------------------------------------+
| bonded_shear               | Bonded elements undergo pure shear                                                        |
+----------------------------+-------------------------------------------------------------------------------------------+
| bonded_tension             | Bonded elements undergo pure tension                                                      |
+----------------------------+-------------------------------------------------------------------------------------------+
| bonded_tension_shear       | Bonded elements undergo both tension and shear                                            |
+----------------------------+-------------------------------------------------------------------------------------------+
| unbonded_compression_force | Unbonded elements undergo pure compression. Bond force is output.                         |
+----------------------------+-------------------------------------------------------------------------------------------+

Initial condition
-----------------

.. code::

   > ../../utils/visualization/make_particle_plot.py -g grid.nc -i particles_in_bonded.nc -o two_particles_bonded_initial_condition.png --removeticks

.. figure:: https://i.imgur.com/9KKKyKY.png
   :width: 50%

   **two_particles_bonded_initial_condition.png**: Two particle initial condition consists of two particles in contact.

Results
-------

.. code::

   > ../../demsi config_bonded_rotate.xml
   > ../../utils/visualization/make_particle_plot.py -g grid.nc -i ./output/particles_out.0001-01-01_00\:01\:00.nc -o two_particles_bonded_rotate.png --removeticks -b

.. figure:: https://i.imgur.com/fQQDJKI.png
   :width: 50%

   **two_particles_bonded_rotate.png**: The bonded rotate test case keeps one element fixed while rotating and translating the other. The bond consists of a midline (*black arrow*) and two sides attached to each element (*green*). The bond undergoes fracture from one end (*red line*). 

Analysis
--------

- **analyze_bond_failure.py**: This prints the relative bond motion (at either end of the bond) and whether the bond has completely failed. This script allows the user to verify the bonds in the above test cases fail at 5000m relative motion.

- **analyze_compressive_force.py**: For the two_bonded_particles_compression_force.xml and two_unbonded_particles_compression_force.xml test cases this script plots the relative element motion and resulting bond force as a function of time.
