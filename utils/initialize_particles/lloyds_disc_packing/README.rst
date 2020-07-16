===============================================================
Initial particle distribution generation with Lloyd's algorithm
===============================================================

Overview
========

These scripts generate tightly pack but not overlapping initial particle distributions for DEMSI. They use a variant of Lloyds algorithm where the generator points of the radical Voronoi tessellation are updated as the centre of the largest inscribed circle of the tessellated polygons during the Lloyd's iteration.

Scripts
=======

* **lloyds_disc_packing.py**: This script generates a double periodic initial particle distribution with a range of radii.

* **particle_distribution_quality.py**: This script evaluates the quality of a particle distribution.

* **expand_particle_distribution.py**: This script duplicates a periodic input particle distribution until a desired domain size is filled.

Example
=======

Generate a periodic particle distribution:

.. code::

   > python lloyds_disc_packing.py -v /Path/To/voro++ --x0 0.0 --x1 1000000.0 --y0 0.0 --y1 1000000.0 --r0 20000.0 --r1 30000.0

The script generates a plot showing the generated particle distribution with radical Voronoi tessellation:

.. figure:: https://i.imgur.com/Ze7wDQl.png
   :width: 80%

   **particles_lloyds.png**: Generated particle distribution with radical Voronoi tessellation


Duplicate the previously generated distribution to fill a desired domain size:

.. code::

   > python expand_particle_distribution.py -i particles_in.nc -o particles_in_expand.nc --dx 11000000 --dy 11000000

Plot the expanded particle distribution:

.. code::

   > ../../visualization/make_particle_plot.py --x0 0 --x1 11000000 --y0 0 --y1 11000000 -i particles_in_expand.nc -o particles_in_expand.png

.. figure:: https://i.imgur.com/Tkcf8yO.png
   :width: 80%

   **particles_in_expand.png**: Duplicated initial particle distribution.
