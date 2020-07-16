Quick Start Guide
=================

Clone the DEMSI repository:

.. code::

   > git clone --recurse-submodules git@github.com:E3SM-Project/DEMSI.git

Go to the main DEMSI directory:

.. code::

   > cd DEMSI

Build the DEMSI model:

.. code::

   > build

With the model built, go to the vortex test case directory:

.. code::

   > cd testcases/vortex

Run the DEMSI model:

.. code::

   > ../../demsi config_fixed_forcing.xml

Create a movie of the vortex test case model output:

.. code::

   > ../../utils/visualization/make_particle_movie.py -g grid -f "./output/particles*" -o vortex.mov