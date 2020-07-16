User Guide
==========

Getting DEMSI
-------------

DEMSI uses github as a code repository and git as its versioning tool. To get DEMSI:

.. code::

   > git clone git@github.com:E3SM-Project/DEMSI.git
   > cd DEMSI
   > git submodule update --init



Building DEMSI
--------------

To build the model, from the upper most DEMSI directory:

.. code::

   > build

The DEMSI build system supports parallel builds. To build DEMSI with 4 threads:

.. code::

   > build -j 4

To clean the build (only the DEMSI source code, not libraries):

.. code::

   > build clean

To compile the DEMSI code in debug mode:

.. code::

   > build -d

In addition to building the DEMSI model, the build system also
prepares many of the test cases.

Initialization of some test cases requires the LAMMPS executable to be
built. To build the LAMMPS executable:

.. code::

   > cd DEMSI/LAMMPS/src
   > make yes-granular
   > make mpi



Running DEMSI
-------------

To run demsi a configuration file, typically called ``configs.xml``, is required. The configuration filename must be given as an argument to the DEMSI command:

.. code::

   > demsi configs.xml

Other input files might be required for a particular configuration of
DEMSI to run:

+---------------------+----------------------------------------------------+
| Filename            | Description                                        |
+=====================+====================================================+
| ``grid.nc``         | Netcdf file defining the model domain              |
+---------------------+----------------------------------------------------+
| ``particles_in.nc`` | Netcdf file defining initial particle distribution |
+---------------------+----------------------------------------------------+
| ``forcing.nc``      | Netcdf file defining model forcing                 |
+---------------------+----------------------------------------------------+

DEMSI contains pre-defined test cases as examples of how to run DEMSI in various configurations.



Running test cases
------------------

The DEMSI repository contains several test cases that demonstrate
various functionalities of DEMSI (see :ref:`test-cases-label`). These are found at
``DEMSI/testcases/``. Test cases contain the following files:

+--------------------------+---------------------------------------------------------------+
| Filename                 | Description                                                   |
+==========================+===============================================================+
| ``CMakeLists.txt``       | Cmake file that allows building of the test case during build |
+--------------------------+---------------------------------------------------------------+
| ``config.xml``           | DEMSI configuration file                                      |
+--------------------------+---------------------------------------------------------------+
| ``config_xxx.xml``       | DEMSI configuration file for sub test `xxx`                   |
+--------------------------+---------------------------------------------------------------+
| ``make_testcase.py``     | Test case script that generates other files needed by DEMSI   |
+--------------------------+---------------------------------------------------------------+
| ``make_testcase_yyy.py`` | Additional test case generation scripts                       |
+--------------------------+---------------------------------------------------------------+

Test cases may contain other test case specific files.

To run a test case first ``cd`` to the test case directory and then create the needed test case files:

.. code::

   > python make_testcase.py

If it hasn't automatically been done, create the appropriate output directory:

.. code::

   > mkdir output

Then run DEMSI:

.. code::

   > ../../demsi config.xml



Getting test case data
----------------------

Several test cases require external binary data to run. Binary data is not stored in the DEMSI git repository so DEMSI includes a tool to download the correct data files from public sources. It is not recommended that data be downloaded fresh with every DEMSI check-out. Instead we recommend DEMSI be checked out once specifically as a data store.

.. code::

   > cd ${DEMSI_DATA_STORE}
   > git clone git@github.com:E3SM-Project/DEMSI.git
   > cd DEMSI/DATA

The ``DATA`` directory contains a python script, ``get_data.py`` that uses ``wget`` to download the appropriate data:

.. code::

   > python get_data.py

With the data downloaded set the ``DEMSI_DATA_DIR`` environment variable so test cases in future check-outs of DEMSI know where to find the test case data.

.. code::

   > export DEMSI_DATA_DIR=${DEMSI_DATA_STORE}/DEMSI/DATA/



Visualization
-------------

DEMSI includes python scripts for easy viewing of model output. These
scripts can be found at ``DEMSI/utils/visualization/``:

+----------------------------+----------------------------------------------------------+
| Filename                   | Description                                              |
+============================+==========================================================+
| ``make_particle_plot.py``  | Plots the elements in a particles netcdf file.           |
+----------------------------+----------------------------------------------------------+
| ``make_particle_movie.py`` | Creates a movie from a series of particles netcdf files. |
+----------------------------+----------------------------------------------------------+
| ``make_grid_plot.py``      | Plots quantities on the DEMSI Eulerian grid.             |
+----------------------------+----------------------------------------------------------+

Run the above files with a ``--help`` argument to see all options available for running them.



Model Documentation
-------------------

The DEMSI model includes this documentation which is generated with the `Sphinx <http://www.sphinx-doc.org/en/master/>`_ tool. Class documentation is generated with `doxygen <http://www.doxygen.org/>`_ and linked to Sphink with the `breathe <https://breathe.readthedocs.io/en/latest/>`_ tool. To generate the documentation:

.. code::

   > cd DEMSI/docs

Run doxygen:

.. code::

   > doxygen

Build the HTML documentation with Sphinx:

.. code::

   > make html

Point a web browser at ``DEMSI/docs/build/html/index.html`` to view.
