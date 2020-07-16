Developer Guide
===============

Testing the model
-----------------

DEMSI includes a ``ctest`` testing system to test the model using the
test cases defined within the model. Several test types are defined in
``DEMSI/utils/testing/``:

+------------+-----------------------------------------------------------------------------------+
| Test name  | Description                                                                       |
+============+===================================================================================+
| Smoke      | Tests that the model completes successfully for a given test case                 |
+------------+-----------------------------------------------------------------------------------+
| Regression | Tests the model produces the same results as a baseline run for a given test case |
+------------+-----------------------------------------------------------------------------------+
| Check      | Runs the model and performs a specified check of the model output                 |
+------------+-----------------------------------------------------------------------------------+

The DEMSI test system is invoked from the main DEMSI directory with:

.. code::

   > ctest

Output logs of the testing are written to ``DEMSI/Testing/Temporary/``

To run an individual test (e.g. test 3):

.. code::

   > ctest -I 3,3

To run a range of tests (e.g. 5 to 8):

.. code::

   > ctest -I 5,8

The DEMSI regression tests require baseline simulations to compare against. To generate baseline results, clone and build DEMSI and run the tests. This generates the required baseline results. To allow other check outs of DEMSI to use these baselines for comparison during regression testing, set the ``DEMSI_TEST_BASELINE`` environment variable to the ``DEMSI`` directory of the baseline simulations.




Profiling the model
-------------------

DEMSI includes Kokkos profiling regions.

To add a region for profiling, wrap the desired code region with these
commands:

.. code::

   Kokkos::Profiling::pushRegion("Region name");
   ... Code to profile ...
   Kokkos::Profiling::popRegion();


To generate output from the profiling regions, a Kokkos tool needs to be installed:

.. code::

   > cd ${TOOL_LOCATION}
   > git clone git@github.com:kokkos/kokkos-tools.git
   > cd kokkos-tools/src/tools/space-time-stack
   > make


Then set the ``KOKKOS_PROFILE_LIBRARY`` environment variable:

.. code::

   > export KOKKOS_PROFILE_LIBRARY=${TOOL_LOCATION}/kokkos-tools/src/tools/space-time-stack/kp_space_time_stack.so


Once DEMSI is run, profiling information should be send to ``stdout``.



Git workflow
------------

Developers should create a fork of the DEMSI and DEMSI-LAMMPS repositories in github.

**Development in DEMSI only:**

Create a feature branch from your DEMSI fork to do the development work in:

.. code::

   > git clone -o <git_username> git@github.com:<git_username>/DEMSI.git
   > git remote add E3SM-Project git@github.com:E3SM-Project/DEMSI.git
   > git fetch --all
   > git checkout -b <development_branch_name> E3SM-Project/master
   > git submodule update --init --recursive

Then perform development and testing of the code. Once the development is complete the changes must be committed and pushed to the developers fork on github:

.. code::

   > git add <files to commit>
   > git commit
   > git push <git_username> <development_branch_name>

Finally, create a pull request for this feature branch in github. The project maintainers will then satisfy this pull request.
