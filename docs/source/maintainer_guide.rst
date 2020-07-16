Maintainer Guide
================

Updating the LAMMPS hash in DEMSI
---------------------------------

Once the hash of DEMSI-LAMMPS has changed this needs to be updated in DEMSI. This is done with creating a feature branch:

.. code::

   > git clone --recursive -o <git_username> git@github.com:<git_username>/DEMSI.git
   > git remote add E3SM-Project git@github.com:E3SM-Project/DEMSI.git
   > git fetch --all
   > git checkout -b <development_branch_name> E3SM-Project/master
   > git submodule update --init --recursive

Update the LAMMPS hash

.. code::

   > git submodule update --remote --merge LAMMPS

Perform any necessary testing, and then commit this change and push back to the repo.

.. code::

   > git add LAMMPS
   > git commit
   > git push <git_username> <development_branch_name>

Then create a PR and satisfy it.
