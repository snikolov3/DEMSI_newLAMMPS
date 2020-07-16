.. raw:: html 

   <p align="center">
     <img width="70%" src="docs/source/DEMSI.svg">
   </p>

DEMSI is a Discrete Element Method (DEM) sea-ice model designed for
global simulations with Earth system models.


Acknowledgements
----------------

Support for DEMSI was provided through the Scientific Discovery
through Advanced Computing (SciDAC) program funded by the US
Department of Energy (DOE), Office of Science, Biological and
Environmental Research, and Advanced Scientific Computing Research
programs.

Model Documentation
-------------------

DEMSI includes documentation which is generated with
the `Sphinx <http://www.sphinx-doc.org/en/master/>`_ tool. Class
documentation is generated with `doxygen <http://www.doxygen.org/>`_
and linked to Sphink with the `breathe
<https://breathe.readthedocs.io/en/latest/>`_ tool. To generate the
documentation:

.. code::

   > cd DEMSI/docs

Run doxygen:

.. code::

   > doxygen

Build the HTML documentation with Sphinx:

.. code::

   > make html

Point a web browser at ``DEMSI/docs/build/html/index.html`` to view.
