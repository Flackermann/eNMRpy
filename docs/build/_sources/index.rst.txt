.. eNMRpy documentation master file, created by
   sphinx-quickstart on Mon Aug  3 13:43:41 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. Welcome to eNMRpy's documentation!

.. _intropage:

============
About eNMRpy
============

eNMRpy is a tool for the import and analysis of electrophoretic NMR (eNMR) Data. It was created for the purpose of automatizing routine analysis procedures, importing meta parameters from measurement folders, while providing a framework for the analysis NMR signal phase shifts in a 2D eNMR experiment.

eNMRpy is based on `nmrglue <https://github.com/jjhelmus/nmrglue>`_ and `lmfit <https://github.com/lmfit/lmfit-py>`_, two exceptional python packages for importing and processing NMR-data, and nonlinear regression.


The submodule Measurement contains all import-related classes while the main module contains the code for the analysis of the respective imported measurement.
If help es needed in order to adapt the module to a specific experimental setup, please consider a `pull-request on GitHub <https://github.com/Flackermann/eNMRpy>`_.


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   /eNMR basics
   /Tutorial
   /pkginfo   


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
