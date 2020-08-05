.. _Tutorialpage:

========
Tutorial
========

Importing and processing the eNMR spectrum
******************************************

loading the data
----------------
Loading the eNMR data takes place by creating an instance (object) of the respective import class with the respective path and expno of the experiment. This is the most critical point when adapting a different experimental setup to eNMRpy. If you are in need for help, please consider a `pull request on GitHub <https://github.com/Flackermann/eNMRpy>`_.

.. code-block::

	>>> import eNMRpy
	>>> m = eNMRpy.Import_eNMR_Measurement(self, path, expno, dependency='U', lineb=.5, d=2.2e-2, cell_resistance=None)
	

Important instance variables
----------------------------

.. table:: important instance variables of the .Pavel()-class *(atm standard import class)*

    +-------------------------+---------------------------------------------------------------------------+
    | variable                |               meaning                                                     |
    +=========================+===========================================================================+
    | self.eNMRraw            | contains the phase data, voltage list, and gradient list                  |
    |                         | data for all eNMR experiments                                             |
    +-------------------------+---------------------------------------------------------------------------+
    | self.d                  | electrode distance in m (cell constant)                                   |
    +-------------------------+---------------------------------------------------------------------------+
    | self.Delta              | observation time in s                                                     |
    +-------------------------+---------------------------------------------------------------------------+
    | self.delta              | magnetic field gradient pulse duration                                    |
    +-------------------------+---------------------------------------------------------------------------+
    | self.gamma              | gyromagnetic ratio in °/Ts                                                |
    +-------------------------+---------------------------------------------------------------------------+
    | self.path               | path used to import the data                                              |
    +-------------------------+---------------------------------------------------------------------------+
    | self.expno              | experiment number according to Bruker's syntax                            |
    +-------------------------+---------------------------------------------------------------------------+
    | self.lineb              | line broadening variable used for processing                              |
    +-------------------------+---------------------------------------------------------------------------+
    | self.ppm                | list of chemical shift for each datapoint                                 |
    +-------------------------+---------------------------------------------------------------------------+
    | self.data               | 2D numpy array containing the complex spectral data                       |
    +-------------------------+---------------------------------------------------------------------------+
    | self.nuc                | investigated nucleus                                                      |
    +-------------------------+---------------------------------------------------------------------------+
    | self.title_page         | imported title page from bruker file -> transformed to be used with print |
    +-------------------------+---------------------------------------------------------------------------+
    | self.dependency         | voltage, gradient, or current dependent measurement --> see docstring     |
    +-------------------------+---------------------------------------------------------------------------+
    | self.cell_resistance    | important to calculate the electric field when const. *I* mode is used    |
    +-------------------------+---------------------------------------------------------------------------+




Processing the spectrum (FFT, linebroadening)
---------------------------------------------

Usually, the imported data represents the free induction decay (FID), recorded for the eNMR measurements.
In order to analyze the eNMR measurement, the FID needs to be transformed into a spectrum via Fast Fourier
Transformation (FFT) via the .proc() method




Plotting the spectrum
---------------------



Phase analysis
**************



automatic phase correction
--------------------------

Approximation of Lorentzian-shaped function
-------------------------------------------

MOSY
----


.. toctree::
   :maxdepth: 2
   :caption: Contents:
