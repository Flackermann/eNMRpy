"""
The package eNMRpy is a tool for the import and analysis of electrophoretic NMR spectroscopy raw data.
The package is mainly based on nmrglue and lmfit other common python imports.

A peer reviewed paper showing the performance and advantages of this package is available under:
Schmidt, Florian, et al. "Spectral deconvolution in electrophoretic NMR to investigate the migration of
neutral molecules in electrolytes." Magnetic Resonance in Chemistry 58.3 (2020): 271-279.

The import classes are grouped in the submodule eNMRpy.Measurement

While the phase analysis via zero order phase correction is included in the import classes,
phase analysis via lorentzian fits is available under eNMRpy.Phasefitting, and via 2D Fast Fourier Transformation
under eNMRpy.MOSY
"""

#from .Phasefitting import SpecModel
#from .MOSY import MOSY
# Pavel Class as the Std Import-Class since it is mostly used in our working group
from .Measurement import Pavel as Import_eNMR_Measurement 
from .Meausrement import Flo as Import_eNMR_Measurement_12Bit
from . import tools
#print('%s imported'%__name__)

