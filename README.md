# Info
This package is a nifty tool for the import and analysis of electrophoretic NMR(eNMR)-spectroscopic data. It is maintained by Florian Ackermann (nairolf.ackermann@gmail.com)

For contributions and questions, please consider the <a href="https://github.com/Flackermann/eNMRpy">GitHub repository</a>.
When using this for scientific pusposes, please cite <a href="https://doi.org/10.1002/mrc.4978">this paper</a>.

For documentation please read <a href="https://enmrpy.readthedocs.io/en/testbranch/"> the read the docs page</a>.


# Installation
Install the latest release simply via <code>$ pip install eNMRpy</code>

# Development
Install the package in editable mode with the test dependencies and run the test suite:

```
pip install -e ".[test]"
pytest
```

or with uv: <code>$ uv run --extra test pytest</code>

# Range of functions covered

- Import of Bruker-based eNMR-Data
  	- 3 different experimental Setups so far
  	- please consider **pull requests** on  <a href="https://github.com/Flackermann/eNMRpy">GitHub</a> to get help for your own experimental setup

- Phase angle analysis
    - Phase correction analysis (old approach)
        - Entropy minimization
        - Spectra matching
        - Comparison of the phase-corrected spectra by overlaying them

    - Phase analysis by fitting (new approach)
        - Lorentz/Voigt peaks
            - superposition of any number of peaks
            - individual fixing of parameters
    
    - Regression
        - Calculation of the respective mobilities from automatically determined experimental parameters
    
    - Comparison of the different results
        (- simple tool for creating graphs)

- Phase analysis via 2D FFT --> Mobility Ordered Spectroscopy (MOSY)
    - States-Haberkorn method
    - Determination of the mobility axis
    - Plotting of slices to compare the results/peaks
    - Automatic normalization of the intensities and detection of the maxima
    - Note: signal loss with increasing voltage leads to an overestimated mobility!
