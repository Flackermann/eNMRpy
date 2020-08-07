# Info
This package is a nifty tool for the import and analysis of electrophoretic NMR(eNMR)-spectroscopic data. It is maintained by Florian Ackermann (nairolf.ackermann@gmail.com)

For contributions and questions, please consider the <a href="https://github.com/Flackermann/eNMRpy">GitHub repository</a>.
When using this for scientific pusposes, please cite <a href="https://doi.org/10.1002/mrc.4978">this paper</a>.

For documentation please read <a href="https://enmrpy.readthedocs.io/en/testbranch/"> the read the docs page</a>.


# Installation
Install the latest release simply via <code>$ pip install eNMRpy</code>

# Range of functions covered

- Import of Bruker-based eNMR-Data
  	- 3 different experimental Setups so far
  	- please consider **pull requests** on  <a href="https://github.com/Flackermann/eNMRpy">GitHub</a> to get help for your own experimental setup

- Phasenwinkelanalyse
    - Phasenkorrektur-Analyse (old approach)
        - Entropieminimierung
        - Spektrenabgleich
        - Vergleich der Phasenkorrigierten Spektren durch Übereinanderlegen

    - Phasenanalyse mittels Fitting (new approach)
        - Lorentz/Voigt-Peaks
            - superposition von beliebig vielen Peaks
            - Individuelles festsetzen von Parametern
    
    - Regressionsrechnung
        - Berechnung der jeweiligen Mobilitäten aus automatisch bestimmten experimentellen Parametern
    
    - Vergleich der verschiedenen Ergebnisse
        (- einfaches Tool zum Erstellen von Graphen)

- Phasenanalyse mittels 2D FFT --> Mobility ordered Spectroscopy (MOSY)
    - States-Haberkorn-Methode
    - Ermittlung der Mobilitätsachse
    - Plotten von Slices zum Vergleich der Ergebnisse/Peaks
    - Automatische normierung der Intensitäten und Auffindung der Maxima.
    - Durch Signalverlust mit steigender Spannung stellt man eine zu hohe Mobilität fest!!!
