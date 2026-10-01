- fehler der Steigung in self.linreg korrigieren --> stimmt nicht mit lin_huber überein

- lin_huber instanzvariablen abschaffen

- add the parameters from a fitmodel to the measurement object to be saved with save_eNMRpy()

- replace Measurement.eNMR_Methods.lin_huber according Scikit-learn documentation to avoid deprecated code

- Gewichtung der Peaks in Fit mit aufnehmen


SQL-Anbindung
=============

- Welche Datenbankstruktur?
	- Nutzer
		-Messungs-ID (wie generieren, um Duplikate zu vermeiden?) --> Aus Metadaten generieren? --> Timestamp?
		-Probenname
		-Datum
		-Substanz
		-Pfad (zum Zeitpunkt der Auswertung)
		-Kern
		-Mobilität - Stabw.
		-Experiment (Fitwerte, Parameter etc)
			-Spannung
			-Phasenwinkel
			-Diffusionszeit
			-Abhängigkeit (U oder G)
			-Gradientendauer
			-Gradientenintensität

