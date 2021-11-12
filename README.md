# CV_Analyzer
## Analysis of Literature Cyclic Voltammograms

Literature data is digitized using svgdigitizer (echemdb/svgdigitizer) and stored in a database.
Important information for each CV, such as electrodes, solvent, electrolyte, scan rate and many more are stored as metadata. \\

Different corrections can be applied.

- reference the potential to another reference electrode
- normalize the areal current density with the scan rate (areal capacitance plot)
- normalize on surface atoms of the respective lattice plane instead of the geometric area.
- integrate the charge in a voltage window (separately for forward and backward scan)

By applying certain filter criteria, literature data can be easily compared and analyzed.



