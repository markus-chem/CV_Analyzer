# CV_Analyzer
### Analysis of Literature Cyclic Voltammograms (CVs)

<p>Literature data is digitized using svgdigitizer (echemdb/svgdigitizer) and stored in a database.
Important information for each CV, such as electrodes, solvent, electrolyte, scan rate and many more are stored as metadata. </p>

<p>Data is stored in datapackages (https://pypi.org/project/datapackage/). By applying certain filter criteria, literature data can be easily compared and analyzed. </p>

<p>Focus on CVs of single crystal metal surfaces in aqueous media.</p>

To the CVs, different corrections can be applied:

- reference the potential to another reference electrode (including pH-dependent Reversible Hydrogen Electrode (RHE), the potential of zero charge (pzc) of the metal, ...)
- normalize the areal current density with the scan rate (areal capacitance plot)
- normalize on surface atoms of the respective lattice plane instead of the geometric area
- integrate the passed charge in a voltage window

### Installation

<b>``pip install -e .``</b>

All important functions are in *CV_analyzer_notebook.ipynb*




