# CV_Analyzer
### Analysis of Literature Cyclic Voltammograms (CVs)

<p>Literature data is digitized using svgdigitizer (echemdb/svgdigitizer) and stored in a database.
Important information for each CV, such as electrodes, solvent, electrolyte, scan rate and many more are stored as metadata. </p>

<p>Data is stored in datapackages (https://pypi.org/project/datapackage/). By applying certain filter criteria, literature data can be easily compared and analyzed. </p>

<p>Focus on CVs of single crystal metal surfaces in aqueous media. CVs can be selected by filtering a database</p>

To the CVs, different corrections can be applied:

- reference the potential of several CVs to a common reference electrode (including pH-dependent Reversible Hydrogen Electrode (RHE), the potential of zero charge (pzc) of the metal, ...)
- normalize the current density with the scan rate (capacitance plot)
- normalize the geometric curren density to surface atoms of the respective lattice plane
- integrate the charge flow within voltage limits

**installation**

``pip install -e .``

Examples in example.ipynb




