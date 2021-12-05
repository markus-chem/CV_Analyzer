from datapackage import Package
from pathlib import Path
import os
import shutil
import yaml
import datetime

import cv_analyzer as cv
from ._literature_params import atomic_density
from .cv_analyzer_class import CV_Analyzer

from svgdigitizer.svgplot import SVGPlot
from svgdigitizer.svg import SVG
from svgdigitizer.electrochemistry.cv import CV

def get_database():
    """
    Get list of .json file directories in the database folder
    
    Inputs: None
    
    Outputs:
    List of directories
    """
    
    datadir = cv.__path__[0]+'/database/'
    files = os.listdir(datadir)
    files.sort()
    files = [os.path.join(datadir, i) for i in files if 'json' in i]
    return files

def filter_db(metal, hkl, component, **kwargs):
    """
    Filter function for the database.
    Searches for datapackages in 'database' folder that match the criteria:
        - metal
        - hkl of surface
        - component (can be a fraction, e.g. 'Br' if component is 'NaBr')
        - include_author (can be a fraction)
        - exclude_author (do not show this author)
        - include_label
        - exclude_label


    If one filter criterium is empty, it is not considered.
    All criteria of different filters have to met a the same time
    (e.g. component and include_author)
    Within one filter the selection is additive (e.g. component=['Pt','Ag'] shows both)

    Inputs:
    **kwargs:
        metal:          list of str, material of working electrode
        hkl:            list of str, lattice plane
        component:      list of str, component of the electrolyte
        include_author:    list, list of first authors last name you wish to include,
                        defaults to []
        exclude_author: list, list of authors last name you wish to exclude, defaults to []
        include_label:  list of labels to be include, default to []
        exclude_label:  list of labels to be exclude, default to []

    Outputs:
    list of selected datapackages
    """

    files = get_database()
    print(f'{len(files)} files loaded')

    include_author  = kwargs.get('include_author', [])
    exclude_author  = kwargs.get('exclude_author', [])
    include_label   = kwargs.get('include_label', [])
    exclude_label   = kwargs.get('exclude_label', [])

    include_author  = [j.lower() for j in include_author]
    exclude_author  = [j.lower() for j in exclude_author]
    include_label   = [j.lower() for j in include_label]
    exclude_label   = [j.lower() for j in exclude_label]

    if not isinstance(metal, list):
        metal = [metal]
    if not isinstance(hkl, list):
        hkl = [hkl]
    hkl = [str(i) for i in hkl]
    if not isinstance(component, list):
        component = [component]

    selxn = set(CV_Analyzer(Package(i), 0) for i in files) # only the first resource is considered
    for i in selxn.copy():  # iterate over copy, set cannot be changed during iteration
        # Apply filter criteria
        if len(metal) > 0:
            if i.metal not in metal:
                try:
                    selxn.remove(i)
                except BaseException:
                    pass
            else:
                pass

        if len(hkl) > 0:
            if i.hkl not in hkl:
                try:
                    selxn.remove(i)
                except BaseException:
                    pass
            else:
                pass
        
        if len(component) > 0:
            # remove element if none of the filter criteria matches
            if any(k in i.electrolyte_name for j, k in enumerate(component)) == False:
                try:
                    selxn.remove(i)
                except BaseException:
                    pass

        if len(include_author) > 0:
            # remove element if none of the filter criteria matches
            if any(k in i.author for j, k in enumerate(include_author)) == False:
                try:
                    selxn.remove(i)
                except BaseException:
                    pass

        if any(i.author is k for j, k in enumerate(exclude_author)):
            try:
                selxn.remove(i)
            except BaseException:
                pass
        else:
            pass
        
        if len(include_label) > 0:
            if any(k in i.label for j, k in enumerate(include_label)) == False:
                try:
                    selxn.remove(i)
                except BaseException:
                    pass

        # if any(i.label.find(k) != -1 for j, k in enumerate(exclude_label)):
        if any(k in i.label for j, k in enumerate(exclude_label)):
            try:
                selxn.remove(i)
            except BaseException:
                pass
        else:
            pass

    if len(selxn) > 0:
        print(f'{len(selxn)} files selected')
    else:
        raise ValueError('No datapackages meet filter criteria.')

    return list(selxn)

def get_exp_CV_data(sel_obj, target_RE='SHE', C_exp=False, atomic=False,
                    c_corr=False):
    """
    This function takes a selection of experimental data and returns the
    corrected applied potential and current (or pseudocapacitance) as array.
    (same corrections as in plot function)

    Inputs:
    sel_obj:        CV_Analyzer class object
    target_RE:      str, name of target reference electrode (s. RE_dict)
    C_exp:          bool, capacitance plot
    atomic:         bool, normalization on atomic surface sites
    c_corr:         bool, Nernst shift for halide concentration to 1M

    Returns:
    U:              np.array, the applied potential as a 1-d array
    j:              np.array, the measured current (or pseudocapacitance) as a 1-d array
    """

    # extra columns for corrections
    sel_obj.CV_df['U_corr'] = sel_obj.CV_df['U']
    sel_obj.CV_df['j_corr'] = sel_obj.CV_df['j']

    # reference correction
    sel_obj.CV_df['U_corr'] = sel_obj.CV_df['U_corr'] - sel_obj.ref_to(target_RE)

    # current to capacitance
    if C_exp:
            sel_obj.CV_df['j_corr'] = sel_obj.CV_df['j_corr'] / sel_obj.scan_rate

    # current in atomic units: µA/cm² / atoms/cm² = µA/atom
    if atomic:
        rho_atomic = atomic_density(sel_obj.metal, sel_obj.hkl)
        sel_obj.CV_df['j_corr'] = sel_obj.CV_df['j_corr'] / rho_atomic

    # apply concentration correction
    if c_corr:
        sel_obj.CV_df['U_corr'] = sel_obj.CV_df['U_corr'] - sel_obj.c_corr()

    U = sel_obj.CV_df['U_corr'].to_numpy()
    j = sel_obj.CV_df['j_corr'].to_numpy()

    return U, j

def create_datapackage(sampling_interval=0.05):
    '''
    Find all .yaml and .svg files in 'data' folder.
    Export them as datapackage into 'database'
    Database will contain .csv .yaml .svg and .json for each CV
    The sampling interval is the parameter of the interpolation
    procedure in the svgdigitizer module.

    Inputs:
    sampling_interval: float

    Outputs:
    None
    '''
    path0 = cv.__path__[0] # base path of the package
    yaml_list = []
    svg_list = []
    name_list = []
    for path in Path(path0,'data').rglob('*.yaml'):
        yaml_list.append(path)
        svg_list.append(path.with_suffix('.svg'))
        name_list.append(str(path.stem))
    n = len(name_list)
    print('Found .yaml files:\n')
    print(*name_list, sep='\n')

    sampling_interval = sampling_interval  # in V

    print('\nExporting ... \n')
    for i in range(n):
        print(name_list[i])
        shutil.copy(svg_list[i], os.getcwd())
        shutil.copy(yaml_list[i], os.getcwd())
        export_datapackage(
            name_list[i] + '.svg',
            name_list[i] + '.yaml',
            sampling_interval)
        for j in ['.csv', '.yaml', '.svg', '.json']:
            shutil.move(name_list[i] + j,
                os.path.join(path0, 'database', name_list[i] + j))

    return None

def export_datapackage(svg, metadata, sampling_interval):

    if metadata:
        with open(metadata, 'r') as f:
            metadata = yaml.load(f, Loader=yaml.SafeLoader)

    cv = CV(SVGPlot(SVG(open(svg, 'rb')),
            sampling_interval=sampling_interval), metadata=metadata)
    csvname = Path(svg).with_suffix('.csv')
    cv.df.to_csv(csvname, index=False)

    def defaultconverter(o):
        if isinstance(o, datetime.datetime):
            return o.__str__()

    p = Package(cv.metadata)
    p.infer(str(csvname))

    import json
    with open(Path(svg).with_suffix('.json'), 'w') as outfile:
        json.dump(p.descriptor, outfile, default=defaultconverter)
