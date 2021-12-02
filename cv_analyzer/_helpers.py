from datapackage import Package
from pathlib import Path
import os
import shutil

import cv_analyzer as cv
from ._literature_params import atomic_density
from .cv_analyzer_class import CV_Analyzer

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
            if i.metadata['electrochemical system']['electrodes']['working electrode']['material'] not in metal:
                try:
                    selxn.remove(i)
                except BaseException:
                    pass
            else:
                pass

        if len(hkl) > 0:
            if i.metadata['electrochemical system']['electrodes']['working electrode']['crystallographic orientation'] not in hkl:
                try:
                    selxn.remove(i)
                except BaseException:
                    pass
            else:
                pass
        
        if len(component) > 0:
            # join to one string and check if actual components and experimental component match
            component_filter = component
            component_exp = ''.join(filter(None,
                [i.metadata['electrochemical system']['electrolyte']['components'][j]['name'] for j in range(4)]))
            # remove element if none of the filter criteria matches
            if any(component_filter[k] in component_exp for k in range(len(component_filter))) == False:
                try:
                    selxn.remove(i)
                except BaseException:
                    pass

        if len(include_author) > 0:
            # join to one string and check for matches
            include_author_filter = include_author
            include_author_exp = i.metadata['source']['bib'].split('_')[0] # author is string before 1st '_'
            # remove element if none of the filter criteria matches
            if any(include_author_filter[k] in include_author_exp for k in range(len(include_author_filter))) == False:
                try:
                    selxn.remove(i)
                except BaseException:
                    pass

        if any(
            i.metadata['source']['bib'].find(exclude_author[j]) != -1
                for j in range(len(exclude_author))):
            try:
                selxn.remove(i)
            except BaseException:
                pass
        else:
            pass
        
        if len(include_label) > 0:
            if any(include_label[k] in i.label for k in range(len(include_label))) == False:
                try:
                    selxn.remove(i)
                except BaseException:
                    pass

        if any(
            i.label.find(exclude_label[j]) != -1 
                for j in range(len(exclude_label))):
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

    Inputs:
    sel_obj:    CV_Analyzer class object
    target_RE:      str, name of target reference electrode (s. RE_dict)
    C_exp:          bool, capacitance plot
    atomic:         bool, normalization on atomic surface sites
    c_corr:        bool, Nernst shift for halide concentration to 1M

    Returns:
    U:  np.array, the applied potential as a 1-d array
    j:  np.array, the measured current (or pseudocapacitance) as a 1-d array
    """

    # extra columns for corrections
    sel_obj.CV_df['U_corr'] = sel_obj.CV_df['U']
    sel_obj.CV_df['j_corr'] = sel_obj.CV_df['j']

    # reference correction
    ref_shift = sel_obj.ref_to(target_RE)
    sel_obj.CV_df['U_corr'] = sel_obj.CV_df['U_corr'] - ref_shift

    # current in capacitance
    if C_exp:
        scu = sel_obj.scan_rate_unit
        sc  = sel_obj.scan_rate
        if 'mV' in scu:
            # scan rate from mV/s to V/s
            sc /= 1000.
        sel_obj.CV_df['j_corr'] = sel_obj.CV_df['j_corr'] / sc

    # current in atomic units: A/m² / atoms/m² = A/atom
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
    n = len(yaml_list)
    print(f'Found .yaml files: \n{yaml_list}')

    sampling_interval = sampling_interval  # in V

    for i in range(n):
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

def export_datapackage(svg, metadata, sampling_interval=1):
    from datapackage import Package
    import shutil
    import copy
    import yaml
    from svgdigitizer.svgplot import SVGPlot
    from svgdigitizer.svg import SVG
    from svgdigitizer.electrochemistry.cv import CV
    from pathlib import Path
    import datetime

    print(metadata, bool(metadata))
    if metadata:
        with open(metadata, "r") as f:
            metadata = yaml.load(f, Loader=yaml.SafeLoader)

    print(metadata, bool(metadata))
    cv = CV(SVGPlot(SVG(open(svg, 'rb')),
            sampling_interval=sampling_interval), metadata=metadata)
    csvname = Path(svg).with_suffix('.csv')
    print(csvname)
    cv.df.to_csv(csvname, index=False)
    print(cv.metadata)

    def defaultconverter(o):
        if isinstance(o, datetime.datetime):
            return o.__str__()

    p = Package(cv.metadata)

    print(p.descriptor.keys())

    p.infer(str(csvname))
    print(p.descriptor['resources'][0]['name'],
          p.descriptor['resources'][0]['path'],
          p.descriptor['resources'][0]['schema'])

    import json
    with open(Path(svg).with_suffix('.json'), "w") as outfile:
        json.dump(p.descriptor, outfile, default=defaultconverter)
