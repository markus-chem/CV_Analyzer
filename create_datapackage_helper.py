import svgdigitizer
from datapackage import Package
import yaml
from svgdigitizer.svgplot import SVGPlot
from svgdigitizer.svg import SVG
from svgdigitizer.electrochemistry.cv import CV
import os
from pathlib import Path
import datetime

def cv_export(svg, sampling_interval, metadata):
    import yaml
    from svgdigitizer.svgplot import SVGPlot
    from svgdigitizer.svg import SVG
    from svgdigitizer.electrochemistry.cv import CV
    print(metadata, bool(metadata))
    if metadata:
        with open(metadata,"r") as f:
            metadata = yaml.load(f, Loader=yaml.SafeLoader)

        
    print(metadata, bool(metadata))
    cv = CV(SVGPlot(SVG(open(svg, 'rb')), sampling_interval=sampling_interval), metadata=metadata)

    from pathlib import Path
    cv.df.to_csv(Path(svg).with_suffix('.csv'), index=False)
    print(cv.metadata)
    import datetime

    def defaultconverter(o):
        if isinstance(o, datetime.datetime):
            return o.__str__()

    import json
    with open(Path(svg).with_suffix('.json'), "w") as outfile:
        json.dump(cv.metadata, outfile, default=defaultconverter)

def export_datapackage(svg, metadata, sampling_interval = 1):
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
        with open(metadata,"r") as f:
            metadata = yaml.load(f, Loader=yaml.SafeLoader)
        
    print(metadata, bool(metadata))
    cv = CV(SVGPlot(SVG(open(svg, 'rb')), sampling_interval=sampling_interval), metadata=metadata)
    csvname = Path(svg).with_suffix('.csv')
    print(csvname)
    cv.df.to_csv( csvname , index=False)
    print(cv.metadata)
    

    def defaultconverter(o):
        if isinstance(o, datetime.datetime):
            return o.__str__()

    p = Package(cv.metadata)

    print(p.descriptor.keys())
    
    p.infer(str(csvname))
    print(p.descriptor['resources'][0]['name'], \
    p.descriptor['resources'][0]['path'], \
    p.descriptor['resources'][0]['schema'])
        
    import json
    with open(Path(svg).with_suffix('.json'), "w") as outfile:
        json.dump(p.descriptor, outfile, default=defaultconverter)