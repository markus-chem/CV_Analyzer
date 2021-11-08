#%%
from datapackage import Package
from CV_Analyzer import CV_Analyzer_class, filter
import glob
import os
import sys
import matplotlib.pyplot as plt

plt.style.use('seaborn-whitegrid')

# list all .json files
files = glob.glob('database/*.json')
print(f'Loaded files: \n{files}')

#%%
# filter for datapackages
metal = ['Ag']
lattice_plane = ['100']
# component = []
component = ['KBr']
author_name = ['Nakamura']
not_this_name = []

selected = filter(files, metal, lattice_plane, component, author_name, not_this_name)
# filt_db = filter_db_class(filter)
# filt._db.function
# %%
# reference electrode
for i in selected:
    p = Package(i)
    entry = CV_Analyzer_class(p, 0) # resource_no 0 as default
    print(f'{p.resource_names[0]} selected')
    fig = entry.reference_electrode(target_RE='SCE') # reference potential to another RE
fig.legend()
# %%
# atomic site normalization
for i in selected:
    p = Package(i)
    entry = CV_Analyzer_class(p, 0) # resource_no 0 as default
    print(f'{p.resource_names[0]} selected')
    fig = entry.site_norm(target_RE='pzc', unit='atomic') # unit = 'atomic' for e-/(s*atom)
fig.legend()
# %%
# charge integration
for i in selected:
    p = Package(i)
    entry = CV_Analyzer_class(p, 0) # resource_no 0 as default
    print(f'{p.resource_names[0]} selected')
    fig, fig2, fig3 = entry.charge_int(lower_lim=-0.4, upper_lim=0.395, target_RE='pzc', unit='atomic') # lower and upper limit in V, unit = 'atomic' for reference current to atom sites
fig.legend()
fig2.legend()
fig3.legend()
# %%
# max and min
for i in selected:
    p = Package(i)
    entry = CV_Analyzer_class(p, 0) # resource_no 0 as default
    print(f'{p.resource_names[0]} selected')
    fig = entry.max_min(lower_lim=0, upper_lim=0.5, target_RE='Ag/AgCl') # lower and upper limit in V
fig.legend()
#%%
# plot original data
for i in selected:
    p = Package(i)
    entry = CV_Analyzer_class(p, 0) # resource_no 0 as default
    print(f'{p.resource_names[0]} selected')
    fig = entry.plot() # plot as reported in literature
fig.legend()
