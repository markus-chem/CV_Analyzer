# %%
import glob, os, shutil
from create_datapackage_helper import export_datapackage, cv_export
from pathlib import Path

# %%
# makes datapackages and saves them as .zip
yaml_list = []
svg_list = []
name_list = []
for path in Path('data').rglob('*.yaml'):
    yaml_list.append(path)
    svg_list.append(path.with_suffix('.svg'))
    name_list.append(str(path.stem))
n = len(yaml_list)
print(f'Found .yaml files: \n{yaml_list}')
# %%
sampling_interval = 0.05  # in V

for i in range(n):
    shutil.copy(svg_list[i], os.getcwd())
    shutil.copy(yaml_list[i], os.getcwd())
    export_datapackage(
        name_list[i] + '.svg',
        name_list[i] + '.yaml',
        sampling_interval)
    for j in ['.csv', '.yaml', '.svg', '.json']:
        shutil.move(
            name_list[i] + j,
            os.path.join(
                'database',
                name_list[i] + j))

# %%
