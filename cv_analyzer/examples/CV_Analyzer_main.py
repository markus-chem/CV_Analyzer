# %%
import glob
import matplotlib.pyplot as plt
from datapackage import Package
from cv_analyzer import filter_db

# apply matplotlib style
plt.style.use('seaborn-whitegrid')

# list all .json files
files = glob.glob('../database/*.json')
print(f'{len(files)} files loaded')

# filter for datapackages
# empty bracktes mean that no filter is applied
metal = ['Ag']
lattice_plane = []
component = ['KBr']
author_name = []
not_this_name = ['Rehim']

selxn = filter_db(
    files,
    metal,
    lattice_plane,
    component,
    author_name,
    not_this_name)

# %%
# generic plot
for i in selxn:
    fig = i.plot(target_RE='pzc', C_exp=False, atomic=False, c_corr=False)
fig.legend()
plt.show()

# # %%
# # charge integration
# for i in selxn:
#     fig, fig2, fig3 = i.charge_int(
#         lower_lim=-0.2, upper_lim=0.2, target_RE='pzc', unit='atomic')
# fig.legend()
# fig2.legend()
# fig3.legend()
# # %%
# # max and min
# for i in selxn:
#     fig = i.max_min(
#         lower_lim=0,
#         upper_lim=0.5,
#         target_RE='Ag/AgCl',
#         capac=True)  # lower and upper limit in V
# fig.legend()


# %%
