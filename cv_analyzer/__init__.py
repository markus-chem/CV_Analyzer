# CV_Analyzer class for datapackages
# for the analysis and the manipulation of digitized CVs

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from datapackage import Package
from ._literature_params import (RE_dict, pzc_dict, lattice_constants_dict,
                                   atomic_density)
from scipy import constants

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

class CV_Analyzer:
    def __init__(self, package, resource_no):
        """
        Initialize with basic information of package defined by the datapackage
        module implemented in ### REFERENCE!!!!
        Inputs:
        package:        Package, ### REFERENCE!!!!
        resource_no:    str, name of the resource which you want to access
                        within the larger datapackage
        """
        self.package = package
        self.CV_df = pd.read_csv(package.resources[0].raw_iter(stream=False))
        self.name = package.resource_names[resource_no]
        self._get_inputs()


    def _get_inputs(self):
        """
        Collects the experimental information from the package information.
        """
        metadata    = self.package.descriptor
        e_chem_sys  = metadata['electrochemical system']

        electrodes  = e_chem_sys['electrodes']
        self.RE     = electrodes['reference electrode']['type']
        work_elec   = electrodes['working electrode']
        self.metal  = work_elec['material']
        self.hkl    = work_elec['crystallographic orientation']

        # Assuming there are only 2 components
        electrolyte = e_chem_sys['electrolyte']
        elec_comps  = electrolyte['components']
        self.electrolyte_name   = [elec_comps[i]['name'] for i in [0, 3]]
        self.c_electrolyte      = [
            elec_comps[i]['concentration']['value'] for i in [0, 3]
        ]
        self.electrolyte_unit   = [
            elec_comps[i]['concentration']['unit'] for i in [0, 3]
        ]
        self.electrolyte        = [
            "{} {} {}".format(self.c_electrolyte[i], self.electrolyte_unit[i],
                              self.electrolyte_name[i]
                             ) for i in range(2)
        ]
        self.T                  = electrolyte['temperature']['value']
        self.pH                 = electrolyte['ph'].get('value', None)

        fig_desc = metadata['figure description']
        self.scan_rate_unit = fig_desc['scan rate']['unit']
        self.scan_rate = fig_desc['scan rate']['value']

        self.label = "{}({}) in {} at {} {}; {}".format(
            self.metal, self.hkl, self.electrolyte, self.scan_rate,
            self.scan_rate_unit, self.name
        )
        self.metadata = metadata
        return None

    def plot_orig(self):
        """
        Recreate the cyclic voltammogram with the original experimental
        parameters.
        """
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_xlabel("V vs. {} (eV)".format(self.RE))
        ax.set_ylabel(r"$\rm j A/m^2")
        ax.plot(self.CV_df['U'], self.CV_df['j'], label=self.label,
                color=colors[0]
        )
        return fig

    def plot(self, target_RE='SHE', C_exp=False, atomic=False, c_corr=False):
        """
        Plots the cyclic voltammogram, however permits changes such as Nernstian
        shifts, or changes to the reference electrode
        """
        settings = f'C_exp={C_exp}, atomic={atomic}, c_corr={c_corr}'

        # extra columns for corrections
        self.CV_df['U_corr'] = self.CV_df['U']
        self.CV_df['j_corr'] = self.CV_df['j']

        # reference correction
        self.CV_df['U_corr'] = self.CV_df['U_corr'] - self.ref_to(target_RE)

        # current to capacitance
        if C_exp:
            scu = self.scan_rate_unit
            if 'mV' in scu:
                sc = self.scan_rate / 1000
                self.CV_df['j_corr'] = self.CV_df['j_corr'] / sc
            else:
                self.CV_df['j_corr'] = self.CV_df['j_corr'] / self.scan_rate

        # current in atomic units: A/m² / atoms/m² = A/atom
        if atomic:
            rho_atomic = atomic_density(self.metal, self.hkl)
            self.CV_df['j_corr'] = self.CV_df['j_corr'] / rho_atomic

        # apply concentration correction
        if c_corr:
            self.CV_df['U_corr'] = self.CV_df['U_corr'] - self.c_corr()

        # plot
        fig = plt.figure(f'CV_settings: {settings}')
        ax  = fig.add_subplot(111)
        ax.plot(self.CV_df['U_corr'], self.CV_df['j_corr'], label=self.label)
        ax.set_xlabel(f'U / V vs. {target_RE}')
        if (C_exp) & (atomic == False):
            ax.set_ylabel('C / F/m$^2$')
        elif (C_exp) & (atomic):
            ax.set_ylabel('C / F/atom')
        elif (C_exp == False) & (atomic == False):
            ax.set_ylabel('j / A/m$^2$')
        elif (C_exp == False) & (atomic == True):
            ax.set_ylabel('j / A/atom')

        return fig

    def charge_int(
            self,
            lower_lim,
            upper_lim,
            target_RE,
            atomic=False,
            c_corr=False):

        # extra columns for corrections
        self.CV_df['U_corr'] = self.CV_df['U']
        self.CV_df['j_corr'] = self.CV_df['j']

        # reference to target_RE
        self.CV_df['U_corr'] = self.CV_df['U'] - \
            self.ref_to(target_RE)  # reference to another RE

        # check voltage limits
        if lower_lim < min(
                self.CV_df['U_corr']) or upper_lim > max(
                self.CV_df['U_corr']):
            raise ValueError(f'Voltage limits out of range for {self.name}')

        # apply concentration correction
        if c_corr:
            self.CV_df['U_corr'] = self.CV_df['U_corr'] - self.c_corr()
        else:
            pass

        # atomic units
        if atomic:
            norm_factor = atomic_density(
                self.metal, self.hkl) * 10**20
            self.CV_df['j_corr'] = self.CV_df['j'] / constants.e / \
                norm_factor  # A/m² / C/e- / atoms/m² = e-/(s*atom)

        # seperate anodic and cathodic scan, apply voltage limits
        forw_scan = self.CV_df.where(
            (self.CV_df['j_corr'] > 0) & (
                self.CV_df['U_corr'] > lower_lim) & (
                self.CV_df['U_corr'] < upper_lim))
        backw_scan = self.CV_df.where(
            (self.CV_df['j_corr'] < 0) & (
                self.CV_df['U_corr'] > lower_lim) & (
                self.CV_df['U_corr'] < upper_lim))

        # voltage for the integration needs to start with 0 V
        int_voltage_forw = forw_scan['U_corr'] - lower_lim
        int_voltage_backw = backw_scan['U_corr'] - lower_lim

        # convert to numpy and remove nan (needed for np.trapz() integration)
        voltage_forw_np = forw_scan['U_corr'].to_numpy()
        voltage_forw_np = voltage_forw_np[~np.isnan(voltage_forw_np)]
        voltage_backw_np = backw_scan['U_corr'].to_numpy()
        voltage_backw_np = voltage_backw_np[~np.isnan(voltage_backw_np)]
        forw_scan_np = forw_scan['j_corr'].to_numpy()
        forw_scan_np = forw_scan_np[~np.isnan(forw_scan_np)]
        backw_scan_np = backw_scan['j_corr'].to_numpy()
        backw_scan_np = backw_scan_np[~np.isnan(backw_scan_np)]
        int_voltage_forw_np = np.array(int_voltage_forw)
        int_voltage_forw_np = int_voltage_forw_np[~np.isnan(
            int_voltage_forw_np)]
        int_voltage_backw_np = np.array(int_voltage_backw)
        int_voltage_backw_np = int_voltage_backw_np[~np.isnan(
            int_voltage_backw_np)]

        # integrate
        Q_forw_int = [np.trapz(forw_scan_np[:i+1], int_voltage_forw_np[:i+1]) /
                      self.scan_rate for i in range(len(forw_scan_np))]
        Q_backw_int = [np.trapz(backw_scan_np[:i+1], int_voltage_backw_np[:i+1]) /
                       self.scan_rate for i in range(len(backw_scan_np))]

        # plot CV with limits
        fig = plt.figure('integrated CV')
        plt.plot(self.CV_df['U_corr'], self.CV_df['j_corr'], label=self.label)
        plt.vlines(
            lower_lim, min(
                self.CV_df['j_corr']), max(
                self.CV_df['j_corr']), color='grey', linestyles='solid')
        plt.vlines(
            upper_lim, min(
                self.CV_df['j_corr']), max(
                self.CV_df['j_corr']), color='grey', linestyles='solid')
        plt.xlabel(f'U / V vs. {target_RE}')
        if atomic:
            plt.ylabel('j / e-/(s*atom)')
        else:
            plt.ylabel('j / µA/cm²')

        # evaluate Q_forw over voltage and plot
        fig2 = plt.figure('charge integration forward')
        plt.plot(voltage_forw_np, Q_forw_int, label=f'{self.label}')
        plt.xlabel(f'U / V vs. {target_RE}')
        if atomic:
            plt.ylabel('Q / e-/atom')
        else:
            plt.ylabel('Q / µC/cm²')

        # evaluate Q_backw over voltage and plot
        fig3 = plt.figure('charge integration backward')
        plt.plot(voltage_backw_np, Q_backw_int, label=f'{self.label}')
        plt.xlabel(f'U / V vs. {target_RE}')
        if atomic:
            plt.ylabel('Q / e-/atom')
        else:
            plt.ylabel('Q / µC/cm²')

        return fig, fig2, fig3

    def max_min(self, lower_lim, upper_lim, target_RE, C_exp=False):

        # extra columns for corrections
        self.CV_df['U_corr'] = self.CV_df['U']
        self.CV_df['j_corr'] = self.CV_df['j']

        # reference correction
        self.CV_df['U_corr'] = self.CV_df['U_corr'] - self.ref_to(target_RE)

        # check voltage limits
        if lower_lim < min(
                self.CV_df['U_corr']) or upper_lim > max(
                self.CV_df['U_corr']):
            raise ValueError(f'Voltage limits out of range for {self.name}')

        # apply voltage limits
        self.CV_df = self.CV_df.where(
                (self.CV_df['U_corr'] > lower_lim) &
                (self.CV_df['U_corr'] < upper_lim))

        # calculate C_expitance if needed
        if C_exp:
            self.CV_df['j_corr'] = self.CV_df['j_corr'] / (self.scan_rate / 1000)
        else:
            pass

        idx = (self.CV_df['j_corr'].idxmax(), self.CV_df['j_corr'].idxmin()) # indices of min and max
        max_ = (self.CV_df.iloc[idx[0]]['U_corr'], self.CV_df.iloc[idx[0]]['j_corr'])  # x and y of max
        min_ = (self.CV_df.iloc[idx[1]]['U_corr'], self.CV_df.iloc[idx[1]]['j_corr'])  # x and y of min

        fig = plt.figure('Max_Min of CV')
        plt.plot(self.CV_df['U_corr'], self.CV_df['j_corr'], label=self.label)
        plt.scatter(max_[0], max_[1], marker='X', color='r')
        plt.scatter(min_[0], min_[1], marker='X', color='r')
        plt.xlabel(f'U / V vs. {target_RE}')
        plt.ylabel('j / µA/cm$^2$')

        print(
            f'Maximum: U = {round(max_[0], 2)}, j = {round(max_[1], 2)} \
            Minimum U = {round(min_[1], 2)}, j = {round(min_[1], 2)} \
            {self.name}')

        return fig


    def ref_to(self, target_RE):
        """
        Switches the reference electrode to one of the options included in the
        RE dictionary.
        Inputs:
        target_RE:  str, name of the target reference electrode
        Output:
        offset:     float, Potential shift due to reference change.
        """
        # check if target_RE exists
        RE_keys = list(RE_dict.keys()) + ['pzc', 'RHE']
        if target_RE not in RE_keys:
            s = "Target reference is unknown, select one of the following: "
            s += ', '.join(RE_keys)
            raise ValueError(s)

        # check for pH dependency
        elif target_RE == 'RHE':
            if self.pH is None:
                raise ValueError('pH value is undefined, RHE not possible.')
            else:
                # offset is later substracted from the voltage
                R, F =  constants.R, constants.value('Faraday constant')
                offset = RE_dict[str(target_RE)] - RE_dict[str(self.RE)]
                # Includes the 59 mV shift per pH unit
                offset += R * self.T / F * self.pH

        # Sets the potential of zero charge as the reference.
        elif target_RE == 'pzc':
            pzc = pzc_dict[self.metal][self.hkl]
            offset =  pzc - RE_dict[str(self.RE)] - RE_dict['U_abs']
        else:
            offset = RE_dict[str(target_RE)] - RE_dict[str(self.RE)]
        return offset


    def c_corr(self):
        """
        Performs a Nernstian shift to shift the concentration to the proper
        potential.
        Returns:
        U_shift: float, the nernstian potential shift that is substracted
        """
        # list of typical halide adsorbates
        ads_set = {'LiCl', 'NaCl', 'KCl', 'RbCl', 'CsCl',
                'LiBr', 'NaBr', 'KBr', 'RbBr', 'CsBr',
                'LiI', 'NaI', 'KI', 'RbI', 'CsI'
                }
        # get list of common elements
        ads = list(ads_set.intersection(self.electrolyte_name))

        if len(ads) == 0:
            raise ValueError('No halide adsorbate found')
        elif len(ads) > 1:
            raise ValueError(f'Cannot handle more than one adsorbate: {ads}')

        idx = self.electrolyte_name.index(str(ads[0]))
        ads_conc = self.c_electrolyte[idx]
        unit = self.electrolyte_unit[idx]

        if str(unit) == 'uM':
            ads_conc = ads_conc / 1.e06
        elif str(unit) == 'mM':
            ads_conc = ads_conc / 1.e03
        elif str(unit) != 'M':
            raise ValueError(f'Unit Error: {unit}')
        R, F = constants.R, constants.value('Faraday constant')
        U_shift = - (R * self.T / F) * np.log(ads_conc)
        return U_shift


def filter_db(metal, hkl, component, **kwargs): #author_name, exclude_author):
    """
    Inputs
    kwargs:
    author_name:    list, list of first authors you wish to include, defaults to
                    []
    exclude_author: list, list of authors you wish to exclude, defaults to []
    """
    from .database import get_database
    files = get_database()
    print(f'{len(files)} files loaded')

    if not isinstance(metal, list):
        metal = [metal]
    if not isinstance(hkl, list):
        hkl = [hkl]
    hkl = [str(i) for i in hkl]
    if not isinstance(component, list):
        component = [component]

    author_name     = kwargs.get("author_name", [])
    exclude_author  = kwargs.get("exclude_author", [])

    selxn = set(CV_Analyzer(Package(i), 0) for i in files)
    for i in selxn.copy():  # iterate over copy, set cannot be changed during iteration
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

        if len(author_name) > 0:
            # join to one string and check for matches
            author_name_filter = author_name
            author_name_exp = i.metadata['source']['bib'].split('_')[0] # author is string before 1st '_'
            # remove element if none of the filter criteria matches
            if any(author_name_filter[k] in author_name_exp for k in range(len(author_name_filter))) == False:
                try:
                    selxn.remove(i)
                except BaseException:
                    pass

        # if len(author_name) > 0:
        #     if any(
        #         i.metadata['source']['bib'].find(
        #             author_name[j]) == -1 for j in range(
        #             len(author_name))):
        #         try:
        #             selxn.remove(i)
        #         except BaseException:
        #             pass
        #     else:
        #         pass

        if any(
            i.metadata['source']['bib'].find(
                exclude_author[j]) != -
                1 for j in range(
                len(exclude_author))):
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
    corrected applied potential and current (or pseudocapacitance)
    Inputs:
    sel_obj:    CV_Analyzer class, the experimental selection
    target_RE:  str, name of the reference electrode to which you want to shift
                the potential.
    C_exp:      bool, decides whether to translate the measured current to the
                pseudocapacitance. This requires knowledge on the scan rate.
    atomic:     bool, norms all the data to number of atoms, instead of on the
                area.
    c_corr:     bool, if True, includes a Nernstian shift depending on the
                electrolyte concentrations.
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
