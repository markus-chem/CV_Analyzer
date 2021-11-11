# CV_Analyzer class for datapackages
# for the analysis and the manipulation of digitized CVs

from literature_parameters import RE_dict, pzc_dict, lattice_constants_dict, atomic_density
from datapackage import Package
from scipy import constants
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


class CV_Analyzer_class:
    def __init__(self, p, resource_no):  # initialize with basic information of package
        self.package = p
        self.CV_df = pd.read_csv(p.resources[0].raw_iter(stream=False))
        self.name = p.resource_names[resource_no]  # the name of the resource
        self.metadata = p.descriptor
        self.RE = p.descriptor['electrochemical system']['electrodes']['reference electrode']['type']
        self.metal = p.descriptor['electrochemical system']['electrodes']['working electrode']['material']
        self.lattice_plane = p.descriptor['electrochemical system'][
            'electrodes']['working electrode']['crystallographic orientation']
        self.electrolyte_name = [
            p.descriptor['electrochemical system']['electrolyte']['components'][i]['name'] for i in [
                0, 3]]  # assumes that there are only 2 components
        self.electrolyte_conc = [p.descriptor['electrochemical system'][
            'electrolyte']['components'][i]['concentration']['value'] for i in [0, 3]]
        self.electrolyte_unit = [p.descriptor['electrochemical system'][
            'electrolyte']['components'][i]['concentration']['unit'] for i in [0, 3]]
        self.electrolyte = [
            f'{self.electrolyte_conc[i]} {self.electrolyte_unit[i]} {self.electrolyte_name[i]}' for i in range(2)]
        self.scan_rate_unit = p.descriptor['figure description']['scan rate']['unit']
        # if self.scan_rate_unit == 'mV / s':
        #     self.scan_rate = p.descriptor['figure description']['scan rate']['value'] / 1000
        self.scan_rate = p.descriptor['figure description']['scan rate']['value']
        self.T = p.descriptor['electrochemical system']['electrolyte']['temperature']['value']
        self.label = f'{self.metal}({self.lattice_plane}) in {self.electrolyte} at {self.scan_rate} {self.scan_rate_unit}; {self.name}'

    def plot_orig(self):
        # plot the non-manipulated CV data
        fig = plt.figure('Original CV')
        plt.plot(self.CV_df['U'], self.CV_df['j'], label=self.label)
        plt.xlabel(f'U / V vs. [original RE]')
        plt.ylabel('j / A/m$^2$')

        return fig

    def plot(
            self,
            target_RE='SHE',
            capac=False,
            atomic=False,
            conc_corr=False):

        settings = f'capac={capac}, atomic={atomic}, conc_corr={conc_corr}'

        # extra columns for corrections
        self.CV_df['U_corr'] = self.CV_df['U']
        self.CV_df['j_corr'] = self.CV_df['j']

        # reference correction
        self.CV_df['U_corr'] = self.CV_df['U_corr'] - self.ref_to(target_RE)

        # current in capacitance
        if capac:
            self.CV_df['j_corr'] = self.CV_df['j_corr'] / (self.scan_rate / 1000) # scan rate from mV/s to V/s
        else:
            pass

        # current in atomic units
        if atomic:
            self.CV_df['j_corr'] = self.CV_df['j_corr'] / atomic_density(
                self.metal, self.lattice_plane)  # A/m² / atoms/m² = A/atom
        else:
            pass

        # apply concentration correction
        if conc_corr:
            self.CV_df['U_corr'] = self.CV_df['U_corr'] - self.conc_corr()
        else:
            pass

        # plot
        fig = plt.figure(f'CV_settings: {settings}')
        plt.plot(self.CV_df['U_corr'], self.CV_df['j_corr'], label=self.label)
        plt.xlabel(f'U / V vs. {target_RE}')
        if (capac) & (atomic == False):
            plt.ylabel('C / F/m$^2$')
        elif (capac) & (atomic):
            plt.ylabel('C / F/atom')
        elif (capac == False) & (atomic == False):
            plt.ylabel('j / A/m$^2$')
        elif (capac == False) & (atomic == True):
            plt.ylabel('j / A/atom')

        return fig

    def charge_int(
            self,
            lower_lim,
            upper_lim,
            target_RE,
            atomic=False,
            conc_corr=False):

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
        if conc_corr:
            self.CV_df['U_corr'] = self.CV_df['U_corr'] - self.conc_corr()
        else:
            pass

        # atomic units
        if atomic:
            # from atoms/Angstrom² to atoms/m²
            norm_factor = atomic_density(
                self.metal, self.lattice_plane) * 10**20
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
                      (self.scan_rate/1000) for i in range(len(forw_scan_np))]
        Q_backw_int = [np.trapz(backw_scan_np[:i+1], int_voltage_backw_np[:i+1]) /
                       (self.scan_rate/1000) for i in range(len(backw_scan_np))]

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
            plt.ylabel('j / A/m²')

        # evaluate Q_forw over voltage and plot
        fig2 = plt.figure('charge integration forward')
        plt.plot(voltage_forw_np, Q_forw_int, label=f'{self.label}')
        plt.xlabel(f'U / V vs. {target_RE}')
        if atomic:
            plt.ylabel('Q / e-/atom')
        else:
            plt.ylabel('Q / C/m²')

        # evaluate Q_backw over voltage and plot
        fig3 = plt.figure('charge integration backward')
        plt.plot(voltage_backw_np, Q_backw_int, label=f'{self.label}')
        plt.xlabel(f'U / V vs. {target_RE}')
        if atomic:
            plt.ylabel('Q / e-/atom')
        else:
            plt.ylabel('Q / C/m²')

        return fig, fig2, fig3

    def max_min(self, lower_lim, upper_lim, target_RE, capac=False):

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

        # calculate capacitance if needed
        if capac:
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
        plt.ylabel('j / A/m$^2$')

        print(
            f'Maximum: U = {round(max_[0], 2)}, j = {round(max_[1], 2)} \
            Minimum U = {round(min_[1], 2)}, j = {round(min_[1], 2)} \
            {self.name}')

        return fig


    def ref_to(self, target_RE):
        # reference to another RE

        # check if target_RE exists
        if target_RE in RE_dict or target_RE == 'pzc' or target_RE == 'RHE':
            pass
        else:
            print('Select one of the following RE')
            print('pzc')
            for key in RE_dict:
                print(key)
            raise ValueError('target reference is unknown')

        # check for pH dependency
        if target_RE == 'RHE':
            pH = self.metadata['electrochemical system']['electrolyte']['ph']['value']
            if pH is None:
                # try to calculate pH?
                raise ValueError('No pH value given. Conversion to RHE not possible.')
            else:
                # offset is later substracted from the voltage
                offset = RE_dict[str(target_RE)] - RE_dict[str(self.RE)] + constants.R * \
                    self.T / constants.value('Faraday constant') * pH  # 59 mV shift per pH unit
        if target_RE == 'pzc':
            offset = pzc_dict[self.metal][self.lattice_plane] - \
                RE_dict[str(self.RE)] - RE_dict['U_abs']
        else:
            offset = RE_dict[str(target_RE)] - RE_dict[str(self.RE)]

        return offset


    def conc_corr(self):
        # list of typical halide adsorbates
        ads_set = {'LiF', 'NaF', 'KF', 'RbF',
                'LiCl', 'NaCl', 'KCl', 'RbCl', 'CsCl',
                'LiBr', 'NaBr', 'KBr', 'RbBr', 'CsBr',
                'LiI', 'NaI', 'KI', 'RbI', 'CsI'
                }
        # get list of common elements
        ads = list(ads_set.intersection(self.electrolyte_name))

        if len(ads) == 0:
            raise ValueError('No halide adsorbate found')
        elif len(ads) > 1:
            raise ValueError(f'Cannot handle more than one adsorbate: {ads}')
        else:
            idx = self.electrolyte_name.index(str(ads[0]))
            ads_conc = self.electrolyte_conc[idx]
            unit = self.electrolyte_unit[idx]

            if str(unit) == 'uM':
                ads_conc = ads_conc / 10**6
            elif str(unit) == 'mM':
                ads_conc = ads_conc / 1000
            elif str(unit) == 'M':
                pass
            else:
                raise ValueError(f'Unit Error: {unit}')
            U_shift = - constants.R * self.T / \
                constants.value('Faraday constant') * np.log(ads_conc)

        return U_shift  # is later substracted


def filter_db(
        files,
        metal,
        lattice_plane,
        component,
        author_name,
        not_this_name):
    selxn = set(CV_Analyzer_class(Package(i), 0) for i in files)
    for i in selxn.copy():  # iterate over copy, set cannot be changed during iteration
        if len(metal) > 0:
            if i.metadata['electrochemical system']['electrodes']['working electrode']['material'] not in metal:
                try:
                    selxn.remove(i)
                except BaseException:
                    pass
            else:
                pass

        if len(lattice_plane) > 0:
            if i.metadata['electrochemical system']['electrodes']['working electrode']['crystallographic orientation'] not in lattice_plane:
                try:
                    selxn.remove(i)
                except BaseException:
                    pass
            else:
                pass

        if len(component) > 0:
            l = list(set(component).intersection(
                [i.metadata['electrochemical system']['electrolyte']['components'][j]['name'] for j in range(4)]))
            if len(l) == 0:
                try:
                    selxn.remove(i)
                except BaseException:
                    pass

        if len(author_name) > 0:
            if any(
                i.metadata['source']['bib'].find(
                    author_name[j]) == -1 for j in range(
                    len(author_name))):
                try:
                    selxn.remove(i)
                except BaseException:
                    pass
            else:
                pass

        if any(
            i.metadata['source']['bib'].find(
                not_this_name[j]) != -
                1 for j in range(
                len(not_this_name))):
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

def get_exp_CV_data(sel_obj, target_RE='SHE',
                    capac=False,
                    atomic=False,
                    conc_corr=False):

    # object is one element of selxn, i.e. a CV_Analyzer object
    # returns U and j as array

    # extra columns for corrections
    sel_obj.CV_df['U_corr'] = sel_obj.CV_df['U']
    sel_obj.CV_df['j_corr'] = sel_obj.CV_df['j']

    # reference correction
    sel_obj.CV_df['U_corr'] = sel_obj.CV_df['U_corr'] - sel_obj.ref_to(target_RE)

    # current in capacitance
    if capac:
        sel_obj.CV_df['j_corr'] = sel_obj.CV_df['j_corr'] / (sel_obj.scan_rate / 1000) # scan rate from mV/s to V/s
    else:
        pass

    # current in atomic units
    if atomic:
        sel_obj.CV_df['j_corr'] = sel_obj.CV_df['j_corr'] / atomic_density(
            sel_obj.metal, sel_obj.lattice_plane)  # A/m² / atoms/m² = A/atom
    else:
        pass

    # apply concentration correction
    if conc_corr:
        sel_obj.CV_df['U_corr'] = sel_obj.CV_df['U_corr'] - sel_obj.conc_corr()
    else:
        pass

    U = sel_obj.CV_df['U_corr'].to_numpy()
    j = sel_obj.CV_df['j_corr'].to_numpy()

    return U, j
