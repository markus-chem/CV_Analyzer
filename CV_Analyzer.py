# CV_Analyzer class for datapackages
# for the analysis and the manipulation of digitized CVs

from literature_parameters import RE_dict, pzc_dict, lattice_constants_dict, atomic_density
from datapackage import Package
from scipy import constants
import matplotlib.pyplot as plt
import pandas as pd
import os
import sys
import numpy as np

class CV_Analyzer_class:
    def __init__(self, p, resource_no): # initialize with basic information of package
        self.package = p
        # self.CV = p.get_resource(p.resource_names[resource_no]).read(keyed=True)
        self.CV_df = pd.read_csv(p.resources[0].raw_iter(stream=False))
        # pd.DataFrame(data=self.CV).apply(pd.to_numeric, downcast='float') # to pandas dataframe, convert from Decimal to Float
        self.name = p.resource_names[resource_no] # the name of the resource
        self.metadata = p.descriptor
        self.RE = p.descriptor['electrochemical system']['electrodes']['reference electrode']['type']
        self.metal = p.descriptor['electrochemical system']['electrodes']['working electrode']['material']
        self.lattice_plane = p.descriptor['electrochemical system']['electrodes']['working electrode']['crystallographic orientation']
        self.electrolyte_name = [p.descriptor['electrochemical system']['electrolyte']['components'][i]['name'] for i in [0,3]] # assumes that there are only 2 components
        self.electrolyte_conc = [p.descriptor['electrochemical system']['electrolyte']['components'][i]['concentration']['value'] for i in [0,3]]
        self.electrolyte_unit = [p.descriptor['electrochemical system']['electrolyte']['components'][i]['concentration']['unit'] for i in [0,3]]
        self.electrolyte = [f'{self.electrolyte_conc[i]} {self.electrolyte_unit[i]} {self.electrolyte_name[i]}' for i in range(2)]
        self.scan_rate = p.descriptor['figure description']['scan rate']['value'] / 1000 # scan rate in V/
        self.label = f'{self.metal}({self.lattice_plane}) in {self.electrolyte}; {self.name}'

    def ref_to(self, target_RE):
        # reference to another RE

        # check if target_RE exists
        if target_RE in RE_dict or target_RE=='pzc' or target_RE=='RHE':
            pass
        else:
            print('Select one of the following RE')
            print('pzc')
            for key in RE_dict:
                print(key)
            sys.exit('target reference is unknown')

        # check for pH dependency
        if target_RE == 'RHE':
            pH = self.metadata['electrochemical system']['electrolyte']['ph']['value']
            if pH == None:
                sys.exit('No pH value given. Conversion to RHE not possible.')
            else:
                offset = RE_dict[str(target_RE)] - RE_dict[str(self.RE)] - 0.059*pH # 59 mV shift per pH unit
        if target_RE == 'pzc':
            offset = pzc_dict[self.metal][self.lattice_plane] - RE_dict[str(self.RE)] - RE_dict['U_abs']
        else:
            offset = RE_dict[str(target_RE)] - RE_dict[str(self.RE)]

        return offset

    def plot(self):
        # plot the non-manipulated CV data
        fig = plt.figure('Original CV')
        plt.plot(self.CV_df['U'], self.CV_df['j'], label=self.label)
        plt.xlabel(f'U / V vs. [original RE]')
        plt.ylabel('j / A/m$^2$')

        return fig

    def reference_electrode(self, target_RE='SHE'):
        # convert to another voltage reference
        # absoulte potentials of other reference electrodes vs. SHE are in reference_electrode_dict

        self.CV_df['U'] = self.CV_df['U'] - self.ref_to(target_RE)
        self.CV_df['C'] = self.CV_df['j'] / self.scan_rate
        print(f'Reference changed from {self.RE} to {target_RE}')

        #plot
        fig = plt.figure(f'CV referenced to {target_RE}')
        plt.plot(self.CV_df['U'], self.CV_df['C'], label=self.label)
        plt.xlabel(f'U / V vs. {target_RE}')
        plt.ylabel('j / A/m$^2$')

        return fig

    def site_norm(self, target_RE, unit='atomic'):
        # normalize current on lattice sites (in atoms/m²)
        self.CV_df['U'] = self.CV_df['U'] - self.ref_to(target_RE) # reference to another RE
        norm_factor = atomic_density(self.metal, self.lattice_plane) * 10**20 # from atoms/Angstrom² to atoms/m²
        self.CV_df['j'] = self.CV_df['j'] / norm_factor # A/m² / atoms/m² = A/atom

        # atomic units for current
        if unit == 'atomic': # convert into e- / (s*atom)
            self.CV_df['j'] = self.CV_df['j'] / constants.value('elementary charge')
            print(f'Current normalized to atomic density of {self.metal} ({self.lattice_plane}) ( e-/(s*atom) ).')
        else:
            print(f'Current normalized to atomic density of {self.metal} ({self.lattice_plane}) (A/atom).')

        # plot
        fig = plt.figure(f'Current normalized to atomic density')
        plt.plot(self.CV_df['U'], self.CV_df['j'], label=self.label)
        plt.xlabel(f'U / V vs. {target_RE}')
        if unit == 'atomic':
            plt.ylabel('j / e-/(s*atom)')
        else:
            plt.ylabel('j / A/atom')
        
        return fig

    def charge_int(self, lower_lim, upper_lim, target_RE, unit='atomic'):
        import scipy.integrate as integrate

        # integrate the charge within the given voltage limits (ref. to target_RE)
        self.CV_df['U'] = self.CV_df['U'] - self.ref_to(target_RE) # reference to another RE

        # check voltage limits
        if lower_lim < min(self.CV_df['U']) or upper_lim > max(self.CV_df['U']):
            sys.exit(f'Voltage limits out of range for {self.name}')

        # atomic units
        if unit == 'atomic':
            norm_factor = atomic_density(self.metal, self.lattice_plane) * 10**20 # from atoms/Angstrom² to atoms/m²
            self.CV_df['j'] = self.CV_df['j'] / constants.value('elementary charge') / norm_factor # A/m² / C/e- / atoms/m² = e-/(s*atom)

        # seperate anodic and cathodic scan, apply voltage limits
        forw_scan = self.CV_df.where((self.CV_df['j'] > 0) & (self.CV_df['U'] > lower_lim) & (self.CV_df['U'] < upper_lim))
        backw_scan = self.CV_df.where((self.CV_df['j'] < 0) & (self.CV_df['U'] > lower_lim) & (self.CV_df['U'] < upper_lim))
        int_voltage_forw = forw_scan['U'] - lower_lim # first point will be U=0V
        int_voltage_backw = backw_scan['U'] - lower_lim

        # integrate
        Q_forw = int_voltage_forw * forw_scan['j'] / self.scan_rate
        Q_forw_sum = Q_forw.sum() # in C/m² or e-/atom (if unit = atomic)
        Q_forw_int = [(forw_scan.iloc[:i]['j'] * int_voltage_forw[:i]).sum() for i in range(len(forw_scan))]
        Q_backw = int_voltage_backw * backw_scan['j'] / self.scan_rate
        Q_backw_sum = Q_backw.sum() # in C/m² or e-/atom (if unit = atomic)
        Q_backw_int = [(backw_scan.iloc[:i]['j'] * int_voltage_backw[:i]).sum() for i in range(len(forw_scan))]

        # prints
        if unit == 'atomic':
            print(f'forward scan: Q={round(Q_forw_sum, 2)} e-/atom for {lower_lim} < U < {upper_lim}')
            print(f'backward scan: Q={round(Q_backw_sum, 2)} e-/atom for {lower_lim} < U < {upper_lim}')
        else:
            print(f'backward scan: Q={round(Q_backw_sum, 2)} C/m² for {lower_lim} < U < {upper_lim}')
            print(f'forward scan: Q={round(Q_forw_sum, 2)} C/m² for {lower_lim} < U < {upper_lim}')

        # plot forward and backward scan
        fig = plt.figure('integrated CV')
        plt.plot(self.CV_df['U'], self.CV_df['j'], label=self.label)
        plt.plot(forw_scan['U'], forw_scan['j'])
        plt.plot(backw_scan['U'], backw_scan['j'])
        plt.vlines(lower_lim, min(self.CV_df['j']), max(self.CV_df['j']), color='grey', linestyles='solid')
        plt.vlines(upper_lim, min(self.CV_df['j']), max(self.CV_df['j']), color='grey', linestyles='solid')
        plt.xlabel(f'U / V vs. {target_RE}')
        if unit == 'atomic':
            plt.ylabel('j / e-/(s*atom)')
        else:
            plt.ylabel('j / A/m²')

        # evaluate Q over voltage
        fig2 = plt.figure('charge integration forward')
        plt.plot(forw_scan['U'], Q_forw_int, label=f'{self.label} forward scan')
        plt.xlabel(f'U / V vs. {target_RE}')
        if unit == 'atomic':
            plt.ylabel('Q / e-/atom')
        else:
            plt.ylabel('Q / C/m²')

        fig3 = plt.figure('charge integration backward')
        plt.plot(backw_scan['U'], Q_backw_int, label=f'{self.label} backward scan')
        plt.xlabel(f'U / V vs. {target_RE}')
        if unit == 'atomic':
            plt.ylabel('Q / e-/atom')
        else:
            plt.ylabel('Q / C/m²')
        
        return fig, fig2, fig3


    # def peaks(self, lower_lim, upper_lim, target_RE):
    #     from scipy.signal import find_peaks
        
    #     self.CV_df['U'] = self.CV_df['U'] - self.ref_to(target_RE) # reference to another RE

    #     # check voltage limits
    #     if lower_lim < min(self.CV_df['U']) or upper_lim > max(self.CV_df['U']):
    #         sys.exit('Voltage limits out of range')

    #     # seperate anodic and cathodic scan, apply voltage limits
    #     forw_scan = self.CV_df.where((self.CV_df['j'] > 0) & (self.CV_df['U'] > lower_lim) & (self.CV_df['U'] < upper_lim))
    #     backw_scan = self.CV_df.where((self.CV_df['j'] < 0) & (self.CV_df['U'] > lower_lim) & (self.CV_df['U'] < upper_lim))

    def max_min(self, lower_lim, upper_lim, target_RE):
        self.CV_df['U'] = self.CV_df['U'] - self.ref_to(target_RE) # reference to another RE

        # check voltage limits
        if lower_lim < min(self.CV_df['U']) or upper_lim > max(self.CV_df['U']):
            sys.exit(f'Voltage limits out of range for {self.name}')
        
        idx = (self.CV_df['j'].idxmax(), self.CV_df['j'].idxmin())
        max_ = self.CV_df.iloc[idx[0]]
        min_ = self.CV_df.iloc[idx[1]]

        fig = plt.figure('Max_Min of CV')
        plt.plot(self.CV_df['U'], self.CV_df['j'], label=self.label)
        plt.scatter(max_[1], max_[2], marker='X', color='r')
        plt.scatter(min_[1], min_[2], marker='X', color='r')
        plt.xlabel(f'U / V vs. {self.RE}')
        plt.ylabel('j / A/m$^2$')

        print(f'Maximum: U = {round(max_[1], 2)}, j = {round(max_[2], 2)} \nMinimum U = {round(min_[1], 2)}, j = {round(min_[2], 2)}')

        return fig

def filter(files, metal, lattice_plane, component, author_name, not_this_name):
    selected = set(files)
    for i in files:
        p = Package(i)
        entry = CV_Analyzer_class(p, 0) # resource_no 0 as default
        
        if len(metal) > 0:
            if entry.metadata['electrochemical system']['electrodes']['working electrode']['material'] not in metal:
                try:
                    selected.remove(i)
                except:
                    pass
            else:
                pass

        if len(lattice_plane) > 0:
            if entry.metadata['electrochemical system']['electrodes']['working electrode']['crystallographic orientation'] not in lattice_plane:
                try:
                    selected.remove(i)
                except:
                    pass
            else:
                pass

        if len(component) > 0:
            l = list(set(component).intersection([entry.metadata['electrochemical system']['electrolyte']['components'][i]['name'] for i in range(4)]))
            if len(l) == 0:
                    try:
                        selected.remove(i)
                    except:
                        pass
            # for j in component:
            #     if j not in [entry.metadata['electrochemical system']['electrolyte']['components'][i]['name'] for i in range(4)]:
            #         try:
            #             selected.remove(i)
            #         except:
            #             pass
            # else:
            #     pass

        if len(author_name) > 0:    
            if any(entry.metadata['source']['bib'].find(author_name[j]) == -1 for j in range(len(author_name))):
                    try:
                        selected.remove(i)
                    except:
                        pass
            else:
                pass
        
        if any(entry.metadata['source']['bib'].find(not_this_name[j]) != -1 for j in range(len(not_this_name))):
                try:
                    selected.remove(i)
                except:
                    pass
        else:
            pass
            

    if len(selected) > 0:
        print(f'{len(selected)} datapackages selected')
        print(selected)
    else:
        sys.exit('No datapackages meet filter criteria.')

    return selected