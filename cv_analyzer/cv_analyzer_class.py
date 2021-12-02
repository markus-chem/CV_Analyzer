import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math
from scipy import constants

from ._literature_params import RE_dict, pzc_dict, atomic_density

class CV_Analyzer:
    def __init__(self, package, resource_no):
        """
        Initialize a datapackage with basic information defined in the metadata
        (.yaml file)
        One CV per datapackage

        Inputs:
        package:        datapackage (https://pypi.org/project/datapackage/)
        resource_no:    integer, resource bumber (typically only 1 .csv file)
                        within the larger datapackage

        Outputs:
        None
        """
        self.package = package
        self.CV_df = pd.read_csv(package.resources[0].raw_iter(stream=False))
        self.CV_df['j'] = self.CV_df['j'] * 100 # change unit from A/m² to µA/cm²
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

        electrolyte = e_chem_sys['electrolyte']
        elec_comps  = electrolyte['components']

        # exclude solvent and empty fields 
        self.electrolyte_name   = [
            elec_comps[i]['name'] \
            for i in range(len(elec_comps)) if \
            (elec_comps[i]['type'] != 'solvent') & (elec_comps[i]['name'] != None)
            ]
        
        self.c_electrolyte      = [
            elec_comps[i]['concentration']['value'] \
            for i in range(len(elec_comps)) if \
            (elec_comps[i]['type'] != 'solvent') & (elec_comps[i]['name'] != None)
            ]
        
        self.electrolyte_unit   = [
            elec_comps[i]['concentration']['unit'] \
            for i in range(len(elec_comps)) if \
            (elec_comps[i]['type'] != 'solvent') & (elec_comps[i]['name'] != None)
            ]

        self.electrolyte        = [
            "{} {} {}".format(self.c_electrolyte[i], self.electrolyte_unit[i],
                                self.electrolyte_name[i])
                                for i in range(len(self.c_electrolyte))
        ]

        self.T                  = electrolyte['temperature']['value']
        self.pH                 = electrolyte['ph']['value']

        fig_desc = metadata['figure description']
        self.scan_rate_unit = fig_desc['scan rate']['unit']
        self.scan_rate = fig_desc['scan rate']['value']

        # change scan rate to V / s if necessary
        if 'mV' in self.scan_rate_unit:
            self.scan_rate = self.scan_rate / 1000
            self.scan_rate_unit = 'V / s'

        # general label for legends in plots
        self.label = "{}({}) in {} at {} {}; {}".format(
            self.metal, self.hkl, self.electrolyte, self.scan_rate,
            self.scan_rate_unit, self.name
        )
        self.metadata = metadata

        return None

    def plot_orig(self):
        """
        Recreate the cyclic voltammogram with the original experimental
        parameters. Possibly different reference electrodes on the same
        axis.

        Inputs:
        None
        
        Outputs:
        fig (matplotlib object)
        """
        fig = plt.figure('Original CV')
        plt.plot(self.CV_df['U'], self.CV_df['j'], label=self.label)
        plt.xlabel(f'U / V vs. [original RE]')
        plt.ylabel('j / µA/cm$^2$')

        return fig

    def plot(self, target_RE='SHE', C_exp=False, atomic=False, c_corr=False):
        """
        Generic plot function for digitized CVs.
        Current and voltage corrections can be applied:
            - reference to another electrode (target_RE)
                For RHE the pH is considered (and calculated from
                components if not given in metadata)
            - normalize on the scan rate and plot capacitance
            - normalize area on atomic sites of the respective hkl
            - apply a Nernstian shift (59 mV/decade) for the halide
                (F-, Cl-, Br-, I-) concentration and normalize to 1M
        
        Corrections can be combined.

        Inputs:
        target_RE:      str, name of target reference electrode (s. RE_dict)
        C_exp:          bool, capacitance plot
        atomic:         bool, normalization on atomic surface sites
        c_corr,         bool, Nernst shift for halide concentration to 1M

        Outputs:
        fig (matplotlib object)
        """

        settings = f'C_exp={C_exp}, atomic={atomic}, c_corr={c_corr}, target_RE={target_RE}'

        # extra columns for corrections
        self.CV_df['U_corr'] = self.CV_df['U']
        self.CV_df['j_corr'] = self.CV_df['j']

        # reference correction
        self.CV_df['U_corr'] = self.CV_df['U_corr'] - self.ref_to(target_RE)

        # current to capacitance
        if C_exp:
                self.CV_df['j_corr'] = self.CV_df['j_corr'] / self.scan_rate

        # current in atomic units: µA/cm² / atoms/cm² = µA/atom
        if atomic:
            rho_atomic = atomic_density(self.metal, self.hkl)
            self.CV_df['j_corr'] = self.CV_df['j_corr'] / rho_atomic

        # apply concentration correction
        if c_corr:
            self.CV_df['U_corr'] = self.CV_df['U_corr'] - self.c_corr()

        # plot
        fig = plt.figure(f'CV_settings: {settings}')
        plt.plot(self.CV_df['U_corr'], self.CV_df['j_corr'], label=self.label)
        plt.xlabel(f'U / V vs. {target_RE}')
        if (C_exp) & (atomic == False):
            plt.ylabel('C / µF/cm$^2$')
        elif (C_exp) & (atomic):
            plt.ylabel('C / µF/atom')
        elif (C_exp == False) & (atomic == False):
            plt.ylabel('j / µA/cm$^2$')
        elif (C_exp == False) & (atomic == True):
            plt.ylabel('j / µA/atom')

        return fig

    def charge_int(
            self,
            lower_lim,
            upper_lim,
            target_RE,
            atomic=False,
            c_corr=False):

        """
        Charge passed within a given potential window
        (scan rate normalized integration of CV)
        Three plots are created:
            - CVs with voltage limits
            - charge integral in forward scan
            - charge integral in backward scan

        Inputs:
        lower_lim:      float, lower voltage limit in V
        upper_lim:      float, upper voltage limit in V
        target_RE:      str, name of target reference electrode (s. RE_dict)
        atomic:         bool, normalization on atomic surface sites
        c_corr:         bool, Nernst shift for halide concentration to 1M

        Outputs:
        fig (matplotlib object)
        fig2 (matplotlib object)
        fig3 (matplotlib object)
        """

        # extra columns for corrections
        self.CV_df['U_corr'] = self.CV_df['U']
        self.CV_df['j_corr'] = self.CV_df['j']

        # reference to target_RE
        self.CV_df['U_corr'] = self.CV_df['U'] - self.ref_to(target_RE)

        # apply concentration correction
        if c_corr:
            self.CV_df['U_corr'] = self.CV_df['U_corr'] - self.c_corr()
        else:
            pass

        # check voltage limits
        if lower_lim < min(
                self.CV_df['U_corr']) or upper_lim > max(self.CV_df['U_corr']):
            raise ValueError(f'Voltage limits out of range for {self.name} !')

        # atomic units
        if atomic:
            norm_factor = atomic_density(self.metal, self.hkl)
            self.CV_df['j_corr'] = self.CV_df['j_corr'] / (constants.e * 10**6) / \
                norm_factor  # µA/cm² / µC/e- / atoms/cm² = e-/(s*atom)

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

    def max_min(self, lower_lim, upper_lim, target_RE, atomic=False, C_exp=False):

        """
        Shows the maximum and minimum current or capacitance within
        a given potential window.

        Inputs:
        lower_lim:      float, lower voltage limit in V
        upper_lim:      float, upper voltage limit in V
        target_RE:      str, name of target reference electrode (s. RE_dict)
        atomic:         bool, normalization on atomic surface sites

        Outputs:
        fig (matplotlib object)
        """

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

        # current to capacitance
        if C_exp:
                self.CV_df['j_corr'] = self.CV_df['j_corr'] / self.scan_rate

        # atomic units
        if atomic:
            norm_factor = atomic_density(
                self.metal, self.hkl)
            self.CV_df['j_corr'] = self.CV_df['j_corr'] / constants.e / \
                norm_factor  # A/m² / C/e- / atoms/m² = e-/(s*atom)

        # find minimum and maximum
        idx = (self.CV_df['j_corr'].idxmax(), self.CV_df['j_corr'].idxmin()) # indices of min and max
        max_ = (self.CV_df.iloc[idx[0]]['U_corr'], self.CV_df.iloc[idx[0]]['j_corr'])  # x and y of max
        min_ = (self.CV_df.iloc[idx[1]]['U_corr'], self.CV_df.iloc[idx[1]]['j_corr'])  # x and y of min

        # plot
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
        Helping function to calculate the voltage offset between
        two reference electrode. Considers references listed in RE_dict.

        Inputs:
        target_RE:  str, name of the target reference electrode

        Output:
        offset:     float, Potential shift that is substracted due to reference change
        """

        # check if target_RE exists
        RE_keys = list(RE_dict.keys()) + ['pzc', 'RHE']
        if target_RE not in RE_keys:
            s = "Target reference is unknown, select one of the following: "
            s += ', '.join(RE_keys)
            raise ValueError(s)
        
        # try to calculate pH based on acid or base components
        if self.pH == None:
            self.pH = self.calc_pH()

        # pH is needed for calculation with RHE
        if 'RHE' in [self.RE, target_RE]:
            if self.pH == None:
                raise ValueError(f'pH value is undefined, RHE not possible ({self.name})')

        # pH shift per decade H+ concentration
        R, F =  constants.R, constants.value('Faraday constant')
        pH_p_dec = R * self.T / F * math.log(10)

        # corrections for the original RE
        if self.RE == 'RHE':
            offset = pH_p_dec * self.pH
        else:
            offset = - RE_dict[str(self.RE)]

        # correction for the target RE
        if target_RE == 'RHE':
            # Add the 59 mV shift per pH unit
            offset -= pH_p_dec * self.pH

        # Sets the potential of zero charge as the reference.
        elif target_RE == 'pzc':
            pzc = pzc_dict[self.metal][self.hkl]
            offset +=  pzc - RE_dict['U_abs']
        else:
            offset += RE_dict[str(target_RE)]

        return offset


    def c_corr(self):
        """
        Helping function to calculate Nernstian shift (59 mV/dec)
        for the halide concentration (F-, Cl-, Br-, I-).
        Normalization on 1M halide concentration.

        Inputs:
        None

        Outputs:
        offset: float, the Nernstian potential shift that is substracted
        """

        # list of typical halide adsorbates
        ads_set = {'LiCl', 'NaCl', 'KCl', 'RbCl', 'CsCl',
                'LiBr', 'NaBr', 'KBr', 'RbBr', 'CsBr',
                'LiI', 'NaI', 'KI', 'RbI', 'CsI'
                }
        # get list of common elements
        ads = list(ads_set.intersection(self.electrolyte_name))

        if len(ads) > 1:
            raise ValueError(f'Cannot handle more than one adsorbate: {ads}')
        elif len(ads) == 0:
            U_shift = 0
        else:
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

    def calc_pH(self):
        """
        Calculates the pH based on the components.
        Only for one pH sensitive component.
        Only for strong acids / bases listed below.
        If pH cannot be calculated, set to pH = 7
        
        Inputs:
        None
        
        Outputs:
        pH: float, pH value
        """

        acids_1H = [
            'HCl', 'HClO4', 
        ]

        acids_2H = [
            'H2SO4'
        ]

        base_1 = [
            'LiOH', 'NaOH', 'KOH', 'CsOH'
        ]

        base_2 = [
            'MgOH2', 'CaOH2'
        ]

        for idx, i in enumerate(self.electrolyte_name):
            # change units to mol/L
            # assuming that only 'M' and 'mM' are used
            # assuming only one pH sensitive species
            if self.electrolyte_unit[idx] == 'mM':
                self.c_electrolyte[idx] = self.c_electrolyte[idx] / 1000
                self.electrolyte_unit[idx] = 'M'
            if i in acids_1H:
                pH = - math.log(self.c_electrolyte[idx]) # concentration in mol/L
            elif i in acids_2H:
                pH = - math.log(2 * self.c_electrolyte[idx])
            elif i in base_1:
                pH = 14 + math.log(self.c_electrolyte[idx])
            elif i in base_2:
                pH = 14 + math.log(2 * self.c_electrolyte[idx])
        
        # set pH=7 if pH not given and not calculable
        try:
            pH
        except:
            pH = 7
            print(f'pH could not be calculated \n assume pH = 7 for {self.name}')
        
        return pH
