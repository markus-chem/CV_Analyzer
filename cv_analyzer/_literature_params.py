# Values are given in V vs. SHE (aqueous medium) and are taken from:
# Gamry web page. Source:
# "Electrochemical Methods: Fundamentals and Applications",
# A J Bard and L R Faulkner, John Wiley & Sons, NY (2000). ISBN 0471405213
# Absolute potential taken from
# "The absolute electrode potential: an explanatory note" (Sergio Trasatti 1986)
# s.a. Wikipedia Absolute Electrode potential ## TODO: link and date.

RE_dict = {
    'SHE':              0.,
    'SCE':              0.241,  # saturated KCl calomel electrode
    'CE_3.5M':          0.25,  # 3.5M KCl
    'CE_1M':            0.28,
    'CE_0.1M':          0.334,
    'Ag/AgCl':          0.197,  # saturated KCl
    'Ag/AgCl_3M':       0.207,  # source: meinsberger-elektroden.de
    'Ag/AgCl_3.5M':     0.205,  # 3.5M KCl
    'Ag/AgCl_0.1M':     0.288,
    'Hg/HgSO4_0.5M':    0.68,  # 0.5M H2SO4
    'Hg/HgSO4_1M':      0.674,  # 1M H2SO4
    'Hg/HgSO4_sat':     0.640,  # saturated K2SO4
    'Hg/HgO_1M_NaOH':   0.140,  # 1M NaOH
    'Hg/HgO/0.1 M NaOH':0.165,  # according to Horswell et al. Langmuir 2004, 20, 10970-10981
    'Hg/HgO_20%_KOH':   0.098,  # 20% KOH
    'Ag/AgSO4_sat':     0.680,  # saturated K2SO4
    'U_abs':            4.44
}


# from SI of 'Grand canonical simulations of electrochemical interfaces in
# implicit solvation models' Nicolas Hoermann
# doi: 10.1063/1.5054580
# point of zero charge in eV for a given metal surface (hkl)
# pc is polycrystalline

pzc_dict = {
    'Pt': {
        '111':      4.71,
        '100[1D]':  4.67,
        '100[2D]':  4.81,
        '110[1D]':  4.66,
        '110[1x1]': 4.56,
        '110[1x2]': 4.55,
        'pc':       4.52
    },
    'Cu': {
        '111':      4.26,
        '111[NaF]': 4.43,
        '100':      4.09,
        '100[NaF]': 4.4,
        '110':      3.85,
        '110[NaF]': 4.37,
        'pc':       4.21
    },
    'Au': {
        '111':                  4.91,
        '111[22xsqrt3]/[1x23]': 4.99,
        '100':                  4.77,
        '100[hex]':             4.99,
        '110':                  4.64,
        '110[1x2]':             4.64,
        'pc':                   4.66
    },
    'Ag': {
        '111':  3.99,
        '100':  3.83,
        '110':  3.69,
        'pc':   3.74
    }
}


# lattice constants from 'PRECISION MEASUREMENTS OF THE LATTICE CONSTANTS OF
# TWELVE COMMON METALS Q/HEELER P. DAVEY', doi: 10.1103/PhysRev.25.753
# in Angstrom

lattice_constants_dict = {
    'fcc':{
        'Cu':   3.597,
        'Ag':   4.08,
        'Pt':   3.912,
        'Au':   4.065
    }
}

# calculated atomic density based on lattice constant
def atomic_density(metal, lattice_plane):
    import math
    lattice_plane = str(lattice_plane)
    a = lattice_constants_dict['fcc'][metal]
    if lattice_plane == '100':
        return (2 / a**2) * 10**16 # in atoms/cm²
    if lattice_plane == '110':
        return (math.sqrt(2) / a**2) * 10**16 # in atoms/cm²
    if lattice_plane == '111':
        return (4 / (math.sqrt(3) * a**2)) * 10**16 # in atoms/cm²
    else:
        raise ValueError(f'{lattice_plane} not listed')
