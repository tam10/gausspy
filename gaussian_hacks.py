__author__ = 'clyde'

import numpy as np
from ase import units

def get_cis_energy(filename):
    """Returns the energy in eV of a CIS calculation"""
    energy = np.nan
    with open(filename) as log_f:
        for line in log_f:
            if 'E(CIS)' in line:
                energy = float(line.split()[4])
    return energy * units.Ha

def get_methods_wo_stability_opt():
    """Returns a list of methods that will throw an error if stable=opt is attempted"""
    a_stable = ['uff', 'huckel', 'mindo3']

    return a_stable