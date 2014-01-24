__author__ = 'clyde'

import numpy as np
from ase import units

def get_cis_energy(filename):
    energy = np.nan
    with open(filename) as log_f:
        for line in log_f:
            if 'E(CIS)' in line:
                energy = float(line.split()[4])
    return energy * units.Ha
