"""
Read/write functions for Gaussian.
Written by:

   Glen R. Jenness
   University of Wisconsin - Madison

See accompanying license files for details.
"""

import os
import numpy as np

import ase.units
import cclib.parser.utils as utils

from ase.atoms import Atoms
from ase.atom import Atom
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io.gaussian_reader import GaussianReader as GR

# http://www.gaussian.com/g_tech/g_ur/k_dft.htm
allowed_methods = ['lsda',  # = 'svwn'
                           'svwn',
                           'svwn5',  # != 'svwn'
                           'blyp',
                           'b3lyp',
                           'bp86',
                           'pbepbe',
                           'pbe1pbe',  # pbe0
                           'm06',
                           'm06hf',
                           'm062x',
                           'tpssh',
                           'tpsstpss',
                           'wb97xd',
                           'uff',
                           'huckel',
                          ]


def read_gaussian_out(filename, index=-1, quantity='atoms', quantities=None):
    """"Interface to GaussianReader and returns various quantities"""
    energy = 0.0

    try:
        data = GR(filename)[index]
    except IndexError:
        return {}

    positions = np.array(data['Positions'])
    numbers = np.array(data['Atomic_numbers'])

    method = data['Method']
    version = data['Version']

    if method.lower()[1:] in allowed_methods:
        method = 'HF'

    atoms = Atoms(positions=positions, numbers=numbers)

    for key, value in data.items():
        if (key in method):
            #covers case where calculation is a scan and hence energy is a list of values not a single one
            try:
                energy = value[-1]
            except TypeError:
                energy = value
    try:
# Re-read in the log file
        f = open(filename, 'r')
        lines = f.readlines()
        f.close()

        forces = list()
        for n, line in enumerate(lines):
            if 'oniom' in method.lower() and 'extrapolated energy' in line:
                energy = float(line.split()[4])
            elif not energy and 'SCF Done' in line:
                energy = float(line.split()[4])

            if ('Forces (Hartrees/Bohr)' in line):
                for j in range(len(atoms)):
                    forces += [[float(lines[n + j + 3].split()[2]),
                                float(lines[n + j + 3].split()[3]),
                                float(lines[n + j + 3].split()[4])]]
        #convert = ase.units.Hartree / ase.units.Bohr
        convert = utils.convertor(1, 'hartree', 'eV') / utils.convertor(1, 'bohr', 'Angstrom')
        forces = np.array(forces) * convert
    except:
        forces = None

    #energy *= ase.units.Hartree  # Convert the energy from a.u. to eV
    energy = utils.convertor(energy, 'hartree', 'eV')

    calc = SinglePointCalculator(atoms, energy=energy, forces=forces)

    calc.path = os.path.split(filename)[0]
    calc.label = os.path.split(filename)[1].replace('.log', '')
    calc.data = data

    atoms.set_calculator(calc)

    data_dict = {'energy': energy, 'forces': forces, 'dipole': data.get('Dipole'), 'atoms': atoms, 'version': version}

    if quantities:
        return [data_dict[q] for q in quantities]
    elif quantity in data_dict:
        return data_dict[quantity]
    else:
        return None

def read_gaussian(filename):
    """Reads a Gaussian input file"""
    f = open(filename, 'r')
    lines = f.readlines()
    f.close()

    atoms = Atoms()
    for n, line in enumerate(lines):
        if ('#' in line):
            i = 0
            while (lines[n + i + 5] != '\n'):
                info = lines[n + i + 5].split()
                symbol = info[0]
                position = [float(info[1]), float(info[2]), float(info[3])]
                atoms += Atom(symbol, position=position)
                i += 1
    return atoms
