__author__ = 'clyde'

import re

from ase import Atoms, Atom
from ase.calculators.gaussian import Gaussian


## currently inactive would be nice to be able to read input (e.g. when we need to use =input only)

list_gaussian_methods = ['B3LYP', 'PM3', 'PM6', 'ONIOM(', 'HF', 'AM1']

def extract_method_basis(route_str):
    return


def extract_calc(inp_file):
    """Constructs an ase atoms object with an attached Gaussian calculator from a gaussian input file
    assumes the input file has format Route.... \n\nTitle...\n\ncoord_sect...\n\nconnect_sect"""

    with open(inp_file, 'r') as fl:
        data = fl.read()

    section_names = ['route', 'title', 'coords', 'connects']
    sections = data.split('\n\n')

    no_sections = len(sections)
    data_dict = dict(zip(section_names[:no_sections], sections))

    label = data_dict.get('title')
    route = data_dict.get('route')
    coords = data_dict.get('coords')
    #connects = data_dict.get('connects')

    if coords:
        charge, mult = coords.split('\n')[0].split()
        coords = '\n'.join(coords.split('\n')[1:])

    atoms = Atoms()
    if coords:
        for line in coords.split('\n'):
            info = line.split()
            is_bq = re.match('bq', info[0], re.I)

            # H -> H, H(ISO=2) -> H
            symbol = info[0].split('(')[0]

            position = [float(info[1]), float(info[2]), float(info[3])]
            #ignoring ghost atoms
            if not is_bq:
                atoms += Atom(symbol.split('(')[0], position=position)

        if not label:
            label = 'Default Gaussian Label'

    route_args = route.split()[1:]
    if 'oniom' in route_args[0]:
        method = route_args[0]
        basis = 'oniom'
        route_args = route_args[1:]
    else:
        method, basis = route_args[0:2]
        route_args = route_args[2:]

    route_kwargs = {'method':method, 'basis':basis}
    for route_arg in route_args:
        if  re.match('IOP', route_arg, re.I):
            route_kwargs.update({'ioplist': route_arg[4:-1].split(',')})
        elif '=' in route_arg:
            route_arg_nm, route_arg_vl = route_arg.split('=')
            route_kwargs.update({route_arg_nm: route_arg_vl})
        else:
            route_kwargs.update({route_arg:route_arg})
        atoms.set_calculator(Gaussian(label=label, multiplicity=mult, **route_kwargs))

    return atoms