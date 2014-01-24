"""
Gaussian calculator for ASE modified by:

    Clyde Fare
    Imperial College London

    originally written by:

    Glen R. Jenness
    University of Wisconsin - Madison

Based off of code written by:

    Glen R. Jenness
    Kuang Yu
    Torsten Kerber, Ecole normale superieure de Lyon (*)
    Paul Fleurat-Lessard, Ecole normale superieure de Lyon (*)
    Martin Krupicka

(*) This work is supported by Award No. UK-C0017, made by King Abdullah
University of Science and Technology (KAUST), Saudi Arabia.

See accompanying license files for details.
"""

import os
import glob
import warnings
import pbs
import gaussian_hacks
import numpy as np
import oniom_utils
from subprocess import Popen, PIPE
from ase.io import read as ase_read
from ase.io.gaussian_reader import GaussianReader as GR
from ase import Atoms
from ase.calculators.general import Calculator

"""
Gaussian has two generic classes of keywords:  link0 and route.
Since both types of keywords have different input styles, we will
distinguish between both types, dividing each type into str's, int's
etc.

For more information on the Link0 commands see:
    http://www.gaussian.com/g_tech/g_ur/k_link0.htm
For more information on the route section keywords, see:
    http://www.gaussian.com/g_tech/g_ur/l_keywords09.htm
"""
link0_str_keys = ['chk',
                  'oldchk',
                  'mem',
                  'rwf',
                  'int',
                  'd2e',
                  'lindaworkers',
                  'kjob',
                  'subst',
                  'save',
                  'nosave',
                 ]

link0_int_keys = ['nprocshared',
                  'nproc',
                 ]

# Multiplicity isn't really a route keyword, but we will put it here anyways
route_int_keys = ['multiplicity',
                  'cachesize',
                  'cbsextrapolate',
                  'constants',
                 ]

route_str_keys = ['method',
                  'functional',
                  'basis',
                  'maxdisk',
                  'cphf',
                  'density',
                  'densityfit',
                  'ept',
                  'field',
                  'pop',
                  'guess',
                  'gvb',
                  'integral',
                  'irc',
                  'ircmax',
                  'name',
                  'nmr',
                  'nodensityfit',
                  'output',
                  'punch',
                  'scf',
                  'symmetry',
                  'td',
                  'units',
                 ]

# This one is a little strange.  Gaussian has several keywords where you just
# specify the keyword, but the keyword itself has several options.
# Ex:  Opt, Opt=QST2, Opt=Conical, etc.
# These keywords are given here.
route_self_keys = ['opt',
                   'force',
                   'freq',
                   'geom',
                   'complex',
                   'fmm',
                   'extrabasis',
                   'genchk',
                   'nosymmetry',
                   'polar',
                   'prop',
                   'pseudo',
                   'restart',
                   'scan',
                   'scrf',
                   'sp',
                   'sparse',
                   'stable',
                   'volume',
                  ]

route_float_keys = ['pressure',
                    'scale',
                    'temperature',
                   ]

route_bool_keys = [
                  ]

oniom_coord_keys =['layers',
             'links',
             'link_connections',
             'layer_mults'
            ]

extra_dict_keys = [#a dictionary of ASE atoms used by qst3/qst2 transition state optimisations
                   'component_configs',
                  ]

extra_list_keys = [ #a list of ASE atoms used by oniom
                    'component_calcs',
                    #lists of integers specifying orbitals to be swapped. This is used in conjuction with guess=alter to populate orbitals in a different order to lowest energy first
                    #this is commonly used in the context of CAS to select an active space
                    'original_alpha_orbitals',
                    'swapped_alpha_orbitals',
                    'original_beta_orbitals',
                    'swapped_beta_orbitals',
                    #list of floats specifying weights to assign state in a state averaged cas calculation
                    'state_weights'
                  ]

extra_keys = ['non_std_route',
              'raw_input',
              'extra_input',
              'gen_fchk'
             ]

job_keys = ['nodes',
            'memory',
            'time',
            'queue',
            'version',
            'qid'
           ]

class Gaussian(Calculator):
    """
    Gaussian calculator
    """
    def __init__(self, label='ase', ioplist=list(), basisfile=None,
                 directory=None, **kwargs):
        Calculator.__init__(self)

# Form a set of dictionaries for each input variable type
        self.link0_int_params = dict()
        self.link0_str_params = dict()
        self.route_str_params = dict()
        self.route_int_params = dict()
        self.route_float_params = dict()
        self.route_bool_params = dict()
        self.route_self_params = dict()
        self.oniom_coord_params = dict()
        self.extra_dict_params = dict()
        self.extra_list_params = dict()
        self.extra_params = dict()
        self.job_params = dict()
        self.redundant_coord_params = []

        for key in link0_int_keys:
            self.link0_int_params[key] = None
        for key in link0_str_keys:
            self.link0_str_params[key] = ""
        for key in route_str_keys:
            self.route_str_params[key] = ""
        for key in route_int_keys:
            self.route_int_params[key] = None
        for key in route_float_keys:
            self.route_float_params[key] = None
        for key in route_bool_keys:
            self.route_bool_params[key] = None
        for key in route_self_keys:
            self.route_self_params[key] = ""
        for key in oniom_coord_keys:
            self.oniom_coord_params[key] = None
        for key in extra_dict_keys:
            self.extra_dict_params[key] = {}
        for key in extra_list_keys:
            self.extra_list_params[key] = []
        for key in extra_keys:
            self.extra_params[key] = None
        for key in job_keys:
            self.job_params[key] = None

        self.set(**kwargs)

        self.atoms = None
        self.positions = None
        self.old_positions = None
        self.old_link0_str_params = None
        self.old_link0_int_params = None
        self.old_route_str_params = None
        self.old_route_int_params = None
        self.old_route_float_params = None
        self.old_route_bool_params = None
        self.old_route_self_params = None
        self.old_basisfile = None
        self.old_label = None
        self.old_ioplist = None

        self.basisfile = basisfile
        self.label = label
        self.ioplist = list(ioplist)[:]
        self.directory = directory
        self.multiplicity = 1
        self.converged = False

        #holds parsed Gaussian data
        #ase gaussian parser making use of machine readable data at the end of the log file
        self._data = {}
        #cclib gaussian parser making use of the body of the log file
        self._max_data = {}
        #molmod fchk parser
        self._fchk_data = None
        #fingerprint to check calculation identity
        self._fingerprint = None

        #hold descriptive info on calculation
        self._status = None
        self._notes = None
        self._calc_complete = None

        #holds runtime info
        self.work_folder = ""

        self.check()

    def check(self):
        """Checks Calculation parameters are valid"""
        allowed_label_chars = '0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ_-/'
        warning_label_chars = '/'

        if any(c not in allowed_label_chars for c in self.label):
            raise(RuntimeError('Invalid label choice: {l}'.format(l=self.label)))
        if any(c in warning_label_chars for c in self.label):
            warnings.warn('Non standard label choice: {l}'.format(l=self.label))

    def set_job(self, nodes=1, memory=400, time=3, queue=None, version='direct_g09'):
        self.job_params['nodes'] = nodes
        self.job_params['memory'] = memory
        self.job_params['time'] = time
        self.job_params['queue'] = queue
        self.job_params['version'] = version
        self.job_params['qid'] = None

        self.link0_int_params['nproc'] = nodes
        self.link0_str_params['mem'] = str(memory) + 'MB'

    #todo
    def set_zmatrix(self, atoms, coord='',repetitions=0, delta=0):
        """Converts input to use a zmatrix instead of cartesian coordinates, optionally includes parameters for scanning across a particular coordinate"""
        self._zmatrix = to_zmatrix(atoms.get_positions())
        return

    def set_orbital_swaps(self, original_a_orbitals, swapped_a_orbitals, original_b_orbitals=None, swapped_b_orbitals=None):
        '''Sets orbitals to be swapped: used in conjuction with guess=alter
        Takes lists of integeres specifying the orbitals to be swapped i.e. we swap original_alpha_orbitals[0] for swapped_alpha_orbitals[0], original_alpha_orbitals[1] for swapped_alpha_orbitals[1], etc.
        In the case of a restricted calculation the orbitals specified in the alpha_orbitals lists are taken for both alpha and beta'''

        self.extra_list_params['original_alpha_orbitals'] = original_a_orbitals
        self.extra_list_params['swapped_alpha_orbitals'] = swapped_a_orbitals

        if original_b_orbitals:
            self.extra_list_params['original_beta_orbitals'] = original_b_orbitals

        if swapped_b_orbitals:
            self.extra_list_params['swapped_beta_orbitals'] = swapped_b_orbitals

    def set_state_weights(self, weights):
        self.extra_list_params['state_weights'] = weights

    def set_redundant(self, coord_type='', atom_nos = None, action_type='', min=None, max=None, diag=None, frames=None, degrees=None):
        """Sets redundant internal coordinate parameters:
        coord_type: cartesian/length/angle/dihedral/bend
        action_type: scan/freeze/activate/add/remove/derivatives/diag
        additional keys: frames/degrees (scan), diag (diag), min/max (all)
        For more details see: http://www.gaussian.com/g_tech/g_ur/k_opt.htm"""

        redundant_coord_params = {}
        if not atom_nos:
            atom_nos = []

        redundant_coord_params['coord_type'] = coord_type.lower()
        redundant_coord_params['action_type'] = action_type.lower()
        redundant_coord_params['atom_nos'] = atom_nos

        if frames:
            redundant_coord_params['frames'] = frames
        if degrees:
            redundant_coord_params['degrees'] = float(degrees)
        if min:
            redundant_coord_params['min'] = min
        if max:
            redundant_coord_params['max'] = max
        if diag:
            redundant_coord_params['diag'] = diag

        self.redundant_coord_params.append(redundant_coord_params)

    def clear_redundant(self):
        self.redundant_coord_params = []

    def set_qst3(self, reactants, ts):
        if 'component_configs' not in self.extra_dict_params:
            self.extra_dict_params['component_configs'] = {}

        self.extra_dict_params['component_configs']['reactants'] = reactants
        self.extra_dict_params['component_configs']['ts'] = ts

    def set_qst2(self, reactants):
        if 'component_configs' not in self.extra_dict_params:
            self.extra_dict_params['component_configs'] = {}

        self.extra_dict_params['component_configs']['reactants'] = reactants

    def set(self, **kwargs):
        """Assigns values to dictionary keys"""
        for key in kwargs:
            if key in self.link0_str_params:
                self.link0_str_params[key] = kwargs[key]
            elif key in self.link0_int_params:
                self.link0_int_params[key] = kwargs[key]
            elif key in self.route_str_params:
                self.route_str_params[key] = kwargs[key]
            elif key in self.route_int_params:
                self.route_int_params[key] = kwargs[key]
            elif key in self.route_float_params:
                self.route_float_params[key] = kwargs[key]
            elif key in self.route_bool_params:
                self.route_bool_params[key] = kwargs[key]
            elif key in self.route_self_params:
                self.route_self_params[key] = kwargs[key]
            elif key in self.oniom_coord_params:
                self.oniom_coord_params[key] = kwargs[key]
            elif key in self.extra_list_params:
                self.extra_list_params[key] = kwargs[key]
            elif key in self.extra_params:
                self.extra_params[key] = kwargs[key]

    def initialize(self, atoms):
        if (self.route_int_params['multiplicity'] is None):
            self.multiplicity = 1
        else:
            self.multiplicity = self.route_int_params['multiplicity']

        # Set some default behavior
        if (self.route_str_params['method'] is None):
            self.route_str_params['method'] = 'hf'

        if (self.route_str_params['method'] == 'hf' and self.route_str_params['basis'] is None):
            self.route_str_params['basis'] = '6-31g*'

        self.converged = None
        self._status = ''

        # Set default job options
        if not self.job_params['nodes']:
            self.set_job()

        #set run time options
        current_fold = os.getcwd()
        try:
            local_fold = os.getcwd().split(os.environ['ASE_HOME'])[1]
        except IndexError:
            raise RuntimeError('Not running from within ASE_HOME')
        self.work_fold = os.environ['GAUSS_SCRATCH'] + local_fold + '/'

        self.link0_str_params['chk'] = self.work_fold + self.label + '.chk'

        if self.link0_str_params['oldchk'] and not os.path.isabs(self.link0_str_params['oldchk']):
            self.link0_str_params['oldchk'] = self.work_fold + self.link0_str_params['oldchk']

        #handle redundant input e.g. if we are specifying component calculations to load in we must also be reading the guess from them
        if self.extra_list_params['component_calcs'] and not self.route_str_params['guess']:
            self.route_str_params['guess'] = 'input'

    def _get_link0(self):
        link0 = ''
        for key, val in self.link0_str_params.items():
            if val:
                link0 += '%%%s=%s\n' % (key, val)

        for key, val in self.link0_int_params.items():
            if val:
                link0 += '%%%s=%i\n' % (key, val)

        return link0

    def _get_route(self):
        """returns the route section of the gaussian input file"""

        #setting extra_params['non_std_route'] allows us to read in non-standard route (note that for oniom calcs we still need to include oniom in self.route_str_params['method'] even though this will be taken care of
        #by our custom route because the section governing atomic coordinates needs to know whether we are performing an oniom calculation or not and checks the route_str_params)
        if self.extra_params['non_std_route']:
            return '# NonStd\n' + self.extra_params['non_std_route']

        #By default we will always use "#p gfprint pop=full" to start.
        if 'oniom' in self.route_str_params['basis'].lower() or not self.route_str_params['basis']:
            route = '#p gfprint %s' % (self.route_str_params['method'])
        else:
            route = '#p gfprint %s/%s' % (self.route_str_params['method'],
                                  self.route_str_params['basis'])

        # Add keywords and IOp options
        # For the 'self' keywords, there are several suboptions available, and if more
        # than 1 is given, then they are wrapped in ()'s and separated by a ','.
        for key, val in self.route_self_params.items():
            if val:
                if (val == key):
                    route += (' ' + val)
                else:
                    if ',' in val:
                        route += ' %s(%s)' % (key, val)
                    else:
                        route += ' %s=%s' % (key, val)

        for key, val in self.route_float_params.items():
            if val:
                route += ' %s=%f' % (key, val)

        for key, val in self.route_int_params.items():
            if val and (key != 'multiplicity'):
                route += ' %s=%i' % (key, val)

        for key, val in self.route_str_params.items():
            if val and (key != 'method') and\
               (key != 'basis'):
                route += ' %s=%s' % (key, val)

        for key, val in self.route_bool_params.items():
            if val:
                route += ' %s=%s' % (key, val)

        if (self.ioplist):
            route += ' IOP(' + ' ,'.join(self.ioplist) + ')'

        return route

    def _get_mol_details(self, atoms):
        """returns the atomic configuration of the gaussian input file"""

        if 'allcheck' in self.route_self_params['geom'].lower():
            return ''

        if 'oniom' in self.route_str_params['method']:
            return self._get_oniom_details(atoms)

        mol_details = ''

        charge = sum(atoms.get_charges())
        mol_details += '%i %i\n' % (charge, self.multiplicity)

        if 'check' in self.route_self_params['geom'].lower():
            return mol_details

        symbols = atoms.get_chemical_symbols()
        coordinates = atoms.get_positions()

        for i in range(len(atoms)):
            mol_details += '%-10s' % symbols[i]
            for j in range(3):
                mol_details += '%20.10f' % coordinates[i, j]
            mol_details += '\n'
        mol_details += '\n'

        return mol_details

    #oniom allows us to read in separate guess wavefunctions for the components of the oniom calculation, we add in these calculations as a 'component_calcs' parameter
    def _get_comp_chks(self):
        comp_calcs = self.extra_list_params['component_calcs']
        if self.route_str_params['guess'] != 'input' or not comp_calcs:
            return ''

        if isinstance(comp_calcs, bool):
            comp_calcs = [self.label + '_low_real', self.label + '_high_model', self.label + '_low_model']

        scratch_dir = self._get_scratch_dir()
        chk_commands = [scratch_dir + '/' + c.calc.label+'.chk' if c else 'generate' for c in self.extra_list_params['component_calcs']]
        return '\n\n'+'\n\n'.join(chk_commands)


    @property
    def log(self):
        current_dir = os.getcwd()
        home_dir = os.environ['ASE_HOME']
        work_dir = os.environ['ASE_SCRATCH']

        try:
            log = work_dir + current_dir.split(home_dir)[1] + '/' + self.label + '.log'
        except IndexError:
            raise RuntimeError('Not running from within ASE_HOME')

        return log

    @property
    def fchk(self):
        current_dir = os.getcwd()
        home_dir = os.environ['ASE_HOME']
        work_dir = os.environ['ASE_SCRATCH']

        try:
            fchk = work_dir + current_dir.split(home_dir)[1] + '/' + self.label + '.fchk'
        except IndexError:
            raise RuntimeError('Not running from within ASE_HOME')

        return fchk

    @property
    def method(self):
        return self.route_str_params['method']

    @method.setter
    def method(self, value):
        self.route_str_params['method'] = value

    @property
    def basis(self):
        return self.route_str_params['basis']

    @basis.setter
    def basis(self, value):
        self.route_str_params['basis'] = value

    @property
    def extra_input(self):
        return self.extra_params['extra_input']

    @extra_input.setter
    def extra_input(self, value):
        self.extra_params['extra_input'] = value

    #todo
    @property
    def active(self):
        return True
#        if self.job_params['qid'] in [e.id for e in pbs.qstat_plain()]:
#            return True
#        else:
#            return False
    #todo
    @property
    def job(self):
        return True
        #try:
        #    return next(e for e in pbs.qstat_plain() if self.job_params['qid'] == e.id)
        #except StopIteration:
        #    return pbs.JobStatus('','')

    def _get_oniom_details(self, atoms):
        if 'check' in self.route_self_params['geom'].lower():
            return '' + self._get_comp_chks()

        oniom_method = self.route_str_params['method'].split('oniom(')[1][0:-1]
        no_layers = len(oniom_method.split(':'))

        if not self.oniom_coord_params['layers'] or len(self.oniom_coord_params['layers']) != no_layers - 1 or no_layers < 2:
            raise RuntimeError('Incorrect specification of layers for Oniom calculation')
        elif all(isinstance(l, Atoms) for l in self.oniom_coord_params['layers']):
            list_atoms = [a for a in atoms]
            try:
                layers = [[list_atoms.index(a) for a in l] for l in self.oniom_coord_params['layers']]
            except ValueError:
                raise RuntimeError('Atoms in a layer do not match any atoms in the total molecule - check the coordinates of your fragment')
        else:
            layers = self.oniom_coord_params['layers']

        if not self.oniom_coord_params['links']:
            links = self.get_links(atoms, layers)
        else:
            links = self.oniom_coord_params['links']

        if not self.oniom_coord_params['link_connections']:
            link_cons = self.get_link_cons(atoms, layers, links)
        else:
            link_cons = self.oniom_coord_params['link_connections']

        if len(link_cons) != len(links):
            raise RuntimeError('Oniom calculation attempted with different number of link atoms and link connections')

        #everything not specified by layers gets assigned to 'L'
        if no_layers == 2:
            layer_chars = ['H']
        elif no_layers == 3:
            layer_chars = ['H', 'M']
        else:
            raise NotImplementedError('Oniom calculations only implemented for 2 or 3 layers')

        layer_charge_strs = []
        atom_coord_strs = []

        #set layer specifid variables
        for i in range(no_layers):
            if i != no_layers -1:
                atom_indexes = layers[i]
            else:
                #last layer is 'L' and includes all atoms not earmarked for a higher layer
                atom_indexes =[ind for ind in range(len(atoms)) if ind not in [e for l in layers for e in l]]

            layer_atoms = atoms[atom_indexes]
            layer_charge = sum(layer_atoms.get_charges())
            layer_multiplicity = self.oniom_coord_params['layer_mults'][i]
            #should check layer_charge is an integer
            layer_charge_strs.append("{c} {m} ".format(c=int(layer_charge), m=layer_multiplicity))

        for i, a in enumerate(atoms):
            a_symbol = a.symbol
            a_layer = 'L'
            a_link_con = None


            for j in range(no_layers):
                #last layer always 'L'
                if j!= no_layers-1 and i in layers[j]:
                    a_layer = layer_chars[j]
                    #first layer always 'H' hence never has link atoms
                if j>0 and i in links[j-1]:
                    link_ind = links[j-1].index(i)
                    a_link_con = link_cons[j-1][link_ind]

            line_strs = ['%-10s' % a_symbol]
            #indicates the atom is free, should really make this a variable
            line_strs.append('%-2s' % 0)
            line_strs += ['%20.10f' % a.position[k] for k in range(3)]
            line_strs.append(' %-2s' % a_layer)

            #link_connections are the index of the higher level atom that the link connects to, if it's the first atom i.e. index 0 this if statement will fail because
            #of the boolean nature of 0 in python so we have to consider that situation specifically
            if a_link_con or a_link_con == 0:
                line_strs.append('%-5s' % 'H')
                #Gaussian uses Fortran style arrays
                line_strs.append('%-1s' % (a_link_con+1))

            atom_coord_strs.append("".join(line_strs))

        layer_coord_strs = "\n".join(atom_coord_strs)

        return " ".join(layer_charge_strs) + "\n" + "".join(layer_coord_strs) + self._get_comp_chks()

    #as these make use of molmods neighbor function care must be taken when building transition structures
    #i.e. if we make a bond very long molmod won't recognise the two atoms as neighbours - this isn't usually a problem as we don't normally cut between atoms that are involved
    # in a transition state
    def get_links(self, atoms, layers):
        from ASE_utils import to_molmod
        mm_atoms = to_molmod(atoms)

        links = []
        for i in range(len(layers)):
            neighbour_as = [mm_atoms.graph.neighbors[layer_a] for layer_a in layers[i]]
            #unwind and remove redundancies
            neighbour_as = set([a for s in neighbour_as for a in s])
            #equivalent to [n for n in neighbour_as if n not in layers[i]] for i =0,
            #if multiple layers ignore neighbours in all previous layers to avoid double counting
            non_link_neighbours = [n for l in layers[:i+1] for n in l]
            link_neighbours = [a for a in neighbour_as if a not in non_link_neighbours]

            if not link_neighbours:
                raise RuntimeError('No link atoms defined in layer {l}'.format(l=i))
            else:
                links.append(link_neighbours)
        return links

    def get_link_cons(self, atoms, layers, links):
        from ASE_utils import to_molmod
        mm_atoms = to_molmod(atoms)
        neighbours = lambda a:mm_atoms.graph.neighbors[a]

        link_cons = []
        for i in range(len(layers)):
            layer_links = []
            layer_atoms = layers[i]
            link_atoms = links[i]

            for link_atom in link_atoms:
                #get neighbours of link_atom that are part of the current layer
                layer_neighbours = [n for n in neighbours(link_atom) if n in layers[i]]
                #check only link_atom only connects to one atom in the current layer
                if len(layer_neighbours) != 1:
                    raise RuntimeError('link atom has invalid number of neighbours')
                else:
                    layer_links.append(layer_neighbours[0])

            link_cons.append(layer_links)
        return link_cons

    def _get_basis(self):
        """returns basis section of a gaussian input file"""
        if not self.route_str_params['basis'].lower() == 'gen' and not self.route_self_params['extrabasis'].lower() == 'extrabasis':
            return ''

        if not self.basisfile:
            raise RuntimeError('Please set basisfile.')
        elif (not os.path.isfile(self.basisfile)):
            raise RuntimeError('Basis file %s does not exist.'\
                               % self.basisfile)
        else:
            with open(self.basisfile, 'r') as f2:
                basis = f2.read().strip()

        return '\n' + basis + '\n\n'

    def _get_pbc(self, atoms):
        """returns periodic boundary conditions section of a gaussian input file"""
        if not atoms.get_pbc().any():
            return ''

        cell = atoms.get_cell()
        line = str()
        for v in cell:
            line += 'TV %20.10f%20.10f%20.10f\n' % (v[0], v[1], v[2])
        return line

    def _get_mod_redundant(self):
        """returns redundant coordinates section of a gaussian input file"""
        if not self.redundant_coord_params:
            return ''

        coord_param_lookup = {'cartesian':'X', 'length':'B', 'angle':'A', 'dihedral':'D', 'bend':'L'}
        action_param_lookup = {'scan':'S', 'freeze':'F', 'activate':'A', 'add':'B', 'kill':'K', 'remove':'R', 'derivatives':'D', 'diag':'H'}

        mod_redundant_strs = []
        for redundant_coord_params in self.redundant_coord_params:
            if not redundant_coord_params['action_type']:
                raise(RuntimeError('No action_type selected'))

            full_coord_type = redundant_coord_params['coord_type']
            full_action_type = redundant_coord_params['action_type']

            coord_type = coord_param_lookup[full_coord_type]
            action_type = action_param_lookup[full_action_type]
            atom_no_str = ' '.join(map(str, redundant_coord_params['atom_nos']))
            m1 = redundant_coord_params.get('min','')
            m2 = redundant_coord_params.get('max','')

            if full_action_type == 'scan':
                p1, p2 = redundant_coord_params['frames'], redundant_coord_params['degrees']
            elif full_action_type == 'diag':
                p1, p2 = redundant_coord_params['diag'], ''
            else:
                p1, p2 = '', ''

            mod_redundant_str =  "{ct} {atms} {at} {param1} {param2} {min} {max}".format(ct=coord_type,
                                  atms=atom_no_str,at=action_type,param1=p1,param2=p2,min=m1,max=m2)
            mod_redundant_strs.append(mod_redundant_str)

        return '\n'.join(mod_redundant_strs) + '\n'

    def _get_orbital_specification(self):
        """currently returns orbitals to be swapped for use with guess=alter"""
        guess_string = self.route_str_params.get('guess', "").lower()
        lines = []

        if not 'alter' in guess_string or not any(
                self.extra_list_params['original_alpha_orbitals'] + self.extra_list_params['swapped_alpha_orbitals'] +
                self.extra_list_params['original_beta_orbitals'] +  self.extra_list_params['swapped_beta_orbitals']):
            return ''

        if not len(self.extra_list_params['original_alpha_orbitals']) == len(self.extra_list_params['swapped_alpha_orbitals']) and \
               len(self.extra_list_params['original_beta_orbitals']) == len(self.extra_list_params['swapped_beta_orbitals']):
            raise(RuntimeError('Invalid specification of orbitals to swap'))

        for o,s in zip(self.extra_list_params['original_alpha_orbitals'], self.extra_list_params['swapped_alpha_orbitals']):
            lines.append('{orig} {swapped}'.format(orig=o, swapped=s))

        lines.append('')

        for o,s in zip(self.extra_list_params['original_beta_orbitals'], self.extra_list_params['swapped_beta_orbitals']):
            lines.append('{orig} {swapped}'.format(orig=o, swapped=s))

        return '\n'.join(lines) + '\n'

    def _get_cas_specification(self):
        method_str = self.route_str_params.get('method','').lower()
        lines = []

        if not 'cas' in method_str or not 'stateaverage' in method_str:
            return ''

        weights = map(str, self.extra_list_params['state_weights'])
        state_weight_line = " ".join(weights) + '\n'
        return state_weight_line

    def _get_extra(self):
        """returns additional input"""

        extra = []

        opt_string = self.route_self_params.get('opt', '').lower()
        extra_inp_str =  self.extra_params.get('extra_input','')

        if 'qst3' in opt_string or 'qst2' in opt_string:
            extra.append(_get_qst_extra)
        if extra_inp_str:
            extra.append('\n' + extra_inp_str)

        return '\n'.join(extra)

    def _get_qst_extra(self):
        """returns additional input for qst2/qst3 frames"""

        opt_string = self.route_self_params.get('opt', "").lower()

        configs = self.extra_dict_params.get('component_configs',{})
        if 'qst3' in opt_string:
            reactants = configs.get('reactants')
            guess_ts = configs.get('ts')
            final_config = 'Reactant configuration\n\n' + self._get_mol_details(reactants) + \
                           'Initial transition state configuration\n\n' + self._get_mol_details(guess_ts) + self._get_mod_redundant() + '\n'
        elif 'qst2' in opt_string:
            reactants = configs.get('reactants')
            final_config = 'Reactant configuration\n\n' + self._get_mol_details(reactants) + self._get_mod_redundant() + '\n'

        else:
            Raise(RuntimeError('Unknown extra qst input specified'))

        return final_config

    #order that gaussian requires various commands here: http://www.gaussian.com/g_tech/g_ur/m_input.htm
    def _get_input(self,atoms):
        """generates input string"""

        if self.extra_params['raw_input']:
            return self.extra_params['raw_input']

        link0 = self._get_link0()
        route = self._get_route()
        mol_details = self._get_mol_details(atoms)
        basis = self._get_basis()
        pbc = self._get_pbc(atoms)
        mod_redundant = self._get_mod_redundant()
        orbital_spec = self._get_orbital_specification()
        cas_spec = self._get_cas_specification()
        extra = self._get_extra()

        if not 'allcheck' in self.route_self_params['geom'].lower():
            title = ' \n\nGaussian input prepared by ASE\n\n'
        else:
            title = '\n\n'

        return link0 + route + title + mol_details + basis + pbc + mod_redundant + orbital_spec + extra + cas_spec + '\n'

    def write_input(self, filename, atoms):
        """Writes the input file"""
        inputfile = open(filename, 'w')
        input_str = self._get_input(atoms)
        inputfile.write(input_str)
        inputfile.close()

        # sends input file to server
        if not self.job_params['version'] == 'direct_g09':
            self.send_to_home(filename)

    def _get_scratch_dir(self):
        scratch = os.environ['GAUSS_SCRATCH']
        local_home = os.environ['ASE_HOME']

        try:
            active_dir = os.getcwd().split(local_home)[1]
        except IndexError:
            raise RuntimeError('Not running from within ASE_HOME')

        scratch_dir = scratch + active_dir
        return scratch_dir

    def _get_home_dir(self):
        home = os.environ['GAUSS_HOME']
        local_home = os.environ['ASE_HOME']

        try:
            active_dir = os.getcwd().split(local_home)[1]
        except IndexError:
            raise RuntimeError('Not running from within ASE_HOME')

        home_dir = home + active_dir
        return home_dir


    def send_to_home(self, filename):
        """copies the input from the working dir to the server's home directory"""
        serv = os.environ['GAUSS_HOST']
        home_dir = self._get_home_dir()

        exitcode= os.system('scp "%s" "%s:%s"' % (filename, serv,home_dir) )

        if (exitcode != 0):
            raise RuntimeError('Unable to send file {f} to server, scp exited with {s}'.format(f=filename, s=exitcode))

    def get_from_scratch(self, filename=None ,frc=False, log=0):
        """if output file not present copies the output from the server's scratch dir to the working dir"""
        #from gaussian_job_manager import server_files_equal

        if not filename:
            filename= self.log

        if not frc and self.calc_complete:
            return

        serv = os.environ['GAUSS_HOST']
        scratch_dir = self._get_scratch_dir()

        #exitcode=0
        #if frc or not server_files_equal(scratch_dir + '/' + filename, filename):
        exitcode = os.system('scp "%s:%s/%s" "."' % (serv,scratch_dir,filename) )

        if log and exitcode != 0:
            raise RuntimeError('Unable to get file {f} from server, scp exited with {s}'.format(f=filename, s=exitcode))

    def gen_fchk(self, frc=False):
        """generates fchk file from chk point file for the molecule specified assumes chk point file exists in the scratch directory"""

        ssh, sftp = pbs.connect_server(ssh=True,sftp=True)
        try:
            sftp.stat(serv_file)
            fchk_exists = True
        except IOError:
            fchk_exists = False
        finally:
            sftp.close()

        if fchk_exists and not frc:
            warnings.warn('.fchk file already generated, aborting')
        elif fchk_exists and frc:
            warnings.warn('.fchk file already generated, overwriting')

        if not fchk_exists or frc:
            i,o,e = ssh.exec_command('/home/gaussian-devel/gaussiandvh13_pgi_118/gdv/formchk {fn}'.format(fn= self._get_scratch_dir() + '/' + c.calc.label+'.chk'))
            formchk_error = e.readlines()
            ssh.close()
            return not bool(formchk_error)


        #this currently reads through the file 24 (3 per call to .read_output()) - which is why it's slow!
#    def read(self, atoms):
#        self.positions = self.read_atoms().get_positions()
#        self.energy_free, self.energy_zero = self.read_energy()
#        self.forces = self.read_forces(atoms)
#        self.dipole = self.read_dipole()
#        self.fermi = self.read_fermi()
#        self.atoms = self.read_atoms().copy()
#        try:
#            self.nbands = self.read_nbands()
#        except NotImplementedError:
#            do_nothing = True
#        except AttributeError:
#            do_nothing = True
#        try:
#            self.stress = self.read_stress()
#        except NotImplementedError:
#            do_nothing = True
#        return

    #reduced to reading through the file 3 times
    def read(self, atoms=None):
        quantities = self.read_multi_output(self.log, quantities=['atoms','energy', 'forces', 'dipole'])

        self.positions = quantities[0].get_positions()
        self.energy_free, self.energy_zero = quantities[1], quantities[1]
        self.forces = quantities[2]
        self.dipole = quantities[3]
        #not implemented
        self.fermi = 0.0
        self.atoms = quantities[0].copy()

        try:
            self.nbands = self.read_nbands()
        except NotImplementedError:
            do_nothing = True
        except AttributeError:
            do_nothing = True
        try:
            self.stress = self.read_stress()
        except NotImplementedError:
            do_nothing = True
        return

    @property
    def data(self):
        if not self._data:
            self.set_data()
        return self._data

    def set_data(self):
        try:
            self._data = GR(self.log)[0]
        except IndexError:
            pass

    @property
    def max_data(self):
        if not self._max_data:
            self.set_max_data()
        return self._max_data

    def set_max_data(self):
        g = GR(self.log,inc_cclib=True)
        try:
            gr_data = g.data[0]
        except IndexError:
            gr_data = {}

        for k in gr_data:
            gr_data[k.title()] = gr_data.pop(k)

        for k in g.cclib_data:
            if not k.islower():
                g.cclib_data[k.lower()] = g.cclib_data.pop(k)

        self._max_data.update(gr_data)
        self._max_data.update(g.cclib_data)


    #kind of annoying that _fchk_data is an object whereas _max_data and _data are dictionaries
    @property
    def fchk_data(self):
        if not self._fchk_data:
            self.set_fchk_data()
        return self._fchk_data

    def set_fchk_data(self):
        from fchk_utils import FCHK
        try:
            self._fchk_data = FCHK(self.fchk)
        except IOError:
            self._fchk_data = None

    @property
    def fingerprint(self):
        if not self._fingerprint:
            self.set_fingerprint()
        return self._fingerprint

    def set_fingerprint(self):
        try:
            with open(self.log) as log_f:
                log_data = log_f.readlines()
        except IOError:
            log_data = []

        try:
            self._fingerprint = log_data[132] + "".join(log_data[-10:])
        except IndexError:
            self._fingerprint = "".join(log_data[-10:])


    @property
    def status(self):
        """extracts the status of the calculation from the log file"""
        if not self._status or self._status == 'Incomplete':
            self.set_status()

        return self._status

    def set_status(self):
        command = "tail -n4 '{fl}'".format(fl=self.log)
        p = Popen(command, stdout=PIPE, stderr=PIPE, shell=True)
        stdout, stderr = p.communicate()

        last_lines = stdout.split('\n')

        if any('Normal termination' in l for l in last_lines):
            self._status =  'Success'
            return
        try:
            i = next(l[0] for l in enumerate(last_lines) if 'Error termination' in l[1])
        except StopIteration:
            self._status = 'Incomplete'
            return

        self._status = 'Error: {e}'.format(e=last_lines[i-1].strip())


    @property
    def notes(self):
        """extracts notes regarding the calculation from the log file"""
        if not self._notes:
            self.set_notes()
        return self._notes

    def set_notes(self):
        command = "grep 'The wavefunction has a' {fl}".format(fl=self.log)
        p = Popen(command, stdout=PIPE, stderr=PIPE, shell=True)
        stdout, stderr = p.communicate()

        if stdout:
            self._notes =  stdout

    @property
    def time(self):
        if not self._time:
            self.set_time()
        return self._time

    def set_time(self):
        command = "tail -n4 '{fl}'".format(fl=self.log)
        p = Popen(command, stdout=PIPE, stderr=PIPE, shell=True)
        stdout, stderr = p.communicate()

        last_lines = stdout.split('\n')
        try:
            time_line = next(l for l in last_lines if l[1:14] == "Job cpu time:")

            days = float(line.split()[3])
            hours = float(line.split()[5])
            minutes = float(line.split()[7])
            seconds = float(line.split()[9])

            self._time =  days*24 + hours + minutes/60 + seconds/3600
        except StopIteration:
            self._time = 0

    @property
    def stable(self):
        """checks whether the calculation is stable"""
        if self.route_self_params.get('stable'):
            return True if 'instability' not in self.notes else False
        else:
            return ''

    #this will mean calls to self.data and self.max_data reparse the data, we need to run this if we fundamentally change the calculation e.g. by restarting it
    def reset_cached_data(self):
        """clears stored parsed data associated with the calculation"""
        self._data = {}
        self._max_data = {}
        self._fchk_data = None
        self._fingerprint = ''
        self._notes = ''
        self._time = 0
        self._status = ''
        self._calc_complete = None

    def read_multi_output(self, filename, quantities):
        import numpy as np
        from ase.io.gaussian import read_gaussian_out
        """reads multiple quantities from the output (avoids repeatedly reading through the same possibly large log file"""

        #if log file not present we always return nan or in the case of atoms the original unmodified atoms
        if not os.path.exists(filename) or not self.calc_complete:
            return [np.nan if q != 'atoms' else self.atoms for q in quantities]

        #gaussian hack
        if 'energy' in quantities and self.method.lower()=='cis':
            energy_ind = quantities.index('energy')
            cis_energy = gaussian_hacks.get_cis_energy(filename)
            quantities.pop(energy_ind)
            quant = read_gaussian_out(filename, quantities=quantities)
            quant.insert(energy_ind, cis_energy)
        else:
            quant = read_gaussian_out(filename, quantities=quantities)

        return quant

    def read_output(self, filename, quantity):
        import numpy as np
        from ase.io.gaussian import read_gaussian_out
        """Reads the output file using GaussianReader"""

        #if log file not present we always return nan
        if not os.path.exists(filename) or not self.calc_complete:
            if quantity != 'atoms':
                return np.nan
            else:
                return self.atoms


        #gaussian hacks
        if quantity == 'energy' and self.method.lower() == 'cis':
            return gaussian_hacks.get_cis_energy(filename)

        #return single quantity
        if (quantity == 'energy'):
            return read_gaussian_out(filename, quantity='energy')
        elif (quantity == 'forces'):
            return read_gaussian_out(filename, quantity='forces')
        elif (quantity == 'dipole'):
            return read_gaussian_out(filename, quantity='dipole')
        elif (quantity == 'version'):
            return read_gaussian_out(filename, quantity='version')
        elif (quantity == 'atoms'):
            return read_gaussian_out(filename, quantity='atoms')
        elif (quantity == 'frequencies'):
            return read_gaussian_out(filename, quantity='frequencies')

    def read_atoms(self):
        """Reads and returns the atoms specified in the log file"""
        atoms = self.read_output(self.log, 'atoms')
        return atoms

    def read_energy(self):
        """Reads and returns the energy"""
        energy = self.read_output(self.log, 'energy')
        return [energy, energy]

    def read_freqs(self, atoms):
        """Reads and returns the frequencies"""
        freqs = self.read_output(self.log, 'frequencies')

    def read_forces(self, atoms):
        """Reads and returns the forces"""
        forces = self.read_output(self.log, 'forces')
        return forces

    def read_dipole(self):
        """Reads and returns the dipole"""
        dipole = self.read_output(self.log, 'dipole')
        return dipole

    def read_fermi(self):
        """No fermi energy, so return 0.0"""
        return 0.0

    def read_stress(self):
        raise NotImplementedError

    def update(self, atoms):
        """Updates and does a check to see if a calculation is required"""
        if self.calculation_required(atoms, ['energy']):
            if (self.atoms is None or
                self.atoms.positions.shape != atoms.positions.shape):
                self.clean()

            if (self.directory):
                curdir = os.getcwd()
                if not os.path.exists(self.directory):
                    os.makedirs(self.directory)
                os.chdir(self.directory)
                self.calculate(atoms)
                os.chdir(curdir)
            else:
                self.calculate(atoms)

    def calculation_required(self, atoms, quantities):
        """Checks if a calculation is required"""
        if (self.positions is None or
           (self.atoms != atoms) or
           (self.link0_str_params != self.old_link0_str_params) or
           (self.link0_int_params != self.old_link0_int_params) or
           (self.route_str_params != self.old_route_str_params) or
           (self.route_int_params != self.old_route_int_params) or
           (self.route_float_params != self.old_route_float_params) or
           (self.route_bool_params != self.old_route_bool_params) or
           (self.route_self_params != self.old_route_self_params) or
           (self.basisfile != self.old_basisfile) or
           (self.label != self.old_label) or
           (self.ioplist != self.old_ioplist)):

            return True
        return False

    def clean(self):
        """Cleans up from a previous run"""
        extensions = ['.chk', '.com', '.log']

        for ext in extensions:
            f = self.label + ext
            try:
                if (self.directory):
                    os.remove(os.path.join(self.directory, f))
                else:
                    os.remove(f)
            except OSError:
                pass

    def get_command(self):
        """Return command string if program installed, otherwise None.  """
        command = None
        if ('GAUSS_EXEDIR' in os.environ) \
                and ('GAUSSIAN_COMMAND' in os.environ):
            command = os.environ['GAUSSIAN_COMMAND']
        return command

    def run(self, test=False):
        """Runs Gaussian"""

        local_home = os.environ['ASE_HOME']

        try:
            active_dir = os.getcwd().split(local_home)[1]
        except IndexError:
            raise RuntimeError('Not running from within ASE_HOME')

        host_dir = os.environ['GAUSS_HOME'] + active_dir + '/'
        scratch_dir = os.environ['GAUSS_SCRATCH'] + active_dir + '/'

        if not self.job_params['version']:
            self.job_params['version'] = 'direct_g09'

        if not self.job_params['queue']:
            queue ='None'
        else:
            queue = self.job_params['queue']

        if self.job_params['version'] == 'direct_g09':
            command = 'module load gaussian; g09 <{inp}> {out}'.format(inp=host_dir + self.label +'.com', out=scratch_dir + self.label +'.log')

        elif self.job_params['version'] == 'local_g09':
            command = 'qsub -l ncpus={n},memory={m}mb,time={t}:00:00 -v inp_f={inp},host_d={fld} -q {q} ~/bin/ase_calcs.py '.format(fld=host_dir, inp=self.label + '.com', q=queue, n=self.job_params['nodes'], m=self.job_params['memory'], t=int(self.job_params['time']))

        elif self.job_params['version'] =='g09':
            command = '~/bin/g09_calcs.py {fld} {inp} {p} {m} {t} {q}'.format(host= os.environ['GAUSS_HOST'], fld=host_dir, inp=self.label + '.com', q=queue, p=self.job_params['nodes'], m=self.job_params['memory'], t=int(self.job_params['time']))
            #os.system("ssh {host} '~/bin/g09_calcs.py {fld} {inp} {p} {m} {t} {q}'".format(host= os.environ['GAUSS_HOST'], fld=host_dir, inp=self.label + '.com', q=queue, p=self.job_params['nodes'], m=self.job_params['memory'], t=int(self.job_params['time'])))

        elif self.job_params['version'] == 'gdv_latest':
            #os.system("ssh {host} '~/bin/gdv_latest_calcs.py {fld} {inp} {p} {m} {t} {q}'".format(host= os.environ['GAUSS_HOST'], fld=host_dir, inp=self.label + '.com', q=queue, p=self.job_params['nodes'], m=self.job_params['memory'], t=int(self.job_params['time'])))
            command = '~/bin/gdv_latest_calcs.py {fld} {inp} {p} {m} {t} {q}'.format(fld=host_dir, inp=self.label + '.com', q=queue, p=self.job_params['nodes'], m=self.job_params['memory'], t=int(self.job_params['time']))

        elif self.job_params['version'] == 'gdvh23':
            #os.system("ssh {host} '~/bin/gdvh23_calcs.py {fld} {inp} {p} {m} {t} {q}'".format(host= os.environ['GAUSS_HOST'], fld=host_dir, inp=self.label + '.com', q=queue, p=self.job_params['nodes'], m=self.job_params['memory'], t=int(self.job_params['time'])))
            command = '~/bin/gdvh23_calcs.py {fld} {inp} {p} {m} {t} {q}'.format(fld=host_dir, inp=self.label + '.com', q=queue, p=self.job_params['nodes'], m=self.job_params['memory'], t=int(self.job_params['time']))

        elif self.job_params['version'] == 'gdvh32':
            #os.system("ssh {host} '~/bin/gdvh23_calcs.py {fld} {inp} {p} {m} {t} {q}'".format(host= os.environ['GAUSS_HOST'], fld=host_dir, inp=self.label + '.com', q=queue, p=self.job_params['nodes'], m=self.job_params['memory'], t=int(self.job_params['time'])))
            command = '~/bin/gdvh32_calcs.py {fld} {inp} {p} {m} {t} {q}'.format(fld=host_dir, inp=self.label + '.com', q=queue, p=self.job_params['nodes'], m=self.job_params['memory'], t=int(self.job_params['time']))


        else:
            raise RuntimeError("Invalid Gaussian version selected")


        if test:
            return command

        elif 'direct' in self.job_params['version'] or 'local' in self.job_params['version']:
            os.system(command)
        else:
            ssh = pbs.connect_server(ssh=True)
            stdin, stdout, stderr = ssh.exec_command(command)
            self.job_params['qid'] = stdout.read().split('.')[0]
            ssh.close()

    def calculate(self, atoms, force=False):
        """initializes calculation and runs Gaussian, force=True will perform a gaussian calculation even if a previous log file already exists"""
        self.initialize(atoms)

        if os.path.isfile(self.log) and self.calc_complete and not force:
           warnings.warn(RuntimeWarning('Prior calculation already performed, not recomputing'))

        elif os.path.isfile(self.log) and not self.calc_complete:
            warnings.warn(RuntimeWarning('Calculation already performed but incomplete, overwriting'))

        elif os.path.isfile(self.log) and self.calc_complete and force:
            warnings.warn(RuntimeWarning('Calculation already performed and complete, forced overwriting'))

        if not self.calc_complete or force:
            self.write_input(self.label + '.com', atoms)
            self.run()

        self.converged = self.read_convergence()
        self.set_results(atoms)

    def get_best_step(self):
        try:
            min_e = min(self.max_data['scfenergies'])
        except StopIteration:
            return -1

        return self.max_data['scfenergies'].index(min_e)

    def get_number_steps(self):
        return len(self.max_data['scfenergies'])

    def get_oniom_components(self, **kwargs):
        return oniom_comp_calcs(self.atoms, **kwargs)

    #remember that calling this changes the label - even if you don't run the calculation you can't call this multiple times with no consequence
    def restart(self, label = "", add_label="", remove_label="", start=False, from_calc=None, no_old_chk= False, auto_opt=False, ioplist=False, **kwargs):
        """restarts a calculation overwrites guess and geom unless specified, otherwise keeps calculation details,

        parameters are:

        label (string) provides a new label
        add_label (string) adds a string to the current label in addition to adding _restart_{n} to it
        remove_label (string) removes a string from the current label in addition to adding _restart_{n} to it
        start (boolean) specifies whether the function call also starts the calculation
        from_calc (ASE_obj with gaussian calc) specifies the calculation we are restarting from
        no_old_chk (boolean) specifies not copying the chk from the previous calculation
        auto_opt (boolean) specifies whether when running an optimisation to intelligently seek out the best step in the previous optimisation to restart from
        ioplist (list) specifies IOPs to add to the restarting calculcation

        additional keywords relating to gaussian can be specified as if we were setting a gaussian calculation object via the Gaussian()

        Note calling this function changes the label of the current calculation object even if you don't run the calculation, thus we you cannot call this multiple times with no consequence"""

        self.set_results(self.atoms)
        self.reset_cached_data()

        if from_calc:
            try:
                restart_chk = from_calc.calc.label + '.chk'
            except AttributeError:
                raise RuntimeError("Invalid restart object require molecule object with attached Gaussian calculator")
        else:
            restart_chk = self.label + '.chk'

        local_home = os.environ['ASE_HOME']
        try:
            active_dir = os.getcwd().split(local_home)[1]
        except IndexError:
            raise RuntimeError('Not running from within ASE_HOME')
        host_dir = os.environ['GAUSS_SCRATCH'] + active_dir + '/'

        restart_chk = host_dir + restart_chk

        #gaussian doesn't let us get the geometry from an oniom checkpoint file whilst reading in the scf guess for the individual fragments from their individual check point files
        #so unless we are specifically not reading from the chk file we set the calculators atoms used to write the geometry to the positions specified by the current calculation (i.e. we achieve the same geometry as if we were reading from the chk file)
        if 'component_calcs' in kwargs and kwargs['component_calcs']:
            if not no_old_chk:
                self.atoms = ase_read(self.log)
            kwargs.update({'guess':'input', 'geom':'', 'oldchk':None})

        elif no_old_chk:
            kwargs.update({'oldchk':None})

        else:
            kwargs.update({'oldchk':restart_chk})

            #by default we read the wavefunction and geometry from the chk point file
            if 'guess' not in kwargs:
                kwargs.update({'guess':'read'})
            if 'geom' not in kwargs and 'layer_mults' not in kwargs and 'multiplicity' not in kwargs:
                kwargs.update({'geom':'allcheck'})
            elif 'geom' not in kwargs:
                kwargs.update({'geom':'check'})

        #if we are restarting a geometry minimum optimisation we restart from the previous lowest energy step
        opt_section1, opt_section2 = kwargs.get('opt', ""), self.route_self_params.get('opt', "")
        is_opt_calc = (opt_section1 and 'ts' not in opt_section1) or (opt_section2 and 'ts' not in opt_section2)
        reading_geom = 'geom' in kwargs and kwargs['geom'] and 'check' in kwargs['geom'].lower()
        step_given = 'geom' in kwargs and kwargs['geom'] and 'step' in kwargs['geom'].lower()

        if reading_geom and is_opt_calc and not step_given and auto_opt:
            best_step = self.get_best_step()
            #if we hit an error trying to extract the best_step or are continuing from the last step we don't need to specify the step (if we try gaussian complains)
            if best_step != -1 and best_step != self.get_number_steps() -1:
                kwargs['geom'] = ', '.join(kwargs['geom'].split(',') + ['step={s}'.format(s=best_step)])
                kwargs['guess'] = None

        self.set(**kwargs)

        if ioplist:
            self.old_ioplist = self.ioplist[:]
            self.ioplist = list(ioplist)

        try:
            file_n = "_restart_".join(self.label.split('_restart_')[0:-1])
            restart_no = int(self.label.split('_restart_')[-1])
        except ValueError:
            file_n = self.label
            restart_no= 0

        #would be better to remove single instance from the right incase the string occurs more than once
        if remove_label:
            if file_n.count(remove_label) > 1:
                raise RuntimeError('Filename contains string to be removed more than once do you really mean to do that?')
            file_n = file_n.replace('_' + remove_label, '')

        label_contents = [file_n, add_label, 'restart', str(restart_no+1)]

        if not label:
            self.label = '_'.join([e for e in label_contents if e])
        else:
            self.label = label

        if start:
            self.start()
        else:
            # starting the calculation automatically reinitializes but if we do not start immediately then in order for other
            # calculation details to be updated (e.g. the link0 information which specifies the .chk file) we reinitialize
            self.initialize(self.atoms)

        #must reset parsed data, otherwise when we call to self.data or self.max_data we will actually be reading from the previous calculation
        self.reset_cached_data()
        return self

    #not really used yet
    #def additional_calc(self, gauss_calc):
    #    if gauss_calc.route_self_params.get('geom') or gauss_calc.route_str_params.get('guess'):
    #        raise RuntimeError('Cannot define geometry or guess wavefunction in additional calculations as they are set by the preceding calculation')
    #    else:
    #        gauss_calc.route_self_params['geom'] = 'allcheck'
    #        gauss_calc.route_str_params['guess'] = 'read'
    #        self.additional_calcs.append(gauss_calc)
    #    return

    def start(self, frc=False):
        self.calculate(self.atoms, force=frc)

    def oniom_stable_start(self, log=0,frc=False):
        import oniom_utils
        temp_mol_obj = self.atoms
        temp_mol_obj.calc = self
        oniom_utils.oniom_stable(temp_mol_obj, log=log, frc=frc, proc=self.job_params['nodes'], mem=self.job_params['memory'])

    @property
    def calc_complete(self):
        #special case as _calc_complete can be True or False if evaluated but is None if it has not been evaluated
        if self._calc_complete is None:
            self.set_calc_complete()
        return self._calc_complete

    def set_calc_complete(self):
        if not os.path.isfile(self.log):
            self._calc_complete = False
            return

        test = 'termination of Gaussian'

        f = open(self.log, 'r')
        lines = f.readlines()
        f.close()

        if not lines:
            self._calc_complete =  False
            return

        for line in lines:
            if (line.rfind(test) > -1):
                complete = True
            else:
                complete = False

        self._calc_complete = complete

    #todo
    def read_convergence(self):
        """Determines if calculations converged"""
        return False
        converged = False

        gauss_dir = os.environ['GAUSS_EXEDIR']
        test = '(Enter ' + gauss_dir + '/l9999.exe)'

        #if calc not complete return false
        if not self.calc_complete:
            return False

        f = open(self.log, 'r')
        lines = f.readlines()
        f.close()

        for line in lines:
            if (line.rfind(test) > -1):
                converged = True
            else:
                converged = False

        return converged

    def set_results(self, atoms):
        """Sets results"""

        #allows us to run the calculator directly without an ASE atoms object
        if not atoms:
            atoms = Atoms('H', [(0,0,0)])

        self.old_positions = atoms.get_positions().copy()
        mod_atoms =  self.read_atoms()
        if mod_atoms:
            self.atoms = mod_atoms.copy()
        else:
            self.atoms = atoms.copy()
        self.old_link0_str_params = self.link0_str_params.copy()
        self.old_link0_int_params = self.link0_int_params.copy()
        self.old_route_str_params = self.route_str_params.copy()
        self.old_route_int_params = self.route_int_params.copy()
        self.old_route_float_params = self.route_float_params.copy()
        self.old_route_bool_params = self.route_bool_params.copy()
        self.old_route_self_params = self.route_self_params.copy()
        self.old_basisfile = self.basisfile
        self.old_label = self.label
        self.old_ioplist = self.ioplist[:]

    def get_version(self):
        return self.read_output(self.log, 'version')

    def get_name(self):
        return self.__class__.__name__
