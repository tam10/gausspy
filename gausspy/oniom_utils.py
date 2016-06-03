__author__ = 'clyde'

from warnings import warn
from gaussian_hacks import get_methods_wo_stability_opt
import itertools
import copy


#def extract_oniom_components():
#    self.restart(method = self.method + '=OnlyInputFiles')


def extract_oniom_inputs(atoms):
    with open(atoms.calc.log, 'r') as log_f:
        contents = log_f.readlines()

    if not any('=OnlyInputFiles' in l for l in contents):
        raise RuntimeError('Can only extract input from ONIOM calculations run with =OnlyInputFiles')

    input_starts = [i for i,e in enumerate(contents) if '#P' in e or '#p' in e][1:-1]
    input_ends = [i for i,e in enumerate(contents) if e.count('-') > 30 and i > min(input_starts) and i+1 not in input_starts]

    if len(input_starts) != len(input_ends):
        raise RuntimeError('Failed to extract input')

    inp_strs = []
    for i in range(len(input_starts)):
        inp_strs.append(''.join(contents[input_starts[i]:input_ends[i]]))

    return inp_strs


def get_route(inp_str):
    return next(l for l in inp_str.split('\n\n') if '#p' in l or '#P' in l)


def add_to_input(inp_str, **kwargs):
    from gaussian import Gaussian
    route = get_route(inp_str)
    command = Gaussian(**kwargs)._get_route().replace('#p', '').replace('#P', '')
    n_route = route + command
    return inp_str.replace(route, n_route)


def switch_input(inp_str, method_lst=None, basis_lst=None):
    """takes an ase atoms object with a gaussian calculator that has raw_input and generates a new raw_input string
    with the method and basis changed to values provided to the function"""
    if not method_lst:
        method_lst = []
    elif not method_lst[0]:
        raise(RuntimeError('I refuse to replace Null'))
    if not basis_lst:
        basis_lst =[]
    elif not basis_lst[0]:
        raise(RuntimeError('I refuse to replace Null'))

    old_route = get_route(inp_str)
    new_route = copy.deepcopy(old_route)

    if method_lst:
        old_method, new_method = method_lst
        new_route = new_route.replace(old_method, new_method)
    if basis_lst:
        old_basis, new_basis = basis_lst
        new_route = new_route.replace(old_basis, new_basis)

    return inp_str.replace(old_route, new_route)


def oniom_comp_calcs(atoms_oniom, **kwargs):
    """extracts component calculations from an atoms_oniom object that has run an =OnlyInputFiles oniom calculation"""
    from gaussian import Gaussian
    # recursion means if we pass in [atoms1s, atoms2, atoms3...] we return [low_reals, high_models, low_models]
    # where low_reals = [low_real1, low_real2, low_real3...] etc.
    if isinstance(atoms_oniom, list):
        return [list(c) for c in zip(*[oniom_comp_calcs(r, **kwargs) for r in atoms_oniom])]

    proc,mem,ver = atoms_oniom.calc.job_params['nodes'], atoms_oniom.calc.job_params['memory'], atoms_oniom.calc.job_params['version']

    init_inp_strs = extract_oniom_inputs(atoms_oniom)
    inp_strs = [add_to_input(s, **kwargs) for s in init_inp_strs]
    old_label = atoms_oniom.calc.label

    components = []
    comp_strs = ['low_real', 'high_model', 'low_model']
    method_strs = get_oniom_calc_methods(atoms_oniom.calc.method.replace('=OnlyInputFiles',''))

    for i in range(len(inp_strs)):
        atoms_comp = copy.deepcopy(atoms_oniom)
        new_label = old_label + '_' + comp_strs[i]

        try:
            method, basis = method_strs[i].split('/')
        except ValueError:
            method, basis = method_strs[i], ''

        #the method/basis variables are not used as we are defining raw_input but it is useful to be able to read them later on
        atoms_comp.set_calculator(Gaussian(label=new_label, raw_input=inp_strs[i], method=method, basis=basis))
        atoms_comp.calc.set_job(nodes=proc, memory=mem, version=ver)

        #the input files gaussian produces have no link information so we have to add that manually
        atoms_comp.calc.initialize(atoms_comp)
        atoms_comp.calc.extra_params.update({'initial_raw_input': inp_strs[i]})
        atoms_comp.calc.extra_params['raw_input'] = atoms_comp.calc._get_link0() + inp_strs[i]

        components.append(atoms_comp)
    return components


#not currently used
def oniom_comp_calcs_v2(atoms, **kwargs):
    import tempfile
    from gaussian_input_parser import extract_calc
    init_inp_strs = extract_oniom_inputs(atoms)
    inp_strs = [add_to_input(s, **kwargs) for s in init_inp_strs]

    temps = []
    for i, inp_str in enumerate(inp_strs):
        temps.append(tempfile.NamedTemporaryFile())
        temps[i].write(inp_str)

    components = [extract_calc(temps[i].name) for i in range(len(init_inp_strs))]

    for temp_fl in temps:
        temp_fl.close()

    return components

def oniom_restart(oniom_atoms, component_calcs):
    low_real,high_model,low_model = component_calcs
    oniom_atoms.restart(guess='Input')
    return


def get_oniom_comp_calcs(oniom_obj, frc=False, **kwargs):
    try:
        calc = oniom_obj.calc
    except AttributeError:
        calc = oniom_obj

    orig_params = calc.job_params['nodes'], calc.job_params['memory'], calc.job_params['time']

    calc.route_str_params['method'] += '=OnlyInputFiles'
    calc.job_params['nodes'], calc.job_params['memory'], calc.job_params['time'] = 1, 100, 1
    calc.start(frc=frc)

    #set the job_parameters back to their original values
    calc.job_params['nodes'], calc.job_params['memory'], calc.job_params['time'] = orig_params
    return oniom_comp_calcs(oniom_obj, **kwargs)


def get_oniom_calc_methods(oniom_str):
    """gets the set of calculation methods used in the individual components of an oniom calculation as specified
    by the gaussian oniom method str, list returned is for the calculations ordered first by decreasing fragment size
    then by decreasing method accuracy so for a two level calculation we will get [L, H, L] for a three level we will get
    [L, M, L, H, M] etc. """
    component_methods = oniom_str.split('oniom(')[1][:-1].split(':')
    no_layers = len(component_methods)
    calc_methods = []
    for i in range(no_layers):
        if i == no_layers-1:
            calc_methods.append(component_methods[-1])
        else:
            calc_methods.append(component_methods[i+1])
            calc_methods.append(component_methods[i])
    return list(reversed(calc_methods))


def oniom_stable(oniom_obj, component_objs=None, log=2, frc=False):
    """performs a stability calculation for an oniom calculation by performing stability calculations on the relevant
    components and then using the subsequent wavefunctions to construct the final stable oniom calculation"""

    if not component_objs:
        component_objs = []

    #can pass in either ASE_molecule objects with a gaussian calculator attached via molecule.calc
    #or can pass in the gaussian calculator itself
    try:
        calc = oniom_obj.calc
    except AttributeError:
        calc = oniom_obj

    #save the no_cpus, memory and label for use later
    proc, mem, ver, orig_label = calc.job_params['nodes'], calc.job_params['memory'], calc.job_params['version'], calc.label

    calc.label += '_init'

    #get_comp_input returns a list of [low_real, high_model, low_model] component calculation objects
    components = get_oniom_comp_calcs(oniom_obj, frc=frc, symmetry='None', stable='opt')

    #get the different methods oniom will use for the various component calculations
    component_methods = get_oniom_calc_methods(calc.method)

    if len(component_methods) != len(components):
        raise RuntimeError('Discrepancy in number of Oniom layers')

    #boolean mask to ignore methods with no stability optimisation implemented in Gaussian
    #force_fields = ['UFF', 'DREIDING', 'AMBER', 'uff', 'dreiding', 'amber']

    force_fields = get_methods_wo_stability_opt()
    force_fields += [f.upper() for f in force_fields]
    mask = []

    for i in range(len(components)):
        if not any(force_field in component_methods[i] for force_field in force_fields):
            mask.append(True)
        else:
            mask.append(False)

    #list of non forcefield component calculations
    filtered_components = list(itertools.compress(components, mask))

    for c in filtered_components:
        #we avoid unnecessary computation by not computing components that are being passed in to the function
        if not any(component.calc.label == c.calc.label for component in component_objs):
            c.calc.start(frc)

    #if stable=opt fails and we are using an ab-initio/dft method, we try using STO-3G if that fails we try guess=mix
    for m in filtered_components:
        orig_inp_str = m.calc.extra_params['initial_raw_input']
        orig_basis = m.calc.basis
        inp_str = copy.deepcopy(orig_inp_str)

        if 'No lower point found' in m.calc.status and orig_basis and orig_basis.upper() != 'STO-3G':
            inp_str = switch_input(inp_str, basis_lst=[orig_basis.upper(), 'STO-3G'])
            m.calc.restart(no_old_chk=True)
            m.calc.extra_params['raw_input'] = m.calc._get_link0() + inp_str
            m.calc.basis = 'STO-3G'
            m.calc.start(frc)

        if 'No lower point found' in m.calc.status:
            inp_str = add_to_input(inp_str, guess='mix')
            m.calc.restart(no_old_chk=True)
            m.calc.extra_params['raw_input'] = m.calc._get_link0() + inp_str
            m.calc.start(frc)

        if m.calc.status == 'Success' and m.calc.basis.upper() == 'STO-3G' and orig_basis.upper() != 'STO-3G':
            restart_inp_str = add_to_input(orig_inp_str, guess='read')
            m.calc.restart(no_old_chk=False)
            m.calc.extra_params['raw_input'] = m.calc._get_link0() + restart_inp_str
            m.calc.basis = orig_basis
            m.calc.start(frc)

    if not all([m.calc.status == 'Success' for m in filtered_components]):
        if log == 2:
            raise RuntimeError('Not all component calculations succeeded!')
        elif log == 1:
            warn(UserWarning('Not all component calculations succeeded!'))

    #list of of components calculations where failed calcs/forcefield calcs are replaced by None entries
    final_components = [components[i] if mask[i] and components[i].calc.status == 'Success' else None for i in range(len(components))]
    method = calc.method.replace('=OnlyInputFiles', '')
    calc.restart(label=orig_label, no_old_chk=True, method=method, component_calcs=final_components)
    calc.set_job(nodes=proc, memory=mem, version=ver)
    calc.start(frc)

    return oniom_obj

def oniom_parser():
    return
    
