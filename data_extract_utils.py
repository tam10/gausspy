__author__ = 'clyde'

import copy
import gzip
import os
import warnings
import re
from pbs_util import pbs


def get_active_dirs():
    scratch, home, local_home = os.environ['GAUSS_SCRATCH'], os.environ['GAUSS_HOME'], os.environ['ASE_HOME']

    try:
        active_dir = os.getcwd().split(local_home)[1]
        home_dir = home + active_dir
        scratch_dir = scratch + active_dir
    except IndexError:
        raise RuntimeError('Not running from within ASE_HOME')

    return home_dir, scratch_dir


#will only work if called on the original calculation, calling on an intermediate restart
#e.g. orig_calc_name_restart_2.log will not return orig_calc_name_restart_5.log even if it is the final calculation
def get_latest_restart_name(file_n, ssh=None):
    """get's the file name of the final restarted calculation from the server"""
    ssh_created = False
    if not ssh:
        ssh = pbs.connect_server(ssh=True)
        ssh_created = True

    i,o,e = ssh.exec_command('ls {fn}_restart_{{?,??}}.log'.format(fn=file_n.replace('.log', '')))
    ls_output = o.readlines()
    clean_ls = [fn.strip() for fn in ls_output]
    clean_ls.sort(key=lambda fn: int(re.search('[0-9]+.log', fn).group().replace('.log' , '')))

    if ssh_created:
        ssh.close()

    try:
        return clean_ls[-1]
    except IndexError:
        return file_n


def load_data(obj_name):
    import cPickle
    """loads gzip pickled ase calculations"""
    with gzip.GzipFile(obj_name + '.pkl.gz', "rb") as data_file:
        list_ase_objs = cPickle.load(data_file)
    return list_ase_objs

def latest_restarts(list_mols):
    home_dir, scratch_dir = get_active_dirs()

    home_files = [mol.calc.label + '.log' for mol in list_mols]
    serv_files = [scratch_dir + '/' + fn for fn in home_files]

    ssh = pbs.connect_server(ssh=True)
    serv_files = [get_latest_restart_name(file_n, ssh) for file_n in serv_files]
    home_files = [sfn.replace(scratch_dir + '/', '') for sfn in serv_files]
    ssh.close()

    new_mols = copy.deepcopy(list_mols)
    for i, mol in enumerate(new_mols):
        mol.calc.label = home_files[i].replace('.log','')

    return new_mols

def import_moldata(data_file):
    """imports saved parseable data from a gziped pickled data_file"""
    try:
        mols_from_file = load_data(data_file)
    except IOError:
        mols_from_file = []
    return mols_from_file


def export_moldata(data_file, list_ase_objs):
    """assumes we have already extracted all parseable data from the calculations and gzip pickles the result"""
    import cPickle

    with gzip.GzipFile(data_file + '.pkl.gz', "wb") as data_file:
        cPickle.dump(list_ase_objs, data_file, -1)

def load_from_server(mol, depth='medium'):
    try:
        mol.calc.get_from_scratch(mol.calc.label + '.log', frc=True, log=True)
    except RuntimeError:
        warnings.warn(RuntimeWarning('Failed to extract {m} from the server'.format(m=mol.calc.label)))
        return mol

    try:
        mol.calc.read(mol)
    except (AttributeError, KeyError):
        pass

    try:
        mol.calc.set_notes()
        mol.calc.set_status()
    except IOError:
        pass

    try:
        mol.calc.set_data()
    except IndexError:
        warnings.warn(RuntimeWarning('Failed to parse {m} with ASE'.format(m=mol.calc.label)))

    mol.calc.set_fingerprint()

    if depth == 'medium':
        try:
            mol.calc.set_max_data()
        except IndexError:
            warnings.warn(RuntimeWarning('Failed to parse {m}.log with cclib'.format(m=mol.calc.label)))

    elif depth == 'heavy':
        mol.calc.get_from_scratch(mol.calc.label + '.fchk', frc=True)
        try:
            mol.calc.set_max_data()
        except IndexError:
            warnings.warn(RuntimeWarning('Failed to parse {m}.log with cclib'.format(m=mol.calc.label)))

        try:
            mol.calc.set_fchk_data()
        except IndexError:
            warnings.warn(RuntimeWarning('Failed to parse {m}.fchk with molmod'.format(m=mol.calc.label)))

    return mol

#todo
#def load_mult_from_server(list_mols):
#    serv = os.environ['GAUSS_HOST']
#    files_names = [mol.calc._get_scratch_dir() + '/' + mol.calc.log for mol in list_mols]

#    exitcode=0
#    #if frc or not server_files_equal(scratch_dir + '/' + filename, filename):
#    #exitcode = os.system('scp "%s:%s/%s" "."' % (serv, scratch_dir, filename) )

#    if exitcode != 0:
#        raise RuntimeError('Unable to get file {f} from server, scp exited with {s}'.format(f=filename, s=exitcode))
#    return

def oniom_components_on_server(mol):
    from gaussian_job_manager import server_file_exists
    home_dir, scratch_dir = get_active_dirs()
    return server_file_exists(scratch_dir + '/' + mol.calc.label + '_init.log')

def unwind(lst_of_lsts):
    return [e for l in lst_of_lsts for e in l]

def clean_local_files(list_mols):
    temp_components = []

    for m in list_mols:
        try:
            coord_calc = copy.deepcopy(m)
            coord_calc.calc.label = m.calc.label + '_init'
            temp_components.append(m.calc.components + [coord_calc])
        except AttributeError:
            pass

    components = unwind(temp_components)

    mol_labels = [m.calc.label for m in list_mols]
    component_labels = [m.calc.label for m in components]
    total_labels = mol_labels + component_labels

    missing_labels = [l for l in total_labels if l + '.log' not in os.listdir(".")]
    if missing_labels:
        warnings.warn(
            RuntimeWarning(
                'Attempting to clean missing files: {m}'.format(m=" ".join(missing_labels))
            )
        )

    remove_labels = [l for l in total_labels if l not in missing_labels]

    for l in remove_labels:
        os.remove(l + '.log')

        try:
            os.remove(l + '.com')
        except OSError:
            pass
        try:
            os.remove(l + '.fchk')
        except OSError:
            pass