
__author__ = 'clyde'

#see https://github.com/jsspencer/job_manager for ideas

"""logic:

perform stable calculation after all calculations.

If we are optimizing the geometry and the final answer is unstable, we restart the optimisation reading the wavefunction
from the stable=opt calc.

If we are optimizing the geometry, we fail to converge but the calculation looks like it's converging we restart the optimization. If we are finding a minimum we restart from the lowest energy step otherwise we restart from the last step?

If we have reached a final geometry following a geometry optimization and the wavefunction is stable we run a frequency calculation.
"""

import os
import time
import copy

from .data_extract_utils import get_last_lines
from ase_extensions import remote
import ConfigParser
import functools

config = ConfigParser.RawConfigParser()
config.read(os.path.expanduser('~/.cc_notebook.ini'))

#, server = config.get('gaussian', 'gauss_host').split('@')


class Job(object):
    def __init__(self, procs=1, memory=400, walltime=1, queue=None):
        self.procs = procs
        self.memory = memory
        self.walltime = walltime
        self.queue = queue

    def exec_command(self, script_fn):
        """creates command to submit the given script to pbs with the given kwargs using the class job options"""
        start = 'qsub -l ncpus={n},memory={m}mb,time={t}:00:00'.format(n=self.procs, m=self.memory, t=self.walltime)
#        args_part = ','.join(['{key}={value}'.format(key=k, value=kwargs[k]) for k in kwargs])
        end = ' -q {q} {script}'.format(q=self.queue, script=script_fn)
        command = start + end
        return command

    def gen_header(self):
        if self.queue:
            header = """#PBS -l ncpus={n}
#PBS -l walltime={t}:00:00
#PBS -l mem={m}mb
#PBS -q {q}
#PBS -j oe
#PBS -V

""".format(n=self.procs, t=self.walltime, m=self.memory, q=self.queue)
        else:
            header = """#PBS -l ncpus={n}
#PBS -l walltime={t}:00:00
#PBS -l mem={m}mb
#PBS -j oe
#PBS -V

""".format(n=self.procs, t=self.walltime, m=self.memory)

        return header


def on_server(nodes=1, memory=400, time=1, queue=''):
    #would be nicer to use nonlocal so that we could just use job all the way through rather than outer_job/job
    #but this is only available in python 3
    if nodes == 1 and memory == 400 and time == 1 and not queue:
        outer_job = None
    else:
        outer_job = Job(nodes, memory + nodes * 150, time, queue)

    def on_server_inner(fn):
        """decorator for functions of ASE_objects"""
        from ase_extensions import ase_utils

        @functools.wraps(fn)
        def server_fn(mol=None, *args, **kwargs):

            job = outer_job

            # if first/mol argument is an ase molecule object with a calculator attached
            # and job details haven't been specified take the job details from the
            # calculation
            if mol and not job:
                try:
                    inner_nodes = mol.calc.job_params['nodes']
                    inner_memory = mol.calc.job_params['memory']
                    inner_time = mol.calc.job_params['time']
                    inner_queue = mol.calc.job_params['queue']
                    job = Job(inner_nodes, inner_memory + inner_nodes * 150, inner_time, inner_queue)
                except AttributeError:
                    pass

            # if mol:   ##unnecessary?
            args = [mol] + list(args)

            return ase_utils.run_on_server((fn, job), *args, **kwargs)
        return server_fn
    return on_server_inner

def send_to_cx1_home(filename):
    """sends the input file in the working dir to the server's home directory"""


    serv = config.get('gaussian', 'gauss_host')
    home = config.get('gaussian', 'gauss_home')
    local_home = config.get('ase', 'ase_home')

    try:
        active_dir = os.getcwd().split(local_home)[1]
    except IndexError:
        raise RuntimeError('Not running from within ASE_HOME')

    home_dir = home + active_dir

    if not server_files_equal(home_dir + '/' + filename, filename):
        exitcode= os.system('scp "%s" "%s:%s"' % (filename, serv,home_dir) )

        if exitcode:
            raise RuntimeError('Unable to send file to server, scp exited with %s' % exitcode)

def get_from_cx1_scratch(filename):
    """copies output from the server's scratch dir to the working dir"""
    #serv = os.environ['GAUSS_HOST']
    #scratch = os.environ['GAUSS_SCRATCH']
    #local_home = os.environ['ASE_HOME']

    serv = config.get('gaussian', 'gauss_host')
    scratch = config.get('gaussian', 'gauss_scratch')
    local_home = config.get('ase', 'ase_home')


    try:
        active_dir = os.getcwd().split(local_home)[1]
    except IndexError:
        raise RuntimeError('Not running from within ASE_HOME')

    scratch_dir = scratch + active_dir

    #only copy if files are not equal
    if not server_files_equal(scratch_dir + '/' + filename, filename):
        exitcode = os.system('scp "%s:%s/%s" "."' % (serv,scratch_dir,filename) )

        if exitcode:
            raise RuntimeError('Unable to get file from server, scp exited with %s' % exitcode)

def server_file_exists(serv_file, sftp=None):
    """checks a given file exists on the server, if an sftp connection is provided that is used (and subsequently left open)"""
    sftp_gen = False

    if not sftp:
        sftp_gen = True
        ssh, sftp = remote.connect_server(ssh=True, sftp=True)
    try:
        sftp.stat(serv_file)
        return True
    except IOError:
        return False
    finally:
        if sftp_gen:
            ssh.close()
            sftp.close()

from subprocess import Popen, PIPE
def server_files_equal(serv_file, local_file):
    """checks whether a local file is the same as the file on a server by comparing last 10 lines of the file"""
    try:
        #{fl} enclosed in quotes to allow for directories with spaces in them
        serv_command = "head -n100 '{fl}'; tail -n10 '{fl}'".format(fl=serv_file)
        local_command = "head -n100 '{fl}'; tail -n10 '{fl}'".format(fl=os.path.realpath(local_file))
        ssh = remote.connect_server(ssh=True)
        stdin, stdout, stderr = ssh.exec_command(serv_command)
        serv_last_lines = stdout.read()
        ssh.close()
    except OSError:
        return False

    try:
        p = Popen(local_command, stdout=PIPE, stderr=PIPE, shell=True)
        stdout, stderr = p.communicate()
        home_last_lines = stdout
    except OSError:
        return False

    if serv_last_lines == home_last_lines:
        return True
    else:
        return False


def server_files_equal_v2(serv_files, local_files):
    """checks whether a list of files on the server are the same as a list of local by comparing line 133 and the last 10 lines of the files,
    returns a list of booleans indicating whether each file is equivalent"""

    if len(serv_files) != len(local_files):
        raise AttributeError('Cannot compare an unequal list of files')

    #{fl} enclosed in quotes to allow for directories with spaces in them
    serv_commands = ["sed -n 133p '{fl}'; tail -n10 '{fl}'; echo server_files_equal_v2_chunk_done;".format(fl=serv_file) for serv_file in serv_files]
    serv_command = "".join(serv_commands)
    ssh = remote.connect_server(ssh=True)
    stdin, stdout, stderr = ssh.exec_command(serv_command)
    #last element is empty space because of the way split works
    serv_lines = stdout.read().split('server_files_equal_v2_chunk_done\n')[0:-1]
    ssh.close()

    local_commands = ["sed -n 133p '{fl}'; tail -n10 '{fl}'; echo server_files_equal_v2_chunk_done;".format(fl=os.path.realpath(local_file)) for local_file in local_files]
    local_command = "".join(local_commands)
    p = Popen(local_command, stdout=PIPE, stderr=PIPE, shell=True)
    stdout, stderr = p.communicate()
    home_lines = stdout.split('server_files_equal_v2_chunk_done\n')[0:-1]

    return [home_lines[i] == serv_lines[i] for i in range(len(home_lines))]


def server_data_unequal(list_mols):
    from ase_extensions.ase_utils import get_active_dirs

    """checks whether calculation data from log files on the server is the same as local calculation data by comparing fingerprints of the files"""
    fingerprints = [mol.calc.fingerprint for mol in list_mols]

    home_dir,scratch_dir = get_active_dirs()

    home_files = [mol.calc.label + '.log' for mol in list_mols]
    serv_files = [scratch_dir + '/' + fn for fn in home_files]


    serv_commands = ["sed -n 133p '{fl}'; tail -n10 '{fl}'; echo server_files_equal_v2_chunk_done;".format(fl=serv_file) for serv_file in serv_files]
    serv_command = "".join(serv_commands)
    ssh = remote.connect_server(ssh=True)
    stdin, stdout, stderr = ssh.exec_command(serv_command)
    #last element is empty space because of the way split works
    serv_fingerprints = stdout.read().split('server_files_equal_v2_chunk_done\n')[0:-1]
    ssh.close()

    return [fingerprints[i] != serv_fingerprints[i] for i in range(len(fingerprints))]

#duplication with gaussian_calculator_obj.status
def extract_status(calc_obj):
    """extracts the status of a gaussian calculation given either an ase molecular object with a gaussian calculator attached or the gaussian log file"""
    try:
        filename = calc_obj.calc.label + '.log'
    except AttributeError:
        if os.path.isfile(calc_obj):
            filename = calc_obj
        elif '.log' in calc_obj:
            return 'Queued'
        else:
            raise RuntimeError('Cannot extract status from {c} - unrecognised object'.format(c=calc_obj))

    command = "tail -n4 '{fl}'".format(fl=filename)
    p = Popen(command, stdout=PIPE, stderr=PIPE, shell=True)
    stdout, stderr = p.communicate()

    last_lines = stdout.split('\n')

    if any('Normal termination' in l for l in last_lines):
        return 'Success'

    try:
        i = next(l[0] for l in enumerate(last_lines) if 'Error termination' in l[1])
    except StopIteration:
        return 'Incomplete'

    return 'Error: {e}'.format(e=last_lines[i-1].strip())

class Job_Manager(object):
    def __init__(self, gaus_calc):
        #starting calculation
        self.calc = gaus_calc
        self.main_calc = copy.deepcopy(self.calc)
        #store of calculations manager performs
        self.calcs = []
        self.restart_details = []

        self.send_to_home = send_to_cx1_home
        self.get_from_scratch = get_from_cx1_scratch

#        self.geom_methods = []
#        self.energy_methods = []

        self.job_complete = False
        self.job_success = False
        self.job_error = False

        self.host_dir = ""
        #self.scratch_dir = ""

        self.set_directories()
        self.id = self.submit_job()
        self.monitor_job(id,master=True)

    def set_directories(self):
        """gets directory on server"""
        #local_home = os.environ['ASE_HOME']
        local_home = config.get('ase', 'ase_home')

        try:
            active_dir = os.getcwd().split(local_home)[1]
        except IndexError:
            raise RuntimeError('Not running from within ASE_HOME')
        #self.host_dir = os.environ['GAUSS_HOME'] + active_dir + '/'
        self.host_dir = config.get('gaussian', 'gauss_home') + active_dir + '/'

    def submit_job(self,verbose=False):
        #store copy of calculation
        self.calcs.append(copy.deepcopy(self.calc))
        #copy the input file to the equivalent location on the server
        self.send_to_home(self.calc.label + '.com')

        #run submission script and collect job info
        ssh = remote.connect_server(ssh=True)
        i,o,e = ssh.exec_command('submit_calc {fld} {inp} {p} {m} {t} {q} {v}'.format(host= config.get('gaussian', 'gauss_host'), fld=self.host_dir, inp=self.calc.label + '.com',
                                                                                      p=self.calc.job_params['nodes'], m=self.calc.job_params['memory'], t=int(self.calc.job_params['time']),
                                                                                      q=self.calc.job_params['queue'], v=self.calc.job_paras['version']))
        qsub_output = o.readlines() + e.readlines()
        ssh.close()

        #sanity check
        if len(qsub_output) == 0:
            raise remote.PBSUtilQSubError("Failed to launch Gaussian %s, qsub gave no stdoutput" % self.calc.label)
        if verbose:
            print '\n%s\n' %  qsub_output

        #extract and return job_id
        pbs_id = remote.parse_qsub_output(qsub_output[0])[0]
        return pbs_id

    #monitor process
    def monitor_job(self, id, master=False):
        #main loop
        while self.is_job_active():
#            if not self.is_job_geom_converging():
#                self.restart_geom()
#            elif not self.is_job_scf_converging():
#                self.restart_energy()
            self.get_job_status()
            time.sleep(5)

        self.get_from_scratch(self.calc.label + '.log')
        self.get_job_status()
        if self.job_success:
            self.check_job_stable()

        if self.job_stable() and self.is_opt_job():
            id = self.run_freq_calc()
            self.monitor_job(id)

        if master:
            return self.get_log()

    def run_freq_calc(self):
        self.calc = self.calc.restart(add_label = 'freq_', freq='freq', geom='allcheck', scf='read',opt=None)
        self.calc.write_input(self.calc.label + '.com', self.calc.atoms)
        self.id = self.submit_job()
        self.restart_details.append({'Option':'freq=freq, opt=None', 'Reason': 'Computing frequencies following geometry optimisation', 'File':self.calc.label +'.log'})

    def run_stable_calc(self):
        self.calc = self.calc.restart(add_label = 'stable', stable='restart',opt=None)
        self.calc.write_input(self.calc.label + '.com', self.calc.atoms)
        self.id = self.submit_job()
        self.restarts_details.append({'Option':'stable=restart', 'Reason': 'Computing calculation stability', 'File':self.calc.label +'.log'})

    def get_log(self):
        log = []
        log.append('Job {id}\n'.format(id=self.id))
        log.append('{n} restarts\n'.format(n=len(self.restart_details)))
        for detail in self.restart_details:
            log.append('Restarted with {o} because {r} restart_file = {f}\n'.format(o=detail['Option'], r=detail['Reason'], f=detail['File']))
        return " ".join(log)

    #creates and returns an ase_object associated with the current calculation
    def get_calc(self):
        from gausspy.io.gaussian import read_gaussian_out
        ase_obj = read_gaussian_out(self.main_calc.label + '.log')
        ase_obj.calc = self.main_calc
        return ase_obj

    #creates returns ase_objects associated with all calculations performed
    def get_all_calcs(self):
        from gausspy.io.gaussian import read_gaussian_out
        list_ase_objs = []

        for calc in self.calcs:
            ase_obj = read_gaussian_out(calc.label +'.log')
            ase_obj.calc=calc
            list_ase_objs.append(ase_obj)

        return list_ase_objs

    def load_data(self,inp_file):
        with open(inp_file,'r') as f:
            return f.readlines()

    def is_job_active(self):
        if self.id in [e.id for e in remote.qstat_plain()]:
            return True
        else:
            return False

    def get_job_status(self):
        if not self.is_job_active():
            last_outp_lines = "".join(get_last_lines(self.calc.label + '.log',5))
            if 'Normal termination of Gaussian' in last_outp_lines:
                self.job_success = True
            elif 'Error' in last_outp_lines:
                self.job_error = True

    def is_opt_job(self):
        inp_data = self.load_data(self.calc.label+'.com')
        if 'opt ' in inp_data[0]:
            return True
        else:
            return False

    def is_qc_job(self):
        inp_data = self.load_data(self.calc.label+'.com')
        if 'qc ' in inp_data[0]:
            return True
        else:
            return False


#    def is_job_geom_converging(self):
#        return True
#
#    def is_job_scf_converging(self):
#        return True

#    def restart_job(self, method, step):
#        return

#    def get_best_geom(self):
#        return

#    def restart_geom(self):
#        best_geom = self.get_best_geom()
#        next_method = next(self.geom_methods)
#        self.restart_job(next_method, best_geom)

#    def restart_energy(self):
#        return
