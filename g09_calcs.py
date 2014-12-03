#! /usr/bin/env python
__author__ = 'clyde'
import os
import sys
import ConfigParser

config = ConfigParser.RawConfigParser()
config.read(os.path.expanduser('~/.cc_notebook.ini'))
user = config.get('pbs', 'user')


def construct_job_script(procs=None, memory=None, time=None, queue='None'):
    """Construct pbs submission script"""
    if not procs:
        procs = 1

    if not memory:
        memory = procs * 1400

    if not time:
        time = 3

    #bash script that we submit to the queue
    if queue == 'None':
        lhead_str = ['#PBS -l ncpus={p}'.format(p=procs),
                     '#PBS -l mem={m}mb'.format(m=memory+(procs*100)),
                     '#PBS -l walltime={t}:00:00'.format(t=time),
                     '#PBS -j oe']
        head_str = "\n".join(lhead_str)

    else:
        lhead_str= ['#PBS -l ncpus={p}'.format(p=procs),
                    '#PBS -l mem={m}mb'.format(m=memory+(procs*100)),
                    '#PBS -l walltime={t}:00:00'.format(t=time),
                    '#PBS -j oe',
                    '#PBS -q {q}'.format(q=queue)]
        head_str = "\n".join(lhead_str) + '\n'

    ljob_str = ["module load gaussian/g09-d01",
                'date=`date "+%d-%m-%y %r"`',
                'echo "$value: $HOME/$FLD$FLNM started on $date1" >> /work/{u}/logfile.txt'.format(u=user),
                'g09 < $HOME/$FLD$FLNM > $WORK/$FLD${FLNM%.*}.log',
                'date2=`date "+%d-%m-%y %r"`',
                'echo "$value: $HOME/$FLD$FLNM finished on $date2" >> /work/{u}/logfile.txt'.format(u=user)]

    return head_str + "\n".join(ljob_str)


def submit_job(fold, inp_fn, procs=None, memory=None, time=None, queue='None'):
    """Submits gaussian job to the pbs queue"""
    if not procs:
        procs = 1

    if not memory:
        memory = procs * 1400

    if not time:
        time = 3

    local_fold = fold.split('/home/' + user + '/')[1]
    work_dir = '/work/' + local_fold

    #bash script that we submit to the queue
    script_str = construct_job_script(procs, memory, time, queue)

    name = inp_fn.split('.')[0]
    job_fn = name + '_job.sh'

    with open(fold+job_fn, 'w') as job_f:
        job_f.write(script_str)

    #addition to top of input file
    linp_str = ['%mem={m}MB'.format(m=memory),
                '%nproc={p}'.format(p=procs),
                '%chk={dir}{file_n}.chk'.format(dir=work_dir, file_n=name)]
    inp_str = "\n".join(linp_str)

    with open(fold+inp_fn, 'r') as inp_f:
        f_contents = inp_f.readlines()

    if '%mem' not in f_contents[0]:
        new_contents = [inp_str] + f_contents

        with open(fold+inp_fn, 'w') as inp_f:
            inp_f.writelines(new_contents)

    #submission to queue
    job_id = os.system("cd {f}; /opt/pbs/default/bin/qsub -v FLD={f1},FLNM={n} {j}".format(f=fold, f1 = local_fold, n=inp_fn, j=fold+job_fn))

    return job_id

if __name__ == '__main__':
    fold = inp_fn = procs = mem = time = queue = None

    args = sys.argv[1:]

    fold = args[0]
    inp_fn = args[1]

    if len(args) > 2:
        procs = args[2]
    if len(args) > 3:
        mem = args[3]
    if len(args) > 4:
        time = args[4]
    if len(args) > 5:
        queue = args[5]

    submit_job(fold, inp_fn, int(procs), int(mem), int(time), queue)
