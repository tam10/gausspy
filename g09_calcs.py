#! /usr/bin/env python
__author__ = 'clyde'
import os
import sys
#from cclib.parser import ccopen

def submit_job(fold, inp_fn, procs, memory, time):

    if not procs:
        procs = 1

    if not memory:
        memory = procs * 1400

    if not time:
        time = 3

    local_fold = fold.split('/home/cjf05/')[1]
    work_dir = '/work/' + local_fold

    #bash script that we submit to the queue
    job_str = """#PBS -l ncpus={p}
#PBS -l mem={m}mb
#PBS -l walltime={t}:00:00
#PBS -j oe

module load gaussian

date=`date "+%d-%m-%y %r"`
echo "$value: $HOME/$FLD$FLNM started on $date1" >> /work/cjf05/logfile.txt

g09 < $HOME/$FLD$FLNM > $WORK/$FLD${{FLNM\%.*}}.log

date2=`date "+%d-%m-%y %r"`
echo "$value: $HOME/$FLD$FLNM finished on $date2" >> /work/cjf05/logfile.txt
""".format(p=procs, t=time,m=memory+(procs*100))

    name = inp_fn.split('.')[0]
    job_fn = name + '_job.sh'

    with open(fold+job_fn, 'w') as job_f:
        job_f.write(job_str)

    #addition to top of input file
    inp_str = r"""%mem={m}MB
%nproc={p}
%chk={dir}{file_n}.chk
""".format(m=memory, p=procs, dir = work_dir, file_n=name)

    with open(fold+inp_fn, 'r') as inp_f:
        f_contents = inp_f.readlines()

    if '%mem' not in f_contents[0]:
        new_contents = [inp_str] + f_contents

        with open(fold+inp_fn, 'w') as inp_f:
            inp_f.writelines(new_contents)

    #submission to queue
    id= os.system("cd {f}; /opt/pbs/default/bin/qsub -v FLD={f1},FLNM={n} {j}".format(f=fold, f1 = local_fold, n=inp_fn, j=fold+job_fn))

    return id

#def load_log(out_fn):
#    out_f = ccopen(out_fn)
#    return out_f.parse() # The following lines are log messages

#def load_logs():
#    gaus_outs = []
#    list_logs = [fn for fn in os.listdir(os.getcwd()) if '.log' in fn]
#
#    for fn in list_logs:
#        gaus_outs.append(load_log(fn))
#    return gaus_outs


def monitor_jobs():
    return

if __name__ == '__main__':
    fold = inp_fn = procs = mem = time = None

    args = sys.argv[1:]

    fold = args[0]
    inp_fn = args[1]

    if len(args) > 2:
        procs = args[2]
    if len(args) > 3:
        mem = args[3]
    if len(args) > 4:
        time= args[4]

    submit_job(fold, inp_fn, int(procs), int(mem), int(time))
