from unittest import TestCase
from ..g09_pbs import construct_job_script
__author__ = 'clyde'


class TestConstruct_job_script(TestCase):
    def setUp(self):
       self.example_job_script = (  '#PBS -l ncpus=1\n'
                                    '#PBS -l mem=500mb\n'
                                    '#PBS -l walltime=3:00:00\n'
                                    '#PBS -j oe\n'
                                    'module load gaussian/g09\n'
                                    'date=`date "+%d-%m-%y %r"`\n'
                                    'echo "$value: $HOME/$FLD$FLNM started on $date1" >> /work/cjf05/logfile.txt\n'
                                    'g09 < $HOME/$FLD$FLNM > $WORK/$FLD${FLNM%.*}.log\n'
                                    'date2=`date "+%d-%m-%y %r"`\n'
                                    'echo "$value: $HOME/$FLD$FLNM finished on $date2" >> /work/cjf05/logfile.txt'  )

    def test_0(self):
        procs = 1
        memory = 400
        time = 3
        queue = 'None'
        version = 'g09'
        test_script = construct_job_script(procs=procs, memory=memory, time=time, queue=queue, version=version)
        print(test_script)
        assert test_script == self.example_job_script

