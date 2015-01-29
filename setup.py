#!/usr/bin/env python

from setuptools import setup
setup(name='gausspy',
      version='0.1',
      packages=['gausspy', 'gausspy.io'],

      entry_points = {
        'console_scripts': [
            'submit_calc = gausspy.g09_pbs:submit_pbs_calc',
        ],
        }
      )
