#!/usr/bin/env python

import os

import sys

import subprocess

path = os.path.abspath(os.path.dirname(__file__))

os.environ['STOCHKIT_ODE'] = path
os.environ['STOCHKIT_HOME'] = os.path.join(path, '../StochKit')
LD_LIBRARY_PATH = ""
DYLD_LIBRARY_PATH = ""

if 'LD_LIBRARY_PATH' in os.environ:
    LD_LIBRARY_PATH = ':' + os.environ['LD_LIBRARY_PATH']

if 'DYLD_LIBRARY_PATH' in os.environ:
    DYLD_LIBRARY_PATH = ':' + os.environ['DYLD_LIBRARY_PATH']

os.environ['LD_LIBRARY_PATH'] = '{0}/../StochKit/libs/boost_1_53_0/stage/lib/'.format(path) + LD_LIBRARY_PATH
os.environ['DYLD_LIBRARY_PATH'] = '{0}/../StochKit/libs/boost_1_53_0/stage/lib/'.format(path) + DYLD_LIBRARY_PATH

process = subprocess.Popen('{0}/bin/stochkit_ode {1}'.format(path, " ".join(sys.argv[1:])).split())

process.wait()
