#!/usr/bin/env python

import fileinput
import numpy
import os
import pickle
import shutil
import struct
import subprocess
import sys
import tempfile

def ode(model, t, i, species = None, parameters = None, rateConstants = None):
    # Build temporary files
    outputdir = tempfile.mkdtemp()
    [mfd, modelfile] = tempfile.mkstemp()

    path = os.path.abspath(os.path.dirname(__file__))

    mfhandle = os.fdopen(mfd, 'w')
    mfhandle.write(model)
    mfhandle.close()

    shutil.rmtree(outputdir)

    args = []

    args.append('-m ' + modelfile)

    args.append('-t ' + str(t))
    args.append('-i ' + str(i))
    if species != None:
        args.append('--species ' + " ".join(map(lambda x : x.strip(), species)))

    if parameters and rateConstants:
        raise Exception("parameters and rateConstants cannot be defined at the same time (mutually exclusive)")

    if parameters or rateConstants:
        args.append('--sensi')

    if parameters != None:
        args.append('--parameters ' + " ".join(map(lambda x : x.strip(), parameters)))
    if rateConstants != None:
        args.append('--rate-constants ' + " ".join(map(lambda x : x.strip(), rateConstants)))

    args.append('--out-dir ' + outputdir)
    args.append('--force')

    print '{0}/bin/stochkit_ode {1}'.format(path, args)

    process = subprocess.Popen('{0}/bin/stochkit_ode {1}'.format(path, " ".join(args)).split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, error = process.communicate()

    if process.returncode != 0:
        raise Exception("""Stochkit_ode failed: 
out : {0}
error : {1}""".format(out,error))
    
    # Collect all the output data
    vhandle = open(outputdir + '/output.txt', 'r')
    values = { 'time' : [], 'trajectories' : {}, 'sensitivities' : {}}
    columnToList = []
    for i, line in enumerate(vhandle):
        if i == 0:
            continue
        elif i == 1:
            names = line.split()
            for name in names:
                if name == 'time':
                    columnToList.append(values['time'])
                elif ':' in name:
                    specie, parameter = name.split(':')
                    if specie not in values['sensitivities']:
                        values['sensitivities'][specie] = {}

                    values['sensitivities'][specie][parameter] = [] # Make a new timeseries for sensitivity
                    columnToList.append(values['sensitivities'][specie][parameter]) # Store a reference here for future use
                else:
                    values['trajectories'][name] = [] # start a new timeseries for this name
                    columnToList.append(values['trajectories'][name]) # Store a reference here for future use
        elif i == 2:
            for storage, value in zip(columnToList, map(float, line.split())):
                storage.append(value)
        elif i == 3:
            continue
        else:
            for storage, value in zip(columnToList, map(float, line.split())):
                storage.append(value)
    vhandle.close()

    shutil.rmtree(outputdir)
    os.remove(modelfile)
    return values
