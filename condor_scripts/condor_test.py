#!/bin/env python
"""Executes arbitrary command on condor. Adds $(ProcId) to arguments.
Edit BASEDIR, DATESTR"""
import sys
import subprocess
from string import Template
import time
import os # path, makedirs

def execute(args):
    """Wrapper function to execute arbitrary shell commands"""
    print '################################'
    print 'args: ', args
    p = subprocess.Popen(args, shell=True, executable='/bin/bash')
    # p = subprocess.call(args, shell=True, executable='/bin/bash')
    p.wait()
    return p
    print '################################'

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-t', '--test', action='store_true', help='Test interactively')
parser.add_argument('queue', type=int)
parser.add_argument('command')
args = parser.parse_args()
print args

command_list = args.command.split()
# BASEDIR = "/afs/cern.ch/user/y/yoshin/work/condor_output"
BASEDIR = "/home/yhshin/data/condor_output"
DATESTR = time.strftime("%Y-%m-%d")
DIRECTORY = os.path.join(BASEDIR, DATESTR) # DIRECTORY = BASEDIR/YY-MM-DD
if not os.path.exists(DIRECTORY):
    os.makedirs(DIRECTORY)
kw_dict = {}
kw_dict['DIRECTORY'] = DIRECTORY
kw_dict['EXECUTABLE'] = command_list[0]
kw_dict['ARGUMENTS']  = " ".join(command_list[1:])
kw_dict['QUEUE']  = "%d" % args.queue
kw_dict['DATE'] = time.strftime("%Y-%m-%d")

# For lxplus or generic condor clusters
jdl_template ="""
universe = vanilla
executable            = ${EXECUTABLE}
arguments             = ${ARGUMENTS} $(ProcId)
output                = ${DIRECTORY}/condor_test.$(ClusterId).$(ProcId).out
error                 = ${DIRECTORY}/condor_test.$(ClusterId).$(ProcId).err
log                   = ${DIRECTORY}/condor_test.$(ClusterId).log
send_credential        = True
queue ${QUEUE}
"""

# For hepcms SL6
jdl_template ="""
Universe = vanilla
Requirements = TARGET.FileSystemDomain == "privnet" && machine != "r510-0-1.privnet"
request_memory = 1024
+IsHighPriorityJob = True
Executable            = ${EXECUTABLE}
Arguments             = ${ARGUMENTS} $(ProcId)
Output                = ${DIRECTORY}/condor_test.$(ClusterId).$(ProcId).out
Error                 = ${DIRECTORY}/condor_test.$(ClusterId).$(ProcId).err
Log                   = ${DIRECTORY}/condor_test.$(ClusterId).log
Queue ${QUEUE}
"""

t = Template(jdl_template)
jdl = t.safe_substitute(kw_dict)

jdlfile = open('test.jdl', 'w')
jdlfile.write(jdl)
jdlfile.close()

if args.test:
    execute(args.command + " 0")
else:
    execute("condor_submit  getenv=True test.jdl")
