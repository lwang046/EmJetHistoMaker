#!/bin/env python
import json
from pprint import pprint
import subprocess

# subprocess.call('condor_history')

with open('condor_history.json') as data_file:
    data = json.load(data_file)

# pprint(data)
# dataf = [d for d in data if (d['ExitCode']!=0 or d['ExitStatus']!=0)]
# dataf = [d['ClusterId'] for d in data if (d['ExitCode']!=0 or d['ExitStatus']!=0)]
dataf = [d['ClusterId'] for d in data]
pprint(dataf)
