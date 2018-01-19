import csv
import time
# modelset = 'cuts/modelset_0105.txt'
modelset = 'cuts/modelset_0105b.txt'
# timestring = time.strftime("%Y-%m%d-%H%M")
timestring = time.strftime("%Y-%m%d")
scriptfile = 'cuts/runSignalCuts-%s.sh' % (timestring)
date = '2018-01-15'
tag = 'acc0'
with open(modelset, 'r') as csvfile, open(scriptfile, 'w+') as ofile:
    ofile.write('DATE=%s; TAG=%s\n' % (date,tag))
    r = csv.reader(csvfile, skipinitialspace=True, delimiter=',')
    for row in r:
        model = row[0]
        cut   = row[1]
        # Run command
        # print model, cut
        command = r'SAMPLE=%s ; CUT=%s; ./condor_scripts/condor_test.py 1 "./main -c configs/config.txt -u $CUT -d ~/data/condor_output/$DATE/histos -l $TAG -s $SAMPLE -n  100"' % (model, cut)
        ofile.write(command+'\n')
print 'head %s' % scriptfile
print 'source %s' % scriptfile

