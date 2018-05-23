import csv
import time
# modelset = 'cuts/modelset_0105_combined.txt'
# modelset = 'cuts/emjet_cutmap_bestcuts_20180215.txt'
modelset = 'cuts/emjet_cutmap_bestcuts_20180402.txt'
# timestring = time.strftime("%Y-%m%d-%H%M")
timestring = time.strftime("%Y-%m%d")
scriptfile = 'cuts/runSignalCuts-%s.sh' % (timestring)
date = '2018-04-26'
tag = 'acc0'
comment = 'nEm >= N'
with open(modelset, 'r') as csvfile, open(scriptfile, 'w+') as ofile:
    ofile.write('# %s\n' % (comment))
    ofile.write('DATE=%s; TAG=%s\n' % (date,tag))
    r = csv.reader(csvfile, skipinitialspace=True, delimiter=',')
    for row in r:
        model = row[0]
        cut   = row[1]
        # Run command
        # print model, cut
        command = r'SAMPLE=%s ; CUT=%s; ./condor_scripts/condor_test.py 1 "./main -c configs/config.txt -u $CUT -d ~/data/condor_output/$DATE/histos -l $TAG -s $SAMPLE -n  1000"' % (model, cut)
        ofile.write(command+'\n')
print 'head %s' % scriptfile
print 'source %s' % scriptfile

