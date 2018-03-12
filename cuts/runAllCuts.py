import csv
import time
modelset = 'cuts/allmodels_withcuts.txt'
# timestring = time.strftime("%Y-%m%d-%H%M")
timestring = time.strftime("%Y-%m%d")
scriptfile = 'cuts/runAllCuts-%s.sh' % (timestring)
date = '2018-03-10'
tag = 'sys0'
comment = 'Fill systematic testing histograms'
with open(modelset, 'r') as csvfile, open(scriptfile, 'w+') as ofile:
    ofile.write('# %s\n' % (comment))
    ofile.write('DATE=%s; TAG=%s\n' % (date,tag))
    r = csv.reader(csvfile, skipinitialspace=True, delimiter=',')
    for row in r:
        for cut in ['cut1', 'cut2', 'cut3', 'cut4', 'cut5', 'cut6', 'cut7', 'cut8',]:
            model = row[0]
            # cut   = row[1]
            # Run command
            # print model, cut
            command = r'SAMPLE=%s ; CUT=%s; ./condor_scripts/condor_test.py 1 "./main -c configs/config.txt -u $CUT -d ~/data/condor_output/$DATE/histos -l $TAG-$CUT -s $SAMPLE -n  100"' % (model, cut)
            ofile.write(command+'\n')
print 'head %s' % scriptfile
print 'source %s' % scriptfile

