import csv
import time
cutset = 'cuts/cutset_combined_1122.txt'
# timestring = time.strftime("%Y-%m%d-%H%M")
timestring = time.strftime("%Y-%m-%d")
scriptfile = 'cuts/runCutTest-%s.sh' % (timestring)
date = timestring
tag = 'cuttest1'
model = 'mass_X_d_1000_mass_pi_d_1_tau_pi_d_5'
with open(cutset, 'r') as csvfile, open(scriptfile, 'w+') as ofile:
    next(csvfile)
    ofile.write('DATE=%s;\n' % (date))
    r = csv.reader(csvfile, skipinitialspace=True, delimiter=',')
    for row in r:
        cut = row[0]
        cuttag = '%s-%s' % (tag, cut)
        # Run command
        # print model, cut
        command = r'SAMPLE=%s ; CUT=%s; ./condor_scripts/condor_test.py 1 "./main -c configs/config.txt -u $CUT -d ~/data/condor_output/$DATE/cuttest/histos -l %s -s $SAMPLE -n  1"' % (model, cut, cuttag)
        ofile.write(command+'\n')
print 'head %s' % scriptfile
print 'source %s' % scriptfile

