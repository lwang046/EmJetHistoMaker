import os
import subprocess
from ROOT import gROOT, TFile, TH1D
from collections import OrderedDict

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-f', '--force', action='store_true')
parser.add_argument('directory')
args = parser.parse_args()

os.chdir(args.directory)

samplesets = OrderedDict() # Dictionary of list of tuples. Tuple contains (label, xsec_in_pb)
samplesets['QCD'] = [
    ( 'QCD_HT100to200'   ,  27990000.00 ), # pb
    # ( 'QCD_HT200to300'   ,   1712000.00 ), # pb FIXME
    ( 'QCD_HT300to500'   ,    347700.00 ), # pb
    ( 'QCD_HT500to700'   ,     32100.00 ), # pb
    ( 'QCD_HT700to1000'  ,      6831.00 ),
    ( 'QCD_HT1000to1500' ,      1207.00 ),
    ( 'QCD_HT1500to2000' ,       119.90 ),
    ( 'QCD_HT2000toInf'  ,        25.24 ),
]

# samplesets['ModelA'] = [( 'ModelA', 14.6 )]

# samplesets['DataBF']   = [
#     ( 'JetHT_B3', 1.0 ),
#     ( 'JetHT_C1', 1.0 ),
#     ( 'JetHT_D1', 1.0 ),
#     ( 'JetHT_E1', 1.0 ),
#     ( 'JetHT_F1', 1.0 ),
# ]

samplesets['DataGH']   = [
    ( 'JetHT_G1', -1.0 ), # Put negative cross section for data
    ( 'JetHT_H2', -1.0 ),
    ( 'JetHT_H3', -1.0 ),
]

# hdirc = OrderedDict()

# hdirc['QCD_HT500to700']   =   32100.00 # pb
# hdirc['QCD_HT700to1000']  =    6831.00
# hdirc['QCD_HT1000to1500'] =    1207.00
# hdirc['QCD_HT1500to2000'] =     119.90
# hdirc['QCD_HT2000toInf']  =      25.24

lumi = 16.13e3 # pb^-1

haddcommand = '/home/yhshin/haddws/haddws'

for sampleset, samplelist in samplesets.iteritems():
    haddwsargs = ''
    print ('# Hadding sample set: %s' % sampleset)
    print ('# ====================================')
    print 'sample, xsec, eventCountPreTrigger, scaling'
    for sample in samplelist:
        samplefile = 'histo-%s.root' % sample[0]
        xsec = sample[1]
        if xsec > 0 :
            filehandle = TFile(samplefile)
            eventCountHist = filehandle.Get("eventCountPreTrigger")
            eventCountPreTrigger = eventCountHist.Integral()
            scaling = xsec / eventCountPreTrigger
            print '%s, %f, %d, %f' % (sample[0], xsec, eventCountPreTrigger, scaling)
        else:
            scaling = 1.0
        haddwsargs += ' %s %f' % (samplefile, scaling)
    print ('# ====================================')
    print(haddcommand + haddwsargs)
    subprocess.call(haddcommand + haddwsargs, shell=True)
    print ('# ------------------------------------')
    mvargs = ' result.root histo-%s.root' % sampleset
    print("mv" + mvargs)
    subprocess.call("mv" + mvargs, shell=True)
    print ('# ====================================')
    print ('')


#
# merge histograms in one HT bin
#
# print 'starting merging different root files in one HT bin to one root file...'
# for idirc in hdirc:
#     subprocess.call("hadd histo_%s.root ./%s/*.root ./%sx/*.root"%(idirc, idirc, idirc), shell=True, executable='/bin/bash')
#     #subprocess.call("hadd histo_%s.root ./%s/*.root"%(idirc, idirc), shell=True, executable='/bin/bash')

# #
# # merge histograms from different HT bins
# #
# print 'starting merging different HT bins to one root file...'
# tstring = ''
# for idirc in hdirc:
#     afile = TFile("histo_%s.root"%idirc)
#     htemp = afile.Get("eventCountProcessed")
#     nevents = htemp.Integral()
#     tscale = lumi * hdirc[idirc]/nevents
#     # hcount = afile.Get("nJets_tag")
#     # ncount = hcount.GetBinContent(3)
#     # ncounterr = hcount.GetBinError(3)
#     # print 'total number of events in %s is %.2f, Emerging events is %.2f+/-%.2f, scale factor is %f'%(idirc, nevents, ncount, ncounterr,tscale)
#     tstring += ' histo_%s.root %f'%(idirc, tscale)

# subprocess.call("/home/yhshin/haddws"+tstring, shell=True)
# subprocess.call("mv result.root histo_QCDHT_merged.root", shell=True)

#
# print out the number of 'Emerging' events after merged
#
# afile = TFile("histo_QCDHT_merged.root")
# hcount = afile.Get("nJets_tag")
# print "QCD weighted count: %.2f+/-%.2f"%(hcount.GetBinContent(3), hcount.GetBinError(3))
