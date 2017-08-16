import ROOT as rt
rt.TH1.SetDefaultSumw2()
import math
from collections import OrderedDict
# from rootutils import getTrees
files= OrderedDict()
# files['WJet'] = "/afs/cern.ch/user/y/yoshin/eos/cms/store/group/phys_exotica/EmergingJets/SelectedJetAnalysis-v0/WJetSkimMuon/SingleMuon/SelectedJetAnalysis/160107_215042/output_merged_WJetSkimMuon.root"
# files['ModelA'] = "/afs/cern.ch/user/y/yoshin/eos/cms/store/group/phys_exotica/EmergingJets/SelectedJetAnalysis-v0/ModelA/EmergingJets_ModelA_TuneCUETP8M1_13TeV_pythia8Mod/SelectedJetAnalysis/160107_220503/output_merged_ModelA.root"
# files['ModelB'] = "/afs/cern.ch/user/y/yoshin/eos/cms/store/group/phys_exotica/EmergingJets/SelectedJetAnalysis-v0/ModelB/EmergingJets_ModelB_TuneCUETP8M1_13TeV_pythia8Mod/SelectedJetAnalysis/160107_220526/output_merged_ModelB.root"
# files['MC'] = "/afs/cern.ch/user/y/yoshin/CMSSW_7_4_12/src/EmergingJetAnalysis/ntuple.root"
# files['ModelA_Analysis-20160314'] = "/afs/cern.ch/user/y/yoshin/eos/cms/store/group/phys_exotica/EmergingJets/Analysis-20160314-v0/ModelA/EmergingJets_ModelA_TuneCUETP8M1_13TeV_pythia8Mod/Analysis-20160314/160314_234710/ntuple_merged_ModelA.root"
# files['Analysis-20160322'] = "/afs/cern.ch/user/y/yoshin/eos/cms/store/group/phys_exotica/EmergingJets/Analysis-20160322-v0/ModelA/EmergingJets_ModelA_TuneCUETP8M1_13TeV_pythia8Mod/Analysis-20160322/160322_171421/ntuple_merged_ModelA.root"
# files['Analysis-20160325'] = "/afs/cern.ch/user/y/yoshin/eos/cms/store/group/phys_exotica/EmergingJets/Analysis-20160325-v0/ModelA/EmergingJets_ModelA_TuneCUETP8M1_13TeV_pythia8Mod/Analysis-20160325/160325_192306/ntuple_merged_ModelA.root"
# files['Analysis-20160615'] = "/afs/cern.ch/user/y/yoshin/eos/cms/store/group/phys_exotica/EmergingJets/Analysis-20160615-v0/ModelA/EmergingJets_ModelA_TuneCUETP8M1_13TeV_pythia8Mod/Analysis-20160615/160615_220449/ntuple_merged_ModelA.root"
# files['Analysis-20160821'] = "/afs/cern.ch/user/y/yoshin/eos/cms/store/group/phys_exotica/EmergingJets/Analysis-20160821-v0/ModelA/EmergingJets_ModelA_TuneCUETP8M1_13TeV_pythia8Mod/Analysis-20160821/160821_201557/ntuple_merged_ModelA.root"
files['Analysis-20170426'] = "/store/user/yoshin/EmJetAnalysis/Analysis-20170426-v2/ModelA/EmergingJets_ModelA_TuneCUETP8M1_13TeV_pythia8Mod/Analysis-20170426/170426_192610/0000/ntuple_1.root"
dirname = 'emJetAnalyzer'
treename = 'emJetTree'
classname = 'BaseClass'

index = 0
for k,v in files.iteritems():
    name = k
    filename = v
    f = rt.TFile(filename)
    tree = f.Get(dirname).Get(treename)
    tree.MakeClass(classname)
    # raw_input("Press Enter to continue...")

# Add appropriate 'using' directives
import fileinput
classheader = classname + '.h'
using_vector = 0
for line in fileinput.input(classheader, inplace=1):
    print line,
    if using_vector==0 and line.startswith('#include "vector"'):
        using_vector = 1
        print 'using std::vector;'

import os
os.rename(classname + '.C', classname + '.cc')

