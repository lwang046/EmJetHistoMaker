# EmJetHistoMaker

Histogram making code for Emerging Jet Analysis

## Instructions

### To checkout latest code:

    git clone git@github.com:yhshin11/EmJetHistoMaker.git
    cd EmJetHistoMaker

### To build binary

EmJetHistomaker does not require CMSSW, but it does require a recent version of GCC. The simplest way to obtain this is to do \(cmsenv\) in a recent CMSSW release, e.g. \(CMSSW_7_6_3\).
For example, if you have an existing \(CMSSW_7_6_3\) release area in your home directory, you can do the following:

    cd ~/CMSSW_7_6_3/src
    cmsenv

Then return to the EmJetHistoMaker directory and run the build script:

    ./build.sh

### To define input samples and files:

The directory `configs` contains text files that define samples and their properties. For example, a config file could look like this:

    # sample group , sample name      , sample Type , MC/DATA , isSignal , xsec[fb] , input file
    ModelA         , ModelA           , SIGNAL      , MC      , true     , 18.45    , /home/yhshin/workspace/EmJetHistoMaker/inputs/ModelA.txt
    ModelB         , ModelB           , SIGNAL      , MC      , true     , 18.45    , /home/yhshin/workspace/EmJetHistoMaker/inputs/ModelB.txt
    QCD            , QCD_HT500to700   , QCD         , MC      , false    , 29370000 , /home/yhshin/workspace/EmJetHistoMaker/inputs/QCD_HT500to700.txt
    QCD            , QCD_HT700to1000  , QCD         , MC      , false    ,  6524000 , /home/yhshin/workspace/EmJetHistoMaker/inputs/QCD_HT700to1000.txt
    QCD            , QCD_HT1000to1500 , QCD         , MC      , false    ,  1064000 , /home/yhshin/workspace/EmJetHistoMaker/inputs/QCD_HT1000to1500.txt
    QCD            , QCD_HT1500to2000 , QCD         , MC      , false    ,   121500 , /home/yhshin/workspace/EmJetHistoMaker/inputs/QCD_HT1500to2000.txt
    QCD            , QCD_HT2000toInf  , QCD         , MC      , false    ,    25450 , /home/yhshin/workspace/EmJetHistoMaker/inputs/QCD_HT2000toInf.txt
    WJetSkim       , WJetSkimMuon     , WJet        , DATA    , false    ,        1 , /home/yhshin/workspace/EmJetHistoMaker/inputs/WJetSkim.txt

This defines two MC samples `ModelA` and `ModelB`, each with 18.45 fb cross section. It also defines a QCD sample group composed of 5 different HT-bin samples, each with different cross section. Finally it defines a `WJetSkimMuon` data sample.

The remaining fields configure various execution switches in the code.

Each sample must be given a list of ROOT files (containing n-tuples created by EmJetAnalyzer) stored in a text file. This is stored in the `inputs` directory. For example, `inputs/QCD_HT2000.txt` could look like this:

    PATH_TO_QCD_HT2000/ntuple_1.root
    PATH_TO_QCD_HT2000/ntuple_10.root
    PATH_TO_QCD_HT2000/ntuple_100.root
    PATH_TO_QCD_HT2000/ntuple_11.root
    PATH_TO_QCD_HT2000/ntuple_12.root
    PATH_TO_QCD_HT2000/ntuple_13.root

**IMPORTANT!** Both the input file path in the config file, and the ROOT file paths in the input file, must be specified in absolute paths.  
**IMPORTANT!** All ROOT files specified in the input file for a given sample are processed to calculate the total generated event count for MC samples. This is used to calculate the normalization.

    
### To run EmJetHistoMaker:

    ./main -c $PWD/configs/config.txt -d <OUTPUTDIR> -s <SAMPLENAME> [-l <LABEL>] -n <NFILES_PER_RUN> <RUN_INDEX>
    
The `-c` argument specifies a full path to the config file.  
The `-d` argument specifies the output directory to store the output ROOT files containing histograms.  
The `-s` argument specifies which sample in the config file to run.  
The `-l` option specifies a label to append to the output files and the output directory.  
The `-n` argument specifies <NFILES_PER_RUN>, i.e. how many files to process per run.  

Finally, a mandatory positional argument `<RUN_INDEX>` must be given to specify which group of `<NFILES_PER_RUN>` to process. For example, `-n 10 0` processes 10 files starting from file index 0, i.e. the first 10 files in the input file. `-n 10 1` processes the next 10 files. And so on.

### Run on condor:

Samples scripts are provided to facilitate running on condor. These scripts can be found in `condor_scripts`.

1. `condor_scripts/condor_test.py <QUEUE> "<COMMAND_TO_SUBMIT>"`

Submits a generic shell command for execution on condor. The `<QUEUE>` argument specifies how many different times to run the command. Each execution is given an argument from `0` to `<QUEUE> - 1`.

Example:

    ./condor_scripts/condor_test.py 100 "./main -c $PWD/configs/config.txt -d ~/data/condor_output/2017-02-10/histos -l v0 -s ModelA -n  1"
    
Submits 100 condor job, each job processing 1 file from `ModelA`. Any extra jobs/files are handled safely, so there is no problem if `ModelA` has less than 100 files.


2. `condor_scripts/haddsamples.py <OUTPUTDIRL>`

Takes an `<OUTPUTDIR>` created by EmJetHistoMaker, merges all the root files in each sub-directory and places them in `<OUTPUTDIR>`.

Example: 

Starting with a directory `OUTPUTDIR=/home/yhshin/data/condor_output/2017-02-10/histos-v0/` that contains the following:

    /home/yhshin/data/condor_output/2017-02-10/histos-v0/ModelA/histo-ModelA-ModelA-v0-0.root
    /home/yhshin/data/condor_output/2017-02-10/histos-v0/ModelA/histo-ModelA-ModelA-v0-1.root
    /home/yhshin/data/condor_output/2017-02-10/histos-v0/ModelA/histo-ModelA-ModelA-v0-2.root
    /home/yhshin/data/condor_output/2017-02-10/histos-v0/ModelA/histo-ModelA-ModelA-v0-3.root
    /home/yhshin/data/condor_output/2017-02-10/histos-v0/ModelA/histo-ModelA-ModelA-v0-4.root
    ...
    /home/yhshin/data/condor_output/2017-02-10/histos-v0/ModelB/histo-ModelB-ModelB-v0-0.root
    /home/yhshin/data/condor_output/2017-02-10/histos-v0/ModelB/histo-ModelB-ModelB-v0-1.root
    /home/yhshin/data/condor_output/2017-02-10/histos-v0/ModelB/histo-ModelB-ModelB-v0-2.root
    /home/yhshin/data/condor_output/2017-02-10/histos-v0/ModelB/histo-ModelB-ModelB-v0-3.root
    /home/yhshin/data/condor_output/2017-02-10/histos-v0/ModelB/histo-ModelB-ModelB-v0-4.root
    
Execute:

    condor_scripts/haddsamples.py ~/data/condor_output/2017-02-10/histos-v0/

This gives you:

    /home/yhshin/data/condor_output/2017-02-10/histos-v0/histo-ModelA.root
    /home/yhshin/data/condor_output/2017-02-10/histos-v0/histo-ModelB.root

