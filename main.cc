#include "EmJetHistoMaker.h"
#include "EmJetFiles.h"
#include "EmJetSample.h"
#include "tclap/CmdLine.h"

#include <iostream>
#include <string>
#include <map>
#include <algorithm>
#include <unistd.h>

using namespace TCLAP;
using std::cout;
using std::cerr;
using std::endl;

int main(int argc, char *argv[])
{
	try {
	// Define the command line object.
	CmdLine cmd("Run EmJetHistoMaker", ' ', "0.9");

	// Define value arguments and add to the command line.
	ValueArg<string> configArg("c","config","Config file",true,"config.txt","string");
	cmd.add( configArg );

	// ValueArg<string> outputArg("o","output","Output file name, e.g. histo_#-version1.root will create files with names histo_[SAMPLE]-version.root where [SAMPLE] is replaced with the mother sample name",true,"histo_#.root","string");
	// cmd.add( outputArg );

	ValueArg<string> outputDirArg("d","directory","Output directory. Do NOT include trailing slash character",true,".","string");
	cmd.add( outputDirArg );

	ValueArg<string> labelArg("l","label","Additional optional label for output file",false,"","string");
	cmd.add( labelArg );

	ValueArg<string> sampleArg("s","sample","Name of sample to run over. Unspecified: Run over all.",false,"","string");
	cmd.add( sampleArg );

	// ValueArg<int> numberOfFilesArg("n","number","Maximum number of files per sample to run over. -1: Run over all.",false,-1,"int");
	// cmd.add( numberOfFilesArg );

	MultiArg<int> inputFileIndexArg("i","input","Index of files to run over for specified sample. Unspecified: Run over all.",false,"int");
	cmd.add( inputFileIndexArg );

	// Define switches and add it to the command line.
	// SwitchArg pileupOnlySwitch("p","pileup","Only output pileup histograms", false);
	// cmd.add( pileupOnlySwitch );

	// Parse the args.
	cmd.parse( argc, argv );

	// Get the value parsed by each arg.
	string iconfig    = configArg.getValue();
	// string ioutput    = outputArg.getValue();
	string ioutputDir = outputDirArg.getValue();
	string ilabel     = labelArg.getValue();
	string isample    = sampleArg.getValue(); // Run over all if ""
  // int inumberOfFiles = numberOfFilesArg.getValue(); // Run over all if -1
  vector<int> iinputFileIndex = inputFileIndexArg.getValue(); // Run over all if empty
	// bool ipileupOnly = pileupOnlySwitch.getValue();

  // // Calculate parameters from command line arguments
  // string prefix, postfix;
  // std::string::size_type pos = ioutput.find('*');
  // if (pos != std::string::npos) {
  //     prefix = ioutput.substr(0, pos);
  //     postfix = ioutput.substr(pos+1);
  // }
  // else {
  //   cerr << "Output file name must contain * character!\n";
  // }

  // Main body of program
  {
    vector<EmJetSample> ejsamples, ejsamples_singlesample;
    ReadSamplesFromConfigFile(iconfig, ejsamples);
    std::cout << "Number of samples read: " << ejsamples.size() << std::endl;
    PrintSample(ejsamples[0]);
    // If input sample is specified, only run over the specified sample
    if (isample != "") {
      for (EmJetSample sample : ejsamples) {
        if (sample.name==isample) ejsamples_singlesample.push_back(sample);
        ejsamples = ejsamples_singlesample;
      }
    } // ejsamples now contains all the samples to run over
    for (EmJetSample sample : ejsamples) {
      std::cout << "--------------------------------\n";
      std::cout << "--------------------------------\n";
      std::cout << "Running over sample: " << sample.name << std::endl;
      // if (iinputFileIndex.empty()) { // Run over all files if file index unspecified
      //   for ( unsigned ifile=0; ifile < sample.files.size(); ifile++ ) iinputFileIndex.push_back(ifile);
      // }
      for (int ifile : iinputFileIndex) {
        std::cout << "--------------------------------\n";
        std::cout << "Running over file: " << ifile << std::endl;
        EmJetHistoMaker hm;
        // Output file name: OUTPUTDIR/histo-SAMPLEGROUP_SAMPLENAME-LABEL-FILEINDEX.root
        hm.OpenOutputFile(ioutputDir+"/histo-"+sample.group+"_"+sample.name+"-"+std::to_string(ifile)+"-"+ilabel+".root");
        int status = hm.SetTree(sample.files[ifile]);
        hm.SetOptions(Sample::WJET, true, 1., 1., false, false);
        if (status==0) hm.LoopOverCurrentTree();
        else { std::cout << "Error! Skipping file\n"; }
        hm.WriteHistograms();
      }
      std::cout << "--------------------------------\n";
    }
  }

  } catch (ArgException &e) {
    cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
  }
  if (0)
  {
    bool pileupOnly = false;
    // Specify any option to turn on pileup only
    if (argc>1) pileupOnly = true;

    std::string prefix = "~/www/2017-01-23/histo-";
    std::string postfix = "-small.root";
    if (pileupOnly) postfix = "_pileupOnly" + postfix;
    for (std::string sample: samples) {
      EmJetHistoMaker hm;
      hm.VERTEXSOURCE = 1;
      std::cout << "Running over sample: " << sample << std::endl;
      auto sample_files = files.equal_range(sample);
      std::cout << "Number of files: " << files.count(sample) << std::endl;
      hm.OpenOutputFile(prefix+sample+postfix);
      for (auto sample_file_pair = sample_files.first; sample_file_pair != sample_files.second; ++sample_file_pair) {
        auto file = sample_file_pair->second;
        std::string filename = file.name;
        int status = hm.SetTree(filename);
        hm.SetOptions(file.sample, file.isData, file.xsec, file.efficiency, file.isSignal, pileupOnly);
        // hm.SetMaxEntries(50000);
        // hm.SetMaxEntries(100);
        if (status==0) {
          std::cout << "Running over file: " << filename << std::endl;
          hm.LoopOverCurrentTree();
        }
        else { std::cout << "Error! Skipping file: " << filename << std::endl; }
      }
      hm.WriteHistograms();
      // int status = hm.Open(file.name, file.sample, file.isData, file.xsec, file.efficiency, file.isSignal);
      // if (status==0) {
      //   // If successfully initialized
      //   hm.Loop(prefix + label + postfix);
      // }
      // else std::cout << "Error! Skipping file: " << file.name << std::endl;
    }
  }
  if (0)
  {
    bool pileupOnly = false;
    // Specify any option to turn on pileup only
    if (argc>1) pileupOnly = true;

    std::string prefix = "~/www/2016-10-14/histo-";
    std::string postfix = "-alltracks.root";
    if (pileupOnly) postfix = "_pileupOnly" + postfix;
    for (std::string sample: samples) {
      EmJetHistoMaker hm;
      hm.VERTEXSOURCE = 2;
      std::cout << "Running over sample: " << sample << std::endl;
      auto sample_files = files.equal_range(sample);
      std::cout << "Number of files: " << files.count(sample) << std::endl;
      hm.OpenOutputFile(prefix+sample+postfix);
      for (auto sample_file_pair = sample_files.first; sample_file_pair != sample_files.second; ++sample_file_pair) {
        auto file = sample_file_pair->second;
        std::string filename = file.name;
        int status = hm.SetTree(filename);
        hm.SetOptions(file.sample, file.isData, file.xsec, file.efficiency, file.isSignal, pileupOnly);
        // hm.SetMaxEntries(50000);
        if (status==0) {
          std::cout << "Running over file: " << filename << std::endl;
          hm.LoopOverCurrentTree();
        }
        else { std::cout << "Error! Skipping file: " << filename << std::endl; }
      }
      hm.WriteHistograms();
      // int status = hm.Open(file.name, file.sample, file.isData, file.xsec, file.efficiency, file.isSignal);
      // if (status==0) {
      //   // If successfully initialized
      //   hm.Loop(prefix + label + postfix);
      // }
      // else std::cout << "Error! Skipping file: " << file.name << std::endl;
    }
  }
}
