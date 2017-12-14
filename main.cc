#include "EmJetHistoMaker.h"
#include "EmJetFiles.h"
#include "EmJetCut.h"
// #include "EmJetSample.h"
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

	ValueArg<string> cutConfigArg("t","cutConfig","Cut config file name. A text file specifying all the cut sets separated by commas.",false,"cuts/cuts.txt","string");
	cmd.add( cutConfigArg );

	ValueArg<string> cutArg("u","cut","Cut name. The cut set to run over, among those specified in the cut config file.",false,"","string");
	cmd.add( cutArg );

	ValueArg<int> numberOfFilesArg("n","number","Number of files per run to process.",false,1,"int");
	cmd.add( numberOfFilesArg );

	// MultiArg<int> inputFileIndexArg("i","input","Index of files to run over for specified sample..",false,"int");
	// cmd.add( inputFileIndexArg );

  UnlabeledValueArg<int>  runArg( "run", "Run number. Use in combination with \"-n\" flag to parallelize. E.g. For arguments -n 10 5, the program runs over files indexed 10*5=50 to 10*5+10-1=59.",
                                  true, 0, "int"  );
  cmd.add( runArg );

	// Define switches and add it to the command line.
	SwitchArg pileupOnlySwitch("p","pileup","Only output pileup histograms", false);
	cmd.add( pileupOnlySwitch );

	// Parse the args.
	cmd.parse( argc, argv );

	// Get the value parsed by each arg.
	string iconfig    = configArg.getValue();
	// string ioutput    = outputArg.getValue();
	string ioutputDir = outputDirArg.getValue();
  std::cout << "Output directory set to : " << ioutputDir << std::endl;
	string ilabel     = labelArg.getValue();
	string labelstring("");
  if (ilabel!="") labelstring = string("-") + ilabel; // If label specified, set labelstring
	string isample    = sampleArg.getValue(); // Run over all if ""
	string icutConfig = cutConfigArg.getValue();
	string icut       = cutArg.getValue();
  int inumberOfFiles = numberOfFilesArg.getValue(); // Run over all if -1
  // vector<int> iinputFileIndex = inputFileIndexArg.getValue(); // Run over all if empty
  int irun = runArg.getValue();
  bool pileupOnly = pileupOnlySwitch.getValue();

  // Main body of program
  {
    // Find cut from cut config file
    EmJetCut cutToRun;
    {
      vector<EmJetCut> cuts;
      bool foundCut = false;
      ReadCutsFromFile(icutConfig, cuts);
      for (cut: cuts) {
        if (cut.name == icut) {
          std::cout << icut << " found in cut config\n";
          PrintCut(cut);
          cutToRun = cut;
          foundCut = true;
          break;
        }
      }
      if (!foundCut) {
        std::cerr << "Error! Cut " << icut << " not found\n";
        return 1;
      }
    }

    vector<EmJetSample> ejsamples, ejsamples_singlesample;
    ReadSamplesFromConfigFile(iconfig, ejsamples);
    std::cout << "Number of samples read: " << ejsamples.size() << std::endl;
    // If input sample is specified, only run over the specified sample
    if (isample != "") {
      std::cout << "Sample specified in command line: " << isample << std::endl;
      for (unsigned i=0; i<ejsamples.size(); i++) {
        auto sample=ejsamples[i];
        if (sample.name==isample) {
          std::cout << "Found matching sample in config file " << std::endl;
          ejsamples_singlesample.push_back(sample);
          PrintSample(sample);
        }
      }
      ejsamples = ejsamples_singlesample;
    } // ejsamples now contains all the samples to run over
    for (EmJetSample sample : ejsamples) {
      std::cout << "--------------------------------\n";
      std::cout << "--------------------------------\n";
      std::cout << "Running over sample: " << sample.name << std::endl;
      string sampledir = ioutputDir + labelstring + "/" + sample.name;
      std::cout << "Creating directory: " << sampledir << std::endl;
      string mkdir_command = "mkdir -p " + sampledir;
      system(mkdir_command.c_str());

      // Calculate total number of events in sample
      long eventCount = 0;
      long eventCountPreTrigger = 0;
      // if (!sample.isData) {
      //   std::cout << "--------------------------------\n";
      //   std::cout << "Calculating event count in files" << std::endl;
      //   for (unsigned ifile=0; ifile < sample.files.size(); ifile++) {
      //     // std::cout << "Calculating event count in file: " << ifile << std::endl;
      //     string filename = sample.files[ifile];
      //     string hname = "eventCountPreTrigger";
      //     TFile *f = new TFile(filename.c_str());
      //     if (f->IsZombie()) {
      //       std::cout << "Zombie file: " << ifile << std::endl;
      //       return 1;
      //     }
      //     TDirectory * dir = (TDirectory*)f->Get((filename+":/"+hname).c_str());
      //     if (!dir) {
      //       std::cout << "File with no directory: " << ifile << std::endl;
      //       return 1;
      //     }
      //     else {
      //       TH1F* eventcounthist = (TH1F*)dir->Get(hname.c_str());
      //       eventCount += eventcounthist->Integral();
      //     }
      //     f->Close();
      //     delete f;
      //   }
      // }
      // std::cout << "--------------------------------\n";
      // std::cout << "Total event count in sample: " << eventCount << std::endl;

      EmJetSample sample_filtered = sample; // Copy of sample with filtered files
      sample_filtered.files.clear();
      OUTPUT(irun);
      OUTPUT(inumberOfFiles);
      unsigned firstfileindex = irun * inumberOfFiles;
      unsigned lastfileindex  = irun * inumberOfFiles + inumberOfFiles - 1;
      for (unsigned ifile=0; ifile < sample.files.size(); ifile++) {
        if (ifile >= firstfileindex && ifile <= lastfileindex) {
          std::cout << "Adding file to list of files to be processed: " << ifile << std::endl;
          sample_filtered.files.push_back(sample.files[ifile]);
        }
      }

      // Process files if there are any to be processed
      if (sample_filtered.files.size()) {
        EmJetHistoMaker hm(sample_filtered);
        hm.OpenOutputFile(sampledir+"/histo-"+sample.group+"-"+sample.name+labelstring+"-"+std::to_string(irun)+".root");
        int histstatus = hm.VerifyHistograms();
        hm.SetOptions(Sample::SIGNAL, eventCount, true, pileupOnly);
        hm.SetCut(cutToRun);
        eventCountPreTrigger = hm.LoopOverCurrentTree(); // Returns sum of eventCountPreTrigger histogram integrals
        hm.FillEventCount(eventCountPreTrigger);
        OUTPUT("WriteHistograms");
        hm.WriteHistograms();
      }
    }
    std::cout << "--------------------------------\n";
  }

  } catch (ArgException &e) {
    cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
  }
}
