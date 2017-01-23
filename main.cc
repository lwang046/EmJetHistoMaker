#include "EmJetHistoMaker.h"
#include "EmJetFiles.h"

#include <iostream>
#include <string>
#include <map>
#include <algorithm>
#include <unistd.h>

int main(int argc, char *argv[])
{
  {
    bool pileupOnly = false;
    // Specify any option to turn on pileup only
    if (argc>1) pileupOnly = true;

    std::string prefix = "~/www/2017-01-17/histo-";
    std::string postfix = "-comparison.root";
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
