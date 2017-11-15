#include "TextFileProcessing.h"

#include <string> // string, stof
#include <vector> // vector
#include <iostream> // cout
#include <fstream> // std::ifstream
#include <sstream> // std::stringstream

using std::string;
using std::vector;
using std::getline;

class EmJetSample
{
 public:
  string group; // Set to the same value for samples that should be combined
  string name; // Unique name
  string type; // Sample type - Used to label samples used for different purposes
  bool isData; // true for data
  bool isSignal; // Set to true for signal MC. Turns on signal histograms
  float xsec; // Cross section in pb
  vector<string> files; // List of files to run over
};

void PrintSample(EmJetSample sample)
{
  std::cout << sample.group << ", " << sample.name << ", " << sample.type << ", " << sample.isData << ", " << sample.isSignal << ", " << sample.xsec << std::endl;
  std::cout << "First 10 files:\n";
  for (unsigned i = 0; i < sample.files.size() && i < 10; i++) {
    std::cout << sample.files[i] << std::endl;
  }
}

void ReadSamplesFromConfigFile(string configFileName, vector<EmJetSample>& samples);
void FieldsToSample(const vector<string>& ifields, EmJetSample& sample);


void
ReadSamplesFromConfigFile(string configFileName, vector<EmJetSample>& samples)
{
  std::cout << "Reading samples from config file: " << configFileName << std::endl;
  vector<string> lines;
  FileToLines(configFileName, lines);
  if (lines.empty()) {
    std::cerr << "Config file empty!\n";
    return;
  }
  // Loop over sample configs
  for (string line : lines) {
    // std::cout << "line: " << line << std::endl;
    vector<string> fields;
    LineToFields(line, fields);
    if (fields.empty()) {
      // std::cerr << "fields empty\n";
      continue;
    }
    else {
      EmJetSample sample;
      FieldsToSample(fields, sample);
      samples.push_back(sample);
      // PrintSample(sample);
    }
  }
}


// Turn vector of text fields into Sample object
void
FieldsToSample(const vector<string>& ifields, EmJetSample& sample)
{
  sample.group = ifields[0];
  sample.name = ifields[1];
  // OUTPUT(sample.name);
  sample.type = ifields[2];

  if      (ifields[3]=="MC"                       ) sample.isData = false;
  else if (ifields[3]=="DATA" or ifields[3]=="Data" ) sample.isData = true;
  else std::cerr << "Sample MC/DATA not specified correctly in config file for sample: " << sample.name << std::endl;

  if      (ifields[4]=="false" ) sample.isSignal = false;
  else if (ifields[4]=="true"  ) sample.isSignal = true;
  else std::cerr << "Sample isSignal not specified correctly in config file for sample: " << sample.name << std::endl;

  sample.xsec = stof(ifields[5]);

  bool found_txt = ( ifields[6].find(".txt") != string::npos );
  bool found_root = ( ifields[6].find(".root") != string::npos );
  if (found_root) {
    std::cout << ifields[6] << std::endl;
    sample.files.push_back(ifields[6]);
  }
  else if (found_txt) {
    vector<string> lines; FileToLines(ifields[6], lines);
    for (string line : lines) {
      sample.files.push_back(line);
    }
  }
  else {
    std::cerr << "Invalid input file" << std::endl;
  }
}

