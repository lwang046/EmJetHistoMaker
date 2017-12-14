#ifndef EMJETCUT_H
#define EMJETCUT_H

#include "TextFileProcessing.h"

#include <string> // string, stof
#include <vector> // vector
#include <iostream> // cout
#include <fstream> // std::ifstream
#include <sstream> // std::stringstream
#include <cassert> // assert

using std::string;
using std::vector;
using std::getline;

class EmJetCut
{
 public:
  vector<float> values; // List of individual cut values
  vector<string> files; // List of files to run over
  string  name;
  float   nEmerging;
  float   ht;
  float   pt0;
  float   pt1;
  float   pt2;
  float   pt3;
  float   met;
  float   pu_dz;
  float   alpha3d_dz;
  float   alpha3d_sig;
  float   alpha3d;
  float   medAbsIp;
  static const int nPar = 12; // Defines total number of parameters (excluding cut name)
};

void PrintCut(EmJetCut cut)
{
  std::cout << "================================================\n";
  std::cout << "Printing cut values\n";
  std::cout << "name , nEmerging , ht , pt0 , pt1 , pt2 , pt3 , met , pu_dz , alpha3d_dz , alpha3d_sig , alpha3d ,  medAbsIp\n";
  std::cout << cut.name << " , ";
  std::cout << cut.values[0];
  for (unsigned i = 0; i < cut.values.size(); i++) {
    std::cout << ", " << cut.values[i];
  }
  std::cout << std::endl;
  std::cout << "================================================\n";
}

void FieldsToCut(const vector<string>& ifields, EmJetCut& cut);

void
ReadCutsFromFile(string cutFileName, vector<EmJetCut>& cuts)
{
  vector<string> lines;
  FileToLines(cutFileName, lines);
  if (lines.empty()) {
    std::cerr << "Cut file empty!\n";
    return;
  }
  for (string line : lines) {
    // std::cout << "line: " << line << std::endl;
    vector<string> fields;
    LineToFields(line, fields);
    if (fields.empty()) {
      // std::cerr << "fields empty\n";
      continue;
    }
    else {
      EmJetCut cut;
      FieldsToCut(fields, cut);
      cuts.push_back(cut);
    }
  }
}


// Turn vector of text fields into Cut object
void
FieldsToCut(const vector<string>& ifields, EmJetCut& cut)
{
  // ifields should have cut.nPar+1 entries, for cut name and one for each cut parameter
  if (ifields.size() != cut.nPar+1) std::cerr << "Number of parameter entries does not match number specified in source code (EmJetCut::nPar)" << std::endl;
  cut.values.clear();
  for (unsigned i=0; i<ifields.size(); i++) {
    string field = ifields[i];
    if (i==0) {
      cut.name = field;
    }
    else {
      cut.values.push_back(std::stof(field));
    }
  }
  assert(cut.values.size() == cut.nPar);
  cut.nEmerging   = cut.values [ 0 ];
  cut.ht          = cut.values [ 1 ];
  cut.pt0         = cut.values [ 2 ];
  cut.pt1         = cut.values [ 3 ];
  cut.pt2         = cut.values [ 4 ];
  cut.pt3         = cut.values [ 5 ];
  cut.met         = cut.values [ 6 ];
  cut.pu_dz       = cut.values [ 7 ];
  cut.alpha3d_dz  = cut.values [ 8 ];
  cut.alpha3d_sig = cut.values [ 9 ];
  cut.alpha3d     = cut.values [10 ];
  cut.medAbsIp    = cut.values [11 ];
}

#endif // EMJETCUT_H
