#ifndef TEXTFILEPROCESSING_H
#define TEXTFILEPROCESSING_H

#include <string> // string, stof
#include <vector> // vector
#include <iostream> // cout
#include <fstream> // std::ifstream
#include <sstream> // std::stringstream

using std::string;
using std::vector;
using std::getline;

void FileToLines(string ifileName, vector<string>& result);
void LineToFields(string iline, vector<string>& result);
void RemoveSpaces(string& input);

// Convert input file to vector of lines
void
FileToLines(string ifileName, vector<string>& result)
{
  std::ifstream lineStream(ifileName);
  string line;
  while(getline(lineStream,line)) {
    result.push_back(line);
  }
  return;
}

// Convert input line from a CSV file to a vector of fields
void
LineToFields(string iline, vector<string>& result)
{
  if (iline.size()==0) {
    return; // Ignore empty lines
  }
  if (iline[0]=='#') {
    // std::cout << "Ignoring comment line" << std::endl;
    return; // Ignore lines that start with '#', return empty vector
  }

  // std::cout << "iline: " << iline << std::endl;
  RemoveSpaces(iline); // Remove spaces
  // std::cout << "iline after space removal: " << iline << std::endl;
  std::stringstream iss(iline);
  string field;
  // Read result separated by comma
  while (getline(iss, field, ',')) {
    result.push_back(field);
  }
  // This checks for a trailing comma with no data after it.
  if (!iss && field.empty()) {
    result.push_back(""); // If there was a trailing comma then add an empty element.
  }
}

void RemoveSpaces(string& input)
{
  input.erase(std::remove(input.begin(),input.end(),' '),input.end());
}

#endif // TEXTFILEPROCESSING_H
