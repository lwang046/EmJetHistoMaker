// Base class for program to loop over single file and fill histograms with arbitrary weights
#include "BaseClass.h"
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TStopwatch.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TMath.h>

#include <cassert>
#include <string>
#include <iostream>
#include <vector>
#include <algorithm>
#include <functional>

#define SET_MEMBER(NAME, TYPE, member) void Set##NAME(TYPE member) {member##_ = member;}
#define SET_MEMBER_DEFAULT(NAME, TYPE, member, DEFAULT) void Set##NAME(TYPE member = DEFAULT) {member##_ = member;}

using std::string;

class HistoMakerBase : protected BaseClass
{
  public:
    HistoMakerBase() {};
    ~HistoMakerBase();
    virtual void OpenOutputFile (string ofilename);
    virtual void InitHistograms() = 0;
    virtual void FillHistograms(long eventnumber) = 0;
    virtual void WriteHistograms();
    virtual int SetTree(string ifilename) = 0;
    SET_MEMBER_DEFAULT(MaxEntries, long, nentries_max, -1);
    long LoopOverCurrentTree ();

  protected:
    int reportEvery_ = 100;
    TFile* ofile_;
    TTree* tree_; // Current tree
    long nentries_; // Number of entries in current tree
    long nentries_max_; // Maximum number of entries to process for current tree (Set to -1 to run over all entries)
    long nentries_to_process_; // Actual number of entries to process for current tree
    TStopwatch timer_total_;
    TStopwatch timer_;
};

// Copied from BaseClass::Loop()
long HistoMakerBase::LoopOverCurrentTree()
{
  if (fChain == 0) {
    std::cout << "Invalid tree!" << std::endl;
    return -1;
  }
  nentries_ = fChain->GetEntriesFast();
  if ( nentries_ == 0 ) {
    std::cout << "No entries!" << std::endl;
    return -1;
  }
  else {
    std::cout << "Number of entries: " << nentries_ << std::endl;
  }

  // Calculate number of entries to process
  if (nentries_max_==-1) { nentries_to_process_ = nentries_                          ; }
  else                   { nentries_to_process_ = std::min(nentries_, nentries_max_) ; }
  std::cout << "Number of entries to process: " << nentries_to_process_ << std::endl;

  Long64_t nbytes = 0, nb = 0;
  timer_total_.Start();

  Int_t fcurrent=-1;
  long eventCount = 0;
  // Loop over all events in TChain
  for (Long64_t jentry = 0; jentry < nentries_to_process_; jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    if (fChain->GetTreeNumber() != fcurrent) {// New root file opened
      fcurrent = fChain->GetTreeNumber();
      TH1F* eventcounthist=(TH1F*)(fChain->GetDirectory()->GetFile()->GetDirectory("eventCountPreTrigger"))->Get("eventCountPreTrigger");
      eventCount += eventcounthist->Integral();
      std::cout << "number of pre-trigger events in current root file is: "<< eventcounthist->Integral()<< std::endl;
    }
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    // Do stuff with one entry of the TTree
    if ( jentry % reportEvery_ == 0 ) {
      if (jentry!=0) std::cout << "Chunk processing time (s): " << timer_.RealTime() << std::endl;
      timer_.Start();
      std::cout << "Running over global entry: " << jentry << std::endl;
    }
    // Fill histograms
    FillHistograms(jentry);
  }
  std::cout << "Total number of processed events is: "<< eventCount << std::endl;
  double total_time_elapsed = timer_total_.RealTime();
  std::cout << "Total processing time (s): " << total_time_elapsed << std::endl;
  std::cout << std::endl;
  return eventCount;
}

void HistoMakerBase::OpenOutputFile(string ofilename)
{
  ofile_ = new TFile(ofilename.c_str(), "RECREATE");
  InitHistograms();
}

void HistoMakerBase::WriteHistograms()
{
  ofile_->Write();
  std::cout << "Wrote to file:\n" << ofile_->GetName() << std::endl;
}

HistoMakerBase::~HistoMakerBase()
{
  delete ofile_;
}
