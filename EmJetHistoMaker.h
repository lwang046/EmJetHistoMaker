#include "HistoMakerBase.h"
#include "EmJetHistos.h"
#include "EmJetSample.h"
#include "LumiReWeightingStandAlone.h"

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TStopwatch.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TMath.h>

#include <cassert>
#include <string>
#include <iostream>
#include <vector>
#include <memory>
#include <algorithm>
#include <functional>

using std::string;
using std::vector;
using std::unique_ptr;

#ifndef OUTPUT
#define OUTPUT(x) std::cout<<#x << ": " << x << std::endl
#endif

const int debug = 0;

class EmJetHistos;
class EmJetSample;
typedef EmJetHistos Histos;
const int nTrackSort=50;
enum class Sample {SIGNAL, QCD, WJET};

class EmJetHistoMaker : public HistoMakerBase
{
 public:
  EmJetHistoMaker();
  EmJetHistoMaker(string ifilename);
  EmJetHistoMaker(EmJetSample isample);
  ~EmJetHistoMaker() {};
  int TRACKSOURCE;
  int VERTEXSOURCE;
  void InitHistograms();
  void FillHistograms    (long eventnumber);
  void FillEventHistograms  (long eventnumber, string tag);
  void FillJetHistograms    (long eventnumber, int ij, string tag);
  void FillTrackHistograms  (long eventnumber, int ij, int itk, string tag);
  void FillPileupHistograms (long eventnumber, string tag);
  void PrintEvent (long eventnumber, string comment);
  int SetTree(string ifilename);
  int SetTree();
  int GetEventCountHistAndClone(string ihistname);
  int GetEventCount(string ihistname);
  void SetOptions(Sample sample=Sample::SIGNAL, bool isData=false, double xsec=1.0, long nevent=1, bool isSignal=false, bool pileupOnly=false);
  void InitLumiReweighting();
 private:
  double CalculateEventWeight(long eventnumber);
  bool SelectJet(int jet_index);
  bool SelectJet_basic(int jet_index);
  bool SelectJet_alphaMax(int jet_index);
  bool SelectJet_ipCut(int jet_index);
  unique_ptr<Histos> histo_;
  unique_ptr<TH1F> histo_nTrueInt_;
  Sample sample_;
  bool isData_;
  double xsec_;
  double nevent_;
  bool isSignal_;
  bool pileupOnly_;
  string file_; // Path to current input file
  unique_ptr<reweight::LumiReWeighting> LumiWeights_;
  // Calculated variables
  double event_ht_;
};

EmJetHistoMaker::EmJetHistoMaker()
{
  SetMaxEntries(-1);
  std::cout << "EmJetHistoMaker::EmJetHistoMaker()" << std::endl;
  TRACKSOURCE = 0;
  InitLumiReweighting();
}

EmJetHistoMaker::EmJetHistoMaker(string ifilename)
{
  SetMaxEntries(-1);
  std::cout << "EmJetHistoMaker::EmJetHistoMaker(\"" << ifilename << "\")" << std::endl;
  file_ = ifilename;
  int status = SetTree(ifilename);
  if (!status==0) std::cerr << "Error when opening file: " << ifilename << std::endl;
  TRACKSOURCE = 0;
  InitLumiReweighting();
}

EmJetHistoMaker::EmJetHistoMaker(EmJetSample isample)
{
  SetMaxEntries(-1);
  std::cout << "EmJetHistoMaker::EmJetHistoMaker(" << isample.name << ")" << std::endl;
  TChain* chain = new TChain("emJetAnalyzer/emJetTree");
  for (string file : isample.files) {
    OUTPUT(file);
    chain->Add(file.c_str(), -1); // -1: the file is connected and the tree header read in memory to get the number of entries.
  }
  tree_ = chain;
  Init(tree_);
}

void EmJetHistoMaker::InitLumiReweighting()
{
  // Initialize lumi reweighting utility
  std::string mcFile   = "~/www/2016-03-23/pileup_mc_2015_25ns_Startup_PoissonOOTPU.root";
  std::string dataFile = "~/www/2016-03-21/pileup-DataSkim-20160302.root";
  std::string mcHist   = "nTrueInt";
  std::string dataHist = "pileup";
  LumiWeights_ = unique_ptr<reweight::LumiReWeighting>(new reweight::LumiReWeighting(mcFile, dataFile, mcHist, dataHist));
}

int EmJetHistoMaker::SetTree(std::string filename)
{
  TFile *f = new TFile(filename.c_str());
  // File read error
  if (f->IsZombie()) {
    std::cout << "File read error for file: " << filename << std::endl;
    return 1;
  }
  TDirectory * dir = (TDirectory*)f->Get((filename+":/emJetAnalyzer").c_str());
  dir->GetObject("emJetTree",tree_);
  // tree_->Print();
  Init(tree_);
  return 0; // No error
}

int EmJetHistoMaker::SetTree()
{
  string filename = file_;
  TFile *f = new TFile(filename.c_str());
  // File read error
  if (f->IsZombie()) {
    std::cout << "File read error for file: " << filename << std::endl;
    return 1;
  }
  TDirectory * dir = (TDirectory*)f->Get((filename+":/emJetAnalyzer").c_str());
  dir->GetObject("emJetTree",tree_);
  // tree_->Print();
  Init(tree_);
  return 0; // No error
}

int EmJetHistoMaker::GetEventCountHistAndClone(string ihistname)
{
  if (isData_) {
    std::cerr << "Error: Attempting to get event count histogram for data\n";
    return -1;
  }
  string filename = file_;
  TFile *f = new TFile(filename.c_str());
  // File read error
  if (f->IsZombie()) {
    std::cout << "File read error for file: " << filename << std::endl;
    return 1;
  }
  TDirectory * dir = (TDirectory*)f->Get((filename+":/"+ihistname).c_str());
  TH1F* eventcounthist = (TH1F*)dir->Get(ihistname.c_str())->Clone();
  eventcounthist->SetDirectory(ofile_);
  eventcounthist->Write();
  return 0; // No error
}

int EmJetHistoMaker::GetEventCount(string ihistname)
{
  if (isData_) {
    std::cerr << "Error: Attempting to get event count histogram for data\n";
    return -1;
  }
  string filename = file_;
  TFile *f = new TFile(filename.c_str());
  // File read error
  if (f->IsZombie()) {
    std::cout << "File read error for file: " << filename << std::endl;
    return -1;
  }
  TDirectory * dir = (TDirectory*)f->Get((filename+":/"+ihistname).c_str());
  TH1F* eventcounthist = (TH1F*)dir->Get(ihistname.c_str());
  if (eventcounthist) return eventcounthist->Integral();
  return -1;
}

void EmJetHistoMaker::InitHistograms()
{
  TH1::SetDefaultSumw2();
  histo_nTrueInt_ = unique_ptr<TH1F>(new TH1F("nTrueInt", "nTrueInt", 100, 0., 100.));
  if (!pileupOnly_) {
    histo_ = unique_ptr<Histos>(new Histos());
  }
}

void EmJetHistoMaker::FillHistograms(long eventnumber)
{
  double w = CalculateEventWeight(eventnumber);
  if (debug==1) std::cout << "Entering FillHistograms" << std::endl;

  FillPileupHistograms(eventnumber, "");
  if (pileupOnly_) return;

  FillEventHistograms(eventnumber, "");
  int nEgamma = 0;
  int nAlphaMax = 0;
  int nEmerging = 0;
  // bool pt_cuts[4];
  // pt_cuts[0] = (*jet_pt)[0] > 400                 ;
  // pt_cuts[1] = (*jet_pt)[1] > 200 && pt_cuts[1-1] ;
  // pt_cuts[2] = (*jet_pt)[2] > 200 && pt_cuts[2-1] ;
  // pt_cuts[3] = (*jet_pt)[3] > 100 && pt_cuts[3-1] ;
  for (unsigned ij = 0; ij < jet_pt->size(); ij++) {
    if (ij>3) break;
    if ( SelectJet_basic(ij) ) nEgamma++;
    if ( SelectJet_basic(ij) && SelectJet_alphaMax(ij) ) nAlphaMax++;
    if ( SelectJet_basic(ij) && SelectJet_alphaMax(ij) && SelectJet_ipCut(ij) ) nEmerging++;
  }
  {
    double ht4 = (*jet_pt)[0] + (*jet_pt)[1] + (*jet_pt)[2] + (*jet_pt)[3];
    const int nCut = 9;
    bool cuts[nCut]; string labels[nCut];
    cuts[0] = true                            ; labels[0]=("nocut");
    cuts[1] = cuts[1-1] && nEgamma>=4         ; labels[1]=("nEgamma>=4        ");
    cuts[2] = cuts[2-1] && ht4 > 1000         ; labels[2]=("ht4 > 1000        ");
    cuts[3] = cuts[3-1] && (*jet_pt)[0] > 400 ; labels[3]=("(*jet_pt)[0] > 400");
    cuts[4] = cuts[4-1] && (*jet_pt)[1] > 200 ; labels[4]=("(*jet_pt)[1] > 200");
    cuts[5] = cuts[5-1] && (*jet_pt)[2] > 200 ; labels[5]=("(*jet_pt)[2] > 200");
    cuts[6] = cuts[6-1] && (*jet_pt)[3] > 100 ; labels[6]=("(*jet_pt)[3] > 100");
    cuts[7] = cuts[7-1] && nAlphaMax < 4      ; labels[7]=("nAlphaMax < 4     ");
    cuts[8] = cuts[8-1] && nEmerging >= 2     ; labels[8]=("nEmerging >= 2    ");
    for (int ic = 0; ic < nCut; ic ++) {
      if ( cuts[ic] ) histo_->hist1d["cutflow2"]->Fill(labels[ic].c_str(), w);
    }
    if (cuts[6]) {
      FillEventHistograms(eventnumber, "__EVTkinematic");
    }
    if (cuts[8]) {
      FillEventHistograms(eventnumber, "__EVTpvpass");
    }
  }
  // Event cut 1
  // FillEventHistograms(eventnumber, "__EVTCUT1");
  // Event cut 2
  // FillEventHistograms(eventnumber, "__EVTCUT2");

  // Fill cut flow histogram for comparison with Sarah's cut flow
  {

    double ht4 = (*jet_pt)[0] + (*jet_pt)[1] + (*jet_pt)[2] + (*jet_pt)[3];
    histo_->hist1d["ht4"]->Fill(ht4, w);
    int nJet = 0;
    for (unsigned ij = 0; ij < jet_pt->size(); ij++) {
      nJet++;
      // if ( TMath::Abs((*jet_eta)[ij]) < 2.0 ) nJet++;
      // if ( (*jet_nhf)[ij] < 0.99 && (*jet_nef)[ij] < 0.99 && (*jet_cef)[ij] < 0.99 ) nJet++;
    }
    int nJetPassingCut = 0;
    for (unsigned ij = 0; ij < jet_pt->size(); ij++) {
      if ( SelectJet(ij) ) nJetPassingCut++;
    }
    // std::cout << event << "\t" << (*jet_pt)[0] << "\t" << (*jet_pt)[1] << "\t" << (*jet_pt)[2] << "\t" << (*jet_pt)[3] << "\t" << ht4 << std::endl;

    const int nCut = 9;
    bool cuts[nCut];
    cuts[0] = true                             ;
    cuts[1] = nJet >= 4                        ;
    cuts[2] = cuts[2-1] && ht4 > 1000          ;
    cuts[3] = cuts[3-1] && (*jet_pt)[0] > 400  ;
    cuts[4] = cuts[4-1] && (*jet_pt)[1] > 200  ;
    cuts[5] = cuts[5-1] && (*jet_pt)[2] > 125  ;
    cuts[6] = cuts[6-1] && (*jet_pt)[3] > 50   ;
    cuts[7] = cuts[7-1] && nJetPassingCut >= 1 ;
    cuts[8] = cuts[8-1] && nJetPassingCut >= 2 ;
    for (int ic = 0; ic < nCut; ic ++) {
      if ( cuts[ic] ) histo_->hist1d["cutflow"]->Fill(ic, w);
    }
  }

}

void EmJetHistoMaker::FillEventHistograms(long eventnumber, string tag)
{
  if (debug==1) std::cout << "Entering FillEventHistograms" << std::endl;

  double w = CalculateEventWeight(eventnumber);

  // Calculate ht
  double ht = 0;
  for (unsigned ij = 0; ij < (*jet_pt).size(); ij++) {
    ht += (*jet_pt)[ij];
  }
  event_ht_ = ht;

  int nJet_egammacut = 0, nJet_alphaMax = 0, nJet_ipcut = 0;
  int nEmerging = 0;
  // Jet loop
  for (unsigned ij = 0; ij < (*jet_pt).size(); ij++) {
    FillJetHistograms(eventnumber, ij, ""+tag);
    if ( (*jet_nDarkPions)[ij] > 0 ) {
      FillJetHistograms(eventnumber, ij, "__sig"+tag);
    }
    if (SelectJet_basic(ij)) {
      FillJetHistograms(eventnumber, ij, "__JTegammacut"+tag);
      nJet_egammacut++;
      if (SelectJet_alphaMax(ij)) {
        FillJetHistograms(eventnumber, ij, "__JTalphaMax"+tag);
        nJet_alphaMax++;
        nEmerging++;
        if (SelectJet_ipCut(ij)) {
          FillJetHistograms(eventnumber, ij, "__JTipcut"+tag);
          nJet_ipcut++;
        }
      }
    }
    // Jet cut 1
    // FillJetHistograms(eventnumber, "__JETCUT1");
    // Jet cut 2
    // FillJetHistograms(eventnumber, "__JETCUT2");
  }

  // Existing quantities (in ntuple)
  histo_->hist1d["nVtx"]->Fill(nVtx, w);

  // Calculated quantities
  histo_->hist1d["ht"]->Fill(ht, w);
  histo_->hist1d[string("jet_N")+"__JTegammacut"+tag]->Fill(nJet_egammacut, w);
  // histo_->hist1d[string("nJet")+"__JTegammacut"+tag]->Fill(nEmerging, w);
  // histo_->hist1d[string("nEmerging")+tag]->Fill(nEmerging, w);
  histo_->hist1d[string("jet_N")+"__JTalphaMax"+tag]->Fill(nJet_alphaMax, w);
  histo_->hist1d[string("jet_N")+"__JTipcut"+tag]->Fill(nJet_ipcut, w);
}

void EmJetHistoMaker::FillJetHistograms(long eventnumber, int ij, string tag)
{
  if (debug==1) std::cout << "Entering FillJetHistograms" << std::endl;
  double w = CalculateEventWeight(eventnumber);

  double ipXYcut = 0.025;
  // Calculate median 2D impact parameter (source=0)
  int nTrack = 0;
  double medianIP = 0.;
  double maxIP = 0.;
  {
    vector<double> vector_ipXY;
    for (unsigned itk=0; itk < (*track_pt)[ij].size(); itk++) {
      if ( (*track_source)[ij][itk] == 0 ) {
        vector_ipXY.push_back( (*track_ipXY)[ij][itk] );
        FillTrackHistograms(eventnumber, ij, itk, tag);
        if ( (*track_ipXY)[ij][itk] < ipXYcut ) {
          FillTrackHistograms(eventnumber, ij, itk, "__TKprompt"+tag);
        }
        else {
          FillTrackHistograms(eventnumber, ij, itk, "__TKdisplaced"+tag);
        }
      }
    }
    std::sort(vector_ipXY.begin(), vector_ipXY.end());
    nTrack = vector_ipXY.size();
    medianIP = 0.;
    if (nTrack>0) {
      if ( nTrack % 2 == 0 ) {
        medianIP = (vector_ipXY[nTrack/2 - 1] + vector_ipXY[nTrack/2]) / 2;
      }
      else {
        medianIP = (vector_ipXY[nTrack/2]);
      }
      maxIP = vector_ipXY[nTrack-1];
    }
  }

  // Calculate median 2D impact parameter (source=0, pt>1)
  int nTrackPostCut = 0;
  double medianIPPostCut = 0.;
  {
    vector<double> vector_ipXY;
    for (unsigned itk=0; itk < (*track_pt)[ij].size(); itk++) {
      if ( (*track_source)[ij][itk] == 0 && (*track_algo)[ij][itk] >= 6 ) {
        vector_ipXY.push_back( (*track_ipXY)[ij][itk] );
        FillTrackHistograms(eventnumber, ij, itk, tag);
      }
    }
    std::sort(vector_ipXY.begin(), vector_ipXY.end());
    int nTrack = vector_ipXY.size();
    double medianIP = 0.;
    if (nTrack>0) {
      if ( nTrack % 2 == 0 ) {
        medianIP = (vector_ipXY[nTrack/2 - 1] + vector_ipXY[nTrack/2]) / 2;
      }
      else {
        medianIP = (vector_ipXY[nTrack/2]);
      }
    }
    nTrackPostCut = nTrack;
    medianIPPostCut = medianIP;
  }

  // Calculate prompt/displaced energy fraction (source=0)
  double prompt_frac = 0;
  double disp_frac = 0;
  {
    double prompt_sum = 0;
    double disp_sum = 0;
    for (unsigned itk=0; itk < (*track_pt)[ij].size(); itk++) {
      if ( (*track_source)[ij][itk] == 0 ) {
        if ( (*track_ipXY)[ij][itk] > 0.1 ) {
          disp_sum += (*track_pt)[ij][itk];
        }
        else {
          prompt_sum += (*track_pt)[ij][itk];
        }
      }
    }
    prompt_frac = prompt_sum / (prompt_sum + disp_sum);
    disp_frac = disp_sum / (prompt_sum + disp_sum);
  }

  // Calculate missing hit track energy fraction (source=0)
  double missInnerHit_frac = 0;
  {
    double missInnerHit_sum = 0;
    double sum = 0;
    for (unsigned itk=0; itk < (*track_pt)[ij].size(); itk++) {
      if ( (*track_source)[ij][itk] == 0 ) {
        if ( (*track_nMissInnerHits)[ij][itk] >= 1 ) {
          missInnerHit_sum += (*track_pt)[ij][itk];
        }
        sum += (*track_pt)[ij][itk];
      }
    }
    missInnerHit_frac = missInnerHit_sum / sum;
  }

  // Existing quantities (in ntuple)
  histo_->hist1d["jet_pt"+tag]->Fill((*jet_pt)[ij], w);
  histo_->hist1d["jet_eta"+tag]->Fill((*jet_eta)[ij], w);
  histo_->hist1d["jet_phi"+tag]->Fill((*jet_phi)[ij], w);
  histo_->hist1d["jet_alphaMax"+tag]->Fill((*jet_alphaMax)[ij], w);
  histo_->hist1d["jet_cef"+tag]->Fill((*jet_cef)[ij], w);
  histo_->hist1d["jet_nef"+tag]->Fill((*jet_nef)[ij], w);
  histo_->hist1d["jet_chf"+tag]->Fill((*jet_chf)[ij], w);
  histo_->hist1d["jet_nhf"+tag]->Fill((*jet_nhf)[ij], w);
  histo_->hist1d["jet_phf"+tag]->Fill((*jet_phf)[ij], w);

  // Calculated quantities
  histo_->hist1d["jet_nTrack"+tag]->Fill(nTrack, w);
  histo_->hist1d["jet_maxIP"+tag]->Fill(maxIP, w);
  histo_->hist1d["jet_medianIP"+tag]->Fill(medianIP, w);
  histo_->hist1d["jet_nTrackPostCut"+tag]->Fill(nTrackPostCut, w);
  histo_->hist1d["jet_medianIPPostCut"+tag]->Fill(medianIPPostCut, w);
  histo_->hist1d["jet_prompt_frac"+tag]->Fill(prompt_frac, w);
  histo_->hist1d["jet_disp_frac"+tag]->Fill(disp_frac, w);
  histo_->hist1d["jet_missInnerHit_frac"+tag]->Fill(missInnerHit_frac, w);

  // 2D histos
  histo_->hist2d[string()+"jet_disp_frac"+"_VS_"+"jet_alphaMax"+tag]->Fill((*jet_alphaMax)[ij], disp_frac, w);
  histo_->hist2d[string()+"jet_disp_frac"+"_VS_"+"jet_pt"+tag]->Fill((*jet_pt)[ij], disp_frac, w);
  histo_->hist2d[string()+"jet_alphaMax"+"_VS_"+"jet_pt"+tag]->Fill((*jet_pt)[ij], (*jet_alphaMax)[ij], w);
  histo_->hist2d[string()+"jet_alphaMax"+"_VS_"+"ht"+tag]->Fill(event_ht_, (*jet_alphaMax)[ij], w);
  histo_->hist2d[string()+"jet_alphaMax"+"_VS_"+"nVtx"+tag]->Fill(nVtx, (*jet_alphaMax)[ij], w);
}

void EmJetHistoMaker::FillTrackHistograms(long eventnumber, int ij, int itk, string tag)
{
  if (debug==1) std::cout << "Entering FillTrackHistograms" << std::endl;
  if (debug==2) OUTPUT(tag);

  double w = CalculateEventWeight(eventnumber);

  // Existing quantities (in ntuple)
  histo_->hist1d["track_pt"+tag]->Fill((*track_pt)[ij][itk], w);
  histo_->hist1d["track_eta"+tag]->Fill((*track_eta)[ij][itk], w);
  histo_->hist1d["track_phi"+tag]->Fill((*track_phi)[ij][itk], w);
  histo_->hist1d["track_quality"+tag]->Fill((*track_quality)[ij][itk], w);
  histo_->hist1d["track_algo"+tag]->Fill((*track_algo)[ij][itk], w);
  histo_->hist1d["track_ipXY"+tag]->Fill((*track_ipXY)[ij][itk], w);
  histo_->hist1d["track_ipXYb"+tag]->Fill((*track_ipXY)[ij][itk], w);
  histo_->hist1d["track_nHits"+tag]->Fill((*track_nHits)[ij][itk], w);
  histo_->hist1d["track_nMissInnerHits"+tag]->Fill((*track_nMissInnerHits)[ij][itk], w);
}

void EmJetHistoMaker::FillPileupHistograms(long eventnumber, string tag)
{
  double weight = CalculateEventWeight(eventnumber);
  histo_nTrueInt_->Fill(nTrueInt, weight);
}

double EmJetHistoMaker::CalculateEventWeight(long eventnumber)
{
  double weight = 0.0;
  if (isData_) weight = 1.0;
  else {
    weight = 1.0;
    double generator_weight = xsec_/nevent_;
    // double generator_weight = xsec_;
    weight *= generator_weight;
    // double pileup_lumi_weight = LumiWeights_->weight(nTrueInt);
    // weight *= pileup_lumi_weight;
  }
  return weight;
}
bool EmJetHistoMaker::SelectJet(int ij)
{
  if (!(ij < 4)) return false;
  int nTrackPassingCut = 0;
  for (unsigned itk = 0; itk < (*track_pt)[ij].size(); itk++) {
    if( (*track_source)[ij][itk] == 0 ) {
      if( (*track_pt)[ij][itk] > 1.0 ) {
        nTrackPassingCut++;
      }
    }
  }

  bool result = true;
  bool cut_nef = (*jet_nef)[ij] < 0.9             ; result = result && cut_nef      ;
  bool cut_cef = (*jet_cef)[ij] < 0.9             ; result = result && cut_cef      ;
  bool cut_alphaMax = (*jet_alphaMax)[ij] < 0.2   ; result = result && cut_alphaMax ;
  bool cut_tracks = nTrackPassingCut > 0          ; result = result && cut_tracks   ;
  return result;
}

bool EmJetHistoMaker::SelectJet_basic(int ij)
{
  // Count number of tracks
  int nTrackPassingCut = 0;
  for (unsigned itk = 0; itk < (*track_pt)[ij].size(); itk++) {
    if( (*track_source)[ij][itk] == 0 ) {
      if( (*track_pt)[ij][itk] > 1.0 ) {
        nTrackPassingCut++;
      }
    }
  }

  bool result = true;
  bool cut_nef = (*jet_nef)[ij] < 0.9             ; result = result && cut_nef      ;
  bool cut_cef = (*jet_cef)[ij] < 0.9             ; result = result && cut_cef      ;
  bool cut_eta = abs((*jet_eta)[ij]) < 2.0        ; result = result && cut_eta      ;
  bool cut_tracks = nTrackPassingCut > 0          ; result = result && cut_tracks   ;
  return result;
}

bool EmJetHistoMaker::SelectJet_alphaMax(int ij)
{
  bool result = true;
  bool cut_alphaMax = (*jet_alphaMax)[ij] < 0.04  ; result = result && cut_alphaMax ;
  return result;
}

bool EmJetHistoMaker::SelectJet_ipCut(int ij)
{
  bool result = true;
  double maxIp = 0.;
  for (unsigned itk = 0; itk < (*track_pt)[ij].size(); itk++) {
    if( (*track_source)[ij][itk] == 0 ) {
      if ( (*track_ipXY)[ij][itk] > maxIp ) maxIp = (*track_ipXY)[ij][itk];
    }
  }
  bool ipCut = maxIp > 0.4;
  result = ipCut;
  return result;
}




void EmJetHistoMaker::PrintEvent (long eventnumber, string comment)
{
  std::cout << "run,lumi,event,comment: " << run << ", " << lumi << ", " << event << ", " << comment << "\n";
}

void
EmJetHistoMaker::SetOptions(Sample sample, bool isData, double xsec, long nevent, bool isSignal, bool pileupOnly)
{
  sample_ = sample;
  isData_ = isData;
  xsec_ = xsec;
  nevent_ = nevent;
  isSignal_ = isSignal;
  pileupOnly_ = pileupOnly;
}
