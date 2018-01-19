#include "HistoMakerBase.h"
#include "EmJetHistos.h"
#include "EmJetSample.h"
#include "EmJetSystematics.h"
#include "EmJetCut.h"
#include "LumiReWeightingStandAlone.h"
#include "TriggerReWeighting.h"

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TStopwatch.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TMath.h>
#include <TParameter.h>

#include <cassert>
#include <string>
#include <iostream>
#include <vector>
#include <memory>
#include <algorithm>
#include <functional>

#include "LHAPDF/LHAPDF.h"

using std::string;
using std::vector;
using std::unique_ptr;

#ifndef OUTPUT
#define OUTPUT(x) std::cout<<#x << ": " << x << std::endl
#endif

const int debug = 5;

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
  void FillEventCount (long eventCountPreTrigger);
  void FillHistograms    (long eventnumber);
  void FillCutflowHistograms (long eventnumber, string tag);
  void FillSystematicHistograms (long eventnumber);
  void FillSystematicTestingHistograms (long eventnumber);
  void FillEventHistograms  (long eventnumber, string tag);
  void FillJetHistograms    (long eventnumber, int ij, string tag);
  void FillTrackHistograms  (long eventnumber, int ij, int itk, string tag);
  void FillPileupHistograms (long eventnumber, string tag);
  void FillHltHistograms    (long eventnumber);
  void PrintEvent (long eventnumber, string comment);
  int SetTree(string ifilename);
  int SetTree();
  int SetCut(EmJetCut cut);
  int GetEventCountHistAndClone(string ihistname);
  int GetEventCount(string ihistname);
  void SetOptions(Sample sample=Sample::SIGNAL, long nevent=1, bool isSignal=false, bool pileupOnly=false);
  int VerifyHistograms();
  void InitReweighting();
  void InitPdfSet(string setname);
  // Computation functions
  float DeltaR(float eta1, float phi1, float eta2, float phi2);
  int FlavourTagging(int ij);
  double GetPt(int ij); // Get JEC-shifted pt for given jet (Direction set by sys_.m_direction_jec)
  double GetAlpha(int ij); // Calculate alpha for given jet
  double GetAlpha2DSig(int ij);
  double GetAlpha3DSigM(int ij);
  double GetLTKFrac(int ij);
  double GetNonPUFrac(int ij);
  double GetFrac2DSig(int ij);
  double GetMedianIP(int ij);
  double GetfabsMedianIP(int ij);
  int GetNTrack(int ij);
  bool GetSignalPartonIndex(long eventnumber);
  bool PartonJetMatching(int ij, int igp);
 private:
  double CalculateEventWeight(long eventnumber);
  double CalculatePileupWeight(long eventnumber);
  double CalculatePdfShift(long eventnumber); // direction: 0 = no shift, +1 = shifted up, -1 = shifted down
  double GetPVTrackFraction(long eventnumber);
  bool CheckPV(long eventnumber); // Check for PV reconstruction failure
  bool SelectJet(int jet_index);
  bool SelectJet_basic(int jet_index);
  bool SelectJet_nonpu(int jet_index);
  bool SelectJet_alphaMax(int jet_index);
  bool SelectJet_emerging(int jet_index);
  bool SelectJet_ipCut(int jet_index);
  bool SelectTrack(int jet_index, int track_index);
  // unique_ptr<TTree> ntree_;
  unique_ptr<Histos> histo_;
  unique_ptr<TH1F> histo_nTrueInt_;
  unique_ptr<TH1F> histo_nTrueInt_puweighted_;
  unique_ptr<TH1F> histo_eventCountPreTrigger_;
  EmJetSystematics sys_; // Class used for systematic uncertainty calculations
  EmJetCut cut_;
  Sample sample_;
  bool isData_;
  double xsec_;
  double nevent_;
  bool isSignal_;
  bool pileupOnly_;
  string file_; // Path to current input file
  unique_ptr<reweight::LumiReWeighting> LumiWeights_;
  unique_ptr<TriggerReWeighting> TriggerWeights_;
  vector<LHAPDF::PDF*> PdfMembers_;
  // Calculated variables
  double event_ht_;
  // Quark index
  int dkqk_index; // Dark quark
  int dnqk_index; // Down quark
  int adkqk_index; // Anti dark quark
  int adnqk_index; // Anti down quark
};

EmJetHistoMaker::EmJetHistoMaker()
{
  SetMaxEntries(-1);
  std::cout << "EmJetHistoMaker::EmJetHistoMaker()" << std::endl;
  TRACKSOURCE = 0;
  InitReweighting();
}

EmJetHistoMaker::EmJetHistoMaker(string ifilename)
{
  SetMaxEntries(-1);
  std::cout << "EmJetHistoMaker::EmJetHistoMaker(\"" << ifilename << "\")" << std::endl;
  file_ = ifilename;
  int status = SetTree(ifilename);
  if (!status==0) std::cerr << "Error when opening file: " << ifilename << std::endl;
  TRACKSOURCE = 0;
  InitReweighting();
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
  // OUTPUT(chain->GetCurrentFile()->IsZombie());
  tree_ = chain;
  Init(tree_);
  TRACKSOURCE = 0;
  InitReweighting();

  // Get options from sample
  {
    isData_ = isample.isData;
    xsec_ = isample.xsec;
  }

  if (!isData_) {
    InitPdfSet("NNPDF23_lo_as_0119_qed");
    // InitPdfSet("NNPDF30_lo_as_0118");
    // InitPdfSet("CT14nlo");
    // InitPdfSet("PDF4LHC15_nlo_mc");
    // InitPdfSet("MMHT2014lo68cl");
  }
}

void EmJetHistoMaker::InitReweighting()
{
  // Initialize lumi reweighting utility
  // std::string mcFile   = "/home/yhshin/data/condor_output/2017-11-08/histos-sys1/histo-QCD.root";
  // std::string mcFile   = "/data/users/yhshin/datajsonfiles_received/Analysis-20171103-v0/pileup_Analysis-20171103-v0_20171109_BF.root";
  // std::string dataFile = "/data/users/yhshin/datajsonfiles_received/Analysis-20171103-v0/pileup_Analysis-20171103-v0_20171109_GH.root";
  std::string mcFile   = "mcPileupPdf-2017-11-16.root";
  std::string dataFile = "/data/users/yhshin/datajsonfiles_received/Analysis-20171103-v0/pileup_Analysis-20171103-v0_20171109_GH.root";
  std::string mcHist   = "nTrueInt";
  std::string dataHist = "pileup";
  LumiWeights_ = unique_ptr<reweight::LumiReWeighting>(new reweight::LumiReWeighting(mcFile, dataFile, mcHist, dataHist));
  const bool printLumiWeights = false;
  if (printLumiWeights) {
    std::cout << "nTrueInt, weight\n";
    for (int i=0; i<=100; i++) {
      double w = LumiWeights_->weight(i);
      std::cout << i << ", " << w << "\n";
    }
  }

  // Initialize trigger reweighting utility
  TriggerWeights_ = unique_ptr<TriggerReWeighting>(new TriggerReWeighting( "emjet_hlt_efficiency.root", "hlt_efficiency_mc", "hlt_efficiency_data"));
  const bool printTriggerWeights = false;
  if (printTriggerWeights) {
    std::cout << "ht, weight\n";
    for (int i=0; i<=100; i++) {
      double ht = i * 25.;
      double w = TriggerWeights_->weight(ht);
      std::cout << ht << ", " << w << "\n";
    }
  }
}

void EmJetHistoMaker::InitPdfSet(string setname)
{
  PdfMembers_ = LHAPDF::mkPDFs(setname);
}

float EmJetHistoMaker::DeltaR(float eta1, float phi1, float eta2, float phi2)
{
  float dR=0.;
  float deta = std::fabs(eta1-eta2);
  float dphi = std::fabs(phi1-phi2);
  if( dphi>3.1415926 ) dphi = 2.*3.1415926-dphi;
  dR=std::sqrt(deta*deta+dphi*dphi);
  return dR;
}

int EmJetHistoMaker::FlavourTagging(int ij)
{
  double maxpT = 0.0;
  int flavour = 10;
  for(unsigned int igp=0; igp<(*gp_pt).size()-1; igp++){
    if( abs((*gp_pdgId)[igp])>5 && abs((*gp_pdgId)[igp])!=21 ) continue;
    if( DeltaR((*jet_eta)[ij], (*jet_phi)[ij], (*gp_eta)[igp], (*gp_phi)[igp])>0.4 ) continue;
    if( (*gp_pt)[igp]> maxpT ){
      maxpT = (*gp_pt)[igp];
      flavour = (*gp_pdgId)[igp];
    }
  }
  if( abs(flavour)==21 ){
    double maxpT2 = 0.0;
    for(unsigned int igp=0; igp<(*gp_pt).size()-1; igp++){
      if(abs((*gp_pdgId)[igp])!=5 || fabs((*gp_pt)[igp])<10.0 ) continue;
      if( DeltaR((*jet_eta)[ij], (*jet_phi)[ij], (*gp_eta)[igp], (*gp_phi)[igp])>0.4 ) continue;
      if( (*gp_pt)[igp] > maxpT2 ){
        maxpT2 = (*gp_pt)[igp];
        flavour = 19;
      }
    }
  }
  return abs(flavour);
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

int EmJetHistoMaker::SetCut(EmJetCut cut)
{
  cut_ = cut;
  sys_.SetCut(cut);
  // Initialize Systematic caculations
  // sys_.SetTesting(1);
  // sys_.SetFilenames("/home/yhshin/data/condor_output/2017-11-08/histos-sys1/histo-JetHT_G1.root", "/home/yhshin/data/condor_output/2017-11-08/histos-sys1/histo-QCD.root");
  sys_.SetModelingFilenames("/home/yhshin/data/condor_output/2017-12-08/histos-sys0/histo-DataGH.root", "/home/yhshin/data/condor_output/2017-12-08/histos-sys1/histo-QCD.root");
  TNamed cutname("cutname", cut_.name);
  cutname.Write();
  if (debug>=5) OUTPUT(sys_.CalculateIpXYShift());
  if (debug>=5) OUTPUT(sys_.Calculate3dSigShift());
  TParameter<double> sys_IpXYShift("sys_IpXyShift", sys_.CalculateIpXYShift(), 'f');
  sys_IpXYShift.Write();
  TParameter<double> sys_3dSigShift("sys_3dSigShift", sys_.Calculate3dSigShift(), 'f');
  sys_3dSigShift.Write();
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
  histo_nTrueInt_            = unique_ptr<TH1F>(new TH1F("nTrueInt"            , "nTrueInt"            , 100 , 0. , 100.));
  histo_nTrueInt_puweighted_ = unique_ptr<TH1F>(new TH1F("nTrueInt_puweighted" , "nTrueInt_puweighted" , 100 , 0. , 100.));
  if (!pileupOnly_) {
    histo_ = unique_ptr<Histos>(new Histos());
  }
}

void EmJetHistoMaker::FillEventCount(long eventCountPreTrigger)
{
  // ntree_ = unique_ptr<TTree>(new TTree("ntree_", "to store list info"));
  //auto b = new TParameter<long> ("eventCountPreTrigger", eventCountPreTrigger, '+');
  // std::cout<< "Total Number of Entries "<< eventCountTotal <<std::endl;
  std::cout<< "Number of Processed Entries "<< eventCountPreTrigger <<std::endl;
  // ntree_->GetUserInfo()->AddLast( a );
  //ntree_->GetUserInfo()->AddLast( b );
  histo_eventCountPreTrigger_ = unique_ptr<TH1F>(new TH1F("eventCountPreTrigger", "eventCountPreTrigger", 2, 0., 2.));
  // for(int i=0; i< eventCountPreTrigger; i++){
  histo_eventCountPreTrigger_->Fill(1., eventCountPreTrigger);
  // }
}

void EmJetHistoMaker::FillHistograms(long eventnumber)
{
  double w = CalculateEventWeight(eventnumber);
  if (debug==1) std::cout << "Entering FillHistograms" << std::endl;
  if (debug==3) std::cout << "Eventweight is: " << w << std::endl;

  FillPileupHistograms(eventnumber, "");
  if (pileupOnly_) return;

  FillHltHistograms(eventnumber);

  // FillEventHistograms(eventnumber, "");

  FillSystematicHistograms(eventnumber);
  FillSystematicTestingHistograms(eventnumber);

  string tag = ""; sys_.SetDirectionJec(0);
  {
    sys_.SetDirectionModeling(0);
    FillCutflowHistograms(eventnumber, ""+tag);
    sys_.SetDirectionModeling(1);
    FillCutflowHistograms(eventnumber, "__ModelingUp"+tag);
    sys_.SetDirectionModeling(0);
  }
  tag = "__PtUp"; sys_.SetDirectionJec(1);
  {
    sys_.SetDirectionModeling(0);
    FillCutflowHistograms(eventnumber, ""+tag);
    sys_.SetDirectionModeling(1);
    FillCutflowHistograms(eventnumber, "__ModelingUp"+tag);
    sys_.SetDirectionModeling(0);
  }
  tag = "__PtDn"; sys_.SetDirectionJec(-1);
  {
    sys_.SetDirectionModeling(0);
    FillCutflowHistograms(eventnumber, ""+tag);
    sys_.SetDirectionModeling(1);
    FillCutflowHistograms(eventnumber, "__ModelingUp"+tag);
    sys_.SetDirectionModeling(0);
  }
}

void EmJetHistoMaker::FillCutflowHistograms (long eventnumber, string tag)
{
  double w = CalculateEventWeight(eventnumber);
  if (debug==1) std::cout << "Entering FillCutflowHistograms" << std::endl;
  // if (debug>=5) std::cout << "Entering FillCutflowHistograms for tag: " << tag << std::endl;

  double pdfshift = 0.;
  if(!isData_) {
    pdfshift = CalculatePdfShift(eventnumber);
    if (debug==6) std::cout << "PdfShift is: " << pdfshift << std::endl;
    histo_->hist1d["pdfshift"]->Fill(pdfshift, w);
  }

  int nBasic = 0;
  int nAlphaMax = 0;
  int nEmerging = 0;
  for (unsigned ij = 0; ij < jet_pt->size(); ij++) {
    if (ij>3) break;
    if ( SelectJet_basic(ij) ) nBasic++;
    if ( SelectJet_basic(ij) && SelectJet_alphaMax(ij) ) nAlphaMax++;
    if ( SelectJet_basic(ij) && SelectJet_alphaMax(ij) && SelectJet_ipCut(ij) ) nEmerging++;
  }
  double ht4 = GetPt(0) + GetPt(1) + GetPt(2) + GetPt(3);
  // Calculate ht
  double ht = 0;
  for (unsigned ij = 0; ij < (*jet_pt).size(); ij++) {
    ht += GetPt(ij);
  }
  // Find basic jets
  vector<int> ijs_basic;
  for (unsigned ij = 0; ij < jet_pt->size(); ij++) {
    if ( SelectJet_basic(ij) ) {
      ijs_basic.push_back(ij);
    }
    if ( ijs_basic.size() == 4 ) break;
  }
  // Find emerging jets
  vector<int> ijs_emerging;
  for (unsigned ij = 0; ij < std::min( int(jet_pt->size()), 4 ); ij++) {
    if ( SelectJet_basic(ij) && SelectJet_emerging(ij) ) {
      ijs_emerging.push_back(ij);
    }
  }
  const int nCut = 16;
  bool cuts[nCut]; string labels[nCut];
  int i=0;
  cuts[0] = true                                               ; labels[0]=( "nocut                      ") ; i++ ;
  cuts[i] = cuts[i-1] && HLT_PFHT900 == 1                      ; labels[i]=( "HLT_PFHT900 == 1           ") ; i++ ;
  cuts[i] = cuts[i-1] && jet_pt->size()>=4                     ; labels[i]=( "jet_pt->size()>=4          ") ; i++ ;
  cuts[i] = cuts[i-1] && ht4 > cut_.ht                         ; labels[i]=( "ht4 > ht_min               ") ; i++ ;
  cuts[i] = cuts[i-1] && met_pt > cut_.met                     ; labels[i]=( "met > met_min              ") ; i++ ;
  cuts[i] = cuts[i-1] && (*pv_index)[0] == 0                   ; labels[i]=( "(*pv_index)[0] == 0        ") ; i++ ;
  cuts[i] = cuts[i-1] && abs((*pv_z)[0]) < 15.0                ; labels[i]=( "abs((*pv_z)[0]) < 15.0     ") ; i++ ;
  cuts[i] = cuts[i-1] && ijs_basic.size() == 4                 ; labels[i]=( "ijs_basic.size() == 4      ") ; i++ ;
  cuts[i] = cuts[i-1] && ijs_basic.back() == 3                 ; labels[i]=( "ijs_basic.back() == 3      ") ; i++ ;
  cuts[i] = cuts[i-1] && GetPt(0) > cut_.pt0                   ; labels[i]=( "GetPt(0) > pt0_min         ") ; i++ ;
  cuts[i] = cuts[i-1] && GetPt(1) > cut_.pt1                   ; labels[i]=( "GetPt(1) > pt1_min         ") ; i++ ;
  cuts[i] = cuts[i-1] && GetPt(2) > cut_.pt2                   ; labels[i]=( "GetPt(2) > pt2_min         ") ; i++ ;
  cuts[i] = cuts[i-1] && GetPt(3) > cut_.pt3                   ; labels[i]=( "GetPt(3) > pt3_min         ") ; i++ ;
  cuts[i] = cuts[i-1] && ijs_emerging.size() >= 1              ; labels[i]=( "nEmerging >= 1             ") ; i++ ;
  cuts[i] = cuts[i-1] && ijs_emerging.size() >= cut_.nEmerging ; labels[i]=( "nEmerging >= nEmerging_min ") ; i++ ;
  cuts[i] = cuts[i-1] && GetPVTrackFraction(eventnumber) > 0.1 ; labels[i]=( "PVTrackFraction > 0.1      ") ; i++ ;
  // cuts[i] = cuts[i-1] && CheckPV(eventnumber)                  ; labels[i]=( "CheckPV(eventnumber)       ") ; i++ ;
  float triggerWeight = TriggerWeights_->weight(ht);
  for (int ic = 0; ic < nCut; ic ++) {
    if (isData_) {
      if ( cuts[ic] ) histo_->hist1d["cutflow"+tag]->Fill(labels[ic].c_str(), 1);
    }
    else {
      if ( cuts[ic] ) histo_->hist1d["cutflow"+tag]->Fill(labels[ic].c_str(), w);
      if ( cuts[ic] ) histo_->hist1d["cutflow__PdfUp"+tag]->Fill(labels[ic].c_str(), w*(1+pdfshift));
      if ( cuts[ic] ) histo_->hist1d["cutflow__PdfDn"+tag]->Fill(labels[ic].c_str(), w*(1-pdfshift));
      if ( cuts[ic] ) histo_->hist1d["cutflow__TriggerUp"+tag]->Fill(labels[ic].c_str(), w*(1+triggerWeight));
    }
  }
  if (cuts[12]) {
    for (unsigned ij = 0; ij < (*jet_pt).size(); ij++) {
      if (ij>=4) break; // :JETCUT:
      if (SelectJet_basic(ij)) {
        histo_->hist1d["test_jet_medianIP"+tag]->Fill(GetMedianIP(ij) , w);
        histo_->hist1d["test_jet_medianAbsIP"+tag]->Fill(GetfabsMedianIP(ij) , w);
      }
    }
  }
  if (cuts[6]) {
    // FillEventHistograms(eventnumber, "__EVTkinematic");
  }
  if (cuts[7]) {
    // FillEventHistograms(eventnumber, "__EVTpvpass");
  }
  if (cuts[nCut-1]) {
    // FillEventHistograms(eventnumber, "__EVTallpass");
  }
}

void EmJetHistoMaker::FillSystematicHistograms (long eventnumber)
{
  double w = CalculateEventWeight(eventnumber);
  if (debug==1) std::cout << "Entering FillSystematicHistograms" << std::endl;

  for (unsigned ij = 0; ij < (*jet_pt).size(); ij++) {
    for (unsigned itk=0; itk < (*track_pt)[ij].size(); itk++) {
      if ( !SelectTrack(ij, itk) ) continue;
      // OUTPUT((*track_ipXY)[ij][itk]);
      // histo_->hist1d_double["sys_track_ipXY"]->Fill( abs((*track_ipXY)[ij][itk]) , w);
      histo_->hist1d_double["sys_log_track_ipXY"]->Fill( TMath::Log10(TMath::Abs((*track_ipXY)[ij][itk])) , w);
      histo_->hist1d_double["sys_log_track_ipXYSig"]->Fill( TMath::Log10(TMath::Abs((*track_ipXYSig)[ij][itk])) , w);
      double tk3dsig2 = TMath::Power(((*pv_z)[0]-(*track_ref_z)[ij][itk])/0.01, 2.0) + TMath::Power((*track_ipXYSig)[ij][itk], 2.0);
      double tk3dsig  = TMath::Sqrt(tk3dsig2);
      histo_->hist1d_double["sys_track_3dSig"]->Fill(tk3dsig);
    }
  }
}

void EmJetHistoMaker::FillSystematicTestingHistograms (long eventnumber)
{
  double w = CalculateEventWeight(eventnumber);
  sys_.SetDirectionModeling(1);
  for (unsigned ij = 0; ij < (*jet_pt).size(); ij++) {
    for (unsigned itk=0; itk < (*track_pt)[ij].size(); itk++) {
      if ( !SelectTrack(ij, itk) ) continue;
      double shifted_tk_ipXY = sys_.GetShiftedIpXY( (*track_ipXY)[ij][itk] ) ;
      // OUTPUT( (*track_ipXY)[ij][itk] );
      // OUTPUT(shifted_tk_ipXY);
      histo_->hist1d["systest_log_track_ipXY"]->Fill( TMath::Log10(TMath::Abs(shifted_tk_ipXY) ) , w);
      double tk3dsig2 = TMath::Power(((*pv_z)[0]-(*track_ref_z)[ij][itk])/0.01, 2.0) + TMath::Power((*track_ipXYSig)[ij][itk], 2.0);
      double tk3dsig  = TMath::Sqrt(tk3dsig2);
      double shifted_tk3dsig = sys_.GetShifted3dSig(tk3dsig);
      histo_->hist1d_double["systest_track_3dSig"]->Fill(shifted_tk3dsig);
    }
  }
  sys_.SetDirectionModeling(0);
}

void EmJetHistoMaker::FillEventHistograms(long eventnumber, string tag)
{
  if (debug==1) std::cout << "Entering FillEventHistograms" << std::endl;

  double w = CalculateEventWeight(eventnumber);
  if (!isData_) {
    double pdfshift = CalculatePdfShift(eventnumber);
  }

  // Calculate ht
  double ht = 0;
  for (unsigned ij = 0; ij < (*jet_pt).size(); ij++) {
    ht += (*jet_pt)[ij];
  }
  event_ht_ = ht;

  int nJet_basic = 0, nJet_alphaMax = 0, nJet_ipcut = 0;
  int nEmerging = 0;
  // Jet loop
  for (unsigned ij = 0; ij < (*jet_pt).size(); ij++) {
    if (ij>=4) break; // :JETCUT:
    FillJetHistograms(eventnumber, ij, ""+tag);
    if ( (*jet_nDarkPions)[ij] > 0 ) {
      FillJetHistograms(eventnumber, ij, "__sig"+tag);
    }
    if (SelectJet_basic(ij)) {
      FillJetHistograms(eventnumber, ij, "__JTbasic"+tag);
      nJet_basic++;
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
  histo_->hist1d[string("jet_N")+"__JTbasic"+tag]->Fill(nJet_basic, w);
  // histo_->hist1d[string("nJet")+"__JTbasic"+tag]->Fill(nEmerging, w);
  // histo_->hist1d[string("nEmerging")+tag]->Fill(nEmerging, w);
  histo_->hist1d[string("jet_N")+"__JTalphaMax"+tag]->Fill(nJet_alphaMax, w);
  histo_->hist1d[string("jet_N")+"__JTipcut"+tag]->Fill(nJet_ipcut, w);
}

void EmJetHistoMaker::FillJetHistograms(long eventnumber, int ij, string tag)
{
  if (debug==1) std::cout << "Entering FillJetHistograms" << std::endl;
  double w = CalculateEventWeight(eventnumber);
  // OUTPUT( (*jet_alpha2)[ij] );
  // OUTPUT( GetAlpha(ij) );

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
  histo_->hist1d["jet_alphaMax_dz1um"   +tag]->Fill((*jet_alphaMax_dz1um)[ij], w);
  histo_->hist1d["jet_alphaMax_dz2um"   +tag]->Fill((*jet_alphaMax_dz2um)[ij], w);
  histo_->hist1d["jet_alphaMax_dz5um"   +tag]->Fill((*jet_alphaMax_dz5um)[ij], w);
  histo_->hist1d["jet_alphaMax_dz10um"  +tag]->Fill((*jet_alphaMax_dz10um)[ij], w);
  histo_->hist1d["jet_alphaMax_dz20um"  +tag]->Fill((*jet_alphaMax_dz20um)[ij], w);
  histo_->hist1d["jet_alphaMax_dz50um"  +tag]->Fill((*jet_alphaMax_dz50um)[ij], w);
  histo_->hist1d["jet_alphaMax_dz100um" +tag]->Fill((*jet_alphaMax_dz100um)[ij], w);
  histo_->hist1d["jet_alphaMax_dz200um" +tag]->Fill((*jet_alphaMax_dz200um)[ij], w);
  histo_->hist1d["jet_alphaMax_dz500um" +tag]->Fill((*jet_alphaMax_dz500um)[ij], w);
  histo_->hist1d["jet_alphaMax_dz1mm"   +tag]->Fill((*jet_alphaMax_dz1mm)[ij], w);
  histo_->hist1d["jet_alphaMax_dz2mm"   +tag]->Fill((*jet_alphaMax_dz2mm)[ij], w);
  histo_->hist1d["jet_alphaMax_dz5mm"   +tag]->Fill((*jet_alphaMax_dz5mm)[ij], w);
  histo_->hist1d["jet_alphaMax_dz1cm"   +tag]->Fill((*jet_alphaMax_dz1cm)[ij], w);
  histo_->hist1d["jet_alphaMax_dz2cm"   +tag]->Fill((*jet_alphaMax_dz2cm)[ij], w);
  histo_->hist1d["jet_alphaMax_dz5cm"   +tag]->Fill((*jet_alphaMax_dz5cm)[ij], w);
  histo_->hist1d["jet_cef"+tag]->Fill((*jet_cef)[ij], w);
  histo_->hist1d["jet_nef"+tag]->Fill((*jet_nef)[ij], w);
  histo_->hist1d["jet_chf"+tag]->Fill((*jet_chf)[ij], w);
  histo_->hist1d["jet_nhf"+tag]->Fill((*jet_nhf)[ij], w);
  // histo_->hist1d["jet_pef"+tag]->Fill((*jet_pef)[ij], w);

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
  histo_->hist2d[string()+"jet_maxIP"+"_VS_"+"jet_alphaMax"+tag]->Fill((*jet_alphaMax)[ij], maxIP, w);
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
  double w = CalculateEventWeight(eventnumber);
  histo_nTrueInt_->Fill(nTrueInt, 1);
  histo_nTrueInt_puweighted_->Fill(nTrueInt, w);
  // histo_->hist1d["pileup_nVtx"]->Fill(nVtx, w);
  // FIXME clean up here
  double ht4 = (*jet_pt)[0] + (*jet_pt)[1] + (*jet_pt)[2] + (*jet_pt)[3];
  histo_->hist1d["norm_ht4"]->Fill(ht4, w);
}

void EmJetHistoMaker::FillHltHistograms(long eventnumber)
{
  double w = CalculateEventWeight(eventnumber);
  if ( (*jet_pt).size() < 4 ) return; //:CUT: Skip events with less than 4 saved jets
  // Calculate ht
  double ht = 0;
  for (unsigned ij = 0; ij < (*jet_pt).size(); ij++) {
    ht += (*jet_pt)[ij];
  }
  // OUTPUT((*jet_pt).size());
  double ht4 = (*jet_pt)[0] + (*jet_pt)[1] + (*jet_pt)[2] + (*jet_pt)[3];

  histo_->hist1d["PFHT900"]->Fill(HLT_PFHT900, w);
  histo_->hist2d["PFHT900_VS_ht"]->Fill(ht, HLT_PFHT900, w);
  histo_->hist2d["PFHT900_VS_ht4"]->Fill(ht4, HLT_PFHT900, w);
  if (HLT_PFHT475) {
    histo_->hist1d["PFHT900__PFHT475"]->Fill(HLT_PFHT900, w);
    histo_->hist2d["PFHT900__PFHT475_VS_ht"]->Fill(ht, HLT_PFHT900, w);
    histo_->hist2d["PFHT900__PFHT475_VS_ht4"]->Fill(ht4, HLT_PFHT900, w);
  }

}

double EmJetHistoMaker::CalculateEventWeight(long eventnumber)
{
  double weight = 1.0;
  if (isData_) {
    // For Data
    weight = 1.0;
  }
  else {
    // For MC
    double pileup_weight = CalculatePileupWeight(eventnumber);
    weight *= pileup_weight;
  }
  return weight;
}

double EmJetHistoMaker::CalculatePileupWeight(long eventnumber)
{
  if(!isData_) {
    double pu_weight = LumiWeights_->weight(nTrueInt);
    return pu_weight;
  }
  else {
    return 1.0;
  }
}

double EmJetHistoMaker::CalculatePdfShift(long eventnumber)
{
  const double Q    = pdf_scalePDF;
  const double id1  = pdf_id1;
  const double id2  = pdf_id2;
  const double x1   = pdf_x1;
  const double x2   = pdf_x2;
  const double pdf1 = pdf_pdf1;
  const double pdf2 = pdf_pdf2;
  const int nWeights = PdfMembers_.size() - 1; // Number of variations not counting central
  // OUTPUT(id1);
  // OUTPUT(id2);
  // OUTPUT(x1);
  // OUTPUT(x2);

  //Loop over members of a given PDF set and get the mean
  float mean = 0;
  for(int ii=1; ii <= nWeights; ++ii){
    auto pdf = PdfMembers_[ii];
    const double xpdf1_new = pdf->xfxQ(id1, x1, Q);
    const double xpdf2_new = pdf->xfxQ(id2, x2, Q);
    // OUTPUT(xpdf1_new);
    // OUTPUT(xpdf2_new);
    double weight = xpdf1_new * xpdf2_new;
    mean += weight;
  }
  mean /= nWeights;
  // OUTPUT(mean);
  // OUTPUT(pdf1*pdf2);

  // loop again for the rms
  float rmssq = 0;
  for(int ii=1; ii <= nWeights; ++ii){
    auto pdf = PdfMembers_[ii];
    const double xpdf1_new = pdf->xfxQ(id1, x1, Q);
    const double xpdf2_new = pdf->xfxQ(id2, x2, Q);
    double weight = xpdf1_new * xpdf2_new;
    rmssq += (weight - mean)*(weight - mean);
  }
  rmssq		/= float((nWeights - 1));
  TH1F* hist = new TH1F("pdf_weights", "pdf_weights", 100, 0, 5);
  for(int ii=1; ii <= nWeights; ++ii){
    auto pdf = PdfMembers_[ii];
    const double xpdf1_new = pdf->xfxQ(id1, x1, Q);
    const double xpdf2_new = pdf->xfxQ(id2, x2, Q);
    double weight = xpdf1_new * xpdf2_new;
    hist->Fill(weight/mean);
  }

  if (mean == 0) {
    OUTPUT("Central weight is zero!");
    OUTPUT(id1);
    OUTPUT(id2);
    OUTPUT(x1);
    OUTPUT(x2);
    OUTPUT(nWeights);
    // If central weight is zero, return zero shift, so we don't fill histograms with non-sensical weights
    // OUTPUT(rmssq);
    // OUTPUT(mean);
    return 0;
  }
  float rms = TMath::Sqrt(rmssq);
  float shift = TMath::Abs(rms/mean);
  // OUTPUT(shift);
  static int first = 1;
  // if (shift > 1 && first == 1) {
  //   OUTPUT("shift > 1");
  //   OUTPUT(id1);
  //   OUTPUT(id2);
  //   OUTPUT(x1);
  //   OUTPUT(x2);
  //   OUTPUT(nWeights);
  //   OUTPUT(rms);
  //   OUTPUT(mean);
  //   for(int ii=1; ii <= nWeights; ++ii){
  //     auto pdf = PdfMembers_[ii];
  //     const double xpdf1_new = pdf->xfxQ(id1, x1, Q);
  //     const double xpdf2_new = pdf->xfxQ(id2, x2, Q);
  //     double weight = xpdf1_new * xpdf2_new;
  //     // OUTPUT(ii);
  //     // OUTPUT(weight);
  //   }
  //   // first = 0;
  // }
  return shift;
}

bool EmJetHistoMaker::CheckPV(long eventnumber)
{
  int nA3DSig0 = 0;
  for(int ij=0; ij<4; ij++){
    if( GetAlpha3DSigM(ij)==0 ) nA3DSig0++;
  }
  if( nA3DSig0>=3 ) return false; // PV reco failed
  else return true;
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
    if( SelectTrack(ij, itk) ) {
      if( (*track_pt)[ij][itk] > 1.0 ) {
        nTrackPassingCut++;
      }
    }
  }

  bool result = true;
  bool cut_nef = (*jet_nef)[ij] < 0.9               ; result = result && cut_nef        ;
  bool cut_cef = (*jet_cef)[ij] < 0.9               ; result = result && cut_cef        ;
  bool cut_eta = abs((*jet_eta)[ij]) < 2.0          ; result = result && cut_eta        ;
  bool cut_tracks = nTrackPassingCut > 0            ; result = result && cut_tracks     ;
  bool cut_ltk_frac = GetLTKFrac(ij) < 0.6          ; result = result && cut_ltk_frac   ;
  return result;
}

bool EmJetHistoMaker::SelectJet_nonpu(int ij)
{
  bool result = true;
  bool cut_nonpu_frac = GetNonPUFrac(ij) > 0.4      ; result = result && cut_nonpu_frac ;
  return result;
}

bool EmJetHistoMaker::SelectJet_alphaMax(int ij)
{
  bool result = true;
  bool cut_alphaMax = (*jet_alphaMax)[ij] < 0.04  ; result = result && cut_alphaMax ;
  return result;
}

bool EmJetHistoMaker::SelectJet_emerging(int ij)
{
  bool result = true;
  bool cut_a3dsigM = GetAlpha3DSigM(ij) < cut_.alpha3d       ; result = result && cut_a3dsigM  ;
  bool cut_theta2D = (*jet_theta2D)[ij] >= 0.0               ; result = result && cut_theta2D  ;
  bool cut_medianIP = GetfabsMedianIP(ij) > cut_.medAbsIp    ; result = result && cut_medianIP ;
  return result;
}

bool EmJetHistoMaker::SelectJet_ipCut(int ij)
{
  bool result = true;
  // Calculate median 2D impact parameter (source=0)
  int nTrack = 0;
  double medianIP = 0.;
  double maxIP = 0.;
  {
    vector<double> vector_ipXY;
    for (unsigned itk=0; itk < (*track_pt)[ij].size(); itk++) {
      if ( (*track_source)[ij][itk] == 0 ) {
        vector_ipXY.push_back( (*track_ipXY)[ij][itk] );
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
  bool ipCut = medianIP > 0.09;
  result = ipCut;
  return result;
}

bool EmJetHistoMaker::SelectTrack(int ij, int itk)
{
  bool result = true;
  bool result1 = (*track_source)[ij][itk] == 0                             ; result = result && result1;
  bool result2 = ( (*track_quality)[ij][itk] & 4 ) > 0                     ; result = result && result2;
  bool result3 = ( fabs((*pv_z)[0]-(*track_ref_z)[ij][itk]) ) < cut_.pu_dz ; result = result && result3;
  return result;
}

void EmJetHistoMaker::PrintEvent (long eventnumber, string comment)
{
  std::cout << "run,lumi,event,comment: " << run << ", " << lumi << ", " << event << ", " << comment << "\n";
}

void
EmJetHistoMaker::SetOptions(Sample sample, long nevent, bool isSignal, bool pileupOnly)
{
  sample_ = sample;
  nevent_ = nevent;
  isSignal_ = isSignal;
  pileupOnly_ = pileupOnly;
}

int
EmJetHistoMaker::VerifyHistograms()
{
  // Verify that histograms are properly initialized
  if (debug>=5) OUTPUT(histo_->hist1d["track_pt"]);
  if ( ! histo_->hist1d["track_pt"] ) {
    return 1;
  }
  return 0;
}

double
EmJetHistoMaker::GetPt(int ij) // Get JEC-shifted pt for given jet (Direction set by sys_.m_direction_jec)
{
  int direction = sys_.GetDirectionJec();
  assert(direction == 0 || direction == -1 || direction == 1);
  if (direction == 0) {
    return  (*jet_pt)[ij];
  } else if (direction == 1) {
    return  (*jet_ptUp)[ij];
  } else if (direction == -1) {
    return  (*jet_ptDown)[ij];
  }
  return -1;
}


double
EmJetHistoMaker::GetAlpha(int ij) // Calculate alpha for given jet
{
  double ptsum_total=0, ptsum=0;
  for (unsigned itk=0; itk < (*track_pt)[ij].size(); itk++) {
    if ( !SelectTrack(ij, itk) ) continue;

    ptsum_total += (*track_pt)[ij][itk];
    if ( (*track_pvWeight)[ij][itk] > 0 ) ptsum += (*track_pt)[ij][itk];
  }

  double alpha = ptsum/ptsum_total;
  return alpha;
}

double EmJetHistoMaker::GetAlpha2DSig(int ij) // Calculate alpha for given jet
{
  double ptsum_total=0, ptsum=0;
  for (unsigned itk=0; itk < (*track_pt)[ij].size(); itk++) {
    if ( !SelectTrack(ij, itk) ) continue;

    ptsum_total += (*track_pt)[ij][itk];
    if ( fabs((*track_ipXYSig)[ij][itk]) < 4.0 ) ptsum += (*track_pt)[ij][itk];// 2D significance matching
  }

  double alpha = (ptsum_total > 0 ? ptsum/ptsum_total : 0.);
  return alpha;
}

double EmJetHistoMaker::GetAlpha3DSigM(int ij)
{
  double ptsum_total=0, ptsum=0;
  for (unsigned itk=0; itk< (*track_pt)[ij].size(); itk++){
    if ( !SelectTrack(ij, itk) ) continue;

    ptsum_total += (*track_pt)[ij][itk];
    double tk3dsig2 = TMath::Power(((*pv_z)[0]-(*track_ref_z)[ij][itk])/cut_.alpha3d_dz, 2.0) + TMath::Power( sys_.GetShifted3dSig((*track_ipXYSig)[ij][itk]), 2.0 );
    double tk3dsig  = TMath::Sqrt(tk3dsig2);
    if( tk3dsig< cut_.alpha3d_sig ) ptsum += (*track_pt)[ij][itk];
  }

  double alpha3dsigm = (ptsum_total>0 ? ptsum/ptsum_total: -1.);
  return alpha3dsigm;
}


double EmJetHistoMaker::GetLTKFrac(int ij)
{
  double maxtkpT = -1;
  for (unsigned itk=0; itk < (*track_pt)[ij].size(); itk++) {
    if ( !SelectTrack(ij, itk) ) continue;
    if ( (*track_pt)[ij][itk] > maxtkpT ) maxtkpT=(*track_pt)[ij][itk];
  }
  double ltkfrac = maxtkpT/(*jet_pt)[ij];
  return ltkfrac;
}

double EmJetHistoMaker::GetNonPUFrac(int ij)
{
  double sumpT = 0., ptsum_total = 0.;
  for (unsigned itk=0; itk< (*track_pt)[ij].size(); itk++){
    if ( !SelectTrack(ij, itk) ) continue;

    ptsum_total += (*track_pt)[ij][itk];
    if ( fabs((*pv_z)[0]-(*track_ref_z)[ij][itk]) < 1.5 ) sumpT += (*track_pt)[ij][itk]; // distance bigger than 1.5cm is considered pile-up tracks
  }

  double nonpilefrac = (ptsum_total > 0 ? sumpT/ptsum_total : 0.);
  return nonpilefrac;
}

double EmJetHistoMaker::GetFrac2DSig(int ij)
{
  int nTrack = 0, nTrackpassing = 0;
  for (unsigned itk=0; itk< (*track_pt)[ij].size(); itk++){
    if ( !SelectTrack(ij, itk) ) continue;

    nTrack++;
    if( fabs((*track_ipXYSig)[ij][itk]) < 2.0 ) nTrackpassing++;
  }

  double frac2DSig = (nTrack > 0 ? (double)nTrackpassing/nTrack : 0.);
  return frac2DSig;
}

double EmJetHistoMaker::GetMedianIP(int ij)
{
  double medip = -10.0;
  vector<double> vector_ipXY;
  for (unsigned itk=0; itk<(*track_pt)[ij].size(); itk++) {
    if ( !SelectTrack(ij, itk) ) continue;
    vector_ipXY.push_back( (*track_ipXY)[ij][itk] );
  }
  std::sort(vector_ipXY.begin(), vector_ipXY.end());
  int nTrack = vector_ipXY.size();
  if ( nTrack>0 ) {
    if ( nTrack%2 ==0 )
      medip = (vector_ipXY[nTrack/2 - 1] + vector_ipXY[nTrack/2]) / 2;
    else
      medip = (vector_ipXY[nTrack/2]);
  }
  return medip;
}

double EmJetHistoMaker::GetfabsMedianIP(int ij)
{
  double medip = -10.0;
  vector<double> vector_ipXY;
  for (unsigned itk=0; itk<(*track_pt)[ij].size(); itk++) {
    if ( !SelectTrack(ij, itk) ) continue;
    double tk_ipXY = sys_.GetShiftedIpXY((*track_ipXY)[ij][itk]) ;
    vector_ipXY.push_back( fabs(tk_ipXY) );
  }

  std::sort(vector_ipXY.begin(), vector_ipXY.end());
  int nTrack = vector_ipXY.size();
  if ( nTrack>0 ) {
    if ( nTrack%2 ==0 )
      medip = (vector_ipXY[nTrack/2 - 1] + vector_ipXY[nTrack/2]) / 2;
    else
      medip = (vector_ipXY[nTrack/2]);
  }
  return medip;
}

int EmJetHistoMaker::GetNTrack(int ij)
{
  int nTrack = 0;
  for( unsigned itk=0; itk<(*track_pt)[ij].size(); itk++){
    if ( !SelectTrack(ij, itk) ) continue;
    nTrack++;
  }
  return nTrack;
}

bool EmJetHistoMaker::GetSignalPartonIndex(long eventnumber)
{
  bool foundallpartons = false;
  for (unsigned igp=0; igp< (*gp_pt).size(); igp++){
    if( (*gp_pdgId)[igp]==4900101 && dkqk_index==-1 && (*gp_pdgId)[igp+1]==1 ){
      dkqk_index = igp;
      dnqk_index = igp+1;
    }
    if( ((*gp_pdgId)[igp]==4900101) && dkqk_index==-1 && (*gp_pdgId)[igp-1]==1 ){
      dkqk_index = igp;
      dnqk_index = igp-1;
    }

    if( (*gp_pdgId)[igp]==-4900101 && adkqk_index==-1 && (*gp_pdgId)[igp+1]==-1 ){
      adkqk_index = igp;
      adnqk_index = igp+1;
    }
    else if( (*gp_pdgId)[igp]==-4900101 && adkqk_index==-1 && (*gp_pdgId)[igp-1]==-1 ){
      adkqk_index = igp;
      adnqk_index = igp-1;
    }
  }

  if( dkqk_index!=-1 && dnqk_index!=-1 && adkqk_index!=-1 && adnqk_index!=-1 ) foundallpartons=true;
  return foundallpartons;
}

bool EmJetHistoMaker::PartonJetMatching(int ij, int igp)
{
  bool matched = false;
  if( DeltaR((*jet_eta)[ij], (*jet_phi)[ij], (*gp_eta)[igp], (*gp_phi)[igp])< 0.4 )  matched=true;
  return matched;
}

double EmJetHistoMaker::GetPVTrackFraction(long eventnumber)
{
  int nTrack_total = 0;
  int nPVTrack = 0;
  for(int ij=0; ij<4; ij++) {
    for(unsigned itk=0; itk<(*track_pt)[ij].size(); itk++) {
      if( (*track_source)[ij][itk]!=0 ) continue;
      if( ( (*track_quality)[ij][itk] & 4 ) == 0 ) continue; // Only process tracks with "highPurity" quality
      if( fabs((*pv_z)[0]-(*track_ref_z)[ij][itk])<0.01 ) nPVTrack++; // count tracks from the PV
      nTrack_total++;
    }
  }
  return (double) nPVTrack/nTrack_total;
}
