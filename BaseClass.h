//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue May  2 17:50:52 2017 by ROOT version 6.06/01
// from TTree emJetTree/emJetTree
// found on file: /store/user/yoshin/EmJetAnalysis/Analysis-20170426-v2/ModelA/EmergingJets_ModelA_TuneCUETP8M1_13TeV_pythia8Mod/Analysis-20170426/170426_192610/0000/ntuple_1.root
//////////////////////////////////////////////////////////

#ifndef BaseClass_h
#define BaseClass_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
using std::vector;
#include "vector"
#include "vector"
#include "vector"

class BaseClass {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           run;
   Int_t           lumi;
   Int_t           event;
   Int_t           bx;
   Int_t           nVtx;
   Int_t           nGoodVtx;
   Int_t           nTrueInt;
   Float_t         met_pt;
   Float_t         met_phi;
   Int_t           nTracks;
   Float_t         alpha_event;
   Float_t         pv_x;
   Float_t         pv_y;
   Float_t         pv_z;
   Float_t         pv_xError;
   Float_t         pv_yError;
   Float_t         pv_zError;
   Float_t         pv_chi2;
   Float_t         pv_ndof;
   Float_t         pv_pt2sum;
   Int_t           pv_nTracks;
   Int_t           pv_indexInColl;
   Int_t           pdf_id1;
   Int_t           pdf_id2;
   Float_t         pdf_x1;
   Float_t         pdf_x2;
   Float_t         pdf_pdf1;
   Float_t         pdf_pdf2;
   Float_t         pdf_scalePDF;
   Bool_t          HLT_PFHT400;
   Bool_t          HLT_PFHT475;
   Bool_t          HLT_PFHT600;
   Bool_t          HLT_PFHT800;
   Bool_t          HLT_PFHT900;
   Bool_t          HLT_HT250;
   Bool_t          HLT_HT350;
   Bool_t          HLT_HT400;
   Bool_t          HLT_HT500;
   vector<int>     *jet_index;
   vector<int>     *jet_source;
   vector<float>   *jet_pt;
   vector<float>   *jet_eta;
   vector<float>   *jet_phi;
   vector<float>   *jet_cef;
   vector<float>   *jet_nef;
   vector<float>   *jet_chf;
   vector<float>   *jet_nhf;
   vector<float>   *jet_pef;
   vector<float>   *jet_mef;
   vector<int>     *jet_missHits;
   vector<int>     *jet_muonHits;
   vector<float>   *jet_alpha;
   vector<float>   *jet_alpha2;
   vector<float>   *jet_alphaMax;
   vector<float>   *jet_alphaMax2;
   vector<float>   *jet_alpha_gen;
   vector<float>   *jet_alphaMax_dz100nm;
   vector<float>   *jet_alphaMax_dz200nm;
   vector<float>   *jet_alphaMax_dz500nm;
   vector<float>   *jet_alphaMax_dz1um;
   vector<float>   *jet_alphaMax_dz2um;
   vector<float>   *jet_alphaMax_dz5um;
   vector<float>   *jet_alphaMax_dz10um;
   vector<float>   *jet_alphaMax_dz20um;
   vector<float>   *jet_alphaMax_dz50um;
   vector<float>   *jet_alphaMax_dz100um;
   vector<float>   *jet_alphaMax_dz200um;
   vector<float>   *jet_alphaMax_dz500um;
   vector<float>   *jet_alphaMax_dz1mm;
   vector<float>   *jet_alphaMax_dz2mm;
   vector<float>   *jet_alphaMax_dz5mm;
   vector<float>   *jet_alphaMax_dz1cm;
   vector<float>   *jet_alphaMax_dz2cm;
   vector<float>   *jet_alphaMax_dz5cm;
   vector<float>   *jet_alphaMax_dz10cm;
   vector<float>   *jet_alphaMax_dz20cm;
   vector<float>   *jet_alphaMax_dz50cm;
   vector<float>   *jet_alphaMax2_dz100nm;
   vector<float>   *jet_alphaMax2_dz200nm;
   vector<float>   *jet_alphaMax2_dz500nm;
   vector<float>   *jet_alphaMax2_dz1um;
   vector<float>   *jet_alphaMax2_dz2um;
   vector<float>   *jet_alphaMax2_dz5um;
   vector<float>   *jet_alphaMax2_dz10um;
   vector<float>   *jet_alphaMax2_dz20um;
   vector<float>   *jet_alphaMax2_dz50um;
   vector<float>   *jet_alphaMax2_dz100um;
   vector<float>   *jet_alphaMax2_dz200um;
   vector<float>   *jet_alphaMax2_dz500um;
   vector<float>   *jet_alphaMax2_dz1mm;
   vector<float>   *jet_alphaMax2_dz2mm;
   vector<float>   *jet_alphaMax2_dz5mm;
   vector<float>   *jet_alphaMax2_dz1cm;
   vector<float>   *jet_alphaMax2_dz2cm;
   vector<float>   *jet_alphaMax2_dz5cm;
   vector<float>   *jet_alphaMax2_dz10cm;
   vector<float>   *jet_alphaMax2_dz20cm;
   vector<float>   *jet_alphaMax2_dz50cm;
   vector<int>     *jet_nDarkPions;
   vector<int>     *jet_nDarkGluons;
   vector<float>   *jet_minDRDarkPion;
   vector<float>   *jet_theta2D;
   vector<vector<int> > *track_index;
   vector<vector<int> > *track_source;
   vector<vector<int> > *track_jet_index;
   vector<vector<int> > *track_vertex_index;
   vector<vector<float> > *track_vertex_weight;
   vector<vector<int> > *track_nHitsInFrontOfVert;
   vector<vector<int> > *track_missHitsAfterVert;
   vector<vector<float> > *track_pt;
   vector<vector<float> > *track_eta;
   vector<vector<float> > *track_phi;
   vector<vector<float> > *track_pca_r;
   vector<vector<float> > *track_pca_eta;
   vector<vector<float> > *track_pca_phi;
   vector<vector<int> > *track_quality;
   vector<vector<int> > *track_algo;
   vector<vector<int> > *track_originalAlgo;
   vector<vector<int> > *track_nHits;
   vector<vector<int> > *track_nMissInnerHits;
   vector<vector<int> > *track_nTrkLayers;
   vector<vector<int> > *track_nMissInnerTrkLayers;
   vector<vector<int> > *track_nMissOuterTrkLayers;
   vector<vector<int> > *track_nMissTrkLayers;
   vector<vector<int> > *track_nPxlLayers;
   vector<vector<int> > *track_nMissInnerPxlLayers;
   vector<vector<int> > *track_nMissOuterPxlLayers;
   vector<vector<int> > *track_nMissPxlLayers;
   vector<vector<float> > *track_ipXY;
   vector<vector<float> > *track_ipZ;
   vector<vector<float> > *track_ipXYSig;
   vector<vector<float> > *track_dRToJetAxis;
   vector<vector<float> > *track_distanceToJet;
   vector<vector<float> > *track_minVertexDz;
   vector<vector<float> > *track_pvWeight;
   vector<vector<int> > *vertex_index;
   vector<vector<int> > *vertex_source;
   vector<vector<int> > *vertex_jet_index;
   vector<vector<float> > *vertex_x;
   vector<vector<float> > *vertex_y;
   vector<vector<float> > *vertex_z;
   vector<vector<float> > *vertex_xError;
   vector<vector<float> > *vertex_yError;
   vector<vector<float> > *vertex_zError;
   vector<vector<float> > *vertex_deltaR;
   vector<vector<float> > *vertex_Lxy;
   vector<vector<float> > *vertex_mass;
   vector<vector<float> > *vertex_chi2;
   vector<vector<float> > *vertex_ndof;
   vector<vector<float> > *vertex_pt2sum;
   vector<int>     *gp_index;
   vector<int>     *gp_status;
   vector<int>     *gp_pdgId;
   vector<int>     *gp_charge;
   vector<float>   *gp_mass;
   vector<float>   *gp_pt;
   vector<float>   *gp_eta;
   vector<float>   *gp_phi;
   vector<float>   *gp_vx;
   vector<float>   *gp_vy;
   vector<float>   *gp_vz;
   vector<float>   *gp_min2Ddist;
   vector<float>   *gp_min2Dsig;
   vector<float>   *gp_min3Ddist;
   vector<float>   *gp_min3Dsig;
   vector<float>   *gp_minDeltaR;
   vector<float>   *gp_matched2Ddist;
   vector<float>   *gp_matched2Dsig;
   vector<float>   *gp_matched3Ddist;
   vector<float>   *gp_matched3Dsig;
   vector<float>   *gp_matchedDeltaR;
   vector<float>   *gp_Lxy;
   vector<int>     *gp_isDark;
   vector<int>     *gp_nDaughters;
   vector<int>     *gp_hasSMDaughter;
   vector<int>     *gp_hasDarkMother;
   vector<int>     *gp_hasDarkPionMother;
   vector<int>     *gp_isTrackable;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_event;   //!
   TBranch        *b_bx;   //!
   TBranch        *b_nVtx;   //!
   TBranch        *b_nGoodVtx;   //!
   TBranch        *b_nTrueInt;   //!
   TBranch        *b_met_pt;   //!
   TBranch        *b_met_phi;   //!
   TBranch        *b_nTracks;   //!
   TBranch        *b_alpha_event;   //!
   TBranch        *b_pv_x;   //!
   TBranch        *b_pv_y;   //!
   TBranch        *b_pv_z;   //!
   TBranch        *b_pv_xError;   //!
   TBranch        *b_pv_yError;   //!
   TBranch        *b_pv_zError;   //!
   TBranch        *b_pv_chi2;   //!
   TBranch        *b_pv_ndof;   //!
   TBranch        *b_pv_pt2sum;   //!
   TBranch        *b_pv_nTracks;   //!
   TBranch        *b_pv_indexInColl;   //!
   TBranch        *b_pdf_id1;   //!
   TBranch        *b_pdf_id2;   //!
   TBranch        *b_pdf_x1;   //!
   TBranch        *b_pdf_x2;   //!
   TBranch        *b_pdf_pdf1;   //!
   TBranch        *b_pdf_pdf2;   //!
   TBranch        *b_pdf_scalePDF;   //!
   TBranch        *b_HLT_PFHT400;   //!
   TBranch        *b_HLT_PFHT475;   //!
   TBranch        *b_HLT_PFHT600;   //!
   TBranch        *b_HLT_PFHT800;   //!
   TBranch        *b_HLT_PFHT900;   //!
   TBranch        *b_HLT_HT250;   //!
   TBranch        *b_HLT_HT350;   //!
   TBranch        *b_HLT_HT400;   //!
   TBranch        *b_HLT_HT500;   //!
   TBranch        *b_jet_index;   //!
   TBranch        *b_jet_source;   //!
   TBranch        *b_jet_pt;   //!
   TBranch        *b_jet_eta;   //!
   TBranch        *b_jet_phi;   //!
   TBranch        *b_jet_cef;   //!
   TBranch        *b_jet_nef;   //!
   TBranch        *b_jet_chf;   //!
   TBranch        *b_jet_nhf;   //!
   TBranch        *b_jet_pef;   //!
   TBranch        *b_jet_mef;   //!
   TBranch        *b_jet_missHits;   //!
   TBranch        *b_jet_muonHits;   //!
   TBranch        *b_jet_alpha;   //!
   TBranch        *b_jet_alpha2;   //!
   TBranch        *b_jet_alphaMax;   //!
   TBranch        *b_jet_alphaMax2;   //!
   TBranch        *b_jet_alpha_gen;   //!
   TBranch        *b_jet_alphaMax_dz100nm;   //!
   TBranch        *b_jet_alphaMax_dz200nm;   //!
   TBranch        *b_jet_alphaMax_dz500nm;   //!
   TBranch        *b_jet_alphaMax_dz1um;   //!
   TBranch        *b_jet_alphaMax_dz2um;   //!
   TBranch        *b_jet_alphaMax_dz5um;   //!
   TBranch        *b_jet_alphaMax_dz10um;   //!
   TBranch        *b_jet_alphaMax_dz20um;   //!
   TBranch        *b_jet_alphaMax_dz50um;   //!
   TBranch        *b_jet_alphaMax_dz100um;   //!
   TBranch        *b_jet_alphaMax_dz200um;   //!
   TBranch        *b_jet_alphaMax_dz500um;   //!
   TBranch        *b_jet_alphaMax_dz1mm;   //!
   TBranch        *b_jet_alphaMax_dz2mm;   //!
   TBranch        *b_jet_alphaMax_dz5mm;   //!
   TBranch        *b_jet_alphaMax_dz1cm;   //!
   TBranch        *b_jet_alphaMax_dz2cm;   //!
   TBranch        *b_jet_alphaMax_dz5cm;   //!
   TBranch        *b_jet_alphaMax_dz10cm;   //!
   TBranch        *b_jet_alphaMax_dz20cm;   //!
   TBranch        *b_jet_alphaMax_dz50cm;   //!
   TBranch        *b_jet_alphaMax2_dz100nm;   //!
   TBranch        *b_jet_alphaMax2_dz200nm;   //!
   TBranch        *b_jet_alphaMax2_dz500nm;   //!
   TBranch        *b_jet_alphaMax2_dz1um;   //!
   TBranch        *b_jet_alphaMax2_dz2um;   //!
   TBranch        *b_jet_alphaMax2_dz5um;   //!
   TBranch        *b_jet_alphaMax2_dz10um;   //!
   TBranch        *b_jet_alphaMax2_dz20um;   //!
   TBranch        *b_jet_alphaMax2_dz50um;   //!
   TBranch        *b_jet_alphaMax2_dz100um;   //!
   TBranch        *b_jet_alphaMax2_dz200um;   //!
   TBranch        *b_jet_alphaMax2_dz500um;   //!
   TBranch        *b_jet_alphaMax2_dz1mm;   //!
   TBranch        *b_jet_alphaMax2_dz2mm;   //!
   TBranch        *b_jet_alphaMax2_dz5mm;   //!
   TBranch        *b_jet_alphaMax2_dz1cm;   //!
   TBranch        *b_jet_alphaMax2_dz2cm;   //!
   TBranch        *b_jet_alphaMax2_dz5cm;   //!
   TBranch        *b_jet_alphaMax2_dz10cm;   //!
   TBranch        *b_jet_alphaMax2_dz20cm;   //!
   TBranch        *b_jet_alphaMax2_dz50cm;   //!
   TBranch        *b_jet_nDarkPions;   //!
   TBranch        *b_jet_nDarkGluons;   //!
   TBranch        *b_jet_minDRDarkPion;   //!
   TBranch        *b_jet_theta2D;   //!
   TBranch        *b_track_index;   //!
   TBranch        *b_track_source;   //!
   TBranch        *b_track_jet_index;   //!
   TBranch        *b_track_vertex_index;   //!
   TBranch        *b_track_vertex_weight;   //!
   TBranch        *b_track_nHitsInFrontOfVert;   //!
   TBranch        *b_track_missHitsAfterVert;   //!
   TBranch        *b_track_pt;   //!
   TBranch        *b_track_eta;   //!
   TBranch        *b_track_phi;   //!
   TBranch        *b_track_pca_r;   //!
   TBranch        *b_track_pca_eta;   //!
   TBranch        *b_track_pca_phi;   //!
   TBranch        *b_track_quality;   //!
   TBranch        *b_track_algo;   //!
   TBranch        *b_track_originalAlgo;   //!
   TBranch        *b_track_nHits;   //!
   TBranch        *b_track_nMissInnerHits;   //!
   TBranch        *b_track_nTrkLayers;   //!
   TBranch        *b_track_nMissInnerTrkLayers;   //!
   TBranch        *b_track_nMissOuterTrkLayers;   //!
   TBranch        *b_track_nMissTrkLayers;   //!
   TBranch        *b_track_nPxlLayers;   //!
   TBranch        *b_track_nMissInnerPxlLayers;   //!
   TBranch        *b_track_nMissOuterPxlLayers;   //!
   TBranch        *b_track_nMissPxlLayers;   //!
   TBranch        *b_track_ipXY;   //!
   TBranch        *b_track_ipZ;   //!
   TBranch        *b_track_ipXYSig;   //!
   TBranch        *b_track_dRToJetAxis;   //!
   TBranch        *b_track_distanceToJet;   //!
   TBranch        *b_track_minVertexDz;   //!
   TBranch        *b_track_pvWeight;   //!
   TBranch        *b_vertex_index;   //!
   TBranch        *b_vertex_source;   //!
   TBranch        *b_vertex_jet_index;   //!
   TBranch        *b_vertex_x;   //!
   TBranch        *b_vertex_y;   //!
   TBranch        *b_vertex_z;   //!
   TBranch        *b_vertex_xError;   //!
   TBranch        *b_vertex_yError;   //!
   TBranch        *b_vertex_zError;   //!
   TBranch        *b_vertex_deltaR;   //!
   TBranch        *b_vertex_Lxy;   //!
   TBranch        *b_vertex_mass;   //!
   TBranch        *b_vertex_chi2;   //!
   TBranch        *b_vertex_ndof;   //!
   TBranch        *b_vertex_pt2sum;   //!
   TBranch        *b_gp_index;   //!
   TBranch        *b_gp_status;   //!
   TBranch        *b_gp_pdgId;   //!
   TBranch        *b_gp_charge;   //!
   TBranch        *b_gp_mass;   //!
   TBranch        *b_gp_pt;   //!
   TBranch        *b_gp_eta;   //!
   TBranch        *b_gp_phi;   //!
   TBranch        *b_gp_vx;   //!
   TBranch        *b_gp_vy;   //!
   TBranch        *b_gp_vz;   //!
   TBranch        *b_gp_min2Ddist;   //!
   TBranch        *b_gp_min2Dsig;   //!
   TBranch        *b_gp_min3Ddist;   //!
   TBranch        *b_gp_min3Dsig;   //!
   TBranch        *b_gp_minDeltaR;   //!
   TBranch        *b_gp_matched2Ddist;   //!
   TBranch        *b_gp_matched2Dsig;   //!
   TBranch        *b_gp_matched3Ddist;   //!
   TBranch        *b_gp_matched3Dsig;   //!
   TBranch        *b_gp_matchedDeltaR;   //!
   TBranch        *b_gp_Lxy;   //!
   TBranch        *b_gp_isDark;   //!
   TBranch        *b_gp_nDaughters;   //!
   TBranch        *b_gp_hasSMDaughter;   //!
   TBranch        *b_gp_hasDarkMother;   //!
   TBranch        *b_gp_hasDarkPionMother;   //!
   TBranch        *b_gp_isTrackable;   //!

   BaseClass(TTree *tree=0);
   virtual ~BaseClass();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef BaseClass_cxx
BaseClass::BaseClass(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/store/user/yoshin/EmJetAnalysis/Analysis-20170426-v2/ModelA/EmergingJets_ModelA_TuneCUETP8M1_13TeV_pythia8Mod/Analysis-20170426/170426_192610/0000/ntuple_1.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/store/user/yoshin/EmJetAnalysis/Analysis-20170426-v2/ModelA/EmergingJets_ModelA_TuneCUETP8M1_13TeV_pythia8Mod/Analysis-20170426/170426_192610/0000/ntuple_1.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("/store/user/yoshin/EmJetAnalysis/Analysis-20170426-v2/ModelA/EmergingJets_ModelA_TuneCUETP8M1_13TeV_pythia8Mod/Analysis-20170426/170426_192610/0000/ntuple_1.root:/emJetAnalyzer");
      dir->GetObject("emJetTree",tree);

   }
   Init(tree);
}

BaseClass::~BaseClass()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t BaseClass::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t BaseClass::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void BaseClass::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   jet_index = 0;
   jet_source = 0;
   jet_pt = 0;
   jet_eta = 0;
   jet_phi = 0;
   jet_cef = 0;
   jet_nef = 0;
   jet_chf = 0;
   jet_nhf = 0;
   jet_pef = 0;
   jet_mef = 0;
   jet_missHits = 0;
   jet_muonHits = 0;
   jet_alpha = 0;
   jet_alpha2 = 0;
   jet_alphaMax = 0;
   jet_alphaMax2 = 0;
   jet_alpha_gen = 0;
   jet_alphaMax_dz100nm = 0;
   jet_alphaMax_dz200nm = 0;
   jet_alphaMax_dz500nm = 0;
   jet_alphaMax_dz1um = 0;
   jet_alphaMax_dz2um = 0;
   jet_alphaMax_dz5um = 0;
   jet_alphaMax_dz10um = 0;
   jet_alphaMax_dz20um = 0;
   jet_alphaMax_dz50um = 0;
   jet_alphaMax_dz100um = 0;
   jet_alphaMax_dz200um = 0;
   jet_alphaMax_dz500um = 0;
   jet_alphaMax_dz1mm = 0;
   jet_alphaMax_dz2mm = 0;
   jet_alphaMax_dz5mm = 0;
   jet_alphaMax_dz1cm = 0;
   jet_alphaMax_dz2cm = 0;
   jet_alphaMax_dz5cm = 0;
   jet_alphaMax_dz10cm = 0;
   jet_alphaMax_dz20cm = 0;
   jet_alphaMax_dz50cm = 0;
   jet_alphaMax2_dz100nm = 0;
   jet_alphaMax2_dz200nm = 0;
   jet_alphaMax2_dz500nm = 0;
   jet_alphaMax2_dz1um = 0;
   jet_alphaMax2_dz2um = 0;
   jet_alphaMax2_dz5um = 0;
   jet_alphaMax2_dz10um = 0;
   jet_alphaMax2_dz20um = 0;
   jet_alphaMax2_dz50um = 0;
   jet_alphaMax2_dz100um = 0;
   jet_alphaMax2_dz200um = 0;
   jet_alphaMax2_dz500um = 0;
   jet_alphaMax2_dz1mm = 0;
   jet_alphaMax2_dz2mm = 0;
   jet_alphaMax2_dz5mm = 0;
   jet_alphaMax2_dz1cm = 0;
   jet_alphaMax2_dz2cm = 0;
   jet_alphaMax2_dz5cm = 0;
   jet_alphaMax2_dz10cm = 0;
   jet_alphaMax2_dz20cm = 0;
   jet_alphaMax2_dz50cm = 0;
   jet_nDarkPions = 0;
   jet_nDarkGluons = 0;
   jet_minDRDarkPion = 0;
   jet_theta2D = 0;
   track_index = 0;
   track_source = 0;
   track_jet_index = 0;
   track_vertex_index = 0;
   track_vertex_weight = 0;
   track_nHitsInFrontOfVert = 0;
   track_missHitsAfterVert = 0;
   track_pt = 0;
   track_eta = 0;
   track_phi = 0;
   track_pca_r = 0;
   track_pca_eta = 0;
   track_pca_phi = 0;
   track_quality = 0;
   track_algo = 0;
   track_originalAlgo = 0;
   track_nHits = 0;
   track_nMissInnerHits = 0;
   track_nTrkLayers = 0;
   track_nMissInnerTrkLayers = 0;
   track_nMissOuterTrkLayers = 0;
   track_nMissTrkLayers = 0;
   track_nPxlLayers = 0;
   track_nMissInnerPxlLayers = 0;
   track_nMissOuterPxlLayers = 0;
   track_nMissPxlLayers = 0;
   track_ipXY = 0;
   track_ipZ = 0;
   track_ipXYSig = 0;
   track_dRToJetAxis = 0;
   track_distanceToJet = 0;
   track_minVertexDz = 0;
   track_pvWeight = 0;
   vertex_index = 0;
   vertex_source = 0;
   vertex_jet_index = 0;
   vertex_x = 0;
   vertex_y = 0;
   vertex_z = 0;
   vertex_xError = 0;
   vertex_yError = 0;
   vertex_zError = 0;
   vertex_deltaR = 0;
   vertex_Lxy = 0;
   vertex_mass = 0;
   vertex_chi2 = 0;
   vertex_ndof = 0;
   vertex_pt2sum = 0;
   gp_index = 0;
   gp_status = 0;
   gp_pdgId = 0;
   gp_charge = 0;
   gp_mass = 0;
   gp_pt = 0;
   gp_eta = 0;
   gp_phi = 0;
   gp_vx = 0;
   gp_vy = 0;
   gp_vz = 0;
   gp_min2Ddist = 0;
   gp_min2Dsig = 0;
   gp_min3Ddist = 0;
   gp_min3Dsig = 0;
   gp_minDeltaR = 0;
   gp_matched2Ddist = 0;
   gp_matched2Dsig = 0;
   gp_matched3Ddist = 0;
   gp_matched3Dsig = 0;
   gp_matchedDeltaR = 0;
   gp_Lxy = 0;
   gp_isDark = 0;
   gp_nDaughters = 0;
   gp_hasSMDaughter = 0;
   gp_hasDarkMother = 0;
   gp_hasDarkPionMother = 0;
   gp_isTrackable = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("bx", &bx, &b_bx);
   fChain->SetBranchAddress("nVtx", &nVtx, &b_nVtx);
   fChain->SetBranchAddress("nGoodVtx", &nGoodVtx, &b_nGoodVtx);
   fChain->SetBranchAddress("nTrueInt", &nTrueInt, &b_nTrueInt);
   fChain->SetBranchAddress("met_pt", &met_pt, &b_met_pt);
   fChain->SetBranchAddress("met_phi", &met_phi, &b_met_phi);
   fChain->SetBranchAddress("nTracks", &nTracks, &b_nTracks);
   fChain->SetBranchAddress("alpha_event", &alpha_event, &b_alpha_event);
   fChain->SetBranchAddress("pv_x", &pv_x, &b_pv_x);
   fChain->SetBranchAddress("pv_y", &pv_y, &b_pv_y);
   fChain->SetBranchAddress("pv_z", &pv_z, &b_pv_z);
   fChain->SetBranchAddress("pv_xError", &pv_xError, &b_pv_xError);
   fChain->SetBranchAddress("pv_yError", &pv_yError, &b_pv_yError);
   fChain->SetBranchAddress("pv_zError", &pv_zError, &b_pv_zError);
   fChain->SetBranchAddress("pv_chi2", &pv_chi2, &b_pv_chi2);
   fChain->SetBranchAddress("pv_ndof", &pv_ndof, &b_pv_ndof);
   fChain->SetBranchAddress("pv_pt2sum", &pv_pt2sum, &b_pv_pt2sum);
   fChain->SetBranchAddress("pv_nTracks", &pv_nTracks, &b_pv_nTracks);
   fChain->SetBranchAddress("pv_indexInColl", &pv_indexInColl, &b_pv_indexInColl);
   fChain->SetBranchAddress("pdf_id1", &pdf_id1, &b_pdf_id1);
   fChain->SetBranchAddress("pdf_id2", &pdf_id2, &b_pdf_id2);
   fChain->SetBranchAddress("pdf_x1", &pdf_x1, &b_pdf_x1);
   fChain->SetBranchAddress("pdf_x2", &pdf_x2, &b_pdf_x2);
   fChain->SetBranchAddress("pdf_pdf1", &pdf_pdf1, &b_pdf_pdf1);
   fChain->SetBranchAddress("pdf_pdf2", &pdf_pdf2, &b_pdf_pdf2);
   fChain->SetBranchAddress("pdf_scalePDF", &pdf_scalePDF, &b_pdf_scalePDF);
   fChain->SetBranchAddress("HLT_PFHT400", &HLT_PFHT400, &b_HLT_PFHT400);
   fChain->SetBranchAddress("HLT_PFHT475", &HLT_PFHT475, &b_HLT_PFHT475);
   fChain->SetBranchAddress("HLT_PFHT600", &HLT_PFHT600, &b_HLT_PFHT600);
   fChain->SetBranchAddress("HLT_PFHT800", &HLT_PFHT800, &b_HLT_PFHT800);
   fChain->SetBranchAddress("HLT_PFHT900", &HLT_PFHT900, &b_HLT_PFHT900);
   fChain->SetBranchAddress("HLT_HT250", &HLT_HT250, &b_HLT_HT250);
   fChain->SetBranchAddress("HLT_HT350", &HLT_HT350, &b_HLT_HT350);
   fChain->SetBranchAddress("HLT_HT400", &HLT_HT400, &b_HLT_HT400);
   fChain->SetBranchAddress("HLT_HT500", &HLT_HT500, &b_HLT_HT500);
   fChain->SetBranchAddress("jet_index", &jet_index, &b_jet_index);
   fChain->SetBranchAddress("jet_source", &jet_source, &b_jet_source);
   fChain->SetBranchAddress("jet_pt", &jet_pt, &b_jet_pt);
   fChain->SetBranchAddress("jet_eta", &jet_eta, &b_jet_eta);
   fChain->SetBranchAddress("jet_phi", &jet_phi, &b_jet_phi);
   fChain->SetBranchAddress("jet_cef", &jet_cef, &b_jet_cef);
   fChain->SetBranchAddress("jet_nef", &jet_nef, &b_jet_nef);
   fChain->SetBranchAddress("jet_chf", &jet_chf, &b_jet_chf);
   fChain->SetBranchAddress("jet_nhf", &jet_nhf, &b_jet_nhf);
   fChain->SetBranchAddress("jet_pef", &jet_pef, &b_jet_pef);
   fChain->SetBranchAddress("jet_mef", &jet_mef, &b_jet_mef);
   fChain->SetBranchAddress("jet_missHits", &jet_missHits, &b_jet_missHits);
   fChain->SetBranchAddress("jet_muonHits", &jet_muonHits, &b_jet_muonHits);
   fChain->SetBranchAddress("jet_alpha", &jet_alpha, &b_jet_alpha);
   fChain->SetBranchAddress("jet_alpha2", &jet_alpha2, &b_jet_alpha2);
   fChain->SetBranchAddress("jet_alphaMax", &jet_alphaMax, &b_jet_alphaMax);
   fChain->SetBranchAddress("jet_alphaMax2", &jet_alphaMax2, &b_jet_alphaMax2);
   fChain->SetBranchAddress("jet_alpha_gen", &jet_alpha_gen, &b_jet_alpha_gen);
   fChain->SetBranchAddress("jet_alphaMax_dz100nm", &jet_alphaMax_dz100nm, &b_jet_alphaMax_dz100nm);
   fChain->SetBranchAddress("jet_alphaMax_dz200nm", &jet_alphaMax_dz200nm, &b_jet_alphaMax_dz200nm);
   fChain->SetBranchAddress("jet_alphaMax_dz500nm", &jet_alphaMax_dz500nm, &b_jet_alphaMax_dz500nm);
   fChain->SetBranchAddress("jet_alphaMax_dz1um", &jet_alphaMax_dz1um, &b_jet_alphaMax_dz1um);
   fChain->SetBranchAddress("jet_alphaMax_dz2um", &jet_alphaMax_dz2um, &b_jet_alphaMax_dz2um);
   fChain->SetBranchAddress("jet_alphaMax_dz5um", &jet_alphaMax_dz5um, &b_jet_alphaMax_dz5um);
   fChain->SetBranchAddress("jet_alphaMax_dz10um", &jet_alphaMax_dz10um, &b_jet_alphaMax_dz10um);
   fChain->SetBranchAddress("jet_alphaMax_dz20um", &jet_alphaMax_dz20um, &b_jet_alphaMax_dz20um);
   fChain->SetBranchAddress("jet_alphaMax_dz50um", &jet_alphaMax_dz50um, &b_jet_alphaMax_dz50um);
   fChain->SetBranchAddress("jet_alphaMax_dz100um", &jet_alphaMax_dz100um, &b_jet_alphaMax_dz100um);
   fChain->SetBranchAddress("jet_alphaMax_dz200um", &jet_alphaMax_dz200um, &b_jet_alphaMax_dz200um);
   fChain->SetBranchAddress("jet_alphaMax_dz500um", &jet_alphaMax_dz500um, &b_jet_alphaMax_dz500um);
   fChain->SetBranchAddress("jet_alphaMax_dz1mm", &jet_alphaMax_dz1mm, &b_jet_alphaMax_dz1mm);
   fChain->SetBranchAddress("jet_alphaMax_dz2mm", &jet_alphaMax_dz2mm, &b_jet_alphaMax_dz2mm);
   fChain->SetBranchAddress("jet_alphaMax_dz5mm", &jet_alphaMax_dz5mm, &b_jet_alphaMax_dz5mm);
   fChain->SetBranchAddress("jet_alphaMax_dz1cm", &jet_alphaMax_dz1cm, &b_jet_alphaMax_dz1cm);
   fChain->SetBranchAddress("jet_alphaMax_dz2cm", &jet_alphaMax_dz2cm, &b_jet_alphaMax_dz2cm);
   fChain->SetBranchAddress("jet_alphaMax_dz5cm", &jet_alphaMax_dz5cm, &b_jet_alphaMax_dz5cm);
   fChain->SetBranchAddress("jet_alphaMax_dz10cm", &jet_alphaMax_dz10cm, &b_jet_alphaMax_dz10cm);
   fChain->SetBranchAddress("jet_alphaMax_dz20cm", &jet_alphaMax_dz20cm, &b_jet_alphaMax_dz20cm);
   fChain->SetBranchAddress("jet_alphaMax_dz50cm", &jet_alphaMax_dz50cm, &b_jet_alphaMax_dz50cm);
   fChain->SetBranchAddress("jet_alphaMax2_dz100nm", &jet_alphaMax2_dz100nm, &b_jet_alphaMax2_dz100nm);
   fChain->SetBranchAddress("jet_alphaMax2_dz200nm", &jet_alphaMax2_dz200nm, &b_jet_alphaMax2_dz200nm);
   fChain->SetBranchAddress("jet_alphaMax2_dz500nm", &jet_alphaMax2_dz500nm, &b_jet_alphaMax2_dz500nm);
   fChain->SetBranchAddress("jet_alphaMax2_dz1um", &jet_alphaMax2_dz1um, &b_jet_alphaMax2_dz1um);
   fChain->SetBranchAddress("jet_alphaMax2_dz2um", &jet_alphaMax2_dz2um, &b_jet_alphaMax2_dz2um);
   fChain->SetBranchAddress("jet_alphaMax2_dz5um", &jet_alphaMax2_dz5um, &b_jet_alphaMax2_dz5um);
   fChain->SetBranchAddress("jet_alphaMax2_dz10um", &jet_alphaMax2_dz10um, &b_jet_alphaMax2_dz10um);
   fChain->SetBranchAddress("jet_alphaMax2_dz20um", &jet_alphaMax2_dz20um, &b_jet_alphaMax2_dz20um);
   fChain->SetBranchAddress("jet_alphaMax2_dz50um", &jet_alphaMax2_dz50um, &b_jet_alphaMax2_dz50um);
   fChain->SetBranchAddress("jet_alphaMax2_dz100um", &jet_alphaMax2_dz100um, &b_jet_alphaMax2_dz100um);
   fChain->SetBranchAddress("jet_alphaMax2_dz200um", &jet_alphaMax2_dz200um, &b_jet_alphaMax2_dz200um);
   fChain->SetBranchAddress("jet_alphaMax2_dz500um", &jet_alphaMax2_dz500um, &b_jet_alphaMax2_dz500um);
   fChain->SetBranchAddress("jet_alphaMax2_dz1mm", &jet_alphaMax2_dz1mm, &b_jet_alphaMax2_dz1mm);
   fChain->SetBranchAddress("jet_alphaMax2_dz2mm", &jet_alphaMax2_dz2mm, &b_jet_alphaMax2_dz2mm);
   fChain->SetBranchAddress("jet_alphaMax2_dz5mm", &jet_alphaMax2_dz5mm, &b_jet_alphaMax2_dz5mm);
   fChain->SetBranchAddress("jet_alphaMax2_dz1cm", &jet_alphaMax2_dz1cm, &b_jet_alphaMax2_dz1cm);
   fChain->SetBranchAddress("jet_alphaMax2_dz2cm", &jet_alphaMax2_dz2cm, &b_jet_alphaMax2_dz2cm);
   fChain->SetBranchAddress("jet_alphaMax2_dz5cm", &jet_alphaMax2_dz5cm, &b_jet_alphaMax2_dz5cm);
   fChain->SetBranchAddress("jet_alphaMax2_dz10cm", &jet_alphaMax2_dz10cm, &b_jet_alphaMax2_dz10cm);
   fChain->SetBranchAddress("jet_alphaMax2_dz20cm", &jet_alphaMax2_dz20cm, &b_jet_alphaMax2_dz20cm);
   fChain->SetBranchAddress("jet_alphaMax2_dz50cm", &jet_alphaMax2_dz50cm, &b_jet_alphaMax2_dz50cm);
   fChain->SetBranchAddress("jet_nDarkPions", &jet_nDarkPions, &b_jet_nDarkPions);
   fChain->SetBranchAddress("jet_nDarkGluons", &jet_nDarkGluons, &b_jet_nDarkGluons);
   fChain->SetBranchAddress("jet_minDRDarkPion", &jet_minDRDarkPion, &b_jet_minDRDarkPion);
   fChain->SetBranchAddress("jet_theta2D", &jet_theta2D, &b_jet_theta2D);
   fChain->SetBranchAddress("track_index", &track_index, &b_track_index);
   fChain->SetBranchAddress("track_source", &track_source, &b_track_source);
   fChain->SetBranchAddress("track_jet_index", &track_jet_index, &b_track_jet_index);
   fChain->SetBranchAddress("track_vertex_index", &track_vertex_index, &b_track_vertex_index);
   fChain->SetBranchAddress("track_vertex_weight", &track_vertex_weight, &b_track_vertex_weight);
   fChain->SetBranchAddress("track_nHitsInFrontOfVert", &track_nHitsInFrontOfVert, &b_track_nHitsInFrontOfVert);
   fChain->SetBranchAddress("track_missHitsAfterVert", &track_missHitsAfterVert, &b_track_missHitsAfterVert);
   fChain->SetBranchAddress("track_pt", &track_pt, &b_track_pt);
   fChain->SetBranchAddress("track_eta", &track_eta, &b_track_eta);
   fChain->SetBranchAddress("track_phi", &track_phi, &b_track_phi);
   fChain->SetBranchAddress("track_pca_r", &track_pca_r, &b_track_pca_r);
   fChain->SetBranchAddress("track_pca_eta", &track_pca_eta, &b_track_pca_eta);
   fChain->SetBranchAddress("track_pca_phi", &track_pca_phi, &b_track_pca_phi);
   fChain->SetBranchAddress("track_quality", &track_quality, &b_track_quality);
   fChain->SetBranchAddress("track_algo", &track_algo, &b_track_algo);
   fChain->SetBranchAddress("track_originalAlgo", &track_originalAlgo, &b_track_originalAlgo);
   fChain->SetBranchAddress("track_nHits", &track_nHits, &b_track_nHits);
   fChain->SetBranchAddress("track_nMissInnerHits", &track_nMissInnerHits, &b_track_nMissInnerHits);
   fChain->SetBranchAddress("track_nTrkLayers", &track_nTrkLayers, &b_track_nTrkLayers);
   fChain->SetBranchAddress("track_nMissInnerTrkLayers", &track_nMissInnerTrkLayers, &b_track_nMissInnerTrkLayers);
   fChain->SetBranchAddress("track_nMissOuterTrkLayers", &track_nMissOuterTrkLayers, &b_track_nMissOuterTrkLayers);
   fChain->SetBranchAddress("track_nMissTrkLayers", &track_nMissTrkLayers, &b_track_nMissTrkLayers);
   fChain->SetBranchAddress("track_nPxlLayers", &track_nPxlLayers, &b_track_nPxlLayers);
   fChain->SetBranchAddress("track_nMissInnerPxlLayers", &track_nMissInnerPxlLayers, &b_track_nMissInnerPxlLayers);
   fChain->SetBranchAddress("track_nMissOuterPxlLayers", &track_nMissOuterPxlLayers, &b_track_nMissOuterPxlLayers);
   fChain->SetBranchAddress("track_nMissPxlLayers", &track_nMissPxlLayers, &b_track_nMissPxlLayers);
   fChain->SetBranchAddress("track_ipXY", &track_ipXY, &b_track_ipXY);
   fChain->SetBranchAddress("track_ipZ", &track_ipZ, &b_track_ipZ);
   fChain->SetBranchAddress("track_ipXYSig", &track_ipXYSig, &b_track_ipXYSig);
   fChain->SetBranchAddress("track_dRToJetAxis", &track_dRToJetAxis, &b_track_dRToJetAxis);
   fChain->SetBranchAddress("track_distanceToJet", &track_distanceToJet, &b_track_distanceToJet);
   fChain->SetBranchAddress("track_minVertexDz", &track_minVertexDz, &b_track_minVertexDz);
   fChain->SetBranchAddress("track_pvWeight", &track_pvWeight, &b_track_pvWeight);
   fChain->SetBranchAddress("vertex_index", &vertex_index, &b_vertex_index);
   fChain->SetBranchAddress("vertex_source", &vertex_source, &b_vertex_source);
   fChain->SetBranchAddress("vertex_jet_index", &vertex_jet_index, &b_vertex_jet_index);
   fChain->SetBranchAddress("vertex_x", &vertex_x, &b_vertex_x);
   fChain->SetBranchAddress("vertex_y", &vertex_y, &b_vertex_y);
   fChain->SetBranchAddress("vertex_z", &vertex_z, &b_vertex_z);
   fChain->SetBranchAddress("vertex_xError", &vertex_xError, &b_vertex_xError);
   fChain->SetBranchAddress("vertex_yError", &vertex_yError, &b_vertex_yError);
   fChain->SetBranchAddress("vertex_zError", &vertex_zError, &b_vertex_zError);
   fChain->SetBranchAddress("vertex_deltaR", &vertex_deltaR, &b_vertex_deltaR);
   fChain->SetBranchAddress("vertex_Lxy", &vertex_Lxy, &b_vertex_Lxy);
   fChain->SetBranchAddress("vertex_mass", &vertex_mass, &b_vertex_mass);
   fChain->SetBranchAddress("vertex_chi2", &vertex_chi2, &b_vertex_chi2);
   fChain->SetBranchAddress("vertex_ndof", &vertex_ndof, &b_vertex_ndof);
   fChain->SetBranchAddress("vertex_pt2sum", &vertex_pt2sum, &b_vertex_pt2sum);
   fChain->SetBranchAddress("gp_index", &gp_index, &b_gp_index);
   fChain->SetBranchAddress("gp_status", &gp_status, &b_gp_status);
   fChain->SetBranchAddress("gp_pdgId", &gp_pdgId, &b_gp_pdgId);
   fChain->SetBranchAddress("gp_charge", &gp_charge, &b_gp_charge);
   fChain->SetBranchAddress("gp_mass", &gp_mass, &b_gp_mass);
   fChain->SetBranchAddress("gp_pt", &gp_pt, &b_gp_pt);
   fChain->SetBranchAddress("gp_eta", &gp_eta, &b_gp_eta);
   fChain->SetBranchAddress("gp_phi", &gp_phi, &b_gp_phi);
   fChain->SetBranchAddress("gp_vx", &gp_vx, &b_gp_vx);
   fChain->SetBranchAddress("gp_vy", &gp_vy, &b_gp_vy);
   fChain->SetBranchAddress("gp_vz", &gp_vz, &b_gp_vz);
   fChain->SetBranchAddress("gp_min2Ddist", &gp_min2Ddist, &b_gp_min2Ddist);
   fChain->SetBranchAddress("gp_min2Dsig", &gp_min2Dsig, &b_gp_min2Dsig);
   fChain->SetBranchAddress("gp_min3Ddist", &gp_min3Ddist, &b_gp_min3Ddist);
   fChain->SetBranchAddress("gp_min3Dsig", &gp_min3Dsig, &b_gp_min3Dsig);
   fChain->SetBranchAddress("gp_minDeltaR", &gp_minDeltaR, &b_gp_minDeltaR);
   fChain->SetBranchAddress("gp_matched2Ddist", &gp_matched2Ddist, &b_gp_matched2Ddist);
   fChain->SetBranchAddress("gp_matched2Dsig", &gp_matched2Dsig, &b_gp_matched2Dsig);
   fChain->SetBranchAddress("gp_matched3Ddist", &gp_matched3Ddist, &b_gp_matched3Ddist);
   fChain->SetBranchAddress("gp_matched3Dsig", &gp_matched3Dsig, &b_gp_matched3Dsig);
   fChain->SetBranchAddress("gp_matchedDeltaR", &gp_matchedDeltaR, &b_gp_matchedDeltaR);
   fChain->SetBranchAddress("gp_Lxy", &gp_Lxy, &b_gp_Lxy);
   fChain->SetBranchAddress("gp_isDark", &gp_isDark, &b_gp_isDark);
   fChain->SetBranchAddress("gp_nDaughters", &gp_nDaughters, &b_gp_nDaughters);
   fChain->SetBranchAddress("gp_hasSMDaughter", &gp_hasSMDaughter, &b_gp_hasSMDaughter);
   fChain->SetBranchAddress("gp_hasDarkMother", &gp_hasDarkMother, &b_gp_hasDarkMother);
   fChain->SetBranchAddress("gp_hasDarkPionMother", &gp_hasDarkPionMother, &b_gp_hasDarkPionMother);
   fChain->SetBranchAddress("gp_isTrackable", &gp_isTrackable, &b_gp_isTrackable);
   Notify();
}

Bool_t BaseClass::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void BaseClass::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t BaseClass::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef BaseClass_cxx
