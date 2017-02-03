//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Aug 22 13:34:46 2016 by ROOT version 5.34/36
// from TTree emJetTree/emJetTree
// found on file: /afs/cern.ch/user/y/yoshin/eos/cms/store/group/phys_exotica/EmergingJets/Analysis-20160821-v0/ModelA/EmergingJets_ModelA_TuneCUETP8M1_13TeV_pythia8Mod/Analysis-20160821/160821_201557/ntuple_merged_ModelA.root
//////////////////////////////////////////////////////////

#ifndef BaseClass_h
#define BaseClass_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>
#include <vector>
#include <iostream>

using std::vector;

// Fixed size dimensions of array or collections stored in the TTree if any.

class BaseClass {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

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
   vector<int>     *jet_index;
   vector<int>     *jet_source;
   vector<float>   *jet_pt;
   vector<float>   *jet_eta;
   vector<float>   *jet_phi;
   vector<float>   *jet_cef;
   vector<float>   *jet_nef;
   vector<float>   *jet_chf;
   vector<float>   *jet_nhf;
   vector<float>   *jet_phf;
   vector<int>     *jet_missHits;
   vector<int>     *jet_muonHits;
   vector<float>   *jet_alphaMax;
   vector<int>     *jet_nDarkPions;
   vector<int>     *jet_nDarkGluons;
   vector<float>   *jet_minDRDarkPion;
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
   TBranch        *b_jet_index;   //!
   TBranch        *b_jet_source;   //!
   TBranch        *b_jet_pt;   //!
   TBranch        *b_jet_eta;   //!
   TBranch        *b_jet_phi;   //!
   TBranch        *b_jet_cef;   //!
   TBranch        *b_jet_nef;   //!
   TBranch        *b_jet_chf;   //!
   TBranch        *b_jet_nhf;   //!
   TBranch        *b_jet_phf;   //!
   TBranch        *b_jet_missHits;   //!
   TBranch        *b_jet_muonHits;   //!
   TBranch        *b_jet_alphaMax;   //!
   TBranch        *b_jet_nDarkPions;   //!
   TBranch        *b_jet_nDarkGluons;   //!
   TBranch        *b_jet_minDRDarkPion;   //!
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
   jet_phf = 0;
   jet_missHits = 0;
   jet_muonHits = 0;
   jet_alphaMax = 0;
   jet_nDarkPions = 0;
   jet_nDarkGluons = 0;
   jet_minDRDarkPion = 0;
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
   fChain->SetBranchAddress("jet_index", &jet_index, &b_jet_index);
   fChain->SetBranchAddress("jet_source", &jet_source, &b_jet_source);
   fChain->SetBranchAddress("jet_pt", &jet_pt, &b_jet_pt);
   fChain->SetBranchAddress("jet_eta", &jet_eta, &b_jet_eta);
   fChain->SetBranchAddress("jet_phi", &jet_phi, &b_jet_phi);
   fChain->SetBranchAddress("jet_cef", &jet_cef, &b_jet_cef);
   fChain->SetBranchAddress("jet_nef", &jet_nef, &b_jet_nef);
   fChain->SetBranchAddress("jet_chf", &jet_chf, &b_jet_chf);
   fChain->SetBranchAddress("jet_nhf", &jet_nhf, &b_jet_nhf);
   fChain->SetBranchAddress("jet_phf", &jet_phf, &b_jet_phf);
   fChain->SetBranchAddress("jet_missHits", &jet_missHits, &b_jet_missHits);
   fChain->SetBranchAddress("jet_muonHits", &jet_muonHits, &b_jet_muonHits);
   fChain->SetBranchAddress("jet_alphaMax", &jet_alphaMax, &b_jet_alphaMax);
   fChain->SetBranchAddress("jet_nDarkPions", &jet_nDarkPions, &b_jet_nDarkPions);
   fChain->SetBranchAddress("jet_nDarkGluons", &jet_nDarkGluons, &b_jet_nDarkGluons);
   fChain->SetBranchAddress("jet_minDRDarkPion", &jet_minDRDarkPion, &b_jet_minDRDarkPion);
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
